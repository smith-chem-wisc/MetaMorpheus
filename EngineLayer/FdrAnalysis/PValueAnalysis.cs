using Microsoft.ML;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace EngineLayer.FdrAnalysis
{
    public static class PValueAnalysis
    {
        public static void ComputePValuesForAllPSMs(List<PeptideSpectralMatch> psms, bool useModel)
        {
            MLContext mlContext = new MLContext();
            if (psms.Count() < 50000 || useModel)
            {
                string modelPath = Path.Combine(GlobalVariables.DataDir, "Data", @"pValueUnitTestTrainedModel.zip");
                ITransformer trainedModel = mlContext.Model.Load(modelPath, out var modelInputSchema);

                //var pipeline = mlContext.Transforms.Concatenate("Features", "Intensity", "ScanPrecursorCharge", "DeltaScore", "Notch", "PsmCount", "ModsCount", "MissedCleavagesCount", "Ambiguity", "LongestFragmentIonSeries").Append(mlContext.BinaryClassification.Trainers.FastTree(labelColumnName: "Label", featureColumnName: "Features"));

                var predictionEngine = mlContext.Model.CreatePredictionEngine<PsmData, TruePositivePrediction>(trainedModel);

                foreach (PeptideSpectralMatch psm in psms)
                {
                    if (psm != null)
                    {
                        var prediction = predictionEngine.Predict(CreateOnePsmDataFromPsm(psm));
                        psm.pValueInfo = prediction.Prediction + "|" + prediction.Probability + "|" + prediction.Score;
                    }
                }
            }
            else
            {
                IDataView trainingData = mlContext.Data.LoadFromEnumerable(GetTrainingSet(psms));
                var pipeline = mlContext.Transforms.Concatenate("Features", "Intensity", "ScanPrecursorCharge", "DeltaScore", "Notch", "PsmCount", "ModsCount", "MissedCleavagesCount", "Ambiguity", "LongestFragmentIonSeries")
                    .Append(mlContext.BinaryClassification.Trainers.FastTree(labelColumnName: "Label", featureColumnName: "Features"));
                var trainedModel = pipeline.Fit(trainingData);
                var predictionEngine = mlContext.Model.CreatePredictionEngine<PsmData, TruePositivePrediction>(trainedModel);

                foreach (PeptideSpectralMatch psm in psms)
                {
                    if (psm != null)
                    {
                        var prediction = predictionEngine.Predict(CreateOnePsmDataFromPsm(psm));
                        psm.pValueInfo = prediction.Prediction + "|" + prediction.Probability + "|" + prediction.Score;
                    }
                }

                //if you want to save a model, you can use this example
                //mlContext.Model.Save(trainedModel, trainingData.Schema, @"C:\Users\User\Downloads\TrainedModel.zip");
            }
        }

        public static IEnumerable<PsmData> GetTrainingSet(List<PeptideSpectralMatch> psms, int? numNeeded = null)
        {
            List<PsmData> trainingSetOfPsms = new List<PsmData>();

            if(numNeeded == null)
            {
                numNeeded = (int)(psms.Count / 100d);
            }

            List<int> originalIndex = new List<int>();
            originalIndex.AddRange(Enumerable.Range(0, psms.Count - 1));
            List<int> randomIndexList = new List<int>();

            Random r = new Random();
            int randomIndex = 0;
            while (originalIndex.Count > 0)
            {
                randomIndex = r.Next(0, originalIndex.Count); //Choose a random object in the list
                randomIndexList.Add(originalIndex[randomIndex]); //add it to the new, random list
                originalIndex.RemoveAt(randomIndex); //remove to avoid duplicates
            }

            int targetCount = 0;
            int decoyCount = 0;
            int theIndex = 0;

            Tuple<int, int> decoyTargetMaxima = PeaksInScoreHistogram(psms.Select(s => s.Score));

            List<PeptideSpectralMatch> trueTrainingPsmsFromTsv = new List<PeptideSpectralMatch>();
            List<PeptideSpectralMatch> falseTrainingPsmsFromTsv = new List<PeptideSpectralMatch>();

            if(decoyTargetMaxima.Item1 < decoyTargetMaxima.Item2)
            {
                while (targetCount < numNeeded.Value || decoyCount < numNeeded.Value)
                {
                    if (psms[randomIndexList[theIndex]].Score > (decoyTargetMaxima.Item2) && !psms[randomIndexList[theIndex]].IsDecoy && targetCount < numNeeded.Value)
                    {
                        trueTrainingPsmsFromTsv.Add(psms[randomIndexList[theIndex]]);
                        targetCount++;
                    }
                    else if (psms[randomIndexList[theIndex]].Score > (decoyTargetMaxima.Item1 - 1) && psms[randomIndexList[theIndex]].Score < (decoyTargetMaxima.Item1 + 1) && decoyCount < numNeeded.Value)
                    {
                        falseTrainingPsmsFromTsv.Add(psms[randomIndexList[theIndex]]);
                        decoyCount++;
                    }
                    if (theIndex < (psms.Count - 2))
                    {
                        theIndex++;
                    }
                    else
                    {
                        targetCount = numNeeded.Value;
                        decoyCount = numNeeded.Value;
                    }
                }
            }
           
            trainingSetOfPsms.AddRange(CreatePsmData(trueTrainingPsmsFromTsv, true));
            trainingSetOfPsms.AddRange(CreatePsmData(falseTrainingPsmsFromTsv, false));

            return trainingSetOfPsms.AsEnumerable();
        }

        public static Tuple<int,int> PeaksInScoreHistogram(IEnumerable<double> scores)
        {
            List<int> scoreHistogram = new List<int>();
            int maxScore = (int)scores.Max();
            for (int i = 0; i <= maxScore; i++)
            {
                scoreHistogram.Add(scores.Where(s => s >= i && s < (i + 1)).Count());
            }

            int lowMax = 0;
            for (int i = 1; i < scoreHistogram.Count; i++)
            {
                if (scoreHistogram[i] >= scoreHistogram[i - 1])
                {
                    lowMax++;
                }
                else
                {
                    break;
                }
            }

            int hiMax = maxScore;
            for (int i = maxScore-1; i >= 0; i--)
            {
                if (scoreHistogram[i] >= scoreHistogram[i + 1])
                {
                    hiMax--;
                }
                else
                {
                    break;
                }
            }

            if(hiMax == lowMax )
            {
                if((lowMax + 3) < maxScore)
                {
                    hiMax = lowMax + 3;
                }
            }
            else
            {
                hiMax = Math.Max(lowMax + 3, (int)(lowMax + (hiMax - lowMax) / 10));
            }

            return new Tuple<int, int>(lowMax,hiMax);
        }
        public static IEnumerable<PsmData> CreatePsmData(List<PeptideSpectralMatch> psms, bool? trueOrFalse = null)
        {
            List<PsmData> pd = new List<PsmData>();
            foreach (PeptideSpectralMatch psm in psms)
            {
                bool label;
                if (trueOrFalse != null)
                {
                    label = trueOrFalse.Value;
                }
                else if (psm.IsDecoy)
                {
                    label = false;
                }
                else
                {
                    label = true;
                }

                //TODO deal with ambiguous peptides
                pd.Add(CreateOnePsmDataFromPsm(psm, label));
            }
            return pd.AsEnumerable();
        }

        public static PsmData CreateOnePsmDataFromPsm(PeptideSpectralMatch psm, bool? trueOrFalse = null)
        {
            float ambiguity = (float)psm.PeptidesToMatchingFragments.Count;//(psm.BaseSequence.Split('|').Count());
            float intensity = (float)(psm.Score - (int)psm.Score);
            float charge = psm.ScanPrecursorCharge;
            float deltaScore = (float)psm.DeltaScore;

            float notch = 0;
            if (psm.Notch.HasValue)
            {
                notch = psm.Notch.Value;
            }

            float psmCount = Convert.ToInt32(psm.PsmCount);
            var firstPeptide = psm.BestMatchingPeptides.Select(p => p.Peptide).First();
            float modCount = firstPeptide.AllModsOneIsNterminus.Values.Count();
            float missedCleavages = firstPeptide.MissedCleavages;
            float longestSeq = psm.FdrInfo.LongestSeriesLength;

            bool label;
            if (trueOrFalse != null)
            {
                label = trueOrFalse.Value;
            }
            else if (psm.IsDecoy)
            {
                label = false;
            }
            else
            {
                label = true;
            }

            return new PsmData()
            {
                Intensity = intensity,
                ScanPrecursorCharge = charge,
                DeltaScore = deltaScore,
                Notch = notch,
                PsmCount = psmCount,
                ModsCount = modCount,
                MissedCleavagesCount = missedCleavages,
                Ambiguity = ambiguity,
                LongestFragmentIonSeries = longestSeq,
                Label = label
            };
        }
    }
}