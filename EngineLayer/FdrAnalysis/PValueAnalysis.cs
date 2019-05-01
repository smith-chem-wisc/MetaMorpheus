using Microsoft.ML;
using System;
using System.Collections.Generic;
using System.Linq;

namespace EngineLayer.FdrAnalysis
{
    public static class PValueAnalysis
    {
        public static void ComputePValuesForAllPSMs(List<PeptideSpectralMatch> psms, string modelPath)
        {
            MLContext mlContext = new MLContext();
            if (modelPath != "")
            {
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

        public static IEnumerable<PsmData> GetTrainingSet(List<PeptideSpectralMatch> psms)
        {
            List<PsmData> trainingSetOfPsms = new List<PsmData>();

            int numNeeded = (int)(psms.Count / 100);
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

            List<PeptideSpectralMatch> trueTrainingPsmsFromTsv = new List<PeptideSpectralMatch>();
            List<PeptideSpectralMatch> falseTrainingPsmsFromTsv = new List<PeptideSpectralMatch>();

            while (targetCount < numNeeded || decoyCount < numNeeded)
            {
                if (psms[randomIndexList[theIndex]].Score > 7 && !psms[randomIndexList[theIndex]].IsDecoy && targetCount < numNeeded)
                {
                    trueTrainingPsmsFromTsv.Add(psms[randomIndexList[theIndex]]);
                    targetCount++;
                }
                else if (psms[randomIndexList[theIndex]].Score > 3 && psms[randomIndexList[theIndex]].Score < 6 && decoyCount < numNeeded)
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
                    targetCount = numNeeded;
                    decoyCount = numNeeded;
                }
            }

            trainingSetOfPsms.AddRange(CreatePsmData(trueTrainingPsmsFromTsv, true));
            trainingSetOfPsms.AddRange(CreatePsmData(falseTrainingPsmsFromTsv, false));

            return trainingSetOfPsms.AsEnumerable();
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