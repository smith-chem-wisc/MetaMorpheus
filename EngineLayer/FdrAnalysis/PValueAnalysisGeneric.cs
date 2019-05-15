using EngineLayer.FdrAnalysis;
using Microsoft.ML;
using Microsoft.ML.Data;
using System;
using System.Collections.Generic;
using System.Linq;
using static Microsoft.ML.DataOperationsCatalog;

namespace EngineLayer
{
    public static class PValueAnalysisGeneric
    {
        public static void ComputePValuesForAllPSMsGeneric(List<PeptideSpectralMatch> psms)
        {
            Dictionary<string, int> accessionAppearances = GetAccessionCounts(psms);

            MLContext mlContext = new MLContext();
            IDataView dataView = mlContext.Data.LoadFromEnumerable(CreatePsmData(psms, accessionAppearances));

            TrainTestData trainTestSplit = mlContext.Data.TrainTestSplit(dataView, testFraction: 0.1);
            IDataView trainingData = trainTestSplit.TrainSet;
            IDataView testData = trainTestSplit.TestSet;

            var trainer = mlContext.BinaryClassification.Trainers.FastTree(labelColumnName: "Label", featureColumnName: "Features");

            var pipeline = mlContext.Transforms.Concatenate("Features", "Intensity", "ScanPrecursorCharge", "DeltaScore", "Notch", "PsmCount", "ModsCount", "MissedCleavagesCount", "Ambiguity", "AccessionAppearances", "LongestFragmentIonSeries")
                .Append(mlContext.BinaryClassification.Trainers.FastTree(labelColumnName: "Label", featureColumnName: "Features"));
            var trainedModel = pipeline.Fit(trainingData);

            var predictionEngine = mlContext.Model.CreatePredictionEngine<PsmData, TruePositivePrediction>(trainedModel);

            foreach (PeptideSpectralMatch psm in psms)
            {
                if (psm != null)
                {
                    var prediction = predictionEngine.Predict(CreateOnePsmDataFromPsm(psm, accessionAppearances));
                    psm.pValueInfo = prediction.Prediction + "|" + prediction.Probability + "|" + prediction.Score;
                }
            }

            var predictions = trainedModel.Transform(testData);

            CalibratedBinaryClassificationMetrics metrics;
            try
            {
                metrics = mlContext.BinaryClassification.Evaluate(data: predictions, labelColumnName: "Label", scoreColumnName: "Score");
            }
            catch
            {

            }

            //TODO: add binary classification metrics to results.tsv
            //PrintBinaryClassificationMetrics(trainer.ToString(), metrics);

            //if you want to save a model, you can use this example
            //mlContext.Model.Save(trainedModel, trainingData.Schema, @"C:\Users\User\Downloads\TrainedModel.zip");
        }

        public static Dictionary<string, int> GetAccessionCounts(List<PeptideSpectralMatch> psms)
        {
            Dictionary<string, int> accessionCountDictionary = new Dictionary<string, int>();

            foreach (PeptideSpectralMatch psm in psms)
            {
                if (psm != null)
                {
                    var firstPeptide = psm.BestMatchingPeptides.Select(p => p.Peptide).First();
                    string accession = firstPeptide.Protein.Accession;
                    if (accession != null)
                    {
                        if (accessionCountDictionary.Keys.Contains(accession))
                        {
                            accessionCountDictionary[accession]++;
                        }
                        else
                        {
                            accessionCountDictionary.Add(accession, 1);
                        }
                    }
                }
            }

            return accessionCountDictionary;
        }


        public static Tuple<int, int> PeaksInScoreHistogram(IEnumerable<double> scores)
        {
            List<int> scoreHistogram = new List<int>();
            int maxScore = (int)scores.Max();
            for (int i = 0; i <= maxScore; i++)
            {
                scoreHistogram.Add(scores.Where(s => s >= i && s < (i + 1)).Count());
            }

            int lowMax = scoreHistogram.IndexOf(scoreHistogram.Max());
            int hiMax = lowMax + 1;

            int sum = 0;
            int total = scoreHistogram.Sum();
            for (int i = 0; i < scoreHistogram.Count; i++)
            {
                sum += scoreHistogram[i];
                double ratio = Convert.ToDouble(sum) / Convert.ToDouble(total);
                if (ratio >= 0.9)
                {
                    hiMax = i;
                    break;
                }
            }

            return new Tuple<int, int>(lowMax, hiMax);
        }

        public static IEnumerable<PsmData> CreatePsmData(List<PeptideSpectralMatch> psms, Dictionary<string, int> accessionAppearances, bool? trueOrFalse = null)
        {
            List<PsmData> pd = new List<PsmData>();
            foreach (PeptideSpectralMatch psm in psms)
            {
                bool label;
                if (trueOrFalse != null)
                {
                    label = trueOrFalse.Value;
                }
                else if (psm.IsDecoy || psm.FdrInfo.QValue > 0.25)
                {
                    label = false;
                    pd.Add(CreateOnePsmDataFromPsm(psm, accessionAppearances, label));
                }
                else if (!psm.IsDecoy && psm.FdrInfo.QValue <= 0.01)
                {
                    label = true;
                    pd.Add(CreateOnePsmDataFromPsm(psm, accessionAppearances, label));
                }
            }
            return pd.AsEnumerable();
        }

        public static PsmData CreateOnePsmDataFromPsm(PeptideSpectralMatch psm, Dictionary<string, int> accessionCounts, bool? trueOrFalse = null)
        {
            //todo: consider adding a count for the number of times a peptides protein accession appears in the list. proteins with more psms should be favored

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

            //todo: for non-specific cleavage, ignore missed cleavages
            float missedCleavages = firstPeptide.MissedCleavages;
            float longestSeq = psm.FdrInfo.LongestSeriesLength;
            string accession = firstPeptide.Protein.Accession;
            float appearances = accessionCounts[accession];
            float score = (float)psm.Score;

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
                AccessionAppearances = appearances,
                Label = label
            };
        }

        public static void PrintBinaryClassificationMetrics(string name, Microsoft.ML.Data.CalibratedBinaryClassificationMetrics metrics)
        {
            Console.WriteLine($"************************************************************");
            Console.WriteLine($"*       Metrics for {name} binary classification model      ");
            Console.WriteLine($"*-----------------------------------------------------------");
            Console.WriteLine($"*       Accuracy: {metrics.Accuracy:P2}");
            Console.WriteLine($"*       Area Under Curve:      {metrics.AreaUnderRocCurve:P2}");
            Console.WriteLine($"*       Area under Precision recall Curve:  {metrics.AreaUnderPrecisionRecallCurve:P2}");
            Console.WriteLine($"*       F1Score:  {metrics.F1Score:P2}");
            Console.WriteLine($"*       LogLoss:  {metrics.LogLoss:#.##}");
            Console.WriteLine($"*       LogLossReduction:  {metrics.LogLossReduction:#.##}");
            Console.WriteLine($"*       PositivePrecision:  {metrics.PositivePrecision:#.##}");
            Console.WriteLine($"*       PositiveRecall:  {metrics.PositiveRecall:#.##}");
            Console.WriteLine($"*       NegativePrecision:  {metrics.NegativePrecision:#.##}");
            Console.WriteLine($"*       NegativeRecall:  {metrics.NegativeRecall:P2}");
            Console.WriteLine($"************************************************************");
        }

    }
}