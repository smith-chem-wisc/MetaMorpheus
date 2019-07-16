using EngineLayer.FdrAnalysis;
using MathNet.Numerics.Statistics;
using Microsoft.ML;
using Microsoft.ML.Data;
using Proteomics.ProteolyticDigestion;
using Proteomics.RetentionTimePrediction;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using static Microsoft.ML.DataOperationsCatalog;

namespace EngineLayer
{
    public static class PEP_Analysis
    {
        private static readonly string filepathSubstitue = "SubstitueFilePath";//this is only for unit test that have psms with no filename.

        public static string ComputePEPValuesForAllPSMsGeneric(List<PeptideSpectralMatch> psms)
        {
            Dictionary<string, Dictionary<int, double>> RetentionTimeHydrophobicityDeviation = new Dictionary<string, Dictionary<int, double>>();
            Dictionary<string, Dictionary<int, double>> RetentionTimeHydrophobicityAverage = new Dictionary<string, Dictionary<int, double>>();

            Dictionary<string, Dictionary<int, double>> RetentionTimeHydrophobicityDeviation_M = new Dictionary<string, Dictionary<int, double>>();
            Dictionary<string, Dictionary<int, double>> RetentionTimeHydrophobicityAverage_M = new Dictionary<string, Dictionary<int, double>>();

            (RetentionTimeHydrophobicityAverage, RetentionTimeHydrophobicityDeviation) = FillHydrophobicityDictionaries(psms, false);
            (RetentionTimeHydrophobicityAverage_M, RetentionTimeHydrophobicityDeviation_M) = FillHydrophobicityDictionaries(psms, true);

            Dictionary<string, int> sequenceToPsmCount = GetSequenceToPSMCount(psms);

            MLContext mlContext = new MLContext();
            IDataView dataView = mlContext.Data.LoadFromEnumerable(CreatePsmData(psms, sequenceToPsmCount, RetentionTimeHydrophobicityAverage, RetentionTimeHydrophobicityDeviation, RetentionTimeHydrophobicityAverage_M, RetentionTimeHydrophobicityDeviation_M));

            //
            // Summary:
            //     Split the dataset into the train set and test set according to the given fraction.
            //     Respects the samplingKeyColumnName if provided.
            //
            // Parameters:
            //   data:
            //     The dataset to split.
            //
            //   testFraction:
            //     The fraction of data to go into the test set.
            //
            //   samplingKeyColumnName:
            //     Name of a column to use for grouping rows. If two examples share the same value
            //     of the samplingKeyColumnName, they are guaranteed to appear in the same subset
            //     (train or test). This can be used to ensure no label leakage from the train to
            //     the test set. If null no row grouping will be performed.
            //
            //   seed:
            //     Seed for the random number generator used to select rows for the train-test split.

            TrainTestData trainTestSplit = mlContext.Data.TrainTestSplit(dataView, testFraction: 0.1);
            IDataView trainingData = trainTestSplit.TrainSet;
            IDataView testData = trainTestSplit.TestSet;

            var trainer = mlContext.BinaryClassification.Trainers.FastTree(labelColumnName: "Label", featureColumnName: "Features");

            string[] customFeatures = GetCustomFeatures();

            var pipeline = mlContext.Transforms.Concatenate("Features", customFeatures)
                .Append(mlContext.BinaryClassification.Trainers.FastTree(labelColumnName: "Label", featureColumnName: "Features"));

            var trainedModel = pipeline.Fit(trainingData);

            var predictionEngine = mlContext.Model.CreatePredictionEngine<PsmData, TruePositivePrediction>(trainedModel);

            string ambiguousScans = "";

            //For Debug
            List<string> someOut = new List<string>();
            someOut.Add("Z|Ambiguity|DeltaScore|Intensity|Label|LongestSeries|MissedCleavages|ModsCount|Notch|PsmCount|PrecursorCharge|call|pepValue|Score|QValue");

            foreach (PeptideSpectralMatch psm in psms)
            {
                if (psm != null)
                {
                    List<int> indiciesOfPeptidesToRemove = new List<int>();
                    List<(int notch, PeptideWithSetModifications pwsm)> bestMatchingPeptidesToRemove = new List<(int notch, PeptideWithSetModifications pwsm)>();
                    List<double> pepValuePredictions = new List<double>();

                    //Here we compute the pepvalue predection for each ambiguous peptide in a PSM. Ambiguous peptides with lower pepvalue predictions are removed from the PSM.
                    foreach (var (Notch, Peptide) in psm.BestMatchingPeptides)
                    {
                        PsmData pd = CreateOnePsmDataEntry(psm, sequenceToPsmCount, RetentionTimeHydrophobicityAverage, RetentionTimeHydrophobicityDeviation, RetentionTimeHydrophobicityAverage_M, RetentionTimeHydrophobicityDeviation_M, Peptide, Notch);
                        double z;
                        if (Peptide.BaseSequence.Equals(Peptide.FullSequence))
                        {
                            z = GetSSRCalcHydrophobicityZScore(psm, Peptide, RetentionTimeHydrophobicityAverage, RetentionTimeHydrophobicityDeviation);
                        }
                        else
                        {
                            z = GetSSRCalcHydrophobicityZScore(psm, Peptide, RetentionTimeHydrophobicityAverage_M, RetentionTimeHydrophobicityDeviation_M);
                        }

                        var pepValuePrediction = predictionEngine.Predict(pd);

                        someOut.Add(z.ToString() + "|" + pd.Ambiguity.ToString() + "|" + pd.DeltaScore.ToString() + "|" + pd.Intensity.ToString() + "|" + pd.Label + "|" + pd.LongestFragmentIonSeries + "|" + pd.MissedCleavagesCount + "|" + pd.ModsCount + "|" + pd.Notch + "|" + pd.PsmCount + "|" + pd.ScanPrecursorCharge + "|" + pepValuePrediction.Prediction + "|" + pepValuePrediction.Probability + "|" + pepValuePrediction.Score);

                        pepValuePredictions.Add(pepValuePrediction.Probability);
                        //A score is available using the variable pepvaluePrediction.Score
                    }

                    double highestPredictedPEPValue = pepValuePredictions.Max();
                    int numberOfPredictions = pepValuePredictions.Count - 1;

                    for (int i = numberOfPredictions; i >= 0; i--)
                    {
                        if (Math.Abs(highestPredictedPEPValue - pepValuePredictions[i]) > 0.000001)
                        {
                            indiciesOfPeptidesToRemove.Add(i);
                            pepValuePredictions.RemoveAt(i);
                            //pValuePredictionStrings.RemoveAt(i);
                        }
                    }

                    int index = 0;

                    foreach (var (Notch, Peptide) in psm.BestMatchingPeptides)
                    {
                        if (indiciesOfPeptidesToRemove.Contains(index))
                        {
                            bestMatchingPeptidesToRemove.Add((Notch, Peptide));
                        }
                        index++;
                    }

                    foreach (var (notch, pwsm) in bestMatchingPeptidesToRemove)
                    {
                        ambiguousScans = ambiguousScans + psm.ScanNumber + "|";
                        psm.RemoveThisAmbiguousPeptide(notch, pwsm);
                    }

                    psm.FdrInfo.PEP = 1 - pepValuePredictions[0]; //they should all be the same at this point so it doesn't matter which you take. First is good.
                }
            }

            //For debug
            File.WriteAllLines(@"C:\Users\Michael Shortreed\Downloads\psmDataVAlues.txt", someOut, System.Text.Encoding.UTF8);

            var predictions = trainedModel.Transform(testData);

            CalibratedBinaryClassificationMetrics metrics;
            try
            {
                metrics = mlContext.BinaryClassification.Evaluate(data: predictions, labelColumnName: "Label", scoreColumnName: "Score");
                return PrintBinaryClassificationMetrics(trainer.ToString(), metrics);
            }
            catch
            {
                return "";
            }

            //if you want to save a model, you can use this example
            //mlContext.Model.Save(trainedModel, trainingData.Schema, @"C:\Users\User\Downloads\TrainedModel.zip");
        }

        /// <summary>
        /// This method can be used to select a custom subset of features for use in the Microsoft ML.NET computation of posterior error probability
        /// Probably should have unique sets for regular, crosslinking, top-down, labelled.
        /// </summary>
        /// <returns></returns>
        private static string[] GetCustomFeatures()
        {
            //non-specific and top-down searches don't used missedCleavages
            return new string[] { "Z", "Intensity", "ScanPrecursorCharge", "DeltaScore", "Notch", "PsmCount", "ModsCount", "MissedCleavagesCount", "Ambiguity", "LongestFragmentIonSeries" };
        }

        private static float GetSSRCalcHydrophobicityZScore(PeptideSpectralMatch psm, PeptideWithSetModifications Peptide, Dictionary<string, Dictionary<int, double>> avg, Dictionary<string, Dictionary<int, double>> dev)
        {
            SSRCalc3 calc = new SSRCalc3("SSRCalc 3.0 (300A)", SSRCalc3.Column.A300);
            double z = double.NaN;

            if (avg.ContainsKey(psm.FullFilePath ?? filepathSubstitue) && dev.ContainsKey(psm.FullFilePath ?? filepathSubstitue))
            {
                int time = (int)Math.Round(psm.ScanRetentionTime, 0);
                if (avg[psm.FullFilePath ?? filepathSubstitue].ContainsKey(time))
                {
                    double predictedHydrophobicity = calc.ScoreSequence(Peptide);
                    z = Math.Abs(avg[psm.FullFilePath ?? filepathSubstitue][time] - predictedHydrophobicity) / dev[psm.FullFilePath ?? filepathSubstitue][time];
                }
            }

            if (double.IsNaN(z) || double.IsInfinity(z))
            {
                z = 100;
            }

            return (float)z;
        }

        private static (Dictionary<string, Dictionary<int, double>>, Dictionary<string, Dictionary<int, double>>) FillHydrophobicityDictionaries(List<PeptideSpectralMatch> psms, bool modified)
        {
            SSRCalc3 calc = new SSRCalc3("SSRCalc 3.0 (300A)", SSRCalc3.Column.A300);
            Dictionary<string, Dictionary<int, double>> localRetentionTimeHydrophobicityDeviation = new Dictionary<string, Dictionary<int, double>>();
            Dictionary<string, Dictionary<int, double>> localRetentionTimeHydrophobicityAverage = new Dictionary<string, Dictionary<int, double>>();

            List<string> filenames = psms.Select(f => f.FullFilePath).ToList();

            filenames = filenames.Distinct().ToList();

            foreach (string filename in filenames)
            {
                Dictionary<int, List<double>> hydrobophobicites = new Dictionary<int, List<double>>();
                Dictionary<int, double> averages = new Dictionary<int, double>();
                Dictionary<int, double> deviations = new Dictionary<int, double>();

                foreach (PeptideSpectralMatch psm in psms.Where(f => (f.FullFilePath == filename || f.FullFilePath == null) && f.FdrInfo.QValue <= 0.01))
                {
                    foreach ((int notch, PeptideWithSetModifications pwsm) in psm.BestMatchingPeptides)
                    {
                        if (pwsm.BaseSequence.Equals(pwsm.FullSequence) && !modified)
                        {
                            double predictedHydrophobicity = calc.ScoreSequence(pwsm);
                            int possibleKey = (int)Math.Round(psm.ScanRetentionTime, 0);
                            if (hydrobophobicites.ContainsKey(possibleKey))
                            {
                                hydrobophobicites[possibleKey].Add(predictedHydrophobicity);
                            }
                            else
                            {
                                hydrobophobicites.Add(possibleKey, new List<double>() { predictedHydrophobicity });
                            }
                        }
                        else if (!pwsm.BaseSequence.Equals(pwsm.FullSequence) && modified)
                        {
                            double predictedHydrophobicity = calc.ScoreSequence(pwsm);
                            int possibleKey = (int)Math.Round(psm.ScanRetentionTime, 0);
                            if (hydrobophobicites.ContainsKey(possibleKey))
                            {
                                hydrobophobicites[possibleKey].Add(predictedHydrophobicity);
                            }
                            else
                            {
                                hydrobophobicites.Add(possibleKey, new List<double>() { predictedHydrophobicity });
                            }
                        }
                    }
                }

                foreach (int key in hydrobophobicites.Keys)
                {
                    averages.Add(key, hydrobophobicites[key].Average());
                    deviations.Add(key, hydrobophobicites[key].StandardDeviation());
                }

                localRetentionTimeHydrophobicityAverage.Add(filename ?? filepathSubstitue, averages);
                localRetentionTimeHydrophobicityDeviation.Add(filename ?? filepathSubstitue, deviations);
            }
            return (localRetentionTimeHydrophobicityAverage, localRetentionTimeHydrophobicityDeviation);
        }

        private static Dictionary<string, int> GetSequenceToPSMCount(List<PeptideSpectralMatch> psms)
        {
            Dictionary<string, int> sequenceToPsmCount = new Dictionary<string, int>();

            List<string> sequences = new List<string>();
            foreach (PeptideSpectralMatch psm in psms)
            {
                var ss = psm.BestMatchingPeptides.Select(b => b.Peptide.FullSequence).ToList();
                sequences.Add(String.Join("|", ss));
            }

            var s = sequences.GroupBy(i => i);

            foreach (var grp in s)
            {
                sequenceToPsmCount.Add(grp.Key, grp.Count());
            }
            return sequenceToPsmCount;
        }

        //This method ignores ambiguity and loads only the first peptide in a series for each PSM
        public static IEnumerable<PsmData> CreatePsmData(List<PeptideSpectralMatch> psms, Dictionary<string, int> sequenceToPsmCount, Dictionary<string, Dictionary<int, double>> avg, Dictionary<string, Dictionary<int, double>> dev, Dictionary<string, Dictionary<int, double>> avg_M, Dictionary<string, Dictionary<int, double>> dev_M, bool? trueOrFalse = null)
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
                    pd.Add(CreateOnePsmDataEntry(psm, sequenceToPsmCount, avg, dev, avg_M, dev_M, null, null, label));
                }
                else if (!psm.IsDecoy && psm.FdrInfo.QValue <= 0.01)
                {
                    label = true;
                    pd.Add(CreateOnePsmDataEntry(psm, sequenceToPsmCount, avg, dev, avg_M, dev_M, null, null, label));
                }
            }
            return pd.AsEnumerable();
        }

        /// <summary>
        ///
        /// </summary>
        /// <param name="psm"></param>
        /// <param name="sequenceToPsmCount"></param>
        /// <param name="selectedPeptide"></param>
        /// <param name="notchToUse"></param>
        /// <param name="trueOrFalse"></param>
        /// <returns></returns>
        public static PsmData CreateOnePsmDataEntry(PeptideSpectralMatch psm, Dictionary<string, int> sequenceToPsmCount, Dictionary<string, Dictionary<int, double>> avg, Dictionary<string, Dictionary<int, double>> dev, Dictionary<string, Dictionary<int, double>> avg_M, Dictionary<string, Dictionary<int, double>> dev_M, PeptideWithSetModifications selectedPeptide, int? notchToUse, bool? trueOrFalse = null)
        {
            float ambiguity = (float)psm.PeptidesToMatchingFragments.Keys.Count;
            float intensity = (float)(psm.Score - (int)psm.Score);
            float charge = psm.ScanPrecursorCharge;
            float deltaScore = (float)psm.DeltaScore;
            float psmCount = sequenceToPsmCount[String.Join("|", psm.BestMatchingPeptides.Select(p => p.Peptide.FullSequence).ToList())];
            int notch = 0;
            if (notchToUse.HasValue)
            {
                notch = notchToUse.Value;
            }
            else if (psm.Notch.HasValue)
            {
                notch = psm.Notch.Value;
            }

            if (selectedPeptide == null)
            {
                selectedPeptide = psm.BestMatchingPeptides.Select(p => p.Peptide).First();
            }

            float modCount = selectedPeptide.AllModsOneIsNterminus.Keys.Count();
            float missedCleavages = selectedPeptide.MissedCleavages;
            float longestSeq = psm.GetLongestIonSeriesBidirectional(selectedPeptide);

            float z;

            if (selectedPeptide.BaseSequence.Equals(selectedPeptide.FullSequence))
            {
                z = GetSSRCalcHydrophobicityZScore(psm, selectedPeptide, avg, dev);
            }
            else
            {
                z = GetSSRCalcHydrophobicityZScore(psm, selectedPeptide, avg_M, dev_M);
            }

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
                Z = z,
                Label = label
            };
        }

        /// <summary>
        /// At the time when the ~10% of the data gets chosen for training, another 10% gets chosen for evaluation. Then after training,
        /// the effectiveness of the model gets evaluated on the test set. The results of that evaluation are converted to text values called
        /// BinarySearchTreeMetrics and this gets written to the results.tsv
        /// </summary>
        /// <param name="name"></param>
        /// <param name="metrics"></param>
        /// <returns></returns>
        public static string PrintBinaryClassificationMetrics(string name, CalibratedBinaryClassificationMetrics metrics)
        {
            StringBuilder s = new StringBuilder();
            s.AppendLine("************************************************************");
            s.AppendLine("*       Metrics for Determination of PEP Using Binary Classification      ");
            s.AppendLine("*-----------------------------------------------------------");
            s.AppendLine("*       Accuracy:  " + metrics.Accuracy.ToString());
            s.AppendLine("*       Area Under Curve:  " + metrics.AreaUnderRocCurve.ToString());
            s.AppendLine("*       Area under Precision recall Curve:  " + metrics.AreaUnderPrecisionRecallCurve.ToString());
            s.AppendLine("*       F1Score:  " + metrics.F1Score.ToString());
            s.AppendLine("*       LogLoss:  " + metrics.LogLoss.ToString());
            s.AppendLine("*       LogLossReduction:  " + metrics.LogLossReduction.ToString());
            s.AppendLine("*       PositivePrecision:  " + metrics.PositivePrecision.ToString());
            s.AppendLine("*       PositiveRecall:  " + metrics.PositiveRecall.ToString());
            s.AppendLine("*       NegativePrecision:  " + metrics.NegativePrecision.ToString());
            s.AppendLine("*       NegativeRecall:  " + metrics.NegativeRecall.ToString());
            s.AppendLine("************************************************************");
            return s.ToString();
        }
    }
}