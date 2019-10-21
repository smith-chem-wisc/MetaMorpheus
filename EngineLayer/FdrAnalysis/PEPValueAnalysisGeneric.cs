using EngineLayer.FdrAnalysis;
using MathNet.Numerics.Statistics;
using Microsoft.ML;
using Microsoft.ML.Data;
using Proteomics.ProteolyticDigestion;
using Proteomics.RetentionTimePrediction;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using static Microsoft.ML.DataOperationsCatalog;

namespace EngineLayer
{
    public static class PEP_Analysis
    {
        private static readonly double AbsoluteProbabilityThatDistinguishesPeptides = 0.05;

        public static string ComputePEPValuesForAllPSMsGeneric(List<PeptideSpectralMatch> psms)
        {
            string searchType = DetermineSearchType(psms);
            string[] trainingVariables = PsmData.trainingInfos[searchType];

            //These two dictionaries contain the average and standard deviations of hydrophobicitys measured in 1 minute increments accross each raw
            //file separately. An individully measured hydrobophicty calculated for a specific PSM sequence is compared to these values by computing
            //the z-score. That z-score is used in the the machine learning.
            //Separate dictionaries are created for peptides with modifications becuase SSRcalc doesn't really do a good job predicting hyrophobicity

            //The first string in the dictionary is the filename
            //The value of the dictionary is another dictionary that profiles the hydrophobicity behavior. Each key is a retention time rounded to the nearest minute. The value Tuple is the average and standard deviation, respectively, of the predicted hydrophobicities of the observed peptides eluting at that rounded retention time.
            Dictionary<string, Dictionary<int, Tuple<double, double>>> fileSpecificTimeDependantHydrophobicityAverageAndDeviation_unmodified = new Dictionary<string, Dictionary<int, Tuple<double, double>>>();
            Dictionary<string, Dictionary<int, Tuple<double, double>>> fileSpecificTimeDependantHydrophobicityAverageAndDeviation_modified = new Dictionary<string, Dictionary<int, Tuple<double, double>>>();

            if (trainingVariables.Contains("HydrophobicityZScore"))
            {
                fileSpecificTimeDependantHydrophobicityAverageAndDeviation_unmodified = ComputeHydrophobicityValues(psms, false);
                fileSpecificTimeDependantHydrophobicityAverageAndDeviation_modified = ComputeHydrophobicityValues(psms, true);
            }

            Dictionary<string, int> sequenceToPsmCount = GetSequenceToPSMCount(psms);
            int chargeStateMode = GetChargeStateMode(psms);

            MLContext mlContext = new MLContext();
            IDataView dataView = mlContext.Data.LoadFromEnumerable(CreatePsmData(psms, sequenceToPsmCount, fileSpecificTimeDependantHydrophobicityAverageAndDeviation_unmodified, fileSpecificTimeDependantHydrophobicityAverageAndDeviation_modified, chargeStateMode, trainingVariables));

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
            //     The seed, '42', is not random but fixed for consistancy. According to the supercomputer Deep Thought the answer to the question of life, the universe and everything was 42 (in Douglas Adam’s Hitchhikers Guide to the Galaxy).

            TrainTestData trainTestSplit = mlContext.Data.TrainTestSplit(dataView, testFraction: 0.1, null, 42);
            IDataView trainingData = trainTestSplit.TrainSet;
            IDataView testData = trainTestSplit.TestSet;

            var trainer = mlContext.BinaryClassification.Trainers.FastTree(labelColumnName: "Label", featureColumnName: "Features");

            var pipeline = mlContext.Transforms.Concatenate("Features", trainingVariables)
                .Append(mlContext.BinaryClassification.Trainers.FastTree(labelColumnName: "Label", featureColumnName: "Features"));

            var trainedModel = pipeline.Fit(trainingData);

            var predictionEngine = mlContext.Model.CreatePredictionEngine<PsmData, TruePositivePrediction>(trainedModel);

            string ambiguousScans = "";
            int ambiguousPeptidesRemovedCount = 0;

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
                        PsmData pd = CreateOnePsmDataEntry(psm, sequenceToPsmCount, fileSpecificTimeDependantHydrophobicityAverageAndDeviation_unmodified, fileSpecificTimeDependantHydrophobicityAverageAndDeviation_modified, chargeStateMode, Peptide, trainingVariables, Notch, !Peptide.Protein.IsDecoy);
                        var pepValuePrediction = predictionEngine.Predict(pd);
                        pepValuePredictions.Add(pepValuePrediction.Probability);
                        //A score is available using the variable pepvaluePrediction.Score
                    }

                    double highestPredictedPEPValue = pepValuePredictions.Max();
                    int numberOfPredictions = pepValuePredictions.Count - 1;

                    for (int i = numberOfPredictions; i >= 0; i--)
                    {
                        if (Math.Abs(highestPredictedPEPValue - pepValuePredictions[i]) > AbsoluteProbabilityThatDistinguishesPeptides)
                        {
                            indiciesOfPeptidesToRemove.Add(i);
                            pepValuePredictions.RemoveAt(i);
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
                        ambiguousPeptidesRemovedCount++;
                    }
                    psm.FdrInfo.PEP = 1 - pepValuePredictions.Max();
                }
            }

            var predictions = trainedModel.Transform(testData);

            CalibratedBinaryClassificationMetrics metrics;
            try
            {
                metrics = mlContext.BinaryClassification.Evaluate(data: predictions, labelColumnName: "Label", scoreColumnName: "Score");
                return PrintBinaryClassificationMetrics(trainer.ToString(), metrics, ambiguousPeptidesRemovedCount);
            }
            catch
            {
                return "";
            }

            //if you want to save a model, you can use this example
            //mlContext.Model.Save(trainedModel, trainingData.Schema, @"C:\Users\User\Downloads\TrainedModel.zip");
        }

        /// <summary>
        /// Here we're getting the most common charge state for precursors that are Targets with q<=0.01.

        private static int GetChargeStateMode(List<PeptideSpectralMatch> psms)
        {
            return psms.Where(p => p.IsDecoy != true && p.FdrInfo.QValue <= 0.01).Select(p => p.ScanPrecursorCharge).GroupBy(n => n).OrderByDescending(g => g.Count()).Select(g => g.Key).FirstOrDefault();
        }

        private static Dictionary<string, Dictionary<int, Tuple<double, double>>> ComputeHydrophobicityValues(List<PeptideSpectralMatch> psms, bool computeHydrophobicitiesforModifiedPeptides)
        {
            SSRCalc3 calc = new SSRCalc3("SSRCalc 3.0 (300A)", SSRCalc3.Column.A300);
            Dictionary<string, Dictionary<int, Tuple<double, double>>> rtHydrophobicityAvgDev = new Dictionary<string, Dictionary<int, Tuple<double, double>>>();

            List<string> filenames = psms.Select(f => f.FullFilePath).ToList();

            filenames = filenames.Distinct().ToList();

            foreach (string filename in filenames)
            {
                Dictionary<int, List<double>> hydrobophobicites = new Dictionary<int, List<double>>();
                Dictionary<int, Tuple<double, double>> averagesCommaStandardDeviations = new Dictionary<int, Tuple<double, double>>();

                foreach (PeptideSpectralMatch psm in psms.Where(f => (f.FullFilePath == filename || f.FullFilePath == null) && f.FdrInfo.QValue <= 0.01 && !f.IsDecoy))
                {
                    List<string> fullSequences = new List<string>();
                    foreach ((int notch, PeptideWithSetModifications pwsm) in psm.BestMatchingPeptides)
                    {
                        if (fullSequences.Contains(pwsm.FullSequence))
                        {
                            continue;
                        }
                        fullSequences.Add(pwsm.FullSequence);
                        double predictedHydrophobicity = calc.ScoreSequence(pwsm);
                        int possibleKey = (int)(2 * Math.Round(psm.ScanRetentionTime / 2d, 0));
                        //First block of if statement is for modified peptides.
                        if (pwsm.AllModsOneIsNterminus.Any() && computeHydrophobicitiesforModifiedPeptides)
                        {
                            if (hydrobophobicites.ContainsKey(possibleKey))
                            {
                                hydrobophobicites[possibleKey].Add(predictedHydrophobicity);
                            }
                            else
                            {
                                hydrobophobicites.Add(possibleKey, new List<double>() { predictedHydrophobicity });
                            }
                        }
                        //this second block of if statment is for unmodified peptides.
                        else if (!pwsm.AllModsOneIsNterminus.Any() && !computeHydrophobicitiesforModifiedPeptides)
                        {
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
                    //TODO consider using inner-quartile range instead of standard deviation
                    averagesCommaStandardDeviations.Add(key, new Tuple<double, double>(hydrobophobicites[key].Average(), hydrobophobicites[key].StandardDeviation()));
                }

                //some standard deviations are too small or too large because of random reasons, so we replace those small numbers of oddballs with reasonable numbers.

                List<double> stDevs = averagesCommaStandardDeviations.Select(x => x.Value.Item2).ToList();

                for (int i = stDevs.Count - 1; i >= 0; i--)
                {
                    if (stDevs[i].Equals(double.NaN))
                    {
                        stDevs.RemoveAt(i);
                    }
                }

                double globalStDev;
                if (stDevs.Count > 1)
                {
                    globalStDev = stDevs.Average();
                }
                else if (stDevs.Count == 1)
                {
                    globalStDev = stDevs.First();
                }
                else
                {
                    globalStDev = 1;
                }

                Dictionary<int, Tuple<double, double>> stuffToChange = new Dictionary<int, Tuple<double, double>>();
                foreach (KeyValuePair<int, Tuple<double, double>> item in averagesCommaStandardDeviations)
                {
                    if (Double.IsNaN(item.Value.Item2) || item.Value.Item2 < 0.5 || (item.Value.Item2 / globalStDev) > 3)
                    {
                        Tuple<double, double> pair = new Tuple<double, double>(averagesCommaStandardDeviations[item.Key].Item1, globalStDev);
                        stuffToChange.Add(item.Key, pair);
                    }
                }

                foreach (int key in stuffToChange.Keys)
                {
                    averagesCommaStandardDeviations[key] = stuffToChange[key];
                }

                rtHydrophobicityAvgDev.Add(filename, averagesCommaStandardDeviations);
            }
            return rtHydrophobicityAvgDev;
        }

        /// <summary>
        /// Assuming that we want to use different sets of training features for different search types, the search type needs to be determined.
        /// In the future, we may have that specified in the GUI and this will be moot.

        private static string DetermineSearchType(List<PeptideSpectralMatch> psms)
        {
            //TODO find a way to determine if crosslink serach type
            if (psms[0].DigestionParams.Protease.Name == "top-down")
            {
                return "topDown";
            }
            else
            {
                return "standard";
            }
        }

        private static float GetSSRCalcHydrophobicityZScore(PeptideSpectralMatch psm, PeptideWithSetModifications Peptide, Dictionary<string, Dictionary<int, Tuple<double, double>>> d)
        {
            //Using SSRCalc3 but probably any number of different calculators could be used instead. One could also use the CE mobility.
            SSRCalc3 calc = new SSRCalc3("SSRCalc 3.0 (300A)", SSRCalc3.Column.A300);
            double hydrophobicityZscore = double.NaN;

            if (d.ContainsKey(psm.FullFilePath))
            {
                int time = (int)(2 * Math.Round(psm.ScanRetentionTime / 2d, 0));
                if (d[psm.FullFilePath].Keys.Contains(time))
                {
                    double predictedHydrophobicity = calc.ScoreSequence(Peptide);
                    hydrophobicityZscore = Math.Abs(d[psm.FullFilePath][time].Item1 - predictedHydrophobicity) / d[psm.FullFilePath][time].Item2;
                }
            }

            double maxHydrophobicityZscore = 10;
            if (double.IsNaN(hydrophobicityZscore) || double.IsInfinity(hydrophobicityZscore) || hydrophobicityZscore > maxHydrophobicityZscore)
            {
                hydrophobicityZscore = maxHydrophobicityZscore;
            }

            return (float)hydrophobicityZscore;
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
        public static IEnumerable<PsmData> CreatePsmData(List<PeptideSpectralMatch> psms, Dictionary<string, int> sequenceToPsmCount, Dictionary<string, Dictionary<int, Tuple<double, double>>> timeDependantHydrophobicityAverageAndDeviation_unmodified, Dictionary<string, Dictionary<int, Tuple<double, double>>> timeDependantHydrophobicityAverageAndDeviation_modified, int chargeStateMode, string[] trainingVariables)
        {
            List<PsmData> pd = new List<PsmData>();
            foreach (PeptideSpectralMatch psm in psms)
            {
                foreach (var (notch, peptideWithSetMods) in psm.BestMatchingPeptides)
                {
                    bool label;
                    if (peptideWithSetMods.Protein.IsDecoy || psm.FdrInfo.QValue > 0.25)
                    {
                        label = false;
                        pd.Add(CreateOnePsmDataEntry(psm, sequenceToPsmCount, timeDependantHydrophobicityAverageAndDeviation_unmodified, timeDependantHydrophobicityAverageAndDeviation_modified, chargeStateMode, peptideWithSetMods, trainingVariables, notch, label));
                    }
                    else if (!peptideWithSetMods.Protein.IsDecoy && psm.FdrInfo.QValue <= 0.01)
                    {
                        label = true;
                        pd.Add(CreateOnePsmDataEntry(psm, sequenceToPsmCount, timeDependantHydrophobicityAverageAndDeviation_unmodified, timeDependantHydrophobicityAverageAndDeviation_modified, chargeStateMode, peptideWithSetMods, trainingVariables, notch, label));
                    }
                }
            }
            return pd.AsEnumerable();
        }

        public static PsmData CreateOnePsmDataEntry(PeptideSpectralMatch psm, Dictionary<string, int> sequenceToPsmCount, Dictionary<string, Dictionary<int, Tuple<double, double>>> timeDependantHydrophobicityAverageAndDeviation_unmodified, Dictionary<string, Dictionary<int, Tuple<double, double>>> timeDependantHydrophobicityAverageAndDeviation_modified, int chargeStateMode, PeptideWithSetModifications selectedPeptide, string[] trainingVariables, int notchToUse, bool label)
        {
            float totalMatchingFragmentCount = 0;
            if (trainingVariables.Contains("TotalMatchingFragmentCount"))
            {
                totalMatchingFragmentCount = (float)Math.Floor(psm.Score);
            }

            float ambiguity = 0;
            if (trainingVariables.Contains("Ambiguity"))
            {
                ambiguity = Math.Min((float)(psm.PeptidesToMatchingFragments.Keys.Count - 1), 10);
            }
            float intensity = 0;
            if (trainingVariables.Contains("Intensity"))
            {
                intensity = (float)(psm.Score - (int)psm.Score);
            }
            float chargeDifference = 0;
            if (trainingVariables.Contains("PrecursorChargeDiffToMode"))
            {
                chargeDifference = -Math.Abs(chargeStateMode - psm.ScanPrecursorCharge);
            }
            float deltaScore = 0;
            if (trainingVariables.Contains("DeltaScore"))
            {
                deltaScore = (float)psm.DeltaScore;
            }
            float psmCount = 1;
            if (trainingVariables.Contains("PsmCount"))
            {
                psmCount = sequenceToPsmCount[String.Join("|", psm.BestMatchingPeptides.Select(p => p.Peptide.FullSequence).ToList())];
                List<int> psmCountList = new List<int> { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 30, 40, 50, 75, 100, 200, 300, 400, 500 };
                int closest = psmCountList.OrderBy(item => Math.Abs(psmCount - item)).First();
                psmCount = closest;
            }

            int notch = 0;
            if (trainingVariables.Contains("Notch"))
            {
                notch = notchToUse;
            }

            if (selectedPeptide == null)
            {
                selectedPeptide = psm.BestMatchingPeptides.Select(p => p.Peptide).First();
            }

            float modCount = 0;
            if (trainingVariables.Contains("ModsCount"))
            {
                modCount = Math.Min((float)selectedPeptide.AllModsOneIsNterminus.Keys.Count(), 10);
            }

            float missedCleavages = 0;
            if (trainingVariables.Contains("MissedCleavagesCount"))
            {
                missedCleavages = selectedPeptide.MissedCleavages;
            }

            float longestSeq = 0;
            if (trainingVariables.Contains("LongestFragmentIonSeries"))
            {
                longestSeq = psm.GetLongestIonSeriesBidirectional(selectedPeptide);
            }

            float hydrophobicityZscore = float.NaN;

            if (selectedPeptide.BaseSequence.Equals(selectedPeptide.FullSequence) && trainingVariables.Contains("HydrophobicityZScore"))
            {
                hydrophobicityZscore = GetSSRCalcHydrophobicityZScore(psm, selectedPeptide, timeDependantHydrophobicityAverageAndDeviation_unmodified);
            }
            else if (trainingVariables.Contains("HydrophobicityZScore"))
            {
                hydrophobicityZscore = GetSSRCalcHydrophobicityZScore(psm, selectedPeptide, timeDependantHydrophobicityAverageAndDeviation_modified);
            }
            bool isVariantPeptide = PeptideIsVariant(selectedPeptide);

            if (psm.IsDecoy)
            {
                label = false;
            }
            else
            {
                label = true;
            }

            psm.PsmData_forPEPandPercolator = new PsmData
            {
                TotalMatchingFragmentCount = totalMatchingFragmentCount,
                Intensity = intensity,
                PrecursorChargeDiffToMode = chargeDifference,
                DeltaScore = deltaScore,
                Notch = notch,
                PsmCount = psmCount,
                ModsCount = modCount,
                MissedCleavagesCount = missedCleavages,
                Ambiguity = ambiguity,
                LongestFragmentIonSeries = longestSeq,
                HydrophobicityZScore = hydrophobicityZscore,
                IsVariantPeptide = Convert.ToSingle(isVariantPeptide),
                Label = label
            };

            return psm.PsmData_forPEPandPercolator;
        }

        private static bool PeptideIsVariant(PeptideWithSetModifications pwsm)
        {
            bool identifiedVariant = false;
            if (pwsm.Protein.AppliedSequenceVariations.Count() > 0)
            {
                foreach (var variant in pwsm.Protein.AppliedSequenceVariations)
                {
                    if (pwsm.IntersectsAndIdentifiesVariation(variant).identifies)
                    {
                        identifiedVariant = true;
                        break;
                    }
                }
            }
            return identifiedVariant;
        }

        /// <summary>
        /// At the time when the ~10% of the data gets chosen for training, another 10% gets chosen for evaluation. Then after training,
        /// the effectiveness of the model gets evaluated on the test set. The results of that evaluation are converted to text values called
        /// BinarySearchTreeMetrics and this gets written to the results.tsv

        public static string PrintBinaryClassificationMetrics(string name, CalibratedBinaryClassificationMetrics metrics, int ambiguousPeptidesRemovedCount)
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
            s.AppendLine("*       Count of Ambiguous Peptides Removed:  " + ambiguousPeptidesRemovedCount.ToString());
            s.AppendLine("************************************************************");
            return s.ToString();
        }
    }
}