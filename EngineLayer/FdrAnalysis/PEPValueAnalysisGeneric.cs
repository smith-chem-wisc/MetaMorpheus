using Chemistry;
using EngineLayer.CrosslinkSearch;
using EngineLayer.FdrAnalysis;
using MathNet.Numerics.Statistics;
using Microsoft.ML;
using Microsoft.ML.Data;
using Proteomics;
using Proteomics.Fragmentation;
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
        private static string SeparationType;

        public static string ComputePEPValuesForAllPSMsGeneric(List<PeptideSpectralMatch> psms, string searchType, string separationType)
        {
            SeparationType = separationType;
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

            Dictionary<string, float> fileSpecificMedianFragmentMassErrors = GetFileSpecificMedianFragmentMassError(psms);

            MLContext mlContext = new MLContext();
            IDataView dataView = mlContext.Data.LoadFromEnumerable(CreatePsmData(searchType, psms, sequenceToPsmCount, fileSpecificTimeDependantHydrophobicityAverageAndDeviation_unmodified, fileSpecificTimeDependantHydrophobicityAverageAndDeviation_modified, fileSpecificMedianFragmentMassErrors, chargeStateMode, trainingVariables));

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
                        PsmData pd = CreateOnePsmDataEntry(searchType, psm, sequenceToPsmCount, fileSpecificTimeDependantHydrophobicityAverageAndDeviation_unmodified, fileSpecificTimeDependantHydrophobicityAverageAndDeviation_modified, fileSpecificMedianFragmentMassErrors, chargeStateMode, Peptide, trainingVariables, Notch, !Peptide.Protein.IsDecoy);
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

        public static Dictionary<string, Dictionary<int, Tuple<double, double>>> ComputeHydrophobicityValues(List<PeptideSpectralMatch> psms, bool computeHydrophobicitiesforModifiedPeptides)
        {
            SSRCalc3 calc = new SSRCalc3("SSRCalc 3.0 (300A)", SSRCalc3.Column.A300);

            //TODO change the tuple so the values have names
            Dictionary<string, Dictionary<int, Tuple<double, double>>> rtHydrophobicityAvgDev = new Dictionary<string, Dictionary<int, Tuple<double, double>>>();

            List<string> filenames = psms.Select(f => f.FullFilePath).ToList();

            filenames = filenames.Distinct().ToList();

            foreach (string filename in filenames)
            {
                Dictionary<int, List<double>> hydrophobicities = new Dictionary<int, List<double>>();
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

                        double predictedHydrophobicity = 0;

                        if (SeparationType == "CZE")
                        {
                            predictedHydrophobicity = CZE.PredictedElectrophoreticMobility(pwsm.BaseSequence, pwsm.MonoisotopicMass);
                        }
                        else
                        {
                            predictedHydrophobicity = calc.ScoreSequence(pwsm);
                        }

                        //here i'm grouping this in 2 minute increments becuase there are cases where you get too few data points to get a good standard deviation an average. This is for stability.
                        int possibleKey = (int)(2 * Math.Round(psm.ScanRetentionTime / 2d, 0));

                        //First block of if statement is for modified peptides.
                        if ((pwsm.AllModsOneIsNterminus.Any() && computeHydrophobicitiesforModifiedPeptides && SeparationType == "HPLC") || (ContainsModificationsThatShiftMobility(pwsm.AllModsOneIsNterminus.Values.AsEnumerable()) && computeHydrophobicitiesforModifiedPeptides && SeparationType == "CZE"))
                        {
                            if (hydrophobicities.ContainsKey(possibleKey))
                            {
                                hydrophobicities[possibleKey].Add(predictedHydrophobicity);
                            }
                            else
                            {
                                hydrophobicities.Add(possibleKey, new List<double>() { predictedHydrophobicity });
                            }
                        }
                        //this second block of if statment is for unmodified peptides.
                        else if ((!pwsm.AllModsOneIsNterminus.Any() && !computeHydrophobicitiesforModifiedPeptides && SeparationType == "HPLC") || (!ContainsModificationsThatShiftMobility(pwsm.AllModsOneIsNterminus.Values.AsEnumerable()) && !computeHydrophobicitiesforModifiedPeptides && SeparationType == "CZE"))
                        {
                            if (hydrophobicities.ContainsKey(possibleKey))
                            {
                                hydrophobicities[possibleKey].Add(predictedHydrophobicity);
                            }
                            else
                            {
                                hydrophobicities.Add(possibleKey, new List<double>() { predictedHydrophobicity });
                            }
                        }
                    }
                }

                List<double> allSquaredHyrophobicityDifferences = new List<double>();

                foreach (int retentionTimeBin in hydrophobicities.Keys)
                {
                    //TODO consider using inner-quartile range instead of standard deviation
                    double averageHydrophobicity = hydrophobicities[retentionTimeBin].Average();
                    averagesCommaStandardDeviations.Add(retentionTimeBin, new Tuple<double, double>(averageHydrophobicity, hydrophobicities[retentionTimeBin].StandardDeviation()));
                    foreach (double hydrophobicity in hydrophobicities[retentionTimeBin])
                    {
                        double difference = Math.Abs(hydrophobicity - averageHydrophobicity);
                        if (!double.IsNaN(difference) && difference > 0)
                        {
                            allSquaredHyrophobicityDifferences.Add(Math.Pow(difference, 2));
                        }
                    }
                }

                //some standard deviations are too small or too large because of random reasons, so we replace those small numbers of oddballs with reasonable numbers.
                double globalStDev = 1;
                if (allSquaredHyrophobicityDifferences.Count() > 1)
                {
                    globalStDev = Math.Sqrt(allSquaredHyrophobicityDifferences.Sum() / (allSquaredHyrophobicityDifferences.Count() - 1));
                }

                Dictionary<int, Tuple<double, double>> stDevsToChange = new Dictionary<int, Tuple<double, double>>();
                foreach (KeyValuePair<int, Tuple<double, double>> item in averagesCommaStandardDeviations)
                {
                    //add stability. not allowing stdevs that are too small or too large at one position relative to the global stdev
                    //here we are finding which stdevs are out of whack.
                    if (Double.IsNaN(item.Value.Item2) || item.Value.Item2 < 0.5 || (item.Value.Item2 / globalStDev) > 3)
                    {
                        Tuple<double, double> pair = new Tuple<double, double>(averagesCommaStandardDeviations[item.Key].Item1, globalStDev);
                        stDevsToChange.Add(item.Key, pair);
                    }
                }
                //here we are replacing the stdevs that are out of whack.
                foreach (int key in stDevsToChange.Keys)
                {
                    averagesCommaStandardDeviations[key] = stDevsToChange[key];
                }

                rtHydrophobicityAvgDev.Add(filename, averagesCommaStandardDeviations);
            }
            return rtHydrophobicityAvgDev;
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
                    double predictedHydrophobicity;

                    if (SeparationType == "CZE")
                    {
                        predictedHydrophobicity = CZE.PredictedElectrophoreticMobility(Peptide.BaseSequence, Peptide.MonoisotopicMass);
                    }
                    else//assume hplc
                    {
                        predictedHydrophobicity = calc.ScoreSequence(Peptide);
                    }

                    hydrophobicityZscore = Math.Abs(d[psm.FullFilePath][time].Item1 - predictedHydrophobicity) / d[psm.FullFilePath][time].Item2;
                }
            }

            double maxHydrophobicityZscore = 10; // each "Z" is one standard deviation. so, maxHydrophobicityZscore 10 is quite large
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
        public static IEnumerable<PsmData> CreatePsmData(string searchType, List<PeptideSpectralMatch> psms, Dictionary<string, int> sequenceToPsmCount, Dictionary<string, Dictionary<int, Tuple<double, double>>> timeDependantHydrophobicityAverageAndDeviation_unmodified, Dictionary<string, Dictionary<int, Tuple<double, double>>> timeDependantHydrophobicityAverageAndDeviation_modified, Dictionary<string, float> fileSpecificMedianFragmentMassErrors, int chargeStateMode, string[] trainingVariables)
        {
            List<PsmData> pd = new List<PsmData>();

            foreach (PeptideSpectralMatch psm in psms)
            {
                if (searchType == "crosslink")
                {
                    CrosslinkSpectralMatch csm = (CrosslinkSpectralMatch)psm;

                    bool label;
                    if (csm.IsDecoy || csm.BetaPeptide.IsDecoy || psm.FdrInfo.QValue > 0.25)
                    {
                        label = false;
                        pd.Add(CreateOnePsmDataEntry(searchType, psm, sequenceToPsmCount, timeDependantHydrophobicityAverageAndDeviation_unmodified, timeDependantHydrophobicityAverageAndDeviation_modified, fileSpecificMedianFragmentMassErrors, chargeStateMode, csm.BestMatchingPeptides.First().Peptide, trainingVariables, 0, label));
                    }
                    else if (!csm.IsDecoy && !csm.BetaPeptide.IsDecoy && psm.FdrInfo.QValue <= 0.01)
                    {
                        label = true;
                        pd.Add(CreateOnePsmDataEntry(searchType, psm, sequenceToPsmCount, timeDependantHydrophobicityAverageAndDeviation_unmodified, timeDependantHydrophobicityAverageAndDeviation_modified, fileSpecificMedianFragmentMassErrors, chargeStateMode, csm.BestMatchingPeptides.First().Peptide, trainingVariables, 0, label));
                    }
                }
                else
                {
                    foreach (var (notch, peptideWithSetMods) in psm.BestMatchingPeptides)
                    {
                        bool label;
                        if (peptideWithSetMods.Protein.IsDecoy || psm.FdrInfo.QValue > 0.25)
                        {
                            label = false;
                            pd.Add(CreateOnePsmDataEntry(searchType, psm, sequenceToPsmCount, timeDependantHydrophobicityAverageAndDeviation_unmodified, timeDependantHydrophobicityAverageAndDeviation_modified, fileSpecificMedianFragmentMassErrors, chargeStateMode, peptideWithSetMods, trainingVariables, notch, label));
                        }
                        else if (!peptideWithSetMods.Protein.IsDecoy && psm.FdrInfo.QValue <= 0.01)
                        {
                            label = true;
                            pd.Add(CreateOnePsmDataEntry(searchType, psm, sequenceToPsmCount, timeDependantHydrophobicityAverageAndDeviation_unmodified, timeDependantHydrophobicityAverageAndDeviation_modified, fileSpecificMedianFragmentMassErrors, chargeStateMode, peptideWithSetMods, trainingVariables, notch, label));
                        }
                    }
                }
            }
            return pd.AsEnumerable();
        }

        public static PsmData CreateOnePsmDataEntry(string searchType, PeptideSpectralMatch psm, Dictionary<string, int> sequenceToPsmCount, Dictionary<string, Dictionary<int, Tuple<double, double>>> timeDependantHydrophobicityAverageAndDeviation_unmodified, Dictionary<string, Dictionary<int, Tuple<double, double>>> timeDependantHydrophobicityAverageAndDeviation_modified, Dictionary<string, float> fileSpecificMedianFragmentMassErrors, int chargeStateMode, PeptideWithSetModifications selectedPeptide, string[] trainingVariables, int notchToUse, bool label)
        {
            float totalMatchingFragmentCount = 0;
            float intensity = 0;
            float chargeDifference = 0;
            float deltaScore = 0;
            float psmCount = 1;
            int notch = 0;
            float ambiguity = 0;
            float modCount = 0;
            float absoluteFragmentMassError = 0;

            float missedCleavages = 0;
            float longestSeq = 0;
            float hydrophobicityZscore = float.NaN;
            bool isVariantPeptide = false;

            //crosslink specific features
            float alphaIntensity = 0;
            float betaIntensity = 0;
            float longestFragmentIonSeries_Alpha = 0;
            float longestFragmentIonSeries_Beta = 0;
            float isDeadEnd = 0;
            float isLoop = 0;
            float isInter = 0;
            float isIntra = 0;

            if (searchType != "crosslink")
            {
                totalMatchingFragmentCount = (float)Math.Floor(psm.Score);
                intensity = (float)(psm.Score - (int)psm.Score);
                chargeDifference = -Math.Abs(chargeStateMode - psm.ScanPrecursorCharge);
                deltaScore = (float)psm.DeltaScore;
                notch = notchToUse;
                modCount = Math.Min((float)selectedPeptide.AllModsOneIsNterminus.Keys.Count(), 10);
                if (psm.PeptidesToMatchingFragments[selectedPeptide]?.Count() > 0)
                {
                    absoluteFragmentMassError = (float)Math.Abs(GetAverageFragmentMassError(psm.PeptidesToMatchingFragments[selectedPeptide]) - fileSpecificMedianFragmentMassErrors[psm.FullFilePath]);
                }

                ambiguity = Math.Min((float)(psm.PeptidesToMatchingFragments.Keys.Count - 1), 10);
                longestSeq = psm.GetLongestIonSeriesBidirectional(psm.PeptidesToMatchingFragments, selectedPeptide);

                //grouping psm counts as follows is done for stability. you get very nice numbers at low psms to get good statistics. But you get a few peptides with high psm counts that could be either targets or decoys and the values swing between extremes. So grouping psms in bundles really adds stability.
                psmCount = sequenceToPsmCount[String.Join("|", psm.BestMatchingPeptides.Select(p => p.Peptide.FullSequence).ToList())];
                List<int> psmCountList = new List<int> { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 30, 40, 50, 75, 100, 200, 300, 400, 500 };
                int closest = psmCountList.OrderBy(item => Math.Abs(psmCount - item)).First();
                psmCount = closest;
                isVariantPeptide = PeptideIsVariant(selectedPeptide);

                if (psm.DigestionParams.Protease.Name != "top-down")
                {
                    missedCleavages = selectedPeptide.MissedCleavages;
                    if (String.IsNullOrEmpty(SeparationType) || SeparationType == "HPLC")
                    {
                        if (selectedPeptide.BaseSequence.Equals(selectedPeptide.FullSequence))
                        {
                            hydrophobicityZscore = GetSSRCalcHydrophobicityZScore(psm, selectedPeptide, timeDependantHydrophobicityAverageAndDeviation_unmodified);
                        }
                        else
                        {
                            hydrophobicityZscore = GetSSRCalcHydrophobicityZScore(psm, selectedPeptide, timeDependantHydrophobicityAverageAndDeviation_modified);
                        }
                    }
                    else if (SeparationType == "CZE")
                    {
                        if (!ContainsModificationsThatShiftMobility(selectedPeptide.AllModsOneIsNterminus.Values.AsEnumerable()))
                        {
                            hydrophobicityZscore = GetSSRCalcHydrophobicityZScore(psm, selectedPeptide, timeDependantHydrophobicityAverageAndDeviation_unmodified);
                        }
                        else
                        {
                            hydrophobicityZscore = GetSSRCalcHydrophobicityZScore(psm, selectedPeptide, timeDependantHydrophobicityAverageAndDeviation_modified);
                        }
                    }
                }
                //this is not for actual crosslinks but for the byproducts of crosslink loop links, deadends, etc.
                if (psm is CrosslinkSpectralMatch)
                {
                    CrosslinkSpectralMatch csm = (CrosslinkSpectralMatch)psm;
                    isDeadEnd = Convert.ToSingle((csm.CrossType == PsmCrossType.DeadEnd) || (csm.CrossType == PsmCrossType.DeadEndH2O) || (csm.CrossType == PsmCrossType.DeadEndNH2) || (csm.CrossType == PsmCrossType.DeadEndTris));
                    isLoop = Convert.ToSingle(csm.CrossType == PsmCrossType.Loop);
                }
            }
            else
            {
                CrosslinkSpectralMatch csm = (CrosslinkSpectralMatch)psm;
                PeptideWithSetModifications selectedAlphaPeptide = csm.BestMatchingPeptides.Select(p => p.Peptide).First();
                PeptideWithSetModifications selectedBetaPeptide = csm.BetaPeptide?.BestMatchingPeptides.Select(p => p.Peptide).First();

                totalMatchingFragmentCount = (float)Math.Round(csm.XLTotalScore, 0);

                //Compute fragment mass error
                int alphaCount = 0;
                float alphaError = 0;
                if (csm.PeptidesToMatchingFragments[selectedAlphaPeptide]?.Count > 0)
                {
                    alphaCount = csm.PeptidesToMatchingFragments[selectedAlphaPeptide].Count;
                    alphaError = (float)Math.Abs(GetAverageFragmentMassError(csm.PeptidesToMatchingFragments[selectedAlphaPeptide]));
                }
                int betaCount = 0;
                float betaError = 0;
                if (csm.BetaPeptide.PeptidesToMatchingFragments[selectedBetaPeptide]?.Count > 0)
                {
                    betaCount = csm.BetaPeptide.PeptidesToMatchingFragments[selectedBetaPeptide].Count;
                    betaError = (float)Math.Abs(GetAverageFragmentMassError(csm.BetaPeptide.PeptidesToMatchingFragments[selectedBetaPeptide]));
                }

                float averageError = 0;
                if ((alphaCount + betaCount) > 0)
                {
                    averageError = (alphaCount * alphaError + betaCount * betaError) / (alphaCount + betaCount);
                }

                absoluteFragmentMassError = averageError - fileSpecificMedianFragmentMassErrors[csm.FullFilePath];
                //End compute fragment mass error

                deltaScore = (float)csm.DeltaScore;
                chargeDifference = -Math.Abs(chargeStateMode - psm.ScanPrecursorCharge);
                alphaIntensity = (float)(csm.Score - (int)csm.Score);
                betaIntensity = csm.BetaPeptide == null ? (float)0 : (float)(csm.BetaPeptide.Score - (int)csm.BetaPeptide.Score);
                longestFragmentIonSeries_Alpha = psm.GetLongestIonSeriesBidirectional(csm.PeptidesToMatchingFragments, selectedAlphaPeptide);
                longestFragmentIonSeries_Beta = selectedBetaPeptide == null ? (float)0 : csm.GetLongestIonSeriesBidirectional(csm.BetaPeptide.PeptidesToMatchingFragments, selectedBetaPeptide);
                isInter = Convert.ToSingle(csm.CrossType == PsmCrossType.Inter);
                isIntra = Convert.ToSingle(csm.CrossType == PsmCrossType.Intra);
            }

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
                AbsoluteAverageFragmentMassErrorFromMedian = absoluteFragmentMassError,
                MissedCleavagesCount = missedCleavages,
                Ambiguity = ambiguity,
                LongestFragmentIonSeries = longestSeq,
                HydrophobicityZScore = hydrophobicityZscore,
                IsVariantPeptide = Convert.ToSingle(isVariantPeptide),

                AlphaIntensity = alphaIntensity,
                BetaIntensity = betaIntensity,
                LongestFragmentIonSeries_Alpha = longestFragmentIonSeries_Alpha,
                LongestFragmentIonSeries_Beta = longestFragmentIonSeries_Beta,
                IsDeadEnd = isDeadEnd,
                IsLoop = isLoop,
                IsInter = isInter,
                IsIntra = isIntra,

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

        public static bool ContainsModificationsThatShiftMobility(IEnumerable<Modification> modifications)
        {
            List<string> shiftingModifications = new List<string> { "Acetylation", "Ammonia loss", "Carbamyl", "Deamidation", "Formylation",
                "N2-acetylarginine", "N6-acetyllysine", "N-acetylalanine", "N-acetylaspartate", "N-acetylcysteine", "N-acetylglutamate", "N-acetylglycine",
                "N-acetylisoleucine", "N-acetylmethionine", "N-acetylproline", "N-acetylserine", "N-acetylthreonine", "N-acetyltyrosine", "N-acetylvaline",
                "Phosphorylation", "Phosphoserine", "Phosphothreonine", "Phosphotyrosine", "Sulfonation" };

            return shiftingModifications.Concat(modifications.Select(m => m.OriginalId).Distinct()).GroupBy(s => s).Where(s => s.Count() > 1).Any();
        }

        public static Dictionary<string, float> GetFileSpecificMedianFragmentMassError(List<PeptideSpectralMatch> psms)
        {
            Dictionary<string, float> fileSpecificMassErrors = new Dictionary<string, float>();
            foreach (string fullFilePath in psms.Select(p => p.FullFilePath).Distinct())
            {
                fileSpecificMassErrors.Add(fullFilePath, GetMedianAverageMassError(psms.Where(p => p.FullFilePath == fullFilePath)));
            }
            return fileSpecificMassErrors;
        }

        public static float GetMedianAverageMassError(IEnumerable<PeptideSpectralMatch> psms)
        {
            List<float> averageMassErrors = new List<float>();
            foreach (PeptideSpectralMatch psm in psms)
            {
                {
                    foreach (KeyValuePair<PeptideWithSetModifications, List<MatchedFragmentIon>> peptide_MFI in psm.PeptidesToMatchingFragments)
                    {
                        if (peptide_MFI.Value != null && peptide_MFI.Value.Count > 0)
                        {
                            averageMassErrors.Add(GetAverageFragmentMassError(peptide_MFI.Value));
                        }
                    }
                }
            }
            return averageMassErrors.Median();
        }

        public static float GetAverageFragmentMassError(IEnumerable<MatchedFragmentIon> matchedIons)
        {
            var matchedIonsGroupedByProductType = matchedIons.GroupBy(i => i.NeutralTheoreticalProduct.ProductType).OrderBy(i => i.Key).ToList();
            List<float> massErrors = new List<float>();
            foreach (var productType in matchedIonsGroupedByProductType)
            {
                var products = productType.OrderBy(p => p.NeutralTheoreticalProduct.TerminusFragment.FragmentNumber)
                    .ToList();

                for (int i = 0; i < products.Count; i++)
                {
                    MatchedFragmentIon ion = products[i];

                    float massError = (float)(ion.Mz.ToMass(ion.Charge) - ion.NeutralTheoreticalProduct.NeutralMass);
                    float ppmMassError = (float)(massError / ion.NeutralTheoreticalProduct.NeutralMass * 1e6);
                    massErrors.Add(ppmMassError);
                }
            }

            return massErrors.Average();
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