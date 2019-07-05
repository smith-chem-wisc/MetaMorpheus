using EngineLayer.FdrAnalysis;
using Microsoft.ML;
using Microsoft.ML.Data;
using Proteomics.ProteolyticDigestion;
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
        public static string ComputePEPValuesForAllPSMsGeneric(List<PeptideSpectralMatch> psms)
        {
            Dictionary<string, int> accessionAppearances = GetAccessionCounts(psms);
            Dictionary<string, int> sequenceToPsmCount = GetSequenceToPSMCount(psms);

            MLContext mlContext = new MLContext();
            IDataView dataView = mlContext.Data.LoadFromEnumerable(CreatePsmData(psms, accessionAppearances, sequenceToPsmCount));

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

            var pipeline = mlContext.Transforms.Concatenate("Features", "Intensity", "ScanPrecursorCharge", "DeltaScore", "Notch", "PsmCount", "ModsCount", "MissedCleavagesCount", "Ambiguity", "AccessionAppearances", "LongestFragmentIonSeries")
                .Append(mlContext.BinaryClassification.Trainers.FastTree(labelColumnName: "Label", featureColumnName: "Features"));

            var trainedModel = pipeline.Fit(trainingData);

            var predictionEngine = mlContext.Model.CreatePredictionEngine<PsmData, TruePositivePrediction>(trainedModel);

            string ambiguousScans = "";

            //For Debug
            List<string> someOut = new List<string>();
            someOut.Add("Accessions|Ambiguity|DeltaScore|Intensity|Label|LongestSeries|MissedCleavages|ModsCount|Notch|PsmCount|PrecursorCharge|call|pepValue|Score|QValue");

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
                        PsmData pd = CreateOnePsmDataFromPsm2(psm, Notch, Peptide, accessionAppearances, sequenceToPsmCount);

                        var pepValuePrediction = predictionEngine.Predict(pd);

                        someOut.Add(pd.AccessionAppearances.ToString() + "|" + pd.Ambiguity.ToString() + "|" + pd.DeltaScore.ToString() + "|" + pd.Intensity.ToString() + "|" + pd.Label + "|" + pd.LongestFragmentIonSeries + "|" + pd.MissedCleavagesCount + "|" + pd.ModsCount + "|" + pd.Notch + "|" + pd.PsmCount + "|" + pd.ScanPrecursorCharge + "|" + pepValuePrediction.Prediction + "|" + pepValuePrediction.Probability + "|" + pepValuePrediction.Score);

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

                    psm.FdrInfo.PEP = 1-pepValuePredictions[0]; //they should all be the same at this point so it doesn't matter which you take. First is good.
                }
            }

            //For debug
            //File.WriteAllLines(@"C:\Users\Michael Shortreed\Downloads\psmDataVAlues.txt", someOut, System.Text.Encoding.UTF8);

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

        /// <summary>
        /// This method goes throught the complete set of peptide spectral matches and counts the number of appearances for each accession number.
        /// This includes all of the ambiguous entries. So, if one peptide spectral match has three peptides with three accessions, all are counted.
        /// </summary>
        /// <param name="psms"></param>
        /// <returns></returns>
        public static Dictionary<string, int> GetAccessionCounts(List<PeptideSpectralMatch> psms)
        {
            Dictionary<string, int> accessionCountDictionary = new Dictionary<string, int>();

            foreach (PeptideSpectralMatch psm in psms)
            {
                if (psm != null)
                {
                    foreach (var (Notch, Peptide) in psm.BestMatchingPeptides)
                    {
                        string accession = Peptide.Protein.Accession;
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
            }

            return accessionCountDictionary;
        }

        //This method ignores ambiguity and loads only the first peptide in a series for each PSM
        public static IEnumerable<PsmData> CreatePsmData(List<PeptideSpectralMatch> psms, Dictionary<string, int> accessionAppearances, Dictionary<string, int> sequenceToPsmCount, bool? trueOrFalse = null)
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
                    pd.Add(CreateOnePsmDataFromPsm(psm, accessionAppearances, sequenceToPsmCount, label));
                }
                else if (!psm.IsDecoy && psm.FdrInfo.QValue <= 0.01)
                {
                    label = true;
                    pd.Add(CreateOnePsmDataFromPsm(psm, accessionAppearances, sequenceToPsmCount, label));
                }
            }
            return pd.AsEnumerable();
        }

        public static PsmData CreateOnePsmDataFromPsm2(PeptideSpectralMatch psm, int notch, PeptideWithSetModifications firstPeptide, Dictionary<string, int> accessionCounts, Dictionary<string, int> sequenceToPsmCount, bool? trueOrFalse = null)
        {
            //dont' think ambiguity is helping so not using currently
            float ambiguity = (float)psm.PeptidesToMatchingFragments.Keys.Count;

            float intensity = (float)(psm.Score - (int)psm.Score);
            float charge = psm.ScanPrecursorCharge;
            float deltaScore = (float)psm.DeltaScore;
            float psmCount = sequenceToPsmCount[String.Join("|", psm.BestMatchingPeptides.Select(p => p.Peptide.FullSequence).ToList())];

            float modCount = firstPeptide.AllModsOneIsNterminus.Keys.Count();

            //todo: for non-specific cleavage, ignore missed cleavages
            float missedCleavages = firstPeptide.MissedCleavages;

            float longestSeq = psm.GetLongestIonSeriesBidirectional(firstPeptide);

            string accession = firstPeptide.Protein.Accession;
            float appearances;
            if (accessionCounts.Keys.Count != 0 && accessionCounts.ContainsKey(accession))
            {
                appearances = (float)accessionCounts[accession];
            }
            else
            {
                appearances = 1;
            }
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

        public static PsmData CreateOnePsmDataFromPsm(PeptideSpectralMatch psm, Dictionary<string, int> accessionCounts, Dictionary<string, int> sequenceToPsmCount, bool? trueOrFalse = null)
        {
            // TODO: some properties like DeltaScore will need to be recalculated if we keep top N peptides per PSM
            // and rerank them because the top-scoring peptide can change

            float ambiguity = (float)psm.PeptidesToMatchingFragments.Count;//(psm.BaseSequence.Split('|').Count());
            float intensity = (float)(psm.Score - (int)psm.Score);
            float charge = psm.ScanPrecursorCharge;
            float deltaScore = (float)psm.DeltaScore;
            float psmCount = sequenceToPsmCount[String.Join("|", psm.BestMatchingPeptides.Select(p => p.Peptide.FullSequence).ToList())];

            var firstPeptide = psm.BestMatchingPeptides.Select(p => p.Peptide).First();
            float modCount = firstPeptide.AllModsOneIsNterminus.Keys.Count();

            float notch = 0;
            if (psm.Notch.HasValue)
            {
                notch = psm.Notch.Value;
            }

            
            
            

            //todo: for non-specific cleavage, ignore missed cleavages
            float missedCleavages = firstPeptide.MissedCleavages;
            float longestSeq = psm.GetLongestIonSeriesBidirectional(firstPeptide);
            string accession = firstPeptide.Protein.Accession;
            float appearances;
            if(accessionCounts.Keys.Count != 0 && accessionCounts.ContainsKey(accession))
            { 
                appearances = (float)accessionCounts[accession];
            }
            else
            {
                appearances = 1;
            }
             
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