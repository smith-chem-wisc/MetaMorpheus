using EngineLayer.FdrAnalysis;
using Microsoft.ML;
using Microsoft.ML.Data;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using static Microsoft.ML.DataOperationsCatalog;

namespace EngineLayer
{
    public class ProteinProbabilityAnalysis
    {
        private const double PeptidePValueCutoff = 0.0;

        public static string ComputeProteinProbabilities(List<ProteinGroup> proteinGroups, bool modPeptidesAreDifferent)
        {
            // create ML context
            MLContext mlContext = new MLContext();
            IDataView dataView = mlContext.Data.LoadFromEnumerable(CreateProteinData(proteinGroups, modPeptidesAreDifferent));

            // split data into train and test
            TrainTestData trainTestSplit = mlContext.Data.TrainTestSplit(dataView, testFraction: 0.1, null, 42);
            IDataView trainingData = trainTestSplit.TrainSet;
            IDataView testData = trainTestSplit.TestSet;

            // get training features
            var features = ProteinData.featuresForTraining;

            // create model
            var trainer = mlContext.BinaryClassification.Trainers.FastTree(labelColumnName: "Label", featureColumnName: "Features");

            var pipeline = mlContext.Transforms.Concatenate("Features", features)
                .Append(mlContext.BinaryClassification.Trainers.FastTree(labelColumnName: "Label", featureColumnName: "Features"));

            // train model
            var trainedModel = pipeline.Fit(trainingData);

            // create prediction engine
            var predictionEngine = mlContext.Model.CreatePredictionEngine<ProteinData, TruePositivePrediction>(trainedModel);

            // predict protein probability for each protein entry
            foreach (ProteinGroup pg in proteinGroups)
            {
                ProteinData pd = CreateProteinDataEntry(pg, modPeptidesAreDifferent, false); // fixme
                var proteinPrediction = predictionEngine.Predict(pd);
                pg.Probability = proteinPrediction.Probability;
            }

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
        }

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

        // build the protein data set that we will feed into the model
        public static IEnumerable<ProteinData> CreateProteinData(List<ProteinGroup> proteinGroups, bool treatModPeptidesDifferent)
        {
            List<ProteinData> proteinDataList = new List<ProteinData>();

            // due to parsimony, each protein group has 1 protein
            foreach (ProteinGroup pg in proteinGroups)
            {
                bool label;
                if (pg.IsDecoy || pg.QValue > 0.05)
                {
                    label = false;
                    ProteinData newProteinData = CreateProteinDataEntry(pg, treatModPeptidesDifferent, label);
                    proteinDataList.Add(newProteinData);
                }
                else if (!pg.IsDecoy && pg.QValue <= 0.01)
                {
                    label = true;
                    ProteinData newProteinData = CreateProteinDataEntry(pg, treatModPeptidesDifferent, label);
                    proteinDataList.Add(newProteinData);
                }
            }

            return proteinDataList.AsEnumerable();
        }

        // each entry contains information about the protein and its features
        public static ProteinData CreateProteinDataEntry(ProteinGroup pg, bool treatModPeptidesDifferent, bool label)
        {
            double averagePEP = 0.0;
            double percentageUnique = 0;
            int totalPeptideCount = 0;
            double sequenceCoverageFraction = 0.0;
            int numProteinAccessions = 0;

            percentageUnique = (double)pg.UniquePeptides.Count / pg.AllPeptides.Count;
            totalPeptideCount = pg.AllPeptides.Count;
            sequenceCoverageFraction = pg.SequenceCoverageFraction.First();
            numProteinAccessions = pg.Proteins.Count;

            List<PeptideSpectralMatch> peptides;
            if (treatModPeptidesDifferent)
            {
                peptides = pg.AllPsmsBelowOnePercentFDR.GroupBy(b => b.FullSequence).Select(b => b.FirstOrDefault()).ToList();
            }
            else
            {
                peptides = pg.AllPsmsBelowOnePercentFDR.GroupBy(b => b.BaseSequence).Select(b => b.FirstOrDefault()).ToList();
            }

            double sumPEP = 0.0;
            foreach (PeptideSpectralMatch peptide in peptides)
            {
                if (peptide.FdrInfo.PEP <= PeptidePValueCutoff) // p value shouldnt be less than 0
                {
                    continue;
                }

                sumPEP += peptide.FdrInfo.PEP;
            }

            averagePEP = sumPEP / peptides.Count;

            return new ProteinData
            {
                AveragePEP = (float)averagePEP,
                PercentageOfUniquePeptides = (float)percentageUnique,
                TotalPeptideCount = totalPeptideCount,
                PercentageSequenceCoverage = (float)sequenceCoverageFraction,
                NumProteinAccessions = numProteinAccessions,
                Label = label
            };
        }

    }

    public class ProteinData
    {
        // todo: identify features to train the model on
        //public static readonly string[] featuresForTraining = new string[] { "AveragePEP", "PercentageOfUniquePeptides", "TotalPeptideCount", "PercentageSequenceCoverage", "NumProteinAccessions" };
        public static readonly string[] featuresForTraining = new string[] { "AveragePEP", "PercentageOfUniquePeptides", "TotalPeptideCount" };

        public float AveragePEP { get; set; }
        public float PercentageOfUniquePeptides { get; set; }
        public float TotalPeptideCount { get; set; }
        public float PercentageSequenceCoverage { get; set; }
        public float NumProteinAccessions { get; set; }
        public bool Label { get; set; }
    }
}
