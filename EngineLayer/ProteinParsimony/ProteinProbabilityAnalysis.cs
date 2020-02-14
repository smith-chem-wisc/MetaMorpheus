using EngineLayer.FdrAnalysis;
using Microsoft.ML;
using System;
using System.Collections.Generic;
using System.Linq;
using static Microsoft.ML.DataOperationsCatalog;

namespace EngineLayer
{
    public class ProteinProbabilityAnalysis
    {
        private const double PeptidePValueCutoff = 0.0;
        private const double ProteinPValueCutoff = 0.7;

        public static void ComputeProteinProbabilities(List<ProteinGroup> proteinGroups, bool modPeptidesAreDifferent)
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
                ProteinData pd = CreateProteinDataEntry(pg, modPeptidesAreDifferent);
                var proteinPrediction = predictionEngine.Predict(pd);
                pg.Probability = proteinPrediction.Probability;
            }
        }

        // build the protein data set that we will feed into the model
        public static IEnumerable<ProteinData> CreateProteinData(List<ProteinGroup> proteinGroups, bool treatModPeptidesDifferent)
        {
            List<ProteinData> proteinDataList = new List<ProteinData>();

            // due to parsimony, each protein group has 1 protein
            foreach (ProteinGroup pg in proteinGroups)
            {
                ProteinData newProteinData = CreateProteinDataEntry(pg, treatModPeptidesDifferent);
                proteinDataList.Add(newProteinData);
            }

            return proteinDataList.AsEnumerable();
        }

        // each entry contains information about the protein and its features
        public static ProteinData CreateProteinDataEntry(ProteinGroup pg, bool treatModPeptidesDifferent)
        {
            double averagePEP = 0.0;
            //bool containUniquePeptides = false;
            double percentageUnique = 0;
            int totalPeptideCount = 0;
            bool label = false;

            percentageUnique = (double)pg.UniquePeptides.Count / pg.AllPeptides.Count;
            //containUniquePeptides = pg.UniquePeptides.Count > 0 ? true : containUniquePeptides;
            totalPeptideCount = pg.AllPeptides.Count;

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
                if (peptide.FdrInfo.PEP <= PeptidePValueCutoff)
                {
                    continue;
                }

                sumPEP += peptide.FdrInfo.PEP;
            }

            averagePEP = sumPEP / peptides.Count;

            // fixme: best way to label training data
            // to get started, compare protein p value with an arbitrary cutoff at 0.7
            label = pg.PValue > ProteinPValueCutoff ? true : label;

            return new ProteinData
            {
                AveragePEP = (float)averagePEP,
                //ContainsUniquePeptides = Convert.ToSingle(containUniquePeptides),
                PercentageOfUniquePeptides = (float)percentageUnique,
                TotalPeptideCount = totalPeptideCount,
                Label = label
            };
        }

    }

    public class ProteinData
    {
        // todo: identify features to train the model on
        public static readonly string[] featuresForTraining = new string[] { "AveragePEP", "PercentageOfUniquePeptides", "TotalPeptideCount" };

        public float AveragePEP { get; set; }
        public float PercentageOfUniquePeptides { get; set; }
        public float TotalPeptideCount { get; set; }
        public bool Label { get; set; }
    }
}
