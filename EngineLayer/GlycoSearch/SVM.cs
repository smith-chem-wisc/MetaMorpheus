using System;
using System.Text;
using System.Collections.Generic;
using System.Linq;
using Microsoft.ML;
using Microsoft.ML.Data;
using Microsoft.ML.Trainers;
using static Microsoft.ML.DataOperationsCatalog;
using System.IO;

namespace EngineLayer.GlycoSearch
{
    public static class SVM
    {
        public static void RunSVM(string fileIn, string filePath)
        {
            List<GlycanIndicator> gis = GlycanIndicator.ReadGsms(fileIn);

            var metric = RunSVM_GlycanIndicator(gis);
            using (StreamWriter output = new StreamWriter(filePath))
            {
                output.Write(PrintMetrics(metric));
            }
        }

        public static BinaryClassificationMetrics RunSVM_GlycanIndicator(List<GlycanIndicator> gis)
        {
            // Create a new context for ML.NET operations. It can be used for
            // exception tracking and logging, as a catalog of available operations
            // and as the source of randomness. Setting the seed to a fixed number
            // in this example to make outputs deterministic.
            var mlContext = new MLContext(seed: 0);

            // Create a list of training data points.
            var dataPoints = DataPoint.GenerateDataPoints(gis);

            // Convert the list of data points to an IDataView object, which is
            // consumable by ML.NET API.
            var allData = mlContext.Data.LoadFromEnumerable(dataPoints);

            //Split data to train and test.
            int randomSeed = 42;
            double fraction = Math.Min(0.25, 1000000.0 / (double)gis.Count);//make training set 25% of the data up to a training set of one million psms
            TrainTestData trainTestSplit = mlContext.Data.TrainTestSplit(allData, testFraction: fraction, null, randomSeed);
            IDataView trainingData = trainTestSplit.TrainSet;
            IDataView testData = trainTestSplit.TestSet;

            // Define trainer options.
            var options = new LinearSvmTrainer.Options
            {
                BatchSize = 10,
                PerformProjection = true,
                NumberOfIterations = 10
            };

            // Define the trainer.
            var pipeline = mlContext.BinaryClassification.Trainers.LinearSvm(options);

            // Train the model.
            var model = pipeline.Fit(trainingData);

            // Run the model on test data set.
            var transformedTestData = model.Transform(testData);

            // Convert IDataView object to a list.
            // var predictions = mlContext.Data.CreateEnumerable<Prediction>(transformedTestData, reuseRowObject: false).ToList();

            // Evaluate the overall metrics.
            var metrics = mlContext.BinaryClassification.EvaluateNonCalibrated(transformedTestData);

            return metrics;
        }

        // Pretty-print BinaryClassificationMetrics objects.
        private static string PrintMetrics(BinaryClassificationMetrics metrics)
        {
            StringBuilder s = new StringBuilder();
            s.AppendLine("************************************************************");
            s.AppendLine("Metrics for Determination of PEP Using Binary Classification      ");
            s.AppendLine("-----------------------------------------------------------");
            s.AppendLine($"Accuracy: {metrics.Accuracy:F2}");
            s.AppendLine($"AUC: {metrics.AreaUnderRocCurve:F2}");
            s.AppendLine($"F1 Score: {metrics.F1Score:F2}");
            s.AppendLine($"Negative Precision: " + $"{metrics.NegativePrecision:F2}");
            s.AppendLine($"Negative Recall: {metrics.NegativeRecall:F2}");
            s.AppendLine($"Positive Precision: " + $"{metrics.PositivePrecision:F2}");
            s.AppendLine(metrics.ConfusionMatrix.GetFormattedConfusionTable());
            s.AppendLine("************************************************************");
            return s.ToString();
        }
    }
}
