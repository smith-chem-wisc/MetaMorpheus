using MassSpectrometry;
using SharpLearning.Common.Interfaces;
using SharpLearning.Containers.Matrices;
using SharpLearning.CrossValidation.TrainingTestSplitters;
using SharpLearning.GradientBoost.Learners;
using SharpLearning.GradientBoost.Models;
using SharpLearning.Metrics.Regression;
using SharpLearning.Optimization;
using SharpLearning.RandomForest.Learners;
using SharpLearning.RandomForest.Models;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace EngineLayer.Calibration
{
    public class CalibrationEngine : MetaMorpheusEngine
    {
        private const double MaximumFracForTraining = 0.70;
        private const double MaximumDatapointsToTrainWith = 20000;
        private const int TrainingIterations = 30;
        private readonly int RandomSeed;

        private readonly MsDataFile MyMsDataFile;
        private readonly DataPointAquisitionResults Datapoints;

        public CalibrationEngine(MsDataFile myMSDataFile, DataPointAquisitionResults datapoints, CommonParameters commonParameters, List<string> nestedIds)
            : base(commonParameters, nestedIds)
        {
            MyMsDataFile = myMSDataFile;
            Datapoints = datapoints;

            // set the random seed based on raw file properties
            if (MyMsDataFile.SourceFile != null && !string.IsNullOrEmpty(MyMsDataFile.SourceFile.CheckSum))
            {
                RandomSeed = MyMsDataFile.SourceFile.CheckSum.GetHashCode();
            }
            else
            {
                RandomSeed = MyMsDataFile.NumSpectra;
            }
        }

        protected override MetaMorpheusEngineResults RunSpecific()
        {
            double ms1fracForTraining = MaximumFracForTraining;
            double ms2fracForTraining = MaximumFracForTraining;

            var myMs1DataPoints = new List<(double[] xValues, double yValue)>();
            var myMs2DataPoints = new List<(double[] xValues, double yValue)>();

            // generate MS1 calibration datapoints
            for (int i = 0; i < Datapoints.Ms1List.Count; i++)
            {
                // x values
                var explanatoryVariables = new double[5];
                explanatoryVariables[0] = Datapoints.Ms1List[i].ExperimentalMz;
                explanatoryVariables[1] = Datapoints.Ms1List[i].RetentionTime;
                explanatoryVariables[2] = Datapoints.Ms1List[i].LogTotalIonCurrent;
                explanatoryVariables[3] = Datapoints.Ms1List[i].LogInjectionTime;
                explanatoryVariables[4] = Datapoints.Ms1List[i].LogIntensity;

                // y value
                double mzError = Datapoints.Ms1List[i].AbsoluteMzError;

                myMs1DataPoints.Add((explanatoryVariables, mzError));
            }

            // generate MS2 calibration datapoints
            for (int i = 0; i < Datapoints.Ms2List.Count; i++)
            {
                // x values
                var explanatoryVariables = new double[5];
                explanatoryVariables[0] = Datapoints.Ms2List[i].ExperimentalMz;
                explanatoryVariables[1] = Datapoints.Ms2List[i].RetentionTime;
                explanatoryVariables[2] = Datapoints.Ms2List[i].LogTotalIonCurrent;
                explanatoryVariables[3] = Datapoints.Ms2List[i].LogInjectionTime;
                explanatoryVariables[4] = Datapoints.Ms2List[i].LogIntensity;

                // y value
                double mzError = Datapoints.Ms2List[i].AbsoluteMzError;

                myMs2DataPoints.Add((explanatoryVariables, mzError));
            }

            if (myMs1DataPoints.Count * MaximumFracForTraining > MaximumDatapointsToTrainWith)
            {
                ms1fracForTraining = MaximumDatapointsToTrainWith / myMs1DataPoints.Count;
            }

            if (myMs2DataPoints.Count * MaximumFracForTraining > MaximumDatapointsToTrainWith)
            {
                ms2fracForTraining = MaximumDatapointsToTrainWith / myMs2DataPoints.Count;
            }

            Status("Generating MS1 calibration function");
            var ms1Model = GetRandomForestModel(myMs1DataPoints, ms1fracForTraining);
            //var ms1Model = GetGradientBoostModel(myMs1DataPoints, ms1fracForTraining);

            Status("Generating MS2 calibration function");
            var ms2Model = GetRandomForestModel(myMs2DataPoints, ms2fracForTraining);
            //var ms2Model = GetGradientBoostModel(myMs2DataPoints, ms2fracForTraining);

            Status("Calibrating spectra");

            CalibrateSpectra(ms1Model, ms2Model);

            return new MetaMorpheusEngineResults(this);
        }

        private void CalibrateSpectra(IPredictorModel<double> ms1predictor, IPredictorModel<double> ms2predictor)
        {
            Parallel.ForEach(Partitioner.Create(1, MyMsDataFile.NumSpectra + 1),
                new ParallelOptions { MaxDegreeOfParallelism = CommonParameters.MaxThreadsToUsePerFile },
                (fff, loopState) =>
            {
                for (int i = fff.Item1; i < fff.Item2; i++)
                {
                    // Stop loop if canceled
                    if (GlobalVariables.StopLoops)
                    {
                        loopState.Stop();
                        return;
                    }

                    var scan = MyMsDataFile.GetOneBasedScan(i);

                    if (scan.MsnOrder == 2)
                    {
                        var precursorScan = MyMsDataFile.GetOneBasedScan(scan.OneBasedPrecursorScanNumber.Value);

                        if (!scan.SelectedIonMonoisotopicGuessIntensity.HasValue && scan.SelectedIonMonoisotopicGuessMz.HasValue)
                        {
                            scan.ComputeMonoisotopicPeakIntensity(precursorScan.MassSpectrum);
                        }

                        double theFunc(MzPeak x) => x.Mz - ms2predictor.Predict(new double[] { x.Mz, scan.RetentionTime, Math.Log(scan.TotalIonCurrent), scan.InjectionTime.HasValue ? Math.Log(scan.InjectionTime.Value) : double.NaN, Math.Log(x.Intensity) });

                        double theFuncForPrecursor(MzPeak x) => x.Mz - ms1predictor.Predict(new double[] { x.Mz, precursorScan.RetentionTime, Math.Log(precursorScan.TotalIonCurrent), precursorScan.InjectionTime.HasValue ? Math.Log(precursorScan.InjectionTime.Value) : double.NaN, Math.Log(x.Intensity) });

                        scan.TransformMzs(theFunc, theFuncForPrecursor);
                    }
                    else
                    {
                        Func<MzPeak, double> theFunc = x => x.Mz - ms1predictor.Predict(new double[] { x.Mz, scan.RetentionTime, Math.Log(scan.TotalIonCurrent), scan.InjectionTime.HasValue ? Math.Log(scan.InjectionTime.Value) : double.NaN, Math.Log(x.Intensity) });
                        scan.MassSpectrum.ReplaceXbyApplyingFunction(theFunc);
                    }
                }
            });
        }

        private RegressionForestModel GetRandomForestModel(List<(double[] xValues, double yValue)> myInputs, double fracForTraining)
        {
            // create a machine learner
            var learner = new RegressionRandomForestLearner();
            var metric = new MeanAbsolutErrorRegressionMetric();

            var splitter = new RandomTrainingTestIndexSplitter<double>(trainingPercentage: fracForTraining, seed: RandomSeed);

            // put x values into a matrix and y values into a 1D array
            var myXValueMatrix = new F64Matrix(myInputs.Count, myInputs.First().xValues.Length);
            for (int i = 0; i < myInputs.Count; i++)
            {
                for (int j = 0; j < myInputs.First().xValues.Length; j++)
                {
                    myXValueMatrix[i, j] = myInputs[i].xValues[j];
                }
            }

            var myYValues = myInputs.Select(p => p.yValue).ToArray();

            // split data into training set and test set
            var splitData = splitter.SplitSet(myXValueMatrix, myYValues);
            var trainingSetX = splitData.TrainingSet.Observations;
            var trainingSetY = splitData.TrainingSet.Targets;

            // parameter ranges for the optimizer
            var parameters = new ParameterBounds[]
            {
                new ParameterBounds(min: 100, max: 200, transform: Transform.Linear),           // trees
                new ParameterBounds(min: 1, max: 5, transform: Transform.Linear),               // min split size
                new ParameterBounds(min: 2000, max: 4000, transform: Transform.Linear),          // max tree depth
                new ParameterBounds(min: 0, max: 2, transform: Transform.Linear),               // featuresPrSplit
                new ParameterBounds(min: 1e-06, max: 1e-05, transform: Transform.Logarithmic),  // min info gain
                new ParameterBounds(min: 0.7, max: 1.5, transform: Transform.Linear)            // subsample ratio
            };

            var validationSplit = new RandomTrainingTestIndexSplitter<double>(trainingPercentage: fracForTraining, seed: RandomSeed)
                .SplitSet(myXValueMatrix, myYValues);

            // define minimization metric
            Func<double[], OptimizerResult> minimize = p =>
            {
                // create the candidate learner using the current optimization parameters
                var candidateLearner = new RegressionRandomForestLearner(
                    trees: (int)p[0],
                    minimumSplitSize: (int)p[1],
                    maximumTreeDepth: (int)p[2],
                    featuresPrSplit: (int)p[3],
                    minimumInformationGain: p[4],
                    subSampleRatio: p[5],
                    seed: RandomSeed,
                    runParallel: false);

                var candidateModel = candidateLearner.Learn(validationSplit.TrainingSet.Observations,
                validationSplit.TrainingSet.Targets);

                var validationPredictions = candidateModel.Predict(validationSplit.TestSet.Observations);
                var candidateError = metric.Error(validationSplit.TestSet.Targets, validationPredictions);

                return new OptimizerResult(p, candidateError);
            };

            // create optimizer
            var optimizer = new RandomSearchOptimizer(parameters, seed: RandomSeed, iterations: TrainingIterations, runParallel: true);

            // find best parameters
            var result = optimizer.OptimizeBest(minimize);
            var best = result.ParameterSet;

            // create the final learner using the best parameters
            // (parameters that resulted in the model with the least error)
            learner = new RegressionRandomForestLearner(
                    trees: (int)best[0],
                    minimumSplitSize: (int)best[1],
                    maximumTreeDepth: (int)best[2],
                    featuresPrSplit: (int)best[3],
                    minimumInformationGain: best[4],
                    subSampleRatio: best[5],
                    seed: RandomSeed,
                    runParallel: true);

            // learn final model with optimized parameters
            var myModel = learner.Learn(trainingSetX, trainingSetY);

            // all done
            return myModel;
        }

        private RegressionGradientBoostModel GetGradientBoostModel(List<(double[] xValues, double yValue)> myInputs, double fracForTraining)
        {
            // create a machine learner
            var learner = new RegressionAbsoluteLossGradientBoostLearner();
            var metric = new MeanAbsolutErrorRegressionMetric();

            var splitter = new RandomTrainingTestIndexSplitter<double>(trainingPercentage: fracForTraining, seed: RandomSeed);

            // put x values into a matrix and y values into a 1D array
            var myXValueMatrix = new F64Matrix(myInputs.Count, myInputs.First().xValues.Length);
            for (int i = 0; i < myInputs.Count; i++)
            {
                for (int j = 0; j < myInputs.First().xValues.Length; j++)
                {
                    myXValueMatrix[i, j] = myInputs[i].xValues[j];
                }
            }

            var myYValues = myInputs.Select(p => p.yValue).ToArray();

            // split data into training set and test set
            var splitData = splitter.SplitSet(myXValueMatrix, myYValues);
            var trainingSetX = splitData.TrainingSet.Observations;
            var trainingSetY = splitData.TrainingSet.Targets;

            // learn an initial model
            var myModel = learner.Learn(trainingSetX, trainingSetY);

            // parameter ranges for the optimizer
            var parameters = new ParameterBounds[]
            {
                new ParameterBounds(min: 100, max: 300, transform: Transform.Linear),           // iterations
                new ParameterBounds(min: 0.05, max: 0.15, transform: Transform.Linear),         // learningrate
                new ParameterBounds(min: 3, max: 10, transform: Transform.Linear),              // max tree depth
                new ParameterBounds(min: 1, max: 3, transform: Transform.Linear),               // min split size
                new ParameterBounds(min: 1e-06, max: 1e-05, transform: Transform.Logarithmic),  // min info gain
                new ParameterBounds(min: 0.8, max: 1.0, transform: Transform.Linear),           // subsample ratio
                new ParameterBounds(min: 0, max: 1, transform: Transform.Linear)                // features per split
            };

            var validationSplit = new RandomTrainingTestIndexSplitter<double>(trainingPercentage: fracForTraining, seed: RandomSeed)
                .SplitSet(myXValueMatrix, myYValues);

            // define minimization metric
            Func<double[], OptimizerResult> minimize = p =>
            {
                // create the candidate learner using the current optimization parameters
                var candidateLearner = new RegressionAbsoluteLossGradientBoostLearner(
                    iterations: (int)p[0],
                    learningRate: p[1],
                    maximumTreeDepth: (int)p[2],
                    minimumSplitSize: (int)p[3],
                    minimumInformationGain: p[4],
                    subSampleRatio: p[5],
                    featuresPrSplit: (int)p[6],
                    runParallel: false);

                var candidateModel = candidateLearner.Learn(validationSplit.TrainingSet.Observations,
                validationSplit.TrainingSet.Targets);

                var validationPredictions = candidateModel.Predict(validationSplit.TestSet.Observations);
                var candidateError = metric.Error(validationSplit.TestSet.Targets, validationPredictions);

                return new OptimizerResult(p, candidateError);
            };

            // create optimizer
            var optimizer = new RandomSearchOptimizer(parameters, seed: RandomSeed, iterations: TrainingIterations, runParallel: true);

            // find best parameters
            var result = optimizer.OptimizeBest(minimize);
            var best = result.ParameterSet;

            // create the final learner using the best parameters
            // (parameters that resulted in the model with the least error)
            learner = new RegressionAbsoluteLossGradientBoostLearner(
                    iterations: (int)best[0],
                    learningRate: best[1],
                    maximumTreeDepth: (int)best[2],
                    minimumSplitSize: (int)best[3],
                    minimumInformationGain: best[4],
                    subSampleRatio: best[5],
                    featuresPrSplit: (int)best[6],
                    runParallel: true);

            // learn final model with optimized parameters
            myModel = learner.Learn(trainingSetX, trainingSetY);

            // all done
            return myModel;
        }
    }
}