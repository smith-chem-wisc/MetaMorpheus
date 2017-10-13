using MassSpectrometry;
using SharpLearning.Common.Interfaces;
using SharpLearning.Containers.Matrices;
using SharpLearning.CrossValidation.TrainingTestSplitters;
using SharpLearning.Metrics.Regression;
using Spectra;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace EngineLayer.Calibration
{
    public class CalibrationEngine : MetaMorpheusEngine
    {
        #region Private Fields

        private const double fracForTraining = 0.75;

        private readonly IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> myMsDataFile;

        private readonly List<ILearner<double>> learners;
        private readonly string varsToUse;
        private readonly DataPointAquisitionResults datapoints;

        #endregion Private Fields

        #region Public Constructors

        public CalibrationEngine(IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> myMSDataFile, DataPointAquisitionResults datapoints, List<ILearner<double>> learners, string varsToUse, List<string> nestedIds) : base(nestedIds)
        {
            this.myMsDataFile = myMSDataFile;
            this.datapoints = datapoints;
            this.learners = learners;
            this.varsToUse = varsToUse;
        }

        #endregion Public Constructors

        #region Protected Methods

        protected override MetaMorpheusEngineResults RunSpecific()
        {
            var splitter = new RandomTrainingTestIndexSplitter<double>(fracForTraining);

            TrainingTestSetSplit ms1splitResult;
            TrainingTestSetSplit ms2splitResult;

            switch (varsToUse)
            {
                case "mz":
                    ms1splitResult = splitter.SplitSet(new F64Matrix(datapoints.Ms1List.SelectMany(b => new[] { b.mz }).ToArray(), datapoints.Ms1List.Count, 1), datapoints.Ms1List.Select(b => b.LabelTh).ToArray());
                    ms2splitResult = splitter.SplitSet(new F64Matrix(datapoints.Ms2List.SelectMany(b => new[] { b.mz }).ToArray(), datapoints.Ms2List.Count, 1), datapoints.Ms2List.Select(b => b.LabelTh).ToArray());
                    break;

                case "mzRtTicInj":
                    ms1splitResult = splitter.SplitSet(new F64Matrix(datapoints.Ms1List.SelectMany(b => new[] { b.mz, b.rt, b.logTotalIonCurrent, b.logInjectionTime }).ToArray(), datapoints.Ms1List.Count, 4), datapoints.Ms1List.Select(b => b.LabelTh).ToArray());
                    ms2splitResult = splitter.SplitSet(new F64Matrix(datapoints.Ms2List.SelectMany(b => new[] { b.mz, b.rt, b.logTotalIonCurrent, b.logInjectionTime }).ToArray(), datapoints.Ms2List.Count, 4), datapoints.Ms2List.Select(b => b.LabelTh).ToArray());
                    break;

                case "mzRtTicInjInt":
                    ms1splitResult = splitter.SplitSet(new F64Matrix(datapoints.Ms1List.SelectMany(b => new[] { b.mz, b.rt, b.logTotalIonCurrent, b.logInjectionTime, b.logIntensity }).ToArray(), datapoints.Ms1List.Count, 5), datapoints.Ms1List.Select(b => b.LabelTh).ToArray());
                    ms2splitResult = splitter.SplitSet(new F64Matrix(datapoints.Ms2List.SelectMany(b => new[] { b.mz, b.rt, b.logTotalIonCurrent, b.logInjectionTime, b.logIntensity }).ToArray(), datapoints.Ms2List.Count, 5), datapoints.Ms2List.Select(b => b.LabelTh).ToArray());
                    break;

                case "mzInt":
                    ms1splitResult = splitter.SplitSet(new F64Matrix(datapoints.Ms1List.SelectMany(b => new[] { b.mz, b.logIntensity }).ToArray(), datapoints.Ms1List.Count, 2), datapoints.Ms1List.Select(b => b.LabelTh).ToArray());
                    ms2splitResult = splitter.SplitSet(new F64Matrix(datapoints.Ms2List.SelectMany(b => new[] { b.mz, b.logIntensity }).ToArray(), datapoints.Ms2List.Count, 2), datapoints.Ms2List.Select(b => b.LabelTh).ToArray());
                    break;

                case "Int":
                    ms1splitResult = splitter.SplitSet(new F64Matrix(datapoints.Ms1List.SelectMany(b => new[] { b.logIntensity }).ToArray(), datapoints.Ms1List.Count, 1), datapoints.Ms1List.Select(b => b.LabelTh).ToArray());
                    ms2splitResult = splitter.SplitSet(new F64Matrix(datapoints.Ms2List.SelectMany(b => new[] { b.logIntensity }).ToArray(), datapoints.Ms2List.Count, 1), datapoints.Ms2List.Select(b => b.LabelTh).ToArray());
                    break;

                default:
                    throw new MetaMorpheusException("unknown vars to use array: " + varsToUse);
            }
            MS1Predictor ms1predictor = new MS1Predictor(DoStuff(ms1splitResult, learners));

            MS2Predictor ms2predictor = new MS2Predictor(DoStuff(ms2splitResult, learners));

            Status("Calibrating Spectra");

            CalibrateSpectra(ms1predictor, ms2predictor);

            return new MetaMorpheusEngineResults(this);
        }

        #endregion Protected Methods

        #region Private Methods

        private static IPredictorModel<double> DoStuff(TrainingTestSetSplit splitResult, List<ILearner<double>> learners)
        {
            var evaluator = new MeanAbsolutErrorRegressionMetric();

            var predictions = new double[splitResult.TestSet.Targets.Length];
            IPredictorModel<double> bestModel = new IdentityCalibrationFunctionPredictorModel();
            var bestError = evaluator.Error(splitResult.TestSet.Targets, predictions);

            foreach (var learner in learners)
            {
                try
                {
                    var model = learner.Learn(splitResult.TrainingSet.Observations, splitResult.TrainingSet.Targets);

                    predictions = new double[splitResult.TestSet.Targets.Length];
                    for (int i = 0; i < splitResult.TestSet.Targets.Length; i++)
                        predictions[i] = model.Predict(splitResult.TestSet.Observations.Row(i));

                    var thisError = evaluator.Error(splitResult.TestSet.Targets, predictions);

                    if (thisError < bestError)
                    {
                        bestError = thisError;
                        bestModel = model;
                    }
                }
                catch
                {
                }
            }

            return bestModel;
        }

        private void CalibrateSpectra(MS1Predictor ms1predictor, MS2Predictor ms2predictor)
        {
            Parallel.ForEach(Partitioner.Create(1, myMsDataFile.NumSpectra + 1), fff =>
              {
                  for (int i = fff.Item1; i < fff.Item2; i++)
                  {
                      var a = myMsDataFile.GetOneBasedScan(i);

                      if (a is IMsDataScanWithPrecursor<IMzSpectrum<IMzPeak>> theScan)
                      {
                          var precursorScan = myMsDataFile.GetOneBasedScan(theScan.OneBasedPrecursorScanNumber.Value);

                          if (!theScan.SelectedIonMonoisotopicGuessIntensity.HasValue)
                              theScan.ComputeMonoisotopicPeakIntensity(precursorScan.MassSpectrum);

                          Func<IPeak, double> theFunc = x => x.X - ms2predictor.Predict(x.X, a.RetentionTime, Math.Log(a.TotalIonCurrent), a.InjectionTime.HasValue ? Math.Log(a.InjectionTime.Value) : double.NaN, Math.Log(x.Y));

                          Func<IPeak, double> theFuncForPrecursor = x => x.X - ms1predictor.Predict(x.X, precursorScan.RetentionTime, Math.Log(precursorScan.TotalIonCurrent), precursorScan.InjectionTime.HasValue ? Math.Log(precursorScan.InjectionTime.Value) : double.NaN, Math.Log(x.Y));

                          theScan.TransformMzs(theFunc, theFuncForPrecursor);
                      }
                      else
                      {
                          Func<IPeak, double> theFunc = x => x.X - ms1predictor.Predict(x.X, a.RetentionTime, Math.Log(a.TotalIonCurrent), a.InjectionTime.HasValue ? Math.Log(a.InjectionTime.Value) : double.NaN, Math.Log(x.Y));
                          a.MassSpectrum.ReplaceXbyApplyingFunction(theFunc);
                      }
                  }
              }
              );
        }
    }

    #endregion Private Methods
}