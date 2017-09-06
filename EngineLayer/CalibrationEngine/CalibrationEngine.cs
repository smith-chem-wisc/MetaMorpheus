using Chemistry;
using MassSpectrometry;
using MathNet.Numerics.Statistics;
using MzLibUtil;
using Proteomics;
using SharpLearning.Common.Interfaces;
using SharpLearning.Containers.Matrices;
using SharpLearning.Metrics.Regression;
using Spectra;
using System;
using System.Collections.Generic;
using System.Linq;

namespace EngineLayer.Calibration
{
    public class CalibrationEngine : MetaMorpheusEngine
    {
        #region Private Fields

        private const double fineResolutionForIsotopeDistCalculation = 0.1;
        private const double fracForTraining = 0.75;

        private readonly int minMS1isotopicPeaksNeededForConfirmedIdentification;
        private readonly int minMS2isotopicPeaksNeededForConfirmedIdentification;
        private readonly int numFragmentsNeededForEveryIdentification;
        private readonly Tolerance mzToleranceForMs1Search;
        private readonly Tolerance mzToleranceForMs2Search;
        private readonly FragmentTypes fragmentTypesForCalibration;
        private readonly Action<List<LabeledMs1DataPoint>, string> ms1ListAction;
        private readonly Action<List<LabeledMs2DataPoint>, string> ms2ListAction;
        private readonly bool doForestCalibration;
        private readonly List<Psm> goodIdentifications;
        private readonly IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> myMsDataFile;
        private readonly Random rnd;
        private readonly int doNotSplitIfUnderThis;
        private readonly int numTrees;
        private int numMs1MassChargeCombinationsConsidered;

        private int numMs1MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks;

        private int numMs2MassChargeCombinationsConsidered;

        private int numMs2MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks;

        private int numFragmentsIdentified;

        #endregion Private Fields

        #region Public Constructors

        public CalibrationEngine(IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> myMSDataFile, Tolerance mzToleranceForMs2Search, List<Psm> goodIdentifications, int minMS1IsotopicPeaksNeededForConfirmedIdentification, int minMS2IsotopicPeaksNeededForConfirmedIdentification, int numFragmentsNeededForEveryIdentification, Tolerance mzToleranceForMs1Search, FragmentTypes fragmentTypesForCalibration, Action<List<LabeledMs1DataPoint>, string> ms1ListAction, Action<List<LabeledMs2DataPoint>, string> ms2ListAction, bool doForestCalibration, Random rnd, List<string> nestedIds, int doNotSplitIfUnderThis, int numTrees) : base(nestedIds)
        {
            this.myMsDataFile = myMSDataFile;
            this.goodIdentifications = goodIdentifications;
            this.minMS1isotopicPeaksNeededForConfirmedIdentification = minMS1IsotopicPeaksNeededForConfirmedIdentification;
            this.minMS2isotopicPeaksNeededForConfirmedIdentification = minMS2IsotopicPeaksNeededForConfirmedIdentification;
            this.numFragmentsNeededForEveryIdentification = numFragmentsNeededForEveryIdentification;
            this.mzToleranceForMs1Search = mzToleranceForMs1Search;
            this.mzToleranceForMs2Search = mzToleranceForMs2Search;
            this.fragmentTypesForCalibration = fragmentTypesForCalibration;
            this.ms1ListAction = ms1ListAction;
            this.ms2ListAction = ms2ListAction;
            this.doForestCalibration = doForestCalibration;
            this.rnd = rnd;
            this.doNotSplitIfUnderThis = doNotSplitIfUnderThis;
            this.numTrees = numTrees;
        }

        #endregion Public Constructors

        #region Protected Methods

        protected override MetaMorpheusEngineResults RunSpecific()
        {
            Status("Calibrating ");
            var trainingPointCounts = new List<int>();
            var result = new CalibrationResults(myMsDataFile, this);
            DataPointAquisitionResults dataPointAcquisitionResult = null;
            for (int linearCalibrationRound = 1; ; linearCalibrationRound++)
            {
                Status("linearCalibrationRound " + linearCalibrationRound);
                dataPointAcquisitionResult = GetDataPoints();
                ms1ListAction(dataPointAcquisitionResult.Ms1List, "beforesc" + linearCalibrationRound.ToString());
                ms2ListAction(dataPointAcquisitionResult.Ms2List, "beforesc" + linearCalibrationRound.ToString());
                result.Add(dataPointAcquisitionResult);
                if (linearCalibrationRound >= 2 && dataPointAcquisitionResult.Count <= trainingPointCounts[linearCalibrationRound - 2])
                    break;
                trainingPointCounts.Add(dataPointAcquisitionResult.Count);
                if (dataPointAcquisitionResult.Ms2List.Count == 0)
                    return new MyErroredResults(this, "No MS2 training points, identification quality is poor. Try to change the Fragment tolerance." + result.ToString());
                if (dataPointAcquisitionResult.Ms1List.Count == 0)
                    return new MyErroredResults(this, "No MS1 training points, identification quality is poor. Try to change the Parent tolerance." + result.ToString());
                Tuple<IPredictorModel<double>, IPredictorModel<double>> combinedCalibration = CalibrateLinear(dataPointAcquisitionResult);
                result.Add(combinedCalibration.Item1, combinedCalibration.Item2);
            }
            if (doForestCalibration)
            {
                trainingPointCounts = new List<int>();
                for (int forestCalibrationRound = 1; ; forestCalibrationRound++)
                {
                    Status("forestCalibrationRound " + forestCalibrationRound);
                    Console.WriteLine("forestCalibrationRound " + forestCalibrationRound);
                    Tuple<IPredictorModel<double>, IPredictorModel<double>> combinedCalibration = CalibrateRF(dataPointAcquisitionResult);
                    result.Add(combinedCalibration.Item1, combinedCalibration.Item2);
                    dataPointAcquisitionResult = GetDataPoints();
                    ms1ListAction(dataPointAcquisitionResult.Ms1List, "afterfc" + forestCalibrationRound.ToString());
                    ms2ListAction(dataPointAcquisitionResult.Ms2List, "afterfc" + forestCalibrationRound.ToString());
                    result.Add(dataPointAcquisitionResult);
                    if (forestCalibrationRound >= 2 && dataPointAcquisitionResult.Count <= trainingPointCounts[forestCalibrationRound - 2])
                        break;
                    trainingPointCounts.Add(dataPointAcquisitionResult.Count);
                    if (dataPointAcquisitionResult.Ms2List.Count == 0)
                        return new MyErroredResults(this, "No MS2 training points, identification quality is poor. Try to change the Fragment tolerance." + result.ToString());
                    if (dataPointAcquisitionResult.Ms1List.Count == 0)
                        return new MyErroredResults(this, "No MS1 training points, identification quality is poor. Try to change the Parent tolerance." + result.ToString());
                }
            }

            return result;
        }

        #endregion Protected Methods

        #region Private Methods

        private Tuple<IPredictorModel<double>, IPredictorModel<double>> CalibrateRF(DataPointAquisitionResults res)
        {
            var metric = new MeanSquaredErrorRegressionMetric();

            Console.WriteLine("  In CalibrateRF");
            var shuffledMs1TrainingPoints = res.Ms1List.OrderBy(item => rnd.Next()).ToList();
            var shuffledMs2TrainingPoints = res.Ms2List.OrderBy(item => rnd.Next()).ToList();

            Console.WriteLine("  All ms2 Points");
            Console.WriteLine(string.Join(Environment.NewLine, shuffledMs2TrainingPoints.OrderBy(b => b.mz).Select(b => "      " + b.mz + " " + b.Label)));

            var trainList1 = shuffledMs1TrainingPoints.Take((int)(shuffledMs1TrainingPoints.Count * fracForTraining)).ToList();
            var testList1 = shuffledMs1TrainingPoints.Skip((int)(shuffledMs1TrainingPoints.Count * fracForTraining)).ToList();
            var trainList2 = shuffledMs2TrainingPoints.Take((int)(shuffledMs2TrainingPoints.Count * fracForTraining)).ToList();
            Warn("  train ms2 Points");
            Console.WriteLine(string.Join(Environment.NewLine, trainList2.OrderBy(b => b.mz).Select(b => "      " + b.mz + " " + b.Label)));

            var testList2 = shuffledMs2TrainingPoints.Skip((int)(shuffledMs2TrainingPoints.Count * fracForTraining)).ToList();
            Console.WriteLine("  testList2 ms2 Points");
            Console.WriteLine(string.Join(Environment.NewLine, testList2.OrderBy(b => b.mz).Select(b => "      " + b.mz + " " + b.Label)));

            IPredictorModel<double> bestMS1predictor;
            double bestMS1MSE;
            IPredictorModel<double> bestMS2predictor;
            double bestMS2MSE;

            {
                var b = new IdentityCalibrationFunction();
                bestMS1predictor = b.Learn(null, null);
                bestMS1MSE = testList1.Select(bc => Math.Pow(bc.Label, 2)).Average();
            }

            {
                var b = new IdentityCalibrationFunction();
                bestMS2predictor = b.Learn(null, null);
                bestMS2MSE = testList2.Select(bc => Math.Pow(bc.Label, 2)).Average();
            }

            Console.WriteLine("  bestMS2MSE identity: " + bestMS2MSE);

            try
            {
                double[] trainList1Concat = trainList1.SelectMany(b => b.Inputs).ToArray();
                F64Matrix trainList1Matrix = new F64Matrix(trainList1Concat, trainList1.Count, trainList1Concat.Length / trainList1.Count);
                double[] trainList1Targets = trainList1.Select(b => b.Label).ToArray();

                double[] testList1Concat = testList1.SelectMany(b => b.Inputs).ToArray();
                F64Matrix testList1Matrix = new F64Matrix(testList1Concat, testList1.Count, testList1Concat.Length / testList1.Count);
                double[] testList1Targets = testList1.Select(b => b.Label).ToArray();

                List<bool[]> boolStuffms1 = new List<bool[]>
            {
                new bool[] {true, true, false, false, false},
                new bool[] {true, true, true, true, true},
            };
                if (trainList1.Count > 0)
                    foreach (var boolStuff in boolStuffms1)
                    {
                        var ms1predictor = new RandomForestCalibrationFunction(numTrees, doNotSplitIfUnderThis, boolStuff);
                        var model = ms1predictor.Learn(trainList1Matrix, trainList1Targets);
                        var predictions = new double[testList1Targets.Length];
                        for (int i = 0; i < testList1Targets.Length; i++)
                            predictions[i] = model.Predict(testList1Matrix.Row(i));
                        var MS1mse = metric.Error(testList1Targets, predictions);
                        if (MS1mse < bestMS1MSE)
                        {
                            bestMS1MSE = MS1mse;
                            bestMS1predictor = model;
                        }
                    }
            }
            catch
            {
            }

            try
            {
                double[] trainList2Concat = trainList2.SelectMany(b => b.Inputs).ToArray();
                F64Matrix trainList2Matrix = new F64Matrix(trainList2Concat, trainList2.Count, trainList2Concat.Length / trainList2.Count);
                double[] trainList2Targets = trainList2.Select(b => b.Label).ToArray();

                double[] testList2Concat = testList2.SelectMany(b => b.Inputs).ToArray();
                F64Matrix testList2Matrix = new F64Matrix(testList2Concat, testList2.Count, testList2Concat.Length / testList2.Count);
                double[] testList2Targets = testList2.Select(b => b.Label).ToArray();

                List<bool[]> boolStuffms2 = new List<bool[]>
            {
                new bool[] {true, false, false, false, false, false},
                new bool[] {true, true,  false, false, false, false},
                new bool[] {true, true,  true,  true,  true,  false},
                new bool[] {true, true,  true,  true,  true,  true },
            };
                if (trainList2.Count > 0)
                    foreach (var boolStuff in boolStuffms2)
                    {
                        var ms2predictor = new RandomForestCalibrationFunction(numTrees, doNotSplitIfUnderThis, boolStuff);
                        var model = ms2predictor.Learn(trainList2Matrix, trainList2Targets);
                        var predictions = new double[testList2Targets.Length];
                        for (int i = 0; i < testList2Targets.Length; i++)
                            predictions[i] = model.Predict(testList2Matrix.Row(i));
                        var MS2mse = metric.Error(testList2Targets, predictions);
                        if (MS2mse < bestMS2MSE)
                        {
                            bestMS2MSE = MS2mse;
                            bestMS2predictor = model;
                        }
                    }
            }
            catch
            {
            }

            Tuple<IPredictorModel<double>, IPredictorModel<double>> bestCf = new Tuple<IPredictorModel<double>, IPredictorModel<double>>(bestMS1predictor, bestMS2predictor);

            Status("Calibrating Spectra");

            CalibrateSpectra(bestCf);

            Console.WriteLine("  Done with CalibrateRF");

            return bestCf;
        }

        private DataPointAquisitionResults GetDataPoints()
        {
            Status("Extracting data points:");
            Console.WriteLine("  In GetDataPoints");
            // The final training point list
            DataPointAquisitionResults res = new DataPointAquisitionResults()
            {
                Ms1List = new List<LabeledMs1DataPoint>(),
                Ms2List = new List<LabeledMs2DataPoint>()
            };

            // Set of peaks, identified by m/z and retention time. If a peak is in here, it means it has been a part of an accepted identification, and should be rejected
            var peaksAddedFromMS1HashSet = new HashSet<Tuple<double, double>>();

            int numIdentifications = goodIdentifications.Count;
            // Loop over identifications
            for (int matchIndex = 0; matchIndex < numIdentifications; matchIndex++)
            {
                Psm identification = goodIdentifications[matchIndex];

                // Progress
                if (numIdentifications < 100 || matchIndex % (numIdentifications / 100) == 0)
                    ReportProgress(new ProgressEventArgs(100 * matchIndex / numIdentifications, "Looking at identifications...", nestedIds));

                // Each identification has an MS2 spectrum attached to it.
                int ms2scanNumber = identification.ScanNumber;
                int peptideCharge = identification.ScanPrecursorCharge;
                if (identification.FullSequence == null)
                    continue;

                var representativeSinglePeptide = identification.CompactPeptides.First().Value.Item2.First();

                // Get the peptide, don't forget to add the modifications!!!!
                var SequenceWithChemicalFormulas = representativeSinglePeptide.SequenceWithChemicalFormulas;
                if (SequenceWithChemicalFormulas == null || representativeSinglePeptide.allModsOneIsNterminus.Any(b => b.Value.neutralLosses.Count != 1 || b.Value.neutralLosses.First() != 0))
                    continue;
                Proteomics.Peptide coolPeptide = new Proteomics.Peptide(SequenceWithChemicalFormulas);

                numMs2MassChargeCombinationsConsidered = 0;
                numMs2MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks = 0;
                numFragmentsIdentified = 0;
                res.Ms2List.AddRange(SearchMS2Spectrum(myMsDataFile.GetOneBasedScan(ms2scanNumber) as IMsDataScanWithPrecursor<IMzSpectrum<IMzPeak>>, coolPeptide, peptideCharge, identification));
                res.numMs2MassChargeCombinationsConsidered += numMs2MassChargeCombinationsConsidered;
                res.numMs2MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks += numMs2MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks;

                // If MS2 has low evidence for peptide, skip and go to next one
                if (numFragmentsIdentified < numFragmentsNeededForEveryIdentification)
                    continue;

                // Calculate isotopic distribution of the full peptide

                var dist = IsotopicDistribution.GetDistribution(coolPeptide.GetChemicalFormula(), fineResolutionForIsotopeDistCalculation, 0.001);

                double[] masses = dist.Masses.ToArray();
                double[] intensities = dist.Intensities.ToArray();

                Array.Sort(intensities, masses, Comparer<double>.Create((x, y) => y.CompareTo(x)));

                numMs1MassChargeCombinationsConsidered = 0;
                numMs1MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks = 0;
                res.Ms1List.AddRange(SearchMS1Spectra(masses, intensities, ms2scanNumber, -1, peaksAddedFromMS1HashSet, peptideCharge, identification));
                res.Ms1List.AddRange(SearchMS1Spectra(masses, intensities, ms2scanNumber, 1, peaksAddedFromMS1HashSet, peptideCharge, identification));
                res.numMs1MassChargeCombinationsConsidered += numMs1MassChargeCombinationsConsidered;
                res.numMs1MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks += numMs1MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks;
            }
            Status("Number of training points: " + res.Count);

            Console.WriteLine("  Done with GetDataPoints");

            return res;
        }

        private Tuple<IPredictorModel<double>, IPredictorModel<double>> CalibrateLinear(DataPointAquisitionResults res)
        {
            var metric = new MeanSquaredErrorRegressionMetric();

            var shuffledMs1TrainingPoints = res.Ms1List.OrderBy(item => rnd.Next()).ToList();
            var shuffledMs2TrainingPoints = res.Ms2List.OrderBy(item => rnd.Next()).ToList();

            var trainList1 = shuffledMs1TrainingPoints.Take((int)(shuffledMs1TrainingPoints.Count * fracForTraining)).ToList();
            var testList1 = shuffledMs1TrainingPoints.Skip((int)(shuffledMs1TrainingPoints.Count * fracForTraining)).ToList();
            var trainList2 = shuffledMs2TrainingPoints.Take((int)(shuffledMs2TrainingPoints.Count * fracForTraining)).ToList();
            var testList2 = shuffledMs2TrainingPoints.Skip((int)(shuffledMs2TrainingPoints.Count * fracForTraining)).ToList();

            IPredictorModel<double> bestMS1predictor;
            double bestMS1MSE;
            IPredictorModel<double> bestMS2predictor;
            double bestMS2MSE;

            {
                var b = new IdentityCalibrationFunction();
                bestMS1predictor = b.Learn(null, null);
                bestMS1MSE = testList1.Select(bc => Math.Pow(bc.Label, 2)).Average();
            }

            {
                var b = new IdentityCalibrationFunction();
                bestMS2predictor = b.Learn(null, null);
                bestMS2MSE = testList2.Select(bc => Math.Pow(bc.Label, 2)).Average();
            }

            var transforms = new List<TransformFunction>
            {
                new TransformFunction(b => new double[] { b[0] }, 1),
                new TransformFunction(b => new double[] { b[1] }, 1),
                new TransformFunction(b => new double[] { Math.Log(b[3]) }, 1),
                new TransformFunction(b => new double[] { Math.Log(b[4]) }, 1),

                new TransformFunction(b => new double[] { b[0], b[1] }, 2),
                new TransformFunction(b => new double[] { b[0], Math.Log(b[3]) }, 2),
                new TransformFunction(b => new double[] { b[0], Math.Log(b[4]) }, 2),
                new TransformFunction(b => new double[] { b[1], Math.Log(b[3]) }, 2),
                new TransformFunction(b => new double[] { b[1], Math.Log(b[4]) }, 2),
                new TransformFunction(b => new double[] { Math.Log(b[3]), Math.Log(b[4]) }, 2),

                new TransformFunction(b => new double[] { b[0], b[1], Math.Log(b[3]) }, 3),
                new TransformFunction(b => new double[] { b[0], b[1], Math.Log(b[4]) }, 3),
                new TransformFunction(b => new double[] { b[0], Math.Log(b[3]), Math.Log(b[4]) }, 3),
                new TransformFunction(b => new double[] { b[1], Math.Log(b[3]), Math.Log(b[4]) }, 3),

                new TransformFunction(b => new double[] { b[0], b[1], Math.Log(b[3]), Math.Log(b[4]) }, 4),
            };

            try
            {
                double[] trainList1Concat = trainList1.SelectMany(b => b.Inputs).ToArray();
                F64Matrix trainList1Matrix = new F64Matrix(trainList1Concat, trainList1.Count, trainList1Concat.Length / trainList1.Count);
                double[] trainList1Targets = trainList1.Select(b => b.Label).ToArray();

                double[] testList1Concat = testList1.SelectMany(b => b.Inputs).ToArray();
                F64Matrix testList1Matrix = new F64Matrix(testList1Concat, testList1.Count, testList1Concat.Length / testList1.Count);
                double[] testList1Targets = testList1.Select(b => b.Label).ToArray();
                {
                    var ms1predictor = new ConstantCalibrationFunction();
                    var model = ms1predictor.Learn(trainList1Matrix, trainList1Targets);
                    var predictions = new double[testList1Targets.Length];
                    for (int i = 0; i < testList1Targets.Length; i++)
                        predictions[i] = model.Predict(testList1Matrix.Row(i));
                    var MS1mse = metric.Error(testList1Targets, predictions);
                    if (MS1mse < bestMS1MSE)
                    {
                        bestMS1MSE = MS1mse;
                        bestMS1predictor = model;
                    }
                }
                foreach (var transform in transforms)
                {
                    try
                    {
                        var ms1predictor = new LinearCalibrationFunctionMathNet(transform);
                        var model = ms1predictor.Learn(trainList1Matrix, trainList1Targets);
                        var predictions = new double[testList1Targets.Length];
                        for (int i = 0; i < testList1Targets.Length; i++)
                            predictions[i] = model.Predict(testList1Matrix.Row(i));
                        var MS1mse = metric.Error(testList1Targets, predictions);
                        if (MS1mse < bestMS1MSE)
                        {
                            bestMS1MSE = MS1mse;
                            bestMS1predictor = model;
                        }
                    }
                    catch
                    {
                    }
                }
            }
            catch
            {
            }
            try
            {
                double[] trainList2Concat = trainList2.SelectMany(b => b.Inputs).ToArray();
                F64Matrix trainList2Matrix = new F64Matrix(trainList2Concat, trainList2.Count, trainList2Concat.Length / trainList2.Count);
                double[] trainList2Targets = trainList2.Select(b => b.Label).ToArray();

                double[] testList2Concat = testList2.SelectMany(b => b.Inputs).ToArray();
                F64Matrix testList2Matrix = new F64Matrix(testList2Concat, testList2.Count, testList2Concat.Length / testList2.Count);
                double[] testList2Targets = testList2.Select(b => b.Label).ToArray();

                {
                    var ms2predictor = new ConstantCalibrationFunction();
                    var model = ms2predictor.Learn(trainList2Matrix, trainList2Targets);
                    var predictions = new double[testList2Targets.Length];
                    for (int i = 0; i < testList2Targets.Length; i++)
                        predictions[i] = model.Predict(testList2Matrix.Row(i));
                    var MS2mse = metric.Error(testList2Targets, predictions);
                    if (MS2mse < bestMS2MSE)
                    {
                        bestMS2MSE = MS2mse;
                        bestMS2predictor = model;
                    }
                }

                foreach (var transform in transforms)
                {
                    try
                    {
                        var ms2predictor = new LinearCalibrationFunctionMathNet(transform);
                        var model = ms2predictor.Learn(trainList2Matrix, trainList2Targets);
                        var predictions = new double[testList2Targets.Length];
                        for (int i = 0; i < testList2Targets.Length; i++)
                            predictions[i] = model.Predict(testList2Matrix.Row(i));
                        var MS2mse = metric.Error(testList2Targets, predictions);
                        if (MS2mse < bestMS2MSE)
                        {
                            bestMS2MSE = MS2mse;
                            bestMS2predictor = model;
                        }
                    }
                    catch
                    {
                    }
                }
            }
            catch
            {
            }

            Tuple<IPredictorModel<double>, IPredictorModel<double>> bestCf = new Tuple<IPredictorModel<double>, IPredictorModel<double>>(bestMS1predictor, bestMS2predictor);

            Status("Calibrating Spectra");

            CalibrateSpectra(bestCf);

            return bestCf;
        }

        private void CalibrateSpectra(Tuple<IPredictorModel<double>, IPredictorModel<double>> bestCf)
        {
            foreach (var a in myMsDataFile)
            {
                if (a is IMsDataScanWithPrecursor<IMzSpectrum<IMzPeak>> theScan)
                {
                    var precursorScan = myMsDataFile.GetOneBasedScan(theScan.OneBasedPrecursorScanNumber);

                    Func<IPeak, double> theFunc = x => x.X - bestCf.Item2.Predict(new double[] { x.X, a.RetentionTime, x.Y, a.TotalIonCurrent, a.InjectionTime ?? double.NaN, theScan.IsolationMz });

                    Func<IPeak, double> theFuncForPrecursor = x => x.X - bestCf.Item1.Predict(new double[] { x.X, precursorScan.RetentionTime, x.Y, precursorScan.TotalIonCurrent, precursorScan.InjectionTime ?? double.NaN });

                    theScan.TransformMzs(theFunc, theFuncForPrecursor);
                }
                else
                {
                    Func<IPeak, double> theFunc = x => x.X - bestCf.Item1.Predict(new double[] { x.X, a.RetentionTime, x.Y, a.TotalIonCurrent, a.InjectionTime ?? double.NaN });
                    a.MassSpectrum.ReplaceXbyApplyingFunction(theFunc);
                }
            }
        }

        private IEnumerable<LabeledMs1DataPoint> SearchMS1Spectra(double[] originalMasses, double[] originalIntensities, int ms2spectrumIndex, int direction, HashSet<Tuple<double, double>> peaksAddedHashSet, int peptideCharge, Psm identification)
        {
            int theIndex;
            if (direction == 1)
                theIndex = ms2spectrumIndex;
            else
                theIndex = ms2spectrumIndex - 1;

            bool addedAscan = true;

            int highestKnownChargeForThisPeptide = peptideCharge;
            while (theIndex >= 1 && theIndex <= myMsDataFile.NumSpectra && addedAscan)
            {
                int countForThisScan = 0;
                if (myMsDataFile.GetOneBasedScan(theIndex).MsnOrder > 1)
                {
                    theIndex += direction;
                    continue;
                }
                addedAscan = false;
                var fullMS1scan = myMsDataFile.GetOneBasedScan(theIndex);
                double ms1RetentionTime = fullMS1scan.RetentionTime;
                var scanWindowRange = fullMS1scan.ScanWindowRange;
                var fullMS1spectrum = fullMS1scan.MassSpectrum;
                if (fullMS1spectrum.Size == 0)
                    break;

                bool startingToAddCharges = false;
                int chargeToLookAt = 1;
                do
                {
                    if (originalMasses[0].ToMz(chargeToLookAt) > scanWindowRange.Maximum)
                    {
                        chargeToLookAt++;
                        continue;
                    }
                    if (originalMasses[0].ToMz(chargeToLookAt) < scanWindowRange.Minimum)
                        break;
                    var trainingPointsToAverage = new List<LabeledMs1DataPoint>();
                    foreach (double a in originalMasses)
                    {
                        double theMZ = a.ToMz(chargeToLookAt);

                        var npwr = fullMS1spectrum.NumPeaksWithinRange(mzToleranceForMs1Search.GetMinimumValue(theMZ), mzToleranceForMs1Search.GetMaximumValue(theMZ));
                        if (npwr == 0)
                        {
                            break;
                        }
                        numMs1MassChargeCombinationsConsidered++;
                        if (npwr > 1)
                        {
                            numMs1MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks++;
                            continue;
                        }

                        int closestPeakIndex = fullMS1spectrum.GetClosestPeakIndex(theMZ);
                        var closestPeakMZ = fullMS1spectrum.XArray[closestPeakIndex];

                        var theTuple = Tuple.Create(closestPeakMZ, ms1RetentionTime);
                        if (!peaksAddedHashSet.Contains(theTuple))
                        {
                            peaksAddedHashSet.Add(theTuple);
                            highestKnownChargeForThisPeptide = Math.Max(highestKnownChargeForThisPeptide, chargeToLookAt);
                            trainingPointsToAverage.Add(new LabeledMs1DataPoint(closestPeakMZ, double.NaN, fullMS1spectrum.YArray[closestPeakIndex], double.NaN, double.NaN, closestPeakMZ - theMZ, null));
                        }
                        else
                            break;
                    }
                    // If started adding and suddnely stopped, go to next one, no need to look at higher charges
                    if (trainingPointsToAverage.Count == 0 && startingToAddCharges == true)
                    {
                        break;
                    }
                    if ((trainingPointsToAverage.Count == 0 || (trainingPointsToAverage.Count == 1 && originalIntensities[0] < 0.65)) && (peptideCharge <= chargeToLookAt))
                    {
                        break;
                    }
                    if ((trainingPointsToAverage.Count == 1 && originalIntensities[0] < 0.65)
                        || trainingPointsToAverage.Count < Math.Min(minMS1isotopicPeaksNeededForConfirmedIdentification, originalIntensities.Count()))
                    {
                    }
                    else
                    {
                        addedAscan = true;
                        startingToAddCharges = true;
                        countForThisScan++;
                        yield return new LabeledMs1DataPoint(trainingPointsToAverage.Select(b => b.mz).Average(),
                                                             fullMS1scan.RetentionTime,
                                                             trainingPointsToAverage.Select(b => b.intensity).Average(),
                                                             fullMS1scan.TotalIonCurrent,
                                                             fullMS1scan.InjectionTime,
                                                             trainingPointsToAverage.Select(b => b.Label).Median(),
                                                             identification);
                    }
                    chargeToLookAt++;
                } while (chargeToLookAt <= highestKnownChargeForThisPeptide + 1);
                theIndex += direction;
            }
        }

        private IEnumerable<LabeledMs2DataPoint> SearchMS2Spectrum(IMsDataScanWithPrecursor<IMzSpectrum<IMzPeak>> ms2DataScan, Proteomics.Peptide peptide, int peptideCharge, Psm identification)
        {
            numFragmentsIdentified = 0;
            // Key: mz value, Value: error
            var addedPeaks = new Dictionary<double, double>();

            var countForThisMS2 = 0;
            var countForThisMS2a = 0;

            var scanWindowRange = ms2DataScan.ScanWindowRange;

            IHasChemicalFormula[] fragmentList = peptide.Fragment(fragmentTypesForCalibration, true).OfType<IHasChemicalFormula>().ToArray();

            foreach (var fragment in fragmentList)
            {
                bool fragmentIdentified = false;
                bool computedIsotopologues = false;
                double[] masses = new double[0];
                double[] intensities = new double[0];
                // First look for monoisotopic masses, do not compute distribution spectrum!

                for (int chargeToLookAt = 1; chargeToLookAt <= peptideCharge; chargeToLookAt++)
                {
                    var monoisotopicMZ = fragment.MonoisotopicMass.ToMz(chargeToLookAt);
                    if (monoisotopicMZ > scanWindowRange.Maximum)
                        continue;
                    if (monoisotopicMZ < scanWindowRange.Minimum)
                        break;
                    var closestPeakMZ = ms2DataScan.MassSpectrum.GetClosestPeakXvalue(monoisotopicMZ);

                    if (mzToleranceForMs2Search.Within(closestPeakMZ, monoisotopicMZ))
                    {
                        if (!computedIsotopologues)
                        {
                            var dist = IsotopicDistribution.GetDistribution(fragment.ThisChemicalFormula, fineResolutionForIsotopeDistCalculation, 0.001);

                            masses = dist.Masses.ToArray();
                            intensities = dist.Intensities.ToArray();

                            Array.Sort(intensities, masses, Comparer<double>.Create((x, y) => y.CompareTo(x)));
                            computedIsotopologues = true;

                            break;
                        }
                    }
                }

                if (computedIsotopologues)
                {
                    bool startingToAdd = false;
                    for (int chargeToLookAt = 1; chargeToLookAt <= peptideCharge; chargeToLookAt++)
                    {
                        if (masses.First().ToMz(chargeToLookAt) > scanWindowRange.Maximum)
                        {
                            continue;
                        }
                        if (masses.Last().ToMz(chargeToLookAt) < scanWindowRange.Minimum)
                        {
                            break;
                        }
                        var trainingPointsToAverage = new List<LabeledMs2DataPoint>();
                        foreach (double a in masses)
                        {
                            double theMZ = a.ToMz(chargeToLookAt);
                            var npwr = ms2DataScan.MassSpectrum.NumPeaksWithinRange(mzToleranceForMs2Search.GetMinimumValue(theMZ), mzToleranceForMs2Search.GetMaximumValue(theMZ));
                            if (npwr == 0)
                            {
                                break;
                            }
                            numMs2MassChargeCombinationsConsidered++;
                            if (npwr > 1)
                            {
                                numMs2MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks++;
                                continue;
                            }
                            int closestPeakIndex = ms2DataScan.MassSpectrum.GetClosestPeakIndex(theMZ);
                            var closestPeakMZ = ms2DataScan.MassSpectrum.XArray[closestPeakIndex];
                            if (!addedPeaks.ContainsKey(closestPeakMZ))
                            {
                                addedPeaks.Add(closestPeakMZ, Math.Abs(closestPeakMZ - theMZ));
                                trainingPointsToAverage.Add(new LabeledMs2DataPoint(closestPeakMZ, double.NaN, ms2DataScan.MassSpectrum.YArray[closestPeakIndex], double.NaN, null, double.NaN, closestPeakMZ - theMZ, null));
                            }
                        }
                        // If started adding and suddnely stopped, go to next one, no need to look at higher charges
                        if (trainingPointsToAverage.Count == 0 && startingToAdd)
                            break;
                        if (trainingPointsToAverage.Count < Math.Min(minMS2isotopicPeaksNeededForConfirmedIdentification, intensities.Count()))
                        {
                        }
                        else
                        {
                            startingToAdd = true;
                            if (!fragmentIdentified)
                            {
                                fragmentIdentified = true;
                                numFragmentsIdentified += 1;
                            }

                            countForThisMS2 += trainingPointsToAverage.Count;
                            countForThisMS2a++;
                            yield return new LabeledMs2DataPoint(trainingPointsToAverage.Select(b => b.mz).Average(),
                                     ms2DataScan.RetentionTime,
                                     trainingPointsToAverage.Select(b => b.intensity).Average(),
                                     ms2DataScan.TotalIonCurrent,
                                     ms2DataScan.InjectionTime,
                                     ms2DataScan.IsolationMz,
                                     trainingPointsToAverage.Select(b => b.Label).Median(),
                                     identification);
                        }
                    }
                }
            }
        }

        #endregion Private Methods
    }
}