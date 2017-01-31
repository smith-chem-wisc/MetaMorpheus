using Chemistry;
using MassSpectrometry;
using MathNet.Numerics.Statistics;
using Proteomics;
using Spectra;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace EngineLayer.Calibration
{
    public class CalibrationEngine : MyEngine
    {

        #region Private Fields

        private const double fineResolutionForIsotopeDistCalculation = 0.1;
        private const double fracForTraining = 0.75;

        private readonly int minMS1isotopicPeaksNeededForConfirmedIdentification;
        private readonly int minMS2isotopicPeaksNeededForConfirmedIdentification;
        private readonly int numFragmentsNeededForEveryIdentification;
        private readonly Tolerance mzToleranceForMs1Search;
        private readonly Tolerance mzToleranceForMs2Search;
        private readonly int randomSeed;
        private readonly FragmentTypes fragmentTypesForCalibration;
        private List<NewPsmWithFdr> identifications;
        private IMsDataFile<IMzSpectrum<MzPeak>> myMsDataFile;

        #endregion Private Fields

        #region Public Constructors

        public CalibrationEngine(IMsDataFile<IMzSpectrum<MzPeak>> myMSDataFile, int randomSeed, Tolerance mzToleranceForMs2Search, List<NewPsmWithFdr> identifications, int minMS1IsotopicPeaksNeededForConfirmedIdentification, int minMS2IsotopicPeaksNeededForConfirmedIdentification, int numFragmentsNeededForEveryIdentification, Tolerance mzToleranceForMS1Search, FragmentTypes fragmentTypesForCalibration) : base(2)
        {
            this.myMsDataFile = myMSDataFile;
            this.randomSeed = randomSeed;
            this.identifications = identifications;
            this.minMS1isotopicPeaksNeededForConfirmedIdentification = minMS1IsotopicPeaksNeededForConfirmedIdentification;
            this.minMS2isotopicPeaksNeededForConfirmedIdentification = minMS2IsotopicPeaksNeededForConfirmedIdentification;
            this.numFragmentsNeededForEveryIdentification = numFragmentsNeededForEveryIdentification;
            this.mzToleranceForMs1Search = mzToleranceForMS1Search;
            this.mzToleranceForMs2Search = mzToleranceForMs2Search;
            this.fragmentTypesForCalibration = fragmentTypesForCalibration;
        }

        #endregion Public Constructors

        #region Protected Methods

        protected override MyResults RunSpecific()
        {
            Status("Calibrating " + Path.GetFileName(myMsDataFile.FilePath));

            var trainingPointCounts = new List<int>();
            var goodResult = new CalibrationResults(myMsDataFile, this);
            List<LabeledDataPoint> pointList;
            for (int calibrationRound = 1; ; calibrationRound++)
            {
                Status("Calibration round " + calibrationRound);
                int numMs1MassChargeCombinationsConsidered = 0;
                int numMs1MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks = 0;
                int numMs2MassChargeCombinationsConsidered = 0;
                int numMs2MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks = 0;
                Status("Getting Training Points");
                pointList = GetDataPoints(ref numMs1MassChargeCombinationsConsidered, ref numMs1MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks, ref numMs2MassChargeCombinationsConsidered, ref numMs2MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks);

                goodResult.Add(numMs1MassChargeCombinationsConsidered, numMs1MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks, pointList.Count, numMs2MassChargeCombinationsConsidered, numMs2MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks);

                if (calibrationRound >= 2 && pointList.Count <= trainingPointCounts[calibrationRound - 2])
                    break;

                trainingPointCounts.Add(pointList.Count);

                var pointList2 = pointList.Where((b) => b.inputs[0] > 0).ToList();
                if (pointList2.Count == 0)
                {
                    return new MyErroredResults(this, "No MS2 training points, identification quality is poor. Try to change the Fragment tolerance." + goodResult.ToString());
                }

                var pointList1 = pointList.Where((b) => b.inputs[0] < 0).ToList();
                if (pointList1.Count == 0)
                {
                    return new MyErroredResults(this, "No MS1 training points, identification quality is poor. Try to change the Parent tolerance." + goodResult.ToString());
                }

                CalibrationFunction combinedCalibration = Calibrate(pointList);

                goodResult.Add(combinedCalibration);

                if (combinedCalibration == null)
                    return new MyErroredResults(this, "Could not calibrate");
            }

            return goodResult;
        }

        #endregion Protected Methods

        #region Private Methods

        private List<LabeledDataPoint> GetDataPoints(ref int numMs1MassChargeCombinationsConsidered, ref int numMs1MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks, ref int numMs2MassChargeCombinationsConsidered, ref int numMs2MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks)
        {
            Status("Extracting data points:");
            // The final training point list
            var trainingPointsToReturn = new List<LabeledDataPoint>();

            // Set of peaks, identified by m/z and retention time. If a peak is in here, it means it has been a part of an accepted identification, and should be rejected
            var peaksAddedFromMS1HashSet = new HashSet<Tuple<double, double>>();

            int numIdentifications = identifications.Count;
            // Loop over identifications
            for (int matchIndex = 0; matchIndex < numIdentifications; matchIndex++)
            {
                var identification = identifications[matchIndex];
                if (identification.qValue > 0.01)
                    break;

                // Progress
                if (numIdentifications < 100 || matchIndex % (numIdentifications / 100) == 0)
                    ReportProgress(new ProgressEventArgs(100 * matchIndex / numIdentifications, "Looking at identifications..."));

                // Skip decoys, they are for sure not there!
                if (identification.IsDecoy)
                    continue;

                // Each identification has an MS2 spectrum attached to it.
                int ms2spectrumIndex = identification.thisPSM.newPsm.scanNumber;

                // Get the peptide, don't forget to add the modifications!!!!
                var SequenceWithChemicalFormulas = identification.thisPSM.SequenceWithChemicalFormulas;
                int peptideCharge = identification.thisPSM.newPsm.scanPrecursorCharge;

                Proteomics.Peptide coolPeptide = new Proteomics.Peptide(SequenceWithChemicalFormulas);

                int numFragmentsIdentified;
                var candidateTrainingPointsForPeptide = SearchMS2Spectrum(myMsDataFile.GetOneBasedScan(ms2spectrumIndex), coolPeptide, peptideCharge, out numFragmentsIdentified, ref numMs2MassChargeCombinationsConsidered, ref numMs2MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks);

                // If MS2 has low evidence for peptide, skip and go to next one
                if (numFragmentsIdentified < numFragmentsNeededForEveryIdentification)
                    continue;

                // Calculate isotopic distribution of the full peptide

                var dist = IsotopicDistribution.GetDistribution(coolPeptide.GetChemicalFormula(), fineResolutionForIsotopeDistCalculation, 0.001);

                double[] masses = new double[dist.masses.Length];
                double[] intensities = new double[dist.intensities.Length];
                for (int i = 0; i < dist.masses.Length; i++)
                {
                    masses[i] = dist.masses[i];
                    intensities[i] = dist.intensities[i];
                }
                Array.Sort(intensities, masses, Comparer<double>.Create((x, y) => y.CompareTo(x)));

                SearchMS1Spectra(masses, intensities, candidateTrainingPointsForPeptide, ms2spectrumIndex, -1, peaksAddedFromMS1HashSet, peptideCharge, ref numMs1MassChargeCombinationsConsidered, ref numMs1MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks);
                SearchMS1Spectra(masses, intensities, candidateTrainingPointsForPeptide, ms2spectrumIndex, 1, peaksAddedFromMS1HashSet, peptideCharge, ref numMs1MassChargeCombinationsConsidered, ref numMs1MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks);

                trainingPointsToReturn.AddRange(candidateTrainingPointsForPeptide);
            }
            Status("Number of training points: " + trainingPointsToReturn.Count);
            return trainingPointsToReturn;
        }

        private CalibrationFunction Calibrate(List<LabeledDataPoint> trainingPoints)
        {
            var rnd = new Random(randomSeed);
            var shuffledTrainingPoints = trainingPoints.OrderBy(item => rnd.Next()).ToArray();

            var trainList = shuffledTrainingPoints.Take((int)(trainingPoints.Count * fracForTraining)).ToList();
            var testList = shuffledTrainingPoints.Skip((int)(trainingPoints.Count * fracForTraining)).ToList();

            var trainList1 = trainList.Where((b) => b.inputs[0] < 0).ToList();
            var trainList2 = trainList.Where((b) => b.inputs[0] > 0).ToList();
            var testList1 = testList.Where((b) => b.inputs[0] < 0).ToList();
            var testList2 = testList.Where((b) => b.inputs[0] > 0).ToList();

            CalibrationFunction bestMS1predictor = new IdentityCalibrationFunction();
            CalibrationFunction bestMS2predictor = new IdentityCalibrationFunction();
            double bestMS1MSE = bestMS1predictor.getMSE(testList1);
            double bestMS2MSE = bestMS2predictor.getMSE(testList2);

            {
                var ms1regressor = new ConstantCalibrationFunction();
                var ms2regressor = new ConstantCalibrationFunction();
                ms1regressor.Train(trainList1);
                ms2regressor.Train(trainList2);
                double MS1mse = ms1regressor.getMSE(testList1);
                double MS2mse = ms2regressor.getMSE(testList2);
                if (MS1mse < bestMS1MSE)
                {
                    bestMS1MSE = MS1mse;
                    bestMS1predictor = ms1regressor;
                }
                if (MS2mse < bestMS2MSE)
                {
                    bestMS2MSE = MS2mse;
                    bestMS2predictor = ms2regressor;
                }
            }
            var transforms = new List<TransformFunction>();

            transforms.Add(new TransformFunction(b => new double[] { b[1] }, 1));
            transforms.Add(new TransformFunction(b => new double[] { b[2] }, 1));
            transforms.Add(new TransformFunction(b => new double[] { Math.Log(b[3]) }, 1));
            transforms.Add(new TransformFunction(b => new double[] { Math.Log(b[4]) }, 1));
            transforms.Add(new TransformFunction(b => new double[] { Math.Log(b[5]) }, 1));

            transforms.Add(new TransformFunction(b => new double[] { b[1], b[2] }, 2));
            transforms.Add(new TransformFunction(b => new double[] { b[1], Math.Log(b[3]) }, 2));
            transforms.Add(new TransformFunction(b => new double[] { b[1], Math.Log(b[4]) }, 2));
            transforms.Add(new TransformFunction(b => new double[] { b[1], Math.Log(b[5]) }, 2));
            transforms.Add(new TransformFunction(b => new double[] { b[2], Math.Log(b[3]) }, 2));
            transforms.Add(new TransformFunction(b => new double[] { b[2], Math.Log(b[4]) }, 2));
            transforms.Add(new TransformFunction(b => new double[] { b[2], Math.Log(b[5]) }, 2));
            transforms.Add(new TransformFunction(b => new double[] { Math.Log(b[3]), Math.Log(b[4]) }, 2));
            transforms.Add(new TransformFunction(b => new double[] { Math.Log(b[3]), Math.Log(b[5]) }, 2));
            transforms.Add(new TransformFunction(b => new double[] { Math.Log(b[4]), Math.Log(b[5]) }, 2));

            transforms.Add(new TransformFunction(b => new double[] { b[1], b[2], Math.Log(b[3]) }, 3));
            transforms.Add(new TransformFunction(b => new double[] { b[1], b[2], Math.Log(b[4]) }, 3));
            transforms.Add(new TransformFunction(b => new double[] { b[1], b[2], Math.Log(b[5]) }, 3));
            transforms.Add(new TransformFunction(b => new double[] { b[1], Math.Log(b[3]), Math.Log(b[4]) }, 3));
            transforms.Add(new TransformFunction(b => new double[] { b[1], Math.Log(b[3]), Math.Log(b[5]) }, 3));
            transforms.Add(new TransformFunction(b => new double[] { b[1], Math.Log(b[4]), Math.Log(b[5]) }, 3));
            transforms.Add(new TransformFunction(b => new double[] { b[2], Math.Log(b[3]), Math.Log(b[4]) }, 3));
            transforms.Add(new TransformFunction(b => new double[] { b[2], Math.Log(b[3]), Math.Log(b[5]) }, 3));
            transforms.Add(new TransformFunction(b => new double[] { b[2], Math.Log(b[4]), Math.Log(b[5]) }, 3));
            transforms.Add(new TransformFunction(b => new double[] { Math.Log(b[3]), Math.Log(b[4]), Math.Log(b[5]) }, 3));

            transforms.Add(new TransformFunction(b => new double[] { b[1], b[2], Math.Log(b[3]), Math.Log(b[4]) }, 4));
            transforms.Add(new TransformFunction(b => new double[] { b[1], b[2], Math.Log(b[3]), Math.Log(b[5]) }, 4));
            transforms.Add(new TransformFunction(b => new double[] { b[1], b[2], Math.Log(b[4]), Math.Log(b[5]) }, 4));
            transforms.Add(new TransformFunction(b => new double[] { b[1], Math.Log(b[3]), Math.Log(b[4]), Math.Log(b[5]) }, 4));
            transforms.Add(new TransformFunction(b => new double[] { b[2], Math.Log(b[3]), Math.Log(b[4]), Math.Log(b[5]) }, 4));

            transforms.Add(new TransformFunction(b => new double[] { b[1], b[2], Math.Log(b[3]), Math.Log(b[4]), Math.Log(b[5]) }, 5));

            try
            {
                foreach (var transform in transforms)
                {
                    var ms1regressorLinear = new LinearCalibrationFunctionMathNet(transform);
                    var ms2regressorLinear = new LinearCalibrationFunctionMathNet(transform);
                    ms1regressorLinear.Train(trainList1);
                    ms2regressorLinear.Train(trainList2);
                    var MS1mse = ms1regressorLinear.getMSE(testList1);
                    var MS2mse = ms2regressorLinear.getMSE(testList2);
                    if (MS1mse < bestMS1MSE)
                    {
                        bestMS1MSE = MS1mse;
                        bestMS1predictor = ms1regressorLinear;
                    }
                    if (MS2mse < bestMS2MSE)
                    {
                        bestMS2MSE = MS2mse;
                        bestMS2predictor = ms2regressorLinear;
                    }
                }
                foreach (var transform in transforms)
                {
                    var ms1regressorQuadratic = new QuadraticCalibrationFunctionMathNet(transform);
                    var ms2regressorQuadratic = new QuadraticCalibrationFunctionMathNet(transform);
                    ms1regressorQuadratic.Train(trainList1);
                    ms2regressorQuadratic.Train(trainList2);
                    var MS1mse = ms1regressorQuadratic.getMSE(testList1);
                    var MS2mse = ms2regressorQuadratic.getMSE(testList2);
                    if (MS1mse < bestMS1MSE)
                    {
                        bestMS1MSE = MS1mse;
                        bestMS1predictor = ms1regressorQuadratic;
                    }
                    if (MS2mse < bestMS2MSE)
                    {
                        bestMS2MSE = MS2mse;
                        bestMS2predictor = ms2regressorQuadratic;
                    }
                }
            }
            catch (ArgumentException)
            {
            }

            CalibrationFunction bestCf = new SeparateCalibrationFunction(bestMS1predictor, bestMS2predictor);

            Status("Calibrating Spectra");

            CalibrateSpectra(bestCf);

            return bestCf;
        }

        private void CalibrateSpectra(CalibrationFunction bestCf)
        {
            foreach (var a in myMsDataFile)
            {
                if (a.MsnOrder == 2)
                {
                    int oneBasedScanNumber;
                    a.TryGetPrecursorOneBasedScanNumber(out oneBasedScanNumber);
                    var precursorScan = myMsDataFile.GetOneBasedScan(oneBasedScanNumber);

                    double precursorMZ;
                    a.TryGetSelectedIonGuessMZ(out precursorMZ);
                    double precursorIntensity;
                    a.TryGetSelectedIonGuessIntensity(out precursorIntensity);
                    double newSelectedMZ = precursorMZ - bestCf.Predict(new double[] { -1, precursorMZ, precursorScan.RetentionTime, precursorIntensity, precursorScan.TotalIonCurrent, precursorScan.InjectionTime });

                    double monoisotopicMZ;
                    a.TryGetSelectedIonGuessMonoisotopicMZ(out monoisotopicMZ);
                    double monoisotopicIntensity;
                    a.TryGetSelectedIonGuessMonoisotopicIntensity(out monoisotopicIntensity);

                    if (double.IsNaN(monoisotopicIntensity))
                        monoisotopicIntensity = precursorScan.MassSpectrum.GetClosestPeak(monoisotopicMZ).Intensity;

                    double newMonoisotopicMZ = monoisotopicMZ - bestCf.Predict(new double[] { -1, monoisotopicMZ, precursorScan.RetentionTime, monoisotopicIntensity, precursorScan.TotalIonCurrent, precursorScan.InjectionTime });

                    double IsolationMZ;
                    a.TryGetIsolationMZ(out IsolationMZ);

                    Func<MzPeak, double> theFunc = x => x.Mz - bestCf.Predict(new double[] { 1, x.Mz, a.RetentionTime, x.Intensity, a.TotalIonCurrent, a.InjectionTime, IsolationMZ, (x.Mz - a.ScanWindowRange.Minimum) / (a.ScanWindowRange.Maximum - a.ScanWindowRange.Minimum) });
                    a.TranformByApplyingFunctionsToSpectraAndReplacingPrecursorMZs(theFunc, newSelectedMZ, newMonoisotopicMZ);
                }
                else
                {
                    Func<MzPeak, double> theFUnc = x => x.Mz - bestCf.Predict(new double[] { -1, x.Mz, a.RetentionTime, x.Intensity, a.TotalIonCurrent, a.InjectionTime });
                    a.TranformByApplyingFunctionsToSpectraAndReplacingPrecursorMZs(theFUnc, double.NaN, double.NaN);
                }
            }
        }

        private int SearchMS1Spectra(double[] originalMasses, double[] originalIntensities, List<LabeledDataPoint> myCandidatePoints, int ms2spectrumIndex, int direction, HashSet<Tuple<double, double>> peaksAddedHashSet, int peptideCharge, ref int numMs1MassChargeCombinationsConsidered, ref int numMs1MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks)
        {
            int goodIndex = -1;
            var scores = new List<int>();
            var theIndex = -1;
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
                var myCandidatePointsForThisMS1scan = new List<LabeledDataPoint>();
                var fullMS1scan = myMsDataFile.GetOneBasedScan(theIndex);
                double ms1RetentionTime = fullMS1scan.RetentionTime;
                var scanWindowRange = fullMS1scan.ScanWindowRange;
                var fullMS1spectrum = fullMS1scan.MassSpectrum;
                if (fullMS1spectrum.Count == 0)
                    break;

                bool startingToAddCharges = false;
                int chargeToLookAt = 1;
                do
                {
                    if (originalMasses[0].ToMassToChargeRatio(chargeToLookAt) > scanWindowRange.Maximum)
                    {
                        chargeToLookAt++;
                        continue;
                    }
                    if (originalMasses[0].ToMassToChargeRatio(chargeToLookAt) < scanWindowRange.Minimum)
                        break;
                    var trainingPointsToAverage = new List<LabeledDataPoint>();
                    foreach (double a in originalMasses)
                    {
                        double theMZ = a.ToMassToChargeRatio(chargeToLookAt);

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

                        var closestPeak = fullMS1spectrum.GetClosestPeak(theMZ);
                        var closestPeakMZ = closestPeak.Mz;

                        var theTuple = Tuple.Create(closestPeakMZ, ms1RetentionTime);
                        if (!peaksAddedHashSet.Contains(theTuple))
                        {
                            peaksAddedHashSet.Add(theTuple);
                            highestKnownChargeForThisPeptide = Math.Max(highestKnownChargeForThisPeptide, chargeToLookAt);
                            trainingPointsToAverage.Add(new LabeledDataPoint(new double[] { closestPeakMZ, closestPeak.Intensity }, closestPeakMZ - theMZ));
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
                    if (trainingPointsToAverage.Count == 1 && originalIntensities[0] < 0.65)
                    {
                    }
                    else if (trainingPointsToAverage.Count < Math.Min(minMS1isotopicPeaksNeededForConfirmedIdentification, originalIntensities.Count()))
                    {
                    }
                    else
                    {
                        addedAscan = true;
                        startingToAddCharges = true;
                        countForThisScan++;
                        double[] inputs = { -1, trainingPointsToAverage.Select(b => b.inputs[0]).Average(), fullMS1scan.RetentionTime, trainingPointsToAverage.Select(b => b.inputs[1]).Average(), fullMS1scan.TotalIonCurrent, fullMS1scan.InjectionTime };
                        var a = new LabeledDataPoint(inputs, trainingPointsToAverage.Select(b => b.output).Median());

                        myCandidatePointsForThisMS1scan.Add(a);
                    }
                    chargeToLookAt++;
                } while (chargeToLookAt <= highestKnownChargeForThisPeptide + 1);

                if (myCandidatePointsForThisMS1scan.Count > 0)
                    goodIndex = theIndex;
                //    SoftwareLockMassRunner.WriteDataToFiles(myCandidatePointsForThisMS1scan, theIndex.ToString());
                myCandidatePoints.AddRange(myCandidatePointsForThisMS1scan);

                scores.Add(countForThisScan);
                theIndex += direction;
            }
            return goodIndex;
        }

        private List<LabeledDataPoint> SearchMS2Spectrum(IMsDataScan<IMzSpectrum<MzPeak>> ms2DataScan, Proteomics.Peptide peptide, int peptideCharge, out int candidateFragmentsIdentified, ref int numMs2MassChargeCombinationsConsidered, ref int numMs2MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks)
        {
            var myCandidatePoints = new List<LabeledDataPoint>();

            // Key: mz value, Value: error
            var addedPeaks = new Dictionary<double, double>();

            double IsolationMZ;
            ms2DataScan.TryGetIsolationMZ(out IsolationMZ);

            var countForThisMS2 = 0;
            var countForThisMS2a = 0;
            var numFragmentsIdentified = 0;

            var scanWindowRange = ms2DataScan.ScanWindowRange;

            Fragment[] fragmentList = peptide.Fragment(fragmentTypesForCalibration, true).ToArray();

            foreach (IHasChemicalFormula fragment in fragmentList)
            {
                bool fragmentIdentified = false;
                bool computedIsotopologues = false;
                double[] masses = new double[0];
                double[] intensities = new double[0];
                // First look for monoisotopic masses, do not compute distribution spectrum!

                for (int chargeToLookAt = 1; chargeToLookAt <= peptideCharge; chargeToLookAt++)
                {
                    var monoisotopicMZ = fragment.MonoisotopicMass.ToMassToChargeRatio(chargeToLookAt);
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

                            masses = new double[dist.masses.Length];
                            intensities = new double[dist.intensities.Length];
                            for (int i = 0; i < dist.masses.Length; i++)
                            {
                                masses[i] = dist.masses[i];
                                intensities[i] = dist.intensities[i];
                            }
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
                        if (masses.First().ToMassToChargeRatio(chargeToLookAt) > scanWindowRange.Maximum)
                        {
                            continue;
                        }
                        if (masses.Last().ToMassToChargeRatio(chargeToLookAt) < scanWindowRange.Minimum)
                        {
                            break;
                        }
                        var trainingPointsToAverage = new List<LabeledDataPoint>();
                        foreach (double a in masses)
                        {
                            double theMZ = a.ToMassToChargeRatio(chargeToLookAt);
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
                            var closestPeak = ms2DataScan.MassSpectrum.GetClosestPeak(theMZ);
                            var closestPeakMZ = closestPeak.Mz;
                            if (!addedPeaks.ContainsKey(closestPeakMZ))
                            {
                                addedPeaks.Add(closestPeakMZ, Math.Abs(closestPeakMZ - theMZ));
                                trainingPointsToAverage.Add(new LabeledDataPoint(new double[] { closestPeakMZ, closestPeak.Intensity }, closestPeakMZ - theMZ));
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

                            double addedMZ = trainingPointsToAverage.Select(b => b.inputs[0]).Average();
                            double relativeMZ = (addedMZ - ms2DataScan.ScanWindowRange.Minimum) / (ms2DataScan.ScanWindowRange.Maximum - ms2DataScan.ScanWindowRange.Minimum);
                            double[] inputs = { 1, addedMZ, ms2DataScan.RetentionTime, trainingPointsToAverage.Select(b => b.inputs[1]).Average(), ms2DataScan.TotalIonCurrent, ms2DataScan.InjectionTime, IsolationMZ, relativeMZ };
                            var a = new LabeledDataPoint(inputs, trainingPointsToAverage.Select(b => b.output).Median());

                            myCandidatePoints.Add(a);
                        }
                    }
                }
            }

            candidateFragmentsIdentified = numFragmentsIdentified;
            return myCandidatePoints;
        }

        #endregion Private Methods

    }
}