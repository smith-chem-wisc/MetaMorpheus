using Chemistry;
using MassSpectrometry;
using MathNet.Numerics.Statistics;
using MzLibUtil;
using Proteomics;
using SharpLearning.Common.Interfaces;
using SharpLearning.Containers.Matrices;
using SharpLearning.CrossValidation.TrainingTestSplitters;
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
        private readonly Action<List<LabeledMs1DataPoint>> ms1ListAction;
        private readonly Action<List<LabeledMs2DataPoint>> ms2ListAction;
        private readonly List<Psm> goodIdentifications;
        private readonly IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> myMsDataFile;
        private int numMs1MassChargeCombinationsConsidered;

        private int numMs1MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks;

        private int numMs2MassChargeCombinationsConsidered;

        private int numMs2MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks;

        private int numFragmentsIdentified;
        private readonly List<ILearner<double>> learners;

        #endregion Private Fields

        #region Public Constructors

        public CalibrationEngine(IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> myMSDataFile, CommonParameters commonParameters, List<Psm> goodIdentifications, Action<List<LabeledMs1DataPoint>> ms1ListAction, Action<List<LabeledMs2DataPoint>> ms2ListAction, CalibrationParameters calibrationParameters, List<ILearner<double>> learners, List<string> nestedIds) : base(nestedIds)
        {
            this.myMsDataFile = myMSDataFile;
            this.goodIdentifications = goodIdentifications;
            this.minMS1isotopicPeaksNeededForConfirmedIdentification = calibrationParameters.minMS1isotopicPeaksNeededForConfirmedIdentification;
            this.minMS2isotopicPeaksNeededForConfirmedIdentification = calibrationParameters.minMS2isotopicPeaksNeededForConfirmedIdentification;
            this.numFragmentsNeededForEveryIdentification = calibrationParameters.numFragmentsNeededForEveryIdentification;
            this.mzToleranceForMs1Search = calibrationParameters.PrecursorMassTolerance;
            this.mzToleranceForMs2Search = commonParameters.ProductMassTolerance;

            fragmentTypesForCalibration = FragmentTypes.None;
            if (commonParameters.BIons)
                fragmentTypesForCalibration = fragmentTypesForCalibration | FragmentTypes.b;
            if (commonParameters.YIons)
                fragmentTypesForCalibration = fragmentTypesForCalibration | FragmentTypes.y;
            if (commonParameters.CIons)
                fragmentTypesForCalibration = fragmentTypesForCalibration | FragmentTypes.c;
            if (commonParameters.ZdotIons)
                fragmentTypesForCalibration = fragmentTypesForCalibration | FragmentTypes.zdot;

            this.ms1ListAction = ms1ListAction;
            this.ms2ListAction = ms2ListAction;
            this.learners = learners;
        }

        #endregion Public Constructors

        #region Protected Methods

        protected override MetaMorpheusEngineResults RunSpecific()
        {
            DataPointAquisitionResults prevResult;
            DataPointAquisitionResults currentResult = new DataPointAquisitionResults
            {
                Ms1List = new List<LabeledMs1DataPoint>(),
                Ms2List = new List<LabeledMs2DataPoint>(),
            };

            int round = 0;
            do
            {
                Console.WriteLine("Round: " + round++);
                prevResult = currentResult;
                var prevms1Info = currentResult.Ms1List.Select(b => b.label).MeanStandardDeviation();
                var prevms2Info = currentResult.Ms2List.Select(b => b.label).MeanStandardDeviation();
                currentResult = GetDataPoints();

                Console.WriteLine("currentResult.Ms1List.Count : " + currentResult.Ms1List.Count + " : currentResult.Ms2List.Count : " + currentResult.Ms2List.Count + " :  currentResult.Count : " + currentResult.Count);

                var ms1Info = currentResult.Ms1List.Select(b => b.label).MeanStandardDeviation();
                var ms2Info = currentResult.Ms2List.Select(b => b.label).MeanStandardDeviation();
                Console.WriteLine("MS1 : " + ms1Info.Item1 + " : " + ms1Info.Item2 + " :  MS2 : " + ms2Info.Item1 + " : " + ms2Info.Item2);
                if (prevms1Info.Item2 < ms1Info.Item2 || prevms2Info.Item2 < ms2Info.Item2 || (prevms1Info.Item2 == ms1Info.Item2 && prevms2Info.Item2 == ms2Info.Item2))
                {
                    Console.WriteLine("done with round!");
                    break;
                }

                ms1ListAction(currentResult.Ms1List);
                ms2ListAction(currentResult.Ms2List);
                Fit(currentResult, learners);
            } while (true);
            return new MetaMorpheusEngineResults(this);
        }

        #endregion Protected Methods

        #region Private Methods

        private DataPointAquisitionResults GetDataPoints()
        {
            Status("Extracting data points:");
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

                //if (matchIndex == 0)
                //    Console.WriteLine("matchIndex=0 basesequence=" + coolPeptide.BaseSequence);
                res.Ms2List.AddRange(SearchMS2Spectrum(myMsDataFile.GetOneBasedScan(ms2scanNumber) as IMsDataScanWithPrecursor<IMzSpectrum<IMzPeak>>, coolPeptide, peptideCharge, identification));
                res.numMs2MassChargeCombinationsConsidered += numMs2MassChargeCombinationsConsidered;
                res.numMs2MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks += numMs2MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks;

                //if (matchIndex == 0)
                //{
                //    Console.WriteLine("numFragmentsIdentified= " + numFragmentsIdentified);
                //    Console.WriteLine("ms2list= " + string.Join(";", res.Ms2List.Select(b => "(" + b.mz + "," + b.Label + ") ")));
                //}

                // If MS2 has low evidence for peptide, skip and go to next one
                if (numFragmentsIdentified < numFragmentsNeededForEveryIdentification)
                    continue;

                //if (matchIndex == 0)
                //    Console.WriteLine("good evidence!");

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

            return res;
        }

        private void Fit(DataPointAquisitionResults res, List<ILearner<double>> learners)
        {
            var splitter = new RandomTrainingTestIndexSplitter<double>(fracForTraining);

            var ms1splitResult = splitter.SplitSet(new F64Matrix(res.Ms1List.SelectMany(b => new[] { b.mz, b.rt, b.logTotalIonCurrent, b.logInjectionTime, b.logIntensity }).ToArray(), res.Ms1List.Count, 5), res.Ms1List.Select(b => b.label).ToArray());

            var ms2splitResult = splitter.SplitSet(new F64Matrix(res.Ms2List.SelectMany(b => new[] { b.mz, b.rt, b.logTotalIonCurrent, b.logInjectionTime, b.logIntensity }).ToArray(), res.Ms2List.Count, 5), res.Ms2List.Select(b => b.label).ToArray());

            MS1predictor ms1predictor = new MS1predictor(DoStuff(ms1splitResult, learners));

            MS2predictor ms2predictor = new MS2predictor(DoStuff(ms2splitResult, learners));

            Status("Calibrating Spectra");

            CalibrateSpectra(ms1predictor, ms2predictor);
        }

        private IPredictorModel<double> DoStuff(TrainingTestSetSplit splitResult, List<ILearner<double>> learners)
        {
            var evaluator = new MeanAbsolutErrorRegressionMetric();

            var predictions = new double[splitResult.TestSet.Targets.Length];
            IPredictorModel<double> bestModel = new IdentityCalibrationFunctionPredictorModel();
            var bestError = evaluator.Error(splitResult.TestSet.Targets, predictions);

            foreach (var learner in learners)
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
            return bestModel;
        }

        private void CalibrateSpectra(MS1predictor ms1predictor, MS2predictor ms2predictor)
        {
            foreach (var a in myMsDataFile)
            {
                if (a is IMsDataScanWithPrecursor<IMzSpectrum<IMzPeak>> theScan)
                {
                    var precursorScan = myMsDataFile.GetOneBasedScan(theScan.OneBasedPrecursorScanNumber.Value);

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

                        var closestPeakIndex = fullMS1spectrum.GetClosestPeakIndex(theMZ);
                        var closestPeakMZ = fullMS1spectrum.XArray[closestPeakIndex];

                        var theTuple = Tuple.Create(closestPeakMZ, ms1RetentionTime);
                        if (!peaksAddedHashSet.Contains(theTuple))
                        {
                            peaksAddedHashSet.Add(theTuple);
                            highestKnownChargeForThisPeptide = Math.Max(highestKnownChargeForThisPeptide, chargeToLookAt);
                            trainingPointsToAverage.Add(new LabeledMs1DataPoint(closestPeakMZ, double.NaN, double.NaN, double.NaN, Math.Log(fullMS1spectrum.YArray[closestPeakIndex]), closestPeakMZ - theMZ, null));
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
                                                             Math.Log(fullMS1scan.TotalIonCurrent),
                                                             fullMS1scan.InjectionTime.HasValue ? Math.Log(fullMS1scan.InjectionTime.Value) : double.NaN,
                                                             trainingPointsToAverage.Select(b => b.logIntensity).Average(),
                                                             trainingPointsToAverage.Select(b => b.label).Median(),
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

            //if (ms2DataScan.OneBasedScanNumber == 16490)
            //    Console.WriteLine("Looking at ms2 scan 16490");
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

                            var closestPeakIndex = ms2DataScan.MassSpectrum.GetClosestPeakIndex(theMZ);
                            var closestPeakMZ = ms2DataScan.MassSpectrum.XArray[closestPeakIndex];

                            if (!addedPeaks.ContainsKey(closestPeakMZ))
                            {
                                addedPeaks.Add(closestPeakMZ, Math.Abs(closestPeakMZ - theMZ));
                                trainingPointsToAverage.Add(new LabeledMs2DataPoint(closestPeakMZ, double.NaN, double.NaN, double.NaN, Math.Log(ms2DataScan.MassSpectrum.YArray[closestPeakIndex]), closestPeakMZ - theMZ, null));
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
                                     Math.Log(ms2DataScan.TotalIonCurrent),
                                     ms2DataScan.InjectionTime.HasValue ? Math.Log(ms2DataScan.InjectionTime.Value) : double.NaN,
                                     trainingPointsToAverage.Select(b => b.logIntensity).Average(),
                                     trainingPointsToAverage.Select(b => b.label).Median(),
                                     identification);
                        }
                    }
                }
            }
        }
    }

    #endregion Private Methods
}