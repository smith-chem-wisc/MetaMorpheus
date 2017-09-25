using Chemistry;
using EngineLayer.Calibration;
using MassSpectrometry;
using MzLibUtil;
using Proteomics;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace EngineLayer
{
    public class DataPointAcquisitionEngine : MetaMorpheusEngine
    {
        #region Private Fields

        private const double fineResolutionForIsotopeDistCalculation = 0.1;

        private readonly List<Psm> goodIdentifications;
        private readonly IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> myMsDataFile;
        private readonly int numFragmentsNeededForEveryIdentification;
        private readonly Tolerance mzToleranceForMs1Search;
        private readonly Tolerance mzToleranceForMs2Search;
        private readonly int minMS1isotopicPeaksNeededForConfirmedIdentification;
        private readonly int minMS2isotopicPeaksNeededForConfirmedIdentification;
        private readonly FragmentTypes fragmentTypesForCalibration;

        #endregion Private Fields

        #region Public Constructors

        public DataPointAcquisitionEngine(
            List<Psm> goodIdentifications,
            IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> myMsDataFile,
            Tolerance mzToleranceForMs1Search,
            Tolerance mzToleranceForMs2Search,
            int numFragmentsNeededForEveryIdentification,
            int minMS1isotopicPeaksNeededForConfirmedIdentification,
            int minMS2isotopicPeaksNeededForConfirmedIdentification,
            FragmentTypes fragmentTypesForCalibration,
            List<string> nestedIds) : base(nestedIds)
        {
            this.goodIdentifications = goodIdentifications;
            this.myMsDataFile = myMsDataFile;
            this.mzToleranceForMs1Search = mzToleranceForMs1Search;
            this.mzToleranceForMs2Search = mzToleranceForMs2Search;
            this.numFragmentsNeededForEveryIdentification = numFragmentsNeededForEveryIdentification;
            this.minMS1isotopicPeaksNeededForConfirmedIdentification = minMS1isotopicPeaksNeededForConfirmedIdentification;
            this.minMS2isotopicPeaksNeededForConfirmedIdentification = minMS2isotopicPeaksNeededForConfirmedIdentification;
            this.fragmentTypesForCalibration = fragmentTypesForCalibration;
        }

        #endregion Public Constructors

        #region Protected Methods

        protected override MetaMorpheusEngineResults RunSpecific()
        {
            Status("Extracting data points:");
            // The final training point list

            int numMs1MassChargeCombinationsConsidered = 0;
            int numMs1MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks = 0;
            int numMs2MassChargeCombinationsConsidered = 0;
            int numMs2MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks = 0;
            List<LabeledMs1DataPoint> Ms1List = new List<LabeledMs1DataPoint>();
            List<LabeledMs2DataPoint> Ms2List = new List<LabeledMs2DataPoint>();

            //// Set of peaks, identified by m/z and retention time. If a peak is in here, it means it has been a part of an accepted identification, and should be rejected
            //var peaksAddedFromMS1HashSet = new HashSet<Tuple<double, double>>();

            int numIdentifications = goodIdentifications.Count;
            // Loop over identifications

            HashSet<string> sequences = new HashSet<string>();

            object lockObj = new object();
            object lockObj2 = new object();
            object lockObj3 = new object();
            Parallel.ForEach(Partitioner.Create(0, numIdentifications), fff =>
            {
                for (int matchIndex = fff.Item1; matchIndex < fff.Item2; matchIndex++)
                {
                    Psm identification = goodIdentifications[matchIndex];

                    //// Progress
                    //if (numIdentifications < 100 || matchIndex % (numIdentifications / 100) == 0)
                    //    ReportProgress(new ProgressEventArgs(100 * matchIndex / numIdentifications, "Looking at identifications...", nestedIds));

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

                    var ms2tuple = SearchMS2Spectrum(myMsDataFile.GetOneBasedScan(ms2scanNumber) as IMsDataScanWithPrecursor<IMzSpectrum<IMzPeak>>, coolPeptide, peptideCharge, identification);

                    // If MS2 has low evidence for peptide, skip and go to next one
                    if (ms2tuple.Item4 < numFragmentsNeededForEveryIdentification)
                        continue;

                    lock (lockObj2)
                    {
                        Ms2List.AddRange(ms2tuple.Item1);
                        numMs2MassChargeCombinationsConsidered += ms2tuple.Item2;
                        numMs2MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks += ms2tuple.Item3;
                        if (sequences.Contains(identification.FullSequence))
                            continue; // Do not search same sequence multiple times in MS1 scans
                        sequences.Add(identification.FullSequence);
                    }

                    // Calculate isotopic distribution of the full peptide
                    var dist = IsotopicDistribution.GetDistribution(coolPeptide.GetChemicalFormula(), fineResolutionForIsotopeDistCalculation, 0.001);

                    double[] theoreticalMasses = dist.Masses.ToArray();
                    double[] theoreticalIntensities = dist.Intensities.ToArray();

                    Array.Sort(theoreticalIntensities, theoreticalMasses, Comparer<double>.Create((x, y) => y.CompareTo(x)));

                    var ms1tupleBack = SearchMS1Spectra(theoreticalMasses, theoreticalIntensities, ms2scanNumber, -1, peptideCharge, identification);

                    var ms1tupleForward = SearchMS1Spectra(theoreticalMasses, theoreticalIntensities, ms2scanNumber, 1, peptideCharge, identification);

                    lock (lockObj)
                    {
                        Ms1List.AddRange(ms1tupleBack.Item1);
                        numMs1MassChargeCombinationsConsidered += ms1tupleBack.Item2;
                        numMs1MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks += ms1tupleBack.Item3;
                        Ms1List.AddRange(ms1tupleForward.Item1);
                        numMs1MassChargeCombinationsConsidered += ms1tupleForward.Item2;
                        numMs1MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks += ms1tupleForward.Item3;
                    }
                }
            });

            return new DataPointAquisitionResults(this,
                Ms1List,
                Ms2List,
                numMs1MassChargeCombinationsConsidered,
                numMs1MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks,
                numMs2MassChargeCombinationsConsidered,
                numMs2MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks
            );
        }

        #endregion Protected Methods

        #region Private Methods

        private (List<LabeledMs1DataPoint>, int, int) SearchMS1Spectra(double[] theoreticalMasses, double[] theoreticalIntensities, int ms2spectrumIndex, int direction, int peptideCharge, Psm identification)
        {
            List<LabeledMs1DataPoint> result = new List<LabeledMs1DataPoint>();
            int numMs1MassChargeCombinationsConsidered = 0;
            int numMs1MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks = 0;

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
                    if (theoreticalMasses[0].ToMz(chargeToLookAt) > scanWindowRange.Maximum)
                    {
                        chargeToLookAt++;
                        continue;
                    }
                    if (theoreticalMasses[0].ToMz(chargeToLookAt) < scanWindowRange.Minimum)
                        break;
                    var trainingPointsToAverage = new List<LabeledMs1DataPoint>();
                    foreach (double a in theoreticalMasses)
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

                        highestKnownChargeForThisPeptide = Math.Max(highestKnownChargeForThisPeptide, chargeToLookAt);
                        trainingPointsToAverage.Add(new LabeledMs1DataPoint(closestPeakMZ, double.NaN, double.NaN, double.NaN, Math.Log(fullMS1spectrum.YArray[closestPeakIndex]), theMZ, null));
                    }
                    // If started adding and suddnely stopped, go to next one, no need to look at higher charges
                    if (trainingPointsToAverage.Count == 0 && startingToAddCharges == true)
                    {
                        break;
                    }
                    if ((trainingPointsToAverage.Count == 0 || (trainingPointsToAverage.Count == 1 && theoreticalIntensities[0] < 0.65)) && (peptideCharge <= chargeToLookAt))
                    {
                        break;
                    }
                    if ((trainingPointsToAverage.Count == 1 && theoreticalIntensities[0] < 0.65)
                        || trainingPointsToAverage.Count < Math.Min(minMS1isotopicPeaksNeededForConfirmedIdentification, theoreticalIntensities.Count()))
                    {
                    }
                    else
                    {
                        addedAscan = true;
                        startingToAddCharges = true;
                        countForThisScan++;
                        result.Add(new LabeledMs1DataPoint(trainingPointsToAverage.Select(b => b.mz).Average(),
                                                             fullMS1scan.RetentionTime,
                                                             Math.Log(fullMS1scan.TotalIonCurrent),
                                                             fullMS1scan.InjectionTime.HasValue ? Math.Log(fullMS1scan.InjectionTime.Value) : double.NaN,
                                                             trainingPointsToAverage.Select(b => b.logIntensity).Average(),
                                                             trainingPointsToAverage.Select(b => b.expectedMZ).Average(),
                                                             identification));
                    }
                    chargeToLookAt++;
                } while (chargeToLookAt <= highestKnownChargeForThisPeptide + 1);
                theIndex += direction;
            }
            return (result, numMs1MassChargeCombinationsConsidered, numMs1MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks);
        }

        private (List<LabeledMs2DataPoint>, int, int, int) SearchMS2Spectrum(IMsDataScanWithPrecursor<IMzSpectrum<IMzPeak>> ms2DataScan, Proteomics.Peptide peptide, int peptideCharge, Psm identification)
        {
            List<LabeledMs2DataPoint> result = new List<LabeledMs2DataPoint>();
            int numMs2MassChargeCombinationsConsidered = 0;
            int numMs2MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks = 0;
            int numFragmentsIdentified = 0;
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

                            var closestPeakIndex = ms2DataScan.MassSpectrum.GetClosestPeakIndex(theMZ);
                            var closestPeakMZ = ms2DataScan.MassSpectrum.XArray[closestPeakIndex];

                            if (!addedPeaks.ContainsKey(closestPeakMZ))
                            {
                                addedPeaks.Add(closestPeakMZ, Math.Abs(closestPeakMZ - theMZ));
                                trainingPointsToAverage.Add(new LabeledMs2DataPoint(closestPeakMZ, double.NaN, double.NaN, double.NaN, Math.Log(ms2DataScan.MassSpectrum.YArray[closestPeakIndex]), theMZ, null));
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
                            result.Add(new LabeledMs2DataPoint(trainingPointsToAverage.Select(b => b.mz).Average(),
                                     ms2DataScan.RetentionTime,
                                     Math.Log(ms2DataScan.TotalIonCurrent),
                                     ms2DataScan.InjectionTime.HasValue ? Math.Log(ms2DataScan.InjectionTime.Value) : double.NaN,
                                     trainingPointsToAverage.Select(b => b.logIntensity).Average(),
                                     trainingPointsToAverage.Select(b => b.expectedMZ).Average(),
                                     identification));
                        }
                    }
                }
            }

            return (result, numMs2MassChargeCombinationsConsidered, numMs2MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks, numFragmentsIdentified);
        }

        #endregion Private Methods
    }
}