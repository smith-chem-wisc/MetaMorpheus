using Chemistry;
using MassSpectrometry;
using MzLibUtil;
using Proteomics;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace EngineLayer.Calibration
{
    public class DataPointAcquisitionEngine : MetaMorpheusEngine
    {
        #region Private Fields

        private const double fineResolutionForIsotopeDistCalculation = 0.1;

        private readonly List<PeptideSpectralMatch> goodIdentifications;
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
            List<PeptideSpectralMatch> goodIdentifications,
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
            List<LabeledDataPoint> Ms1List = new List<LabeledDataPoint>();
            List<LabeledDataPoint> Ms2List = new List<LabeledDataPoint>();

            int numIdentifications = goodIdentifications.Count;

            // Loop over identifications

            HashSet<string> sequences = new HashSet<string>();

            object lockObj = new object();
            object lockObj2 = new object();
            Parallel.ForEach(Partitioner.Create(0, numIdentifications), fff =>
            {
                for (int matchIndex = fff.Item1; matchIndex < fff.Item2; matchIndex++)
                {
                    PeptideSpectralMatch identification = goodIdentifications[matchIndex];

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

                    var ms2tuple = SearchMS2Spectrum(myMsDataFile.GetOneBasedScan(ms2scanNumber) as IMsDataScanWithPrecursor<IMzSpectrum<IMzPeak>>, identification);

                    lock (lockObj2)
                    {
                        Ms2List.AddRange(ms2tuple);
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

            // datapoints are ordered because they were acquired in a parallized search and we want repeatable results
            return new DataPointAquisitionResults(this,
                Ms1List.OrderBy(p => p.rt).ThenBy(p => p.experimentalMz).ToList(),
                Ms2List.OrderBy(p => p.rt).ThenBy(p => p.experimentalMz).ToList(),
                numMs1MassChargeCombinationsConsidered,
                numMs1MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks,
                numMs2MassChargeCombinationsConsidered,
                numMs2MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks
            );
        }

        #endregion Protected Methods

        #region Private Methods

        private (List<LabeledDataPoint>, int, int) SearchMS1Spectra(double[] theoreticalMasses, double[] theoreticalIntensities, int ms2spectrumIndex, int direction, int peptideCharge, PeptideSpectralMatch identification)
        {
            List<LabeledDataPoint> result = new List<LabeledDataPoint>();
            int numMs1MassChargeCombinationsConsidered = 0;
            int numMs1MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks = 0;

            int theIndex;
            theIndex = direction == 1 ? ms2spectrumIndex : ms2spectrumIndex - 1;

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
                    var trainingPointsToAverage = new List<LabeledDataPoint>();
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
                        var closestPeakMZ = fullMS1spectrum.XArray[closestPeakIndex.Value];

                        highestKnownChargeForThisPeptide = Math.Max(highestKnownChargeForThisPeptide, chargeToLookAt);
                        trainingPointsToAverage.Add(new LabeledDataPoint(closestPeakMZ, double.NaN, double.NaN, double.NaN, Math.Log(fullMS1spectrum.YArray[closestPeakIndex.Value]), theMZ, null));
                    }
                    // If started adding and suddnely stopped, go to next one, no need to look at higher charges
                    if (trainingPointsToAverage.Count == 0 && startingToAddCharges)
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
                        result.Add(new LabeledDataPoint(trainingPointsToAverage.Select(b => b.experimentalMz).Average(),
                                                             fullMS1scan.RetentionTime,
                                                             Math.Log(fullMS1scan.TotalIonCurrent),
                                                             fullMS1scan.InjectionTime.HasValue ? Math.Log(fullMS1scan.InjectionTime.Value) : double.NaN,
                                                             trainingPointsToAverage.Select(b => b.logIntensity).Average(),
                                                             trainingPointsToAverage.Select(b => b.theoreticalMz).Average(),
                                                             identification));
                    }
                    chargeToLookAt++;
                } while (chargeToLookAt <= highestKnownChargeForThisPeptide + 1);
                theIndex += direction;
            }
            return (result, numMs1MassChargeCombinationsConsidered, numMs1MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks);
        }

        private static List<LabeledDataPoint> SearchMS2Spectrum(IMsDataScanWithPrecursor<IMzSpectrum<IMzPeak>> ms2DataScan, PeptideSpectralMatch identification)
        {
            List<LabeledDataPoint> result = new List<LabeledDataPoint>();

            if (ms2DataScan.MassSpectrum.Size == 0)
                return result;
            
            foreach(var productType in identification.MatchedIonDictOnlyMatches)
            {
                for (int i = 0; i < productType.Value.Length; i++)
                {
                    double theorMz = productType.Value[i].ToMz(1);
                    int ind = ms2DataScan.MassSpectrum.GetClosestPeakIndex(theorMz).Value;

                    double exptPeakMz = ms2DataScan.MassSpectrum.XArray[ind];
                    double exptPeakIntensity = ms2DataScan.MassSpectrum.YArray[ind];
                    double injTime = ms2DataScan.InjectionTime ?? double.NaN;

                    result.Add(
                        new LabeledDataPoint(
                            exptPeakMz, 
                            ms2DataScan.RetentionTime, 
                            Math.Log(ms2DataScan.TotalIonCurrent),
                            Math.Log(injTime),
                            Math.Log(exptPeakIntensity),
                            theorMz,
                            identification));
                }
            }

            return result;
        }

        #endregion Private Methods
    }
}