using Chemistry;
using MassSpectrometry;
using MzLibUtil;
using Proteomics.AminoAcidPolymer;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace EngineLayer.Calibration
{
    public class DataPointAcquisitionEngine : MetaMorpheusEngine
    {
        private const double FineResolutionForIsotopeDistCalculation = 0.1;

        private readonly List<PeptideSpectralMatch> GoodIdentifications;
        private readonly MsDataFile MyMsDataFile;
        private readonly Tolerance MzToleranceForMs1Search;
        private readonly int MinMS1isotopicPeaksNeededForConfirmedIdentification;

        public DataPointAcquisitionEngine(
            List<PeptideSpectralMatch> goodIdentifications,
            MsDataFile myMsDataFile,
            Tolerance mzToleranceForMs1Search,
            int minMS1isotopicPeaksNeededForConfirmedIdentification,
            CommonParameters commonParameters,
            List<string> nestedIds) : base(commonParameters, nestedIds)
        {
            GoodIdentifications = goodIdentifications;
            MyMsDataFile = myMsDataFile;
            MzToleranceForMs1Search = mzToleranceForMs1Search;
            MinMS1isotopicPeaksNeededForConfirmedIdentification = minMS1isotopicPeaksNeededForConfirmedIdentification;
        }

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

            object lockObj = new object();
            object lockObj2 = new object();

            int maxThreadsPerFile = commonParameters.MaxThreadsToUsePerFile;
            int[] threads = Enumerable.Range(0, maxThreadsPerFile).ToArray();
            Parallel.ForEach(threads, (matchIndex) =>
            {
                for (; matchIndex < GoodIdentifications.Count; matchIndex += maxThreadsPerFile)
                {
                    // Stop loop if canceled
                    if (GlobalVariables.StopLoops) { return; }

                    PeptideSpectralMatch identification = GoodIdentifications[matchIndex];

                    // Each identification has an MS2 spectrum attached to it.
                    int ms2scanNumber = identification.ScanNumber;
                    int peptideCharge = identification.ScanPrecursorCharge;
                    if (identification.FullSequence == null || identification.BestMatchingPeptides.Any(p => p.Peptide.AllModsOneIsNterminus.Any(m => m.Value.ChemicalFormula == null)))
                        continue;

                    var representativeSinglePeptide = identification.BestMatchingPeptides.First().Peptide;

                    // Get the peptide, don't forget to add the modifications!!!!
                    var SequenceWithChemicalFormulas = representativeSinglePeptide.SequenceWithChemicalFormulas;
                    if (SequenceWithChemicalFormulas == null || representativeSinglePeptide.AllModsOneIsNterminus.Any(b => b.Value.NeutralLosses != null))
                        continue;

                    Peptide coolPeptide = new Peptide(SequenceWithChemicalFormulas);

                    var ms2tuple = SearchMS2Spectrum(MyMsDataFile.GetOneBasedScan(ms2scanNumber), identification);

                    lock (lockObj2)
                    {
                        Ms2List.AddRange(ms2tuple);
                    }

                    // Calculate isotopic distribution of the full peptide
                    var dist = IsotopicDistribution.GetDistribution(coolPeptide.GetChemicalFormula(), FineResolutionForIsotopeDistCalculation, 0.001);

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
                GoodIdentifications,
                Ms1List.OrderBy(p => p.ScanNumber).ThenBy(p => p.ExperimentalMz).ToList(),
                Ms2List.OrderBy(p => p.ScanNumber).ThenBy(p => p.ExperimentalMz).ToList(),
                numMs1MassChargeCombinationsConsidered,
                numMs1MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks,
                numMs2MassChargeCombinationsConsidered,
                numMs2MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks
            );
        }

        private (List<LabeledDataPoint>, int, int) SearchMS1Spectra(double[] theoreticalMasses, double[] theoreticalIntensities, int ms2spectrumIndex, int direction, int peptideCharge, PeptideSpectralMatch identification)
        {
            List<LabeledDataPoint> result = new List<LabeledDataPoint>();
            int numMs1MassChargeCombinationsConsidered = 0;
            int numMs1MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks = 0;

            int theIndex;
            theIndex = direction == 1 ? ms2spectrumIndex : ms2spectrumIndex - 1;

            bool addedAscan = true;

            int highestKnownChargeForThisPeptide = peptideCharge;
            while (theIndex >= 1 && theIndex <= MyMsDataFile.NumSpectra && addedAscan)
            {
                int countForThisScan = 0;
                if (MyMsDataFile.GetOneBasedScan(theIndex).MsnOrder > 1)
                {
                    theIndex += direction;
                    continue;
                }
                addedAscan = false;
                var fullMS1scan = MyMsDataFile.GetOneBasedScan(theIndex);
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

                        var npwr = fullMS1spectrum.NumPeaksWithinRange(MzToleranceForMs1Search.GetMinimumValue(theMZ), MzToleranceForMs1Search.GetMaximumValue(theMZ));
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
                        trainingPointsToAverage.Add(new LabeledDataPoint(closestPeakMZ, -1, double.NaN, double.NaN, Math.Log(fullMS1spectrum.YArray[closestPeakIndex.Value]), theMZ, null));
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
                        || trainingPointsToAverage.Count < Math.Min(MinMS1isotopicPeaksNeededForConfirmedIdentification, theoreticalIntensities.Count()))
                    {
                    }
                    else
                    {
                        addedAscan = true;
                        startingToAddCharges = true;
                        countForThisScan++;
                        result.Add(new LabeledDataPoint(trainingPointsToAverage.Select(b => b.ExperimentalMz).Average(),
                                                             fullMS1scan.OneBasedScanNumber,
                                                             Math.Log(fullMS1scan.TotalIonCurrent),
                                                             fullMS1scan.InjectionTime.HasValue ? Math.Log(fullMS1scan.InjectionTime.Value) : double.NaN,
                                                             trainingPointsToAverage.Select(b => b.LogIntensity).Average(),
                                                             trainingPointsToAverage.Select(b => b.TheoreticalMz).Average(),
                                                             identification));
                    }
                    chargeToLookAt++;
                } while (chargeToLookAt <= highestKnownChargeForThisPeptide + 1);
                theIndex += direction;
            }
            return (result, numMs1MassChargeCombinationsConsidered, numMs1MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks);
        }

        private static List<LabeledDataPoint> SearchMS2Spectrum(MsDataScan ms2DataScan, PeptideSpectralMatch identification)
        {
            List<LabeledDataPoint> result = new List<LabeledDataPoint>();

            if (ms2DataScan.MassSpectrum.Size == 0)
                return result;

            foreach (var matchedIon in identification.MatchedFragmentIons)
            {
                double exptPeakMz = matchedIon.Mz;
                double exptPeakIntensity = matchedIon.Intensity;
                double injTime = ms2DataScan.InjectionTime ?? double.NaN;

                result.Add(
                    new LabeledDataPoint(
                        exptPeakMz,
                        ms2DataScan.OneBasedScanNumber,
                        Math.Log(ms2DataScan.TotalIonCurrent),
                        Math.Log(injTime),
                        Math.Log(exptPeakIntensity),
                        matchedIon.NeutralTheoreticalProduct.NeutralMass.ToMz(matchedIon.Charge),
                        identification));
            }
            return result;
        }
    }
}