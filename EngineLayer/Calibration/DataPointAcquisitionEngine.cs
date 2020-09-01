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
        private readonly List<Ms2ScanWithSpecificMass> GoodScans;
        private readonly MsDataFile MyMsDataFile;
        private readonly Tolerance PrecursorMassTolerance;
        private readonly Tolerance ProductMassTolerance;
        private readonly int MinMS1isotopicPeaksNeededForConfirmedIdentification;

        public DataPointAcquisitionEngine(
            List<PeptideSpectralMatch> goodIdentifications,
            List<Ms2ScanWithSpecificMass> goodScans,
            MsDataFile myMsDataFile,
            Tolerance precursorMassTolerance,
            Tolerance productMassTolerance,
            int minMS1isotopicPeaksNeededForConfirmedIdentification,
            CommonParameters commonParameters,
            List<(string FileName, CommonParameters Parameters)> fileSpecificParameters,
            List<string> nestedIds) : base(commonParameters, fileSpecificParameters, nestedIds)
        {
            GoodIdentifications = goodIdentifications;
            GoodScans = goodScans;
            MyMsDataFile = myMsDataFile;
            PrecursorMassTolerance = precursorMassTolerance;
            ProductMassTolerance = productMassTolerance;
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

            int maxThreadsPerFile = CommonParameters.MaxThreadsToUsePerFile;
            int[] threads = Enumerable.Range(0, maxThreadsPerFile).ToArray();
            Parallel.ForEach(threads, (matchIndex) =>
            {
                //foreach psm
                for (; matchIndex < GoodIdentifications.Count; matchIndex += maxThreadsPerFile)
                {
                    // Stop loop if canceled
                    if (GlobalVariables.StopLoops) { return; }

                    PeptideSpectralMatch identification = GoodIdentifications[matchIndex];

                    // Each identification has an MS2 spectrum attached to it.
                    int ms2scanNumber = identification.ScanNumber;
                    int peptideCharge = identification.ScanPrecursorCharge;
                    //skip if ambiguous
                    if (identification.FullSequence == null || identification.BestMatchingPeptides.Any(p => p.Peptide.AllModsOneIsNterminus.Any(m => m.Value.ChemicalFormula == null)))
                        continue;

                    var representativeSinglePeptide = identification.BestMatchingPeptides.First().Peptide;

                    // Get the peptide, don't forget to add the modifications!!!!
                    var SequenceWithChemicalFormulas = representativeSinglePeptide.SequenceWithChemicalFormulas;
                    if (SequenceWithChemicalFormulas == null || representativeSinglePeptide.AllModsOneIsNterminus.Any(b => b.Value.NeutralLosses != null))
                        continue;

                    Peptide coolPeptide = new Peptide(SequenceWithChemicalFormulas);

                    List<LabeledDataPoint> ms2tuple = SearchMS2Spectrum(GoodScans[matchIndex], identification, ProductMassTolerance);

                    lock (lockObj2)
                    {
                        Ms2List.AddRange(ms2tuple);
                    }

                    // Calculate theoretical isotopic distribution of the full peptide
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

            //looking for the index of the ms1 scan following (1) or preceeding (-1) the ms2 scan
            int theIndex;
            theIndex = direction == 1 ? ms2spectrumIndex + 1 : identification.PrecursorScanNumber ?? ms2spectrumIndex;

            bool addedAscan = true;

            int highestKnownChargeForThisPeptide = peptideCharge;
            while (theIndex >= 1 && theIndex <= MyMsDataFile.NumSpectra && addedAscan) //as long as we're finding the peptide in ms1 scans
            {
                //find an ms1 scan
                if (MyMsDataFile.GetOneBasedScan(theIndex).MsnOrder != 1)
                {
                    theIndex += direction;
                    continue;
                }

                int countForThisScan = 0;
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

                        var npwr = fullMS1spectrum.NumPeaksWithinRange(PrecursorMassTolerance.GetMinimumValue(theMZ), PrecursorMassTolerance.GetMaximumValue(theMZ));
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
                        trainingPointsToAverage.Add(new LabeledDataPoint(closestPeakMZ, -1, double.NaN, double.NaN, Math.Log(fullMS1spectrum.YArray[closestPeakIndex]), theMZ, null));
                    }
                    // If started adding and suddenly stopped, go to next one, no need to look at higher charges
                    if (trainingPointsToAverage.Count == 0 && startingToAddCharges)
                    {
                        break;
                    }
                    if ((trainingPointsToAverage.Count == 0 || (trainingPointsToAverage.Count == 1 && theoreticalIntensities[0] < 0.65)) && (peptideCharge <= chargeToLookAt))
                    {
                        break;
                    }
                    if (!((trainingPointsToAverage.Count == 1 && theoreticalIntensities[0] < 0.65)
                        || trainingPointsToAverage.Count < Math.Min(MinMS1isotopicPeaksNeededForConfirmedIdentification, theoreticalIntensities.Count())))
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

        private static List<LabeledDataPoint> SearchMS2Spectrum(Ms2ScanWithSpecificMass ms2DataScan, PeptideSpectralMatch identification, Tolerance ms2Tolerance)
        {
            List<LabeledDataPoint> result = new List<LabeledDataPoint>();

            if (ms2DataScan.TheScan.MassSpectrum.Size == 0)
                return result;


            //deconvolute
            var isotopicEnvelopes = ms2DataScan.ExperimentalFragments;
            //foreach theoretical fragment that was found to match, use those points
            foreach (var matchedIon in identification.MatchedFragmentIons)
            {
                double theoreticalMass = matchedIon.NeutralTheoreticalProduct.NeutralMass;
                //get envelopes that match
                var envelopesThatMatch = isotopicEnvelopes.Where(x => ms2Tolerance.Within(x.MonoisotopicMass, theoreticalMass)).OrderBy(x => Math.Abs(x.MonoisotopicMass - theoreticalMass)).ToList();
                //only allow one envelope per charge state
                bool[] chargeStateFound = new bool[envelopesThatMatch.Max(x => x.Charge) + 1];

                foreach (var envelopeThatMatched in envelopesThatMatch)
                {
                    //if we haven't seen this charge state already
                    if (!chargeStateFound[envelopeThatMatched.Charge])
                    {
                        chargeStateFound[envelopeThatMatched.Charge] = true;

                        double exptPeakMz = envelopeThatMatched.MonoisotopicMass.ToMz(envelopeThatMatched.Charge);
                        double exptPeakIntensity = envelopeThatMatched.TotalIntensity;
                        double injTime = ms2DataScan.TheScan.InjectionTime ?? double.NaN;

                        result.Add(
                            new LabeledDataPoint(
                                exptPeakMz,
                                ms2DataScan.OneBasedScanNumber,
                                Math.Log(ms2DataScan.TotalIonCurrent),
                                Math.Log(injTime),
                                Math.Log(exptPeakIntensity),
                                matchedIon.NeutralTheoreticalProduct.NeutralMass.ToMz(envelopeThatMatched.Charge),
                                identification));
                    }
                }
            }
            return result;
        }
    }
}