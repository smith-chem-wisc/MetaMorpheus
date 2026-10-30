using Chemistry;
using MassSpectrometry;
using MzLibUtil;
using Omics;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using Transcriptomics;
using Transcriptomics.Digestion;

namespace EngineLayer.Calibration
{
    public class DataPointAcquisitionEngine : MetaMorpheusEngine
    {
        private const double FineResolutionForIsotopeDistCalculation = 0.1;

        private readonly List<SpectralMatch> GoodIdentifications;
        private readonly List<Ms2ScanWithSpecificMass> GoodScans;
        private readonly MsDataFile MyMsDataFile;
        private readonly Tolerance PrecursorMassTolerance;
        private readonly Tolerance ProductMassTolerance;
        private readonly int MinMS1isotopicPeaksNeededForConfirmedIdentification;

        public DataPointAcquisitionEngine(
            List<SpectralMatch> goodIdentifications,
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

        private readonly object _ms1Lock = new();
        private readonly object _ms2Lock = new();

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

            int maxThreadsPerFile = CommonParameters.MaxThreadsToUsePerFile;
            int[] threads = Enumerable.Range(0, maxThreadsPerFile).ToArray();
            Parallel.ForEach(threads, (matchIndex) =>
            {
                //foreach psm
                for (; matchIndex < GoodIdentifications.Count; matchIndex += maxThreadsPerFile)
                {
                    // Stop loop if canceled
                    if (GlobalVariables.StopLoops) { return; }

                    SpectralMatch identification = GoodIdentifications[matchIndex];

                    // Each identification has an MS2 spectrum attached to it.
                    int ms2scanNumber = identification.ScanNumber;
                    int peptideCharge = identification.ScanPrecursorCharge;
                    //skip if ambiguous
                    if (identification.FullSequence == null)
                        continue;

                    var representativeSinglePeptide = identification.BestMatchingBioPolymersWithSetMods.First().SpecificBioPolymer;

                    // Fragment masses need no chemical formula, so every match gives MS2 points, whatever its modifications.
                    List<LabeledDataPoint> ms2tuple = SearchMS2Spectrum(GoodScans[matchIndex], identification, ProductMassTolerance);

                    lock (_ms2Lock)
                    {
                        Ms2List.AddRange(ms2tuple);
                    }

                    // Calculate theoretical isotopic distribution of the full peptide
                    var (theoreticalMasses, theoreticalIntensities) = TheoreticalIsotopicDistribution(representativeSinglePeptide);
                    if (theoreticalMasses.Length == 0)
                        continue;

                    Array.Sort(theoreticalIntensities, theoreticalMasses, Comparer<double>.Create((x, y) => y.CompareTo(x)));

                    var ms1tupleBack = SearchMS1Spectra(theoreticalMasses, theoreticalIntensities, ms2scanNumber, -1, peptideCharge, identification);

                    var ms1tupleForward = SearchMS1Spectra(theoreticalMasses, theoreticalIntensities, ms2scanNumber, 1, peptideCharge, identification);

                    lock (_ms1Lock)
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
                numMs2MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks,
                CommonParameters
            );
        }

        /// <summary>
        /// The theoretical isotopic envelope of a match's full sequence, for finding its precursor in MS1 scans.
        /// </summary>
        /// <remarks>
        /// When every modification has a chemical formula, this is the distribution of the full formula, as it always was.
        /// A modification with a mass and no formula (the isobaric labels in Mods\TMT.txt, for one) leaves the full formula
        /// null. Its envelope is then built from the formula that is known (the residues, termini and every modification
        /// that has a formula), with averagine standing in for the rest of the monoisotopic mass and the envelope anchored
        /// at that exact mass (mzLib's AveragineFormula, as FlashLFQ does). Before, such a match gave no calibration point
        /// at all, so a search with a fixed isobaric label calibrated on its few unlabelled matches only (issue #1875).
        /// </remarks>
        public static (double[] Masses, double[] Intensities) TheoreticalIsotopicDistribution(IBioPolymerWithSetMods bioPolymer)
        {
            if (bioPolymer.ThisChemicalFormula is { } fullFormula)
            {
                var dist = IsotopicDistribution.GetDistribution(fullFormula, FineResolutionForIsotopeDistCalculation, 0.001);
                return (dist.Masses.ToArray(), dist.Intensities.ToArray());
            }

            var averagine = bioPolymer is OligoWithSetMods
                ? new OxyriboAveragine().GetAverageChemicalFormula()
                : new Averagine().GetAverageChemicalFormula();
            return AveragineFormula.GetAnchoredDistribution(KnownChemicalFormula(bioPolymer), bioPolymer.MonoisotopicMass,
                averagine, FineResolutionForIsotopeDistCalculation, 0.001);
        }

        /// <summary>
        /// The formula of the parts of <paramref name="bioPolymer"/> that have one: its residues and termini, and every
        /// modification with a chemical formula. A type this does not know contributes nothing, so averagine fills all of it.
        /// </summary>
        public static ChemicalFormula KnownChemicalFormula(IBioPolymerWithSetMods bioPolymer)
        {
            ChemicalFormula formula = bioPolymer switch
            {
                PeptideWithSetModifications => new Proteomics.AminoAcidPolymer.Peptide(bioPolymer.BaseSequence).GetChemicalFormula(),
                OligoWithSetMods oligo => new RNA(oligo.BaseSequence, fivePrimeTerm: oligo.FivePrimeTerminus, threePrimeTerm: oligo.ThreePrimeTerminus).GetChemicalFormula(),
                _ => new ChemicalFormula(),
            };
            foreach (var mod in bioPolymer.AllModsOneIsNterminus.Values)
            {
                if (mod.ChemicalFormula != null)
                    formula.Add(mod.ChemicalFormula);
            }
            return formula;
        }

        private (List<LabeledDataPoint>, int, int) SearchMS1Spectra(double[] theoreticalMasses, double[] theoreticalIntensities, int ms2spectrumIndex, int direction, int peptideCharge, SpectralMatch identification)
        {
            List<LabeledDataPoint> result = new List<LabeledDataPoint>();
            int numMs1MassChargeCombinationsConsidered = 0;
            int numMs1MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks = 0;

            //looking for the index of the ms1 scan following (1) or preceeding (-1) the ms2 scan
            int theIndex;
            theIndex = direction == 1 ? ms2spectrumIndex + 1 : identification.PrecursorScanNumber ?? ms2spectrumIndex;

            bool addedAscan = true;
            bool positiveMode = identification.ScanPrecursorCharge > 0;

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
                int chargeToLookAt = positiveMode ? 1 : -1; 
                int chargeLimit = positiveMode ? highestKnownChargeForThisPeptide + 1 : highestKnownChargeForThisPeptide - 1;

                do
                {
                    double mz = theoreticalMasses[0].ToMz(chargeToLookAt);

                    // If m/z is above the scan window, try next charge
                    if (mz > scanWindowRange.Maximum)
                    {
                        chargeToLookAt += positiveMode ? 1 : -1;
                        continue;
                    }

                    // If m/z is below the scan window, break early, more charge will only make the mz smaller. 
                    if (mz < scanWindowRange.Minimum)
                    {
                        break;
                    }

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

                        highestKnownChargeForThisPeptide = positiveMode
                            ? Math.Max(highestKnownChargeForThisPeptide, chargeToLookAt)
                            : Math.Min(highestKnownChargeForThisPeptide, chargeToLookAt);
                        trainingPointsToAverage.Add(new LabeledDataPoint(closestPeakMZ, -1, double.NaN, double.NaN, Math.Log(fullMS1spectrum.YArray[closestPeakIndex]), theMZ, null));
                    }
                    // If started adding and suddenly stopped, go to next one, no need to look at higher charges
                    if (trainingPointsToAverage.Count == 0 && startingToAddCharges)
                    {
                        break;
                    }
                    if ((trainingPointsToAverage.Count == 0 || (trainingPointsToAverage.Count == 1 && theoreticalIntensities[0] < 0.65)) && (positiveMode ? peptideCharge <= chargeToLookAt : peptideCharge >= chargeToLookAt))
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
                    chargeToLookAt += positiveMode ? 1 : -1;
                } while (positiveMode ? chargeToLookAt <= chargeLimit : chargeToLookAt >= chargeLimit);
                theIndex += direction;
            }
            return (result, numMs1MassChargeCombinationsConsidered, numMs1MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks);
        }

        private static List<LabeledDataPoint> SearchMS2Spectrum(Ms2ScanWithSpecificMass ms2DataScan, SpectralMatch identification, Tolerance ms2Tolerance)
        {
            List<LabeledDataPoint> result = new List<LabeledDataPoint>();

            if (ms2DataScan.TheScan.MassSpectrum.Size == 0)
                return result;

            // LowCID deconvolutes no fragments (Ms2ScanWithSpecificMass), and CalibrationEngine does not
            // calibrate MS2 scans in that mode, so there are no MS2 datapoints to collect.
            if (ms2DataScan.ExperimentalFragments == null)
                return result;


            //deconvolute
            var isotopicEnvelopes = ms2DataScan.ExperimentalFragments;
            //foreach theoretical fragment that was found to match, use those points
            foreach (var matchedIon in identification.MatchedFragmentIons)
            {
                double theoreticalMass = matchedIon.NeutralTheoreticalProduct.NeutralMass;
                //get envelopes that match
                var envelopesThatMatch = isotopicEnvelopes.Where(x => ms2Tolerance.Within(x.MonoisotopicMass, theoreticalMass)).OrderBy(x => Math.Abs(x.MonoisotopicMass - theoreticalMass)).ToList();
                
                if(envelopesThatMatch.Count == 0)
                    continue;
                //only allow one envelope per charge state
                bool[] chargeStateFound = new bool[envelopesThatMatch.Max(x => Math.Abs(x.Charge)) + 1];

                foreach (var envelopeThatMatched in envelopesThatMatch)
                {
                    int charge = Math.Abs(envelopeThatMatched.Charge);
                    //if we haven't seen this charge state already
                    if (!chargeStateFound[charge])
                    {
                        chargeStateFound[charge] = true;

                        double exptPeakMz = envelopeThatMatched.MonoisotopicMass.ToMz(charge);
                        double exptPeakIntensity = envelopeThatMatched.TotalIntensity;
                        double injTime = ms2DataScan.TheScan.InjectionTime ?? double.NaN;

                        result.Add(
                            new LabeledDataPoint(
                                exptPeakMz,
                                ms2DataScan.OneBasedScanNumber,
                                Math.Log(ms2DataScan.TotalIonCurrent),
                                Math.Log(injTime),
                        Math.Log(exptPeakIntensity),
                                matchedIon.NeutralTheoreticalProduct.NeutralMass.ToMz(charge),
                                identification));
                    }
                }
            }
            return result;
        }
    }
}