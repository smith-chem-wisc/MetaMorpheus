using Chemistry;
using MassSpectrometry;
using MzLibUtil;
using Omics.Modifications;
using System;
using System.Collections.Generic;
using System.Linq;

namespace EngineLayer
{
    public class Ms2ScanWithSpecificMass
    {
        public Ms2ScanWithSpecificMass(MsDataScan mzLibScan, double precursorMonoisotopicPeakMz, int precursorCharge, string fullFilePath, CommonParameters commonParam,
            IsotopicEnvelope[] neutralExperimentalFragments = null, double? precursorIntensity = null, int? envelopePeakCount = null, double? precursorFractionalIntensity = null,
            double? precursorMostAbundantMass = null,
            double precursorDeconvolutionScore = 0)
        {
            PrecursorMonoisotopicPeakMz = precursorMonoisotopicPeakMz;
            PrecursorCharge = precursorCharge;
            PrecursorMass = PrecursorMonoisotopicPeakMz.ToMass(precursorCharge);
            PrecursorMostAbundantMass = precursorMostAbundantMass ?? 0;
            PrecursorIntensity = precursorIntensity ?? 1;
            PrecursorEnvelopePeakCount = envelopePeakCount ?? 1;
            PrecursorFractionalIntensity = precursorFractionalIntensity ?? -1;
            PrecursorDeconvolutionScore = precursorDeconvolutionScore;
            FullFilePath = fullFilePath;
            ChildScans = new List<Ms2ScanWithSpecificMass>();
            NativeId = mzLibScan.NativeId;

            TheScan = mzLibScan;

            ScanMetadata = new ScanMetadata(
                OneBasedScanNumber: mzLibScan.OneBasedScanNumber,
                OneBasedPrecursorScanNumber: mzLibScan.OneBasedPrecursorScanNumber,
                RetentionTime: mzLibScan.RetentionTime,
                NumPeaks: mzLibScan.MassSpectrum.Size,
                TotalIonCurrent: mzLibScan.TotalIonCurrent,
                NativeId: mzLibScan.NativeId,
                FullFilePath: fullFilePath,
                PrecursorCharge: precursorCharge,
                PrecursorMonoisotopicPeakMz: precursorMonoisotopicPeakMz,
                PrecursorMass: PrecursorMass,
                PrecursorIntensity: PrecursorIntensity,
                PrecursorEnvelopePeakCount: PrecursorEnvelopePeakCount,
                PrecursorFractionalIntensity: PrecursorFractionalIntensity,
                OneOverK0: mzLibScan is TimsDataScan tims ? tims.OneOverK0 : null);

            if (commonParam.DissociationType != DissociationType.LowCID)
            {
                ExperimentalFragments = neutralExperimentalFragments ?? GetNeutralExperimentalFragments(mzLibScan, commonParam);
            }
            if (ExperimentalFragments != null && ExperimentalFragments.Any())
            {
                DeconvolutedMonoisotopicMasses = ExperimentalFragments.Select(p => p.MonoisotopicMass).ToArray();
            }
            else
            {
                DeconvolutedMonoisotopicMasses = new double[0];
            }
            MassesAreAscendingNumbers = IsAscendingWithoutNaN(DeconvolutedMonoisotopicMasses);
        }

        public MsDataScan TheScan { get; }

        /// <summary>
        /// Lightweight, immutable snapshot of scan and precursor metadata.
        /// Designed to be passed to SpectralMatch so the heavyweight scan objects
        /// can be released from memory after scoring.
        /// </summary>
        public ScanMetadata ScanMetadata { get; private set; }

        /// <summary>
        /// Re-reads the peak count after the spectrum has been altered in place. XCorr pre-processing
        /// runs after these objects are built, so without this the recorded count is the one from
        /// before pre-processing rather than the one the search actually scored against.
        /// </summary>
        public void RefreshPeakCount() =>
            ScanMetadata = ScanMetadata with { NumPeaks = TheScan.MassSpectrum.Size };

        public double PrecursorMonoisotopicPeakMz { get; }

        /// <summary>
        /// The observed monoisotopic precursor mass, in every search. This property does not change
        /// meaning with the search type.
        /// </summary>
        public double PrecursorMass { get; }

        /// <summary>
        /// The observed neutral mass of the most abundant (tallest) isotopologue of the precursor envelope,
        /// or 0 when no envelope was deconvoluted for this precursor (a scan-header precursor, or a neutral
        /// mass read from a pre-deconvoluted file). Like <see cref="PrecursorMass"/> this is a plain
        /// observation, populated in every search and never redefined by one: it is what the detector saw,
        /// not what the current search chose to match on.
        /// <para>
        /// Which of the two a search matches candidates against is decided by the
        /// <see cref="EngineLayer.MassDiffAcceptor"/> — see
        /// <see cref="PrecursorMassExtensions.GetPrecursorMassForSearch(Ms2ScanWithSpecificMass, MassDiffAcceptor)"/>.
        /// </para>
        /// </summary>
        public double PrecursorMostAbundantMass { get; }
        public int PrecursorCharge { get; }
        public double PrecursorIntensity { get; }
        public int PrecursorEnvelopePeakCount { get; }
        public double PrecursorFractionalIntensity { get; }
        /// <summary>
        /// Method-agnostic envelope-quality score in [0, 1] from mzLib's DeconvolutionScorer.
        /// 0 indicates either a maximally low-quality envelope or that no envelope was
        /// deconvoluted for this scan (e.g. scan-header-only precursor path).
        /// </summary>
        public double PrecursorDeconvolutionScore { get; }
        public string FullFilePath { get; }
        public IsotopicEnvelope[] ExperimentalFragments { get; private set; }
        public List<Ms2ScanWithSpecificMass> ChildScans { get; set; } // MS2/MS3 scans that are children of this MS2 scan
        private double[] DeconvolutedMonoisotopicMasses;

        /// <summary>
        /// True when <see cref="DeconvolutedMonoisotopicMasses"/> is non-empty, holds no NaN and never decreases, which is
        /// what <see cref="GetNeutralExperimentalFragments"/> produces. Only then does a plain lower-bound search agree with
        /// Array.BinarySearch (which sorts NaN first); any other array keeps using Array.BinarySearch.
        /// </summary>
        private readonly bool MassesAreAscendingNumbers;
        public string NativeId { get; }
        public int OneBasedScanNumber => TheScan.OneBasedScanNumber;
        public int? OneBasedPrecursorScanNumber => TheScan.OneBasedPrecursorScanNumber;
        public double RetentionTime => TheScan.RetentionTime;
        public int NumPeaks => TheScan.MassSpectrum.Size;
        public double TotalIonCurrent => TheScan.TotalIonCurrent;
        /// <summary>
        /// An array containing the intensities of the reporter ions for isobaric mass tags. 
        /// If multiplex quantification wasn't performed, this will be null
        /// </summary>
        public double[]? IsobaricMassTagReporterIonIntensities { get; private set; }

        public static IsotopicEnvelope[] GetNeutralExperimentalFragments(MsDataScan scan, CommonParameters commonParam)
        {
            var neutralExperimentalFragmentMasses =
                Deconvoluter.Deconvolute(scan, commonParam.ProductDeconvolutionParameters, scan.MassSpectrum.Range).ToList();

            if (!commonParam.AssumeOrphanPeaksAreZ1Fragments || scan.MassSpectrum is NeutralMassSpectrum)
                return neutralExperimentalFragmentMasses.OrderBy(p => p.MonoisotopicMass).ToArray();

            HashSet<double> alreadyClaimedMzs = new HashSet<double>(neutralExperimentalFragmentMasses
                .SelectMany(p => p.Peaks.Select(v => v.mz.RoundedDouble()!.Value)));

            int charge = scan.Polarity == Polarity.Positive ? 1 : -1;
            for (int i = 0; i < scan.MassSpectrum.XArray.Length; i++)
            {
                double mz = scan.MassSpectrum.XArray[i];
                double intensity = scan.MassSpectrum.YArray[i];

                if (!alreadyClaimedMzs.Contains(mz.RoundedDouble()!.Value))
                {
                    neutralExperimentalFragmentMasses.Add(new IsotopicEnvelope(
                        new List<(double mz, double intensity)> { (mz, intensity) },
                        mz.ToMass(charge), charge, intensity, 0));
                }
            }

            return neutralExperimentalFragmentMasses.OrderBy(p => p.MonoisotopicMass).ToArray();
        }

        /// <summary>
        /// Writes the reporter ion intensities into the IsobaricMassTagReporterIonIntensities property
        /// If the scan has
        /// </summary>
        /// <param name="massTag"></param>
        public void SetIsobaricMassTagReporterIonIntensities(IsobaricMassTag massTag)
        {
            IsobaricMassTagReporterIonIntensities = massTag.GetReporterIonIntensities(ReporterIonScan.MassSpectrum);
        }

        /// <summary>
        /// Writes the reporter ion intensities into the IsobaricMassTagReporterIonIntensities property, reading them from
        /// the scans as the spectra file stores them. The spectra this scan and its child scans carry may have been trimmed
        /// at load or binned for a LowCID search since, and either can remove a reporter ion. The scan read is the same one
        /// the overload above reads: the most intense child scan when there is one, otherwise this scan.
        /// </summary>
        /// <param name="storedScanFor">Each loaded scan, by reference, paired with the same scan as the file stores it.
        /// It is keyed by the scan object, not the scan number, because an MGF can repeat a SCANS= value. A scan with no
        /// pair is read from its loaded spectrum, as the overload above reads it.</param>
        public void SetIsobaricMassTagReporterIonIntensities(IsobaricMassTag massTag, IReadOnlyDictionary<MsDataScan, MsDataScan> storedScanFor)
        {
            MsDataScan loaded = ReporterIonScan;
            MsDataScan source = storedScanFor.TryGetValue(loaded, out var stored) ? stored : loaded;
            IsobaricMassTagReporterIonIntensities = massTag.GetReporterIonIntensities(source.MassSpectrum);
        }

        /// <summary>
        /// The scan reporter ions are read from, by both overloads above: the most intense child scan when there is one,
        /// otherwise this scan.
        /// </summary>
        public MsDataScan ReporterIonScan => UseChildScansForIsobaricQuant(out var mostIntenseChildScan)
            ? mostIntenseChildScan.TheScan
            : TheScan;

        /// <summary>
        /// Helper method to determine if child scans contain isobaric mass tag (e.g., TMT) reporter ions
        /// In most cases, if MS3 scans exist, there will only be one MS3 scan per MS2 scan, and the ChildScans list will contain that one MS3 scan
        /// If there are methods that generate multiple child scans, then this method will return the most intense child scan
        /// However, if we ever have multiple child scans with isobaric mass tag reporter ions, we should revisit this logic
        /// </summary>
        /// <param name="mostIntenseChildScan">The child scan that has the highest ion current </param>
        private bool UseChildScansForIsobaricQuant(out Ms2ScanWithSpecificMass mostIntenseChildScan)
        {
            mostIntenseChildScan = null;
            if (ChildScans.IsNullOrEmpty()) return false;
            mostIntenseChildScan = ChildScans.MaxBy(s => s.TotalIonCurrent);
            if(mostIntenseChildScan != null)
            {
                return true;
            }
            return false;
        }
       
        public IsotopicEnvelope GetClosestExperimentalIsotopicEnvelope(double theoreticalNeutralMass)
        {
            if (DeconvolutedMonoisotopicMasses.Length == 0)
            {
                return null;
            }
            return ExperimentalFragments[GetClosestFragmentMass(theoreticalNeutralMass)];
        }

        /// <summary>
        /// The index into <see cref="ExperimentalFragments"/> of the envelope <see cref="GetClosestExperimentalIsotopicEnvelope"/>
        /// would return, or -1 when the scan has no deconvoluted fragments (where that method returns null).
        /// </summary>
        public int GetClosestExperimentalIsotopicEnvelopeIndex(double theoreticalNeutralMass)
        {
            if (DeconvolutedMonoisotopicMasses.Length == 0)
            {
                return -1;
            }
            return GetClosestFragmentMass(theoreticalNeutralMass);
        }

        /// <summary>
        /// The monoisotopic mass of <see cref="ExperimentalFragments"/>[<paramref name="index"/>], read from the scan's own
        /// sorted mass array. Lets a caller reject a fragment match on mass error without dereferencing the envelope -
        /// in a search nearly every theoretical product is rejected, and touching the envelope to do it costs a cache miss.
        /// </summary>
        public double GetDeconvolutedMonoisotopicMass(int index) => DeconvolutedMonoisotopicMasses[index];

        public int GetClosestFragmentMass(double mass)
        {
            double[] masses = DeconvolutedMonoisotopicMasses;

            // A search calls this for every theoretical product of every candidate, and Array.BinarySearch<double> is a
            // call through a generic comparer each time. For a mass that is not in the array, BinarySearch returns the
            // complement of the first index holding a larger mass, which is exactly the lower bound found here with plain
            // comparisons. An exact hit is left to BinarySearch, because with repeated masses (real spectra have them) it
            // may land on any of the equal entries; so is a NaN mass, and any array that is not ascending numbers.
            if (MassesAreAscendingNumbers && !double.IsNaN(mass))
            {
                int lo = 0;
                int hi = masses.Length;
                while (lo < hi)
                {
                    int mid = (int)((uint)(lo + hi) >> 1);
                    if (masses[mid] < mass)
                    {
                        lo = mid + 1;
                    }
                    else
                    {
                        hi = mid;
                    }
                }
                if (lo == masses.Length || masses[lo] != mass)
                {
                    return ClosestOfNeighbours(masses, mass, lo);
                }
            }

            int index = Array.BinarySearch(masses, mass);
            if (index >= 0)
            {
                return index;
            }
            return ClosestOfNeighbours(masses, mass, ~index);
        }

        /// <summary>
        /// For each of the first <paramref name="count"/> entries of <paramref name="ascendingMasses"/>, the index
        /// <see cref="GetClosestFragmentMass"/> returns for it (-1 when the scan has no deconvoluted fragments), written to the
        /// same position of <paramref name="closestIndices"/>. When the masses ascend, as a peptide's sorted product masses do,
        /// the lower bounds are found in one forward walk of the scan's mass array: a sequential read the processor can fetch
        /// ahead, where a binary search per product waits on one cache miss after another. The lower bound of a mass does not
        /// depend on how it was found, and exact hits, NaN, masses out of order and arrays that are not ascending numbers all
        /// go to <see cref="GetClosestFragmentMass"/> itself, so every index is the one it would return.
        /// </summary>
        public void GetClosestFragmentMasses(double[] ascendingMasses, int count, int[] closestIndices)
        {
            double[] masses = DeconvolutedMonoisotopicMasses;
            if (!MassesAreAscendingNumbers)
            {
                for (int k = 0; k < count; k++)
                {
                    closestIndices[k] = GetClosestFragmentMass(ascendingMasses[k]); // -1 for an empty array
                }
                return;
            }

            int lo = 0;
            double previous = double.NegativeInfinity;
            for (int k = 0; k < count; k++)
            {
                double mass = ascendingMasses[k];
                if (double.IsNaN(mass))
                {
                    closestIndices[k] = GetClosestFragmentMass(mass);
                    continue;
                }
                if (mass < previous)
                {
                    lo = 0; // out of order: walk again from the start rather than assume
                }
                previous = mass;

                while (lo < masses.Length && masses[lo] < mass)
                {
                    lo++;
                }
                closestIndices[k] = lo < masses.Length && masses[lo] == mass
                    ? GetClosestFragmentMass(mass)
                    : ClosestOfNeighbours(masses, mass, lo);
            }
        }

        /// <summary>
        /// Given the index of the first mass larger than <paramref name="mass"/> (masses.Length if none), the index of the
        /// closer of it and the mass before it; a tie goes to the one before.
        /// </summary>
        private static int ClosestOfNeighbours(double[] masses, double mass, int index)
        {
            if (index == masses.Length)
            {
                return index - 1;
            }
            if (index == 0 || mass - masses[index - 1] > masses[index] - mass)
            {
                return index;
            }

            return index - 1;
        }

        private static bool IsAscendingWithoutNaN(double[] masses)
        {
            if (masses.Length == 0 || double.IsNaN(masses[0]))
            {
                return false;
            }
            for (int i = 1; i < masses.Length; i++)
            {
                // written so that a NaN fails it too
                if (!(masses[i] >= masses[i - 1]))
                {
                    return false;
                }
            }
            return true;
        }

        //look for IsotopicEnvelopes which are in the range of acceptable mass 
        public IsotopicEnvelope[] GetClosestExperimentalIsotopicEnvelopeList(double minimumMass, double maxMass)
        {

            if (DeconvolutedMonoisotopicMasses.Length == 0)
            {
                return null;
            }

            //if no mass is in this range, then return null
            if (DeconvolutedMonoisotopicMasses[0] > maxMass || DeconvolutedMonoisotopicMasses.Last() < minimumMass)
            {
                return null;
            }

            int startIndex = GetClosestFragmentMass(minimumMass);
            int endIndex = GetClosestFragmentMass(maxMass);

            //the index we get from GetClosestFragmentMass is the closest mass, while the acceptable mass we need is between minimumMass and maxMass
            //so the startIndex mass is supposed to be larger than minimumMass, if not , then startIndex increases by 1;
            //the endIndex mass is supposed to be smaller than maxMass, if not , then endIndex decreases by 1;
            if (DeconvolutedMonoisotopicMasses[startIndex]<minimumMass)
            {
                startIndex = startIndex+1;
            }
            if(DeconvolutedMonoisotopicMasses[endIndex] > maxMass)
            {
                endIndex = endIndex - 1;
            }
            int length = endIndex - startIndex + 1;

            if (length < 1)
            {
                return null;
            }
            IsotopicEnvelope[] isotopicEnvelopes = ExperimentalFragments.Skip(startIndex).Take(length).ToArray();
            return isotopicEnvelopes;
        }
    }
}