using MassSpectrometry;
using MassSpectrometry.MzSpectra;
using MzLibUtil;
using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using Chemistry;
using FlashLFQ;
using IsotopicEnvelope = MassSpectrometry.IsotopicEnvelope;
using ThermoFisher.CommonCore.Data.Business;

namespace EngineLayer.ClassicSearch
{
    public class MiniClassicSearchEngine
    {
        private readonly SpectralLibrary SpectralLibrary;
        private readonly MassDiffAcceptor SearchMode;
        private readonly Ms2ScanWithSpecificMass[] MS2ByRetentionTime;
        private readonly MsDataScan[] _ms2DataScansByRetentionTime;
        private readonly Dictionary<int, List<double>> _oneBasedScanIndexToDeconvolutedPrecursorMz;
        private readonly double[] _arrayOfRTs;
        private CommonParameters myCommonParameters;
        private CommonParameters FileSpecificParameters;
        private Deconvoluter _classicDeconvoluter;
        private MsDataFile _msDataFile;
        private readonly string _fullFilePath;
        
        /// <summary>
        /// PeakRtApex +/- RetentionTimeWindowHalfWidth sets the retention time range for recovered spectra.
        /// </summary>
        public double RetentionTimeWindowHalfWidth { get; }

        // Each instance of MCSE will be specific to one file.
        public MiniClassicSearchEngine(
            Ms2ScanWithSpecificMass[] arrayOfRTSortedMS2Scans,
            MassDiffAcceptor searchMode,
            CommonParameters commonParameters,
            SpectralLibrary spectralLibrary,
            CommonParameters fileSpecificParameters,
            MsDataFile dataFile,
            string fullFilePath,
            double retentionTimeWindowHalfWidth = 1.0)
        {
            _msDataFile = dataFile;
            _ms2DataScansByRetentionTime = _msDataFile
                .GetAllScansList()
                .Where(s => s.MsnOrder == 2)
                .OrderBy(s => s.RetentionTime)
                .ToArray();
            _arrayOfRTs = _ms2DataScansByRetentionTime.Select(s => s.RetentionTime).ToArray();
            _fullFilePath = fullFilePath;

            SearchMode = searchMode;
            SpectralLibrary = spectralLibrary;
            myCommonParameters = commonParameters;
            FileSpecificParameters = fileSpecificParameters ?? commonParameters;
            RetentionTimeWindowHalfWidth = retentionTimeWindowHalfWidth;
            
            DeconvolutionParameters deconParameters = new ClassicDeconvolutionParameters(
                minCharge: 1,
                maxCharge: 10,
                deconPpm: fileSpecificParameters.PrecursorMassTolerance.Value,
                intensityRatio: fileSpecificParameters.DeconvolutionIntensityRatio);
            _classicDeconvoluter = new Deconvoluter(DeconvolutionTypes.ClassicDeconvolution, deconParameters);
            
        }

        /// <summary>
        /// Searches all Ms2 scans in a 2 minute window around a peak apex for possible PSMs.
        /// Psms are scored against a peptideWithSetModifications that acted as a donor in MBR.
        /// Calculates traditional and spectral contrast angle scores
        /// </summary>
        /// <param name="donorPwsm"> Ms2 scans in window are searched for matches to this donor peptide</param>
        /// <param name="peakApexRT"> The center of the 2 minute window where the search occurs</param>
        /// <returns></returns>
        public IEnumerable<PeptideSpectralMatch> SearchAroundPeak(PeptideWithSetModifications donorPwsm,
            double peakApexRT, int peakCharge, ChromatographicPeak peak = null)
        {
            var targetFragmentsForEachDissociationType = new Dictionary<DissociationType, List<Product>>();

            // check if we're supposed to autodetect dissociation type from the scan header or not
            if (FileSpecificParameters.DissociationType == DissociationType.Autodetect)
            {
                foreach (var item in GlobalVariables.AllSupportedDissociationTypes.Where(p => p.Value != DissociationType.Autodetect))
                {
                    targetFragmentsForEachDissociationType.Add(item.Value, new List<Product>());
                }
            }
            else
            {
                targetFragmentsForEachDissociationType.Add(FileSpecificParameters.DissociationType, new List<Product>());
            }

            // score each scan that has an acceptable precursor mass
            IEnumerable<RecoveredMs2ScanWithSpecificMass> acceptableScans = RecoverScans(donorPwsm.MonoisotopicMass, peak);
            List<PeptideSpectralMatch> acceptablePsms = new();

            foreach (RecoveredMs2ScanWithSpecificMass scan in acceptableScans)
            {
                var dissociationType = FileSpecificParameters.DissociationType == DissociationType.Autodetect & scan.TheScan.DissociationType != null ?
                    scan.TheScan.DissociationType.Value : FileSpecificParameters.DissociationType;

                if (!targetFragmentsForEachDissociationType.TryGetValue(dissociationType, out var peptideTheorProducts))
                {
                    //TODO: print some kind of warning here. the scan header dissociation type was unknown
                    continue;
                }

                // check if we've already generated theoretical fragments for this peptide+dissociation type
                if (peptideTheorProducts.Count == 0)
                {
                    donorPwsm.Fragment(dissociationType, FileSpecificParameters.DigestionParams.FragmentationTerminus, peptideTheorProducts);
                    targetFragmentsForEachDissociationType[dissociationType] = peptideTheorProducts;
                }

                // match theoretical target ions to spectrum
                List<MatchedFragmentIon> matchedIons = MetaMorpheusEngine.MatchFragmentIons(scan, peptideTheorProducts, FileSpecificParameters);

                // calculate the peptide's score
                double thisScore = MetaMorpheusEngine.CalculatePeptideScore(scan.TheScan, matchedIons);

                // Add psm to list
                acceptablePsms.Add(new PeptideSpectralMatch(donorPwsm, notch: -1, thisScore, scanIndex: -1 , scan, FileSpecificParameters, matchedIons));
            }

            foreach (PeptideSpectralMatch psm in acceptablePsms)
            {
                psm.ResolveAllAmbiguities();
            }

            CalculateSpectralAngles(SpectralLibrary, acceptablePsms, donorPwsm, peak.Apex.ChargeState, myCommonParameters, FileSpecificParameters);

            return acceptablePsms;
        }

        private IEnumerable<RecoveredMs2ScanWithSpecificMass> RecoverScans(double donorPeptideMonoisotopicMass, ChromatographicPeak peak)
        {
            MsDataScan[] arrayOfSortedDataScans = GetDataScansInWindow(peak.Apex.IndexedPeak.RetentionTime);
            double donorPeptideMz = donorPeptideMonoisotopicMass.ToMz(peak.Apex.ChargeState);
            double[] myScanIsolationWindowMzStart = arrayOfSortedDataScans.Select(p => p.IsolationRange.Minimum).ToArray();
            int scanIndex = GetFirstDataScanWithIsolationOverOrEqual(donorPeptideMz, myScanIsolationWindowMzStart);

            while (scanIndex < arrayOfSortedDataScans.Length)
            {
                var scan = arrayOfSortedDataScans[scanIndex];
                if (scan.IsolationRange.Contains(donorPeptideMz))
                {

                    yield return DeconvoluteAcceptorScan(scan, donorPeptideMz, peak.Apex.ChargeState);
                    scanIndex++;
                }
                else
                {
                    break;
                }
            }
        }

        // TODO: Write method to deconvolute and store only these msDataScans
        private MsDataScan[] GetDataScansInWindow(double apexRT)
        {
            int startIndex = Array.BinarySearch(_arrayOfRTs, apexRT - RetentionTimeWindowHalfWidth);
            if (startIndex < 0)
                startIndex = ~startIndex;
            int endIndex = Array.BinarySearch(_arrayOfRTs, apexRT + RetentionTimeWindowHalfWidth);
            if (endIndex < 0)
                endIndex = ~endIndex;

            return _ms2DataScansByRetentionTime[startIndex..endIndex].OrderBy(b => b.IsolationRange.Minimum).ToArray();
        }

        private RecoveredMs2ScanWithSpecificMass DeconvoluteAcceptorScan(MsDataScan ms2scan, double donorPeptideMz, int peakCharge)
        {
            if (ms2scan.OneBasedPrecursorScanNumber == null)
            {
                return null;
            }
            MsDataScan precursorSpectrum = _msDataFile.GetOneBasedScan(ms2scan.OneBasedPrecursorScanNumber.Value);
            int closestPrecursorPeakIndex = precursorSpectrum.MassSpectrum.GetClosestPeakIndex(donorPeptideMz);
            MzPeak closestPrecursorPeak = new MzPeak(precursorSpectrum.MassSpectrum.XArray[closestPrecursorPeakIndex],
                precursorSpectrum.MassSpectrum.YArray[closestPrecursorPeakIndex]);
            try
            {
                ms2scan.RefineSelectedMzAndIntensity(precursorSpectrum.MassSpectrum);
            }
            catch (MzLibException ex)
            {
                return null;
            }

            // Not sure if this is necessary
            if (ms2scan.SelectedIonMonoisotopicGuessMz.HasValue)
            {
                ms2scan.ComputeMonoisotopicPeakIntensity(precursorSpectrum.MassSpectrum);
            }

            IsotopicEnvelope closestPrecursorEnvelope = null;
            double closestDeconvolutedPrecursorMz = 0;
            double minMzDistance = ms2scan.IsolationRange.Maximum - donorPeptideMz;
            #pragma warning disable CS0612 // Deconvolution method is obsolete
            foreach (IsotopicEnvelope envelope in ms2scan.GetIsolatedMassesAndCharges(
                         precursorSpectrum.MassSpectrum, 
                         minAssumedChargeState: peakCharge,
                         maxAssumedChargeState: peakCharge,
                         FileSpecificParameters.DeconvolutionMassTolerance.Value,
                         FileSpecificParameters.DeconvolutionIntensityRatio))
            #pragma warning restore CS0612
            {
                double monoPeakMz = envelope.MonoisotopicMass.ToMz(envelope.Charge);
                double thisMzDistance = Math.Abs(monoPeakMz - donorPeptideMz);
                if (thisMzDistance < minMzDistance)
                {
                    minMzDistance = thisMzDistance;
                    closestDeconvolutedPrecursorMz = monoPeakMz;
                    closestPrecursorEnvelope = envelope;
                }
            }

            return new RecoveredMs2ScanWithSpecificMass(
                ms2scan, closestDeconvolutedPrecursorMz, peakCharge, _fullFilePath, FileSpecificParameters,
                closestPrecursorPeak, closestPrecursorEnvelope, neutralExperimentalFragments: null); // neutralFragments are generated in constructor
        }

        /// <summary>
        /// Finds MS2 scans falling within the relevant time window
        /// </summary>
        /// <param name="mbrPeak"> Acceptor peak </param>
        /// <param name="arrayOfRTs"> </param>
        /// <param name="arrayOfMs2ScansSortedByRT"> </param>
        /// <returns> An array of MS2 scans falling within a 2 minute retention time window of the Acceptor peak apex.
        ///           This array is sorted by precursor mass. </returns>
        private Ms2ScanWithSpecificMass[] GetScansInWindow(double apexRT)
        {
            int startIndex = Array.BinarySearch(_arrayOfRTs, apexRT - RetentionTimeWindowHalfWidth);
            if (startIndex < 0)
                startIndex = ~startIndex;
            int endIndex = Array.BinarySearch(_arrayOfRTs, apexRT + RetentionTimeWindowHalfWidth);
            if (endIndex < 0)
                endIndex = ~endIndex;

            return MS2ByRetentionTime[startIndex..endIndex].OrderBy(b => b.PrecursorMass).ToArray();
        }

        private int GetFirstScanWithMassOverOrEqual(double minimum, double[] precursorMasses)
        {
            int index = Array.BinarySearch(precursorMasses, minimum);
            if (index < 0)
            {
                index = ~index;
            }

            // index of the first element that is larger than value
            return index;
        }

        private int GetFirstDataScanWithIsolationOverOrEqual(double minimum, double[] isolationWindowMinimums)
        {
            int index = Array.BinarySearch(isolationWindowMinimums, minimum);
            if (index < 0)
            {
                index = ~index;
            }

            return index;
        }

        public static void CalculateSpectralAngles(SpectralLibrary spectralLibrary, List<PeptideSpectralMatch> psms, PeptideWithSetModifications donorPeptide, int chargeState,
             CommonParameters commonParameters, CommonParameters fileSpecificParameters)
        {
            if (spectralLibrary != null)
            {
                if(!spectralLibrary.TryGetSpectrum(donorPeptide.FullSequence, chargeState, out var librarySpectrum))
                {
                    throw new MetaMorpheusException("Donor peptide library spectrum not found");
                }

                foreach (var psm in psms)
                {
                    SpectralSimilarity s = new SpectralSimilarity(psm.MsDataScan.MassSpectrum, librarySpectrum.XArray, librarySpectrum.YArray,
                        SpectralSimilarity.SpectrumNormalizationScheme.squareRootSpectrumSum, fileSpecificParameters.ProductMassTolerance.Value, false);
                    psm.SpectralAngle = s.SpectralContrastAngle() ?? -2;
                }
            }
        }

        private IEnumerable<ScanWithIndexAndNotchInfo> GetAcceptableScans(double peptideMonoisotopicMass, double apexRT, MassDiffAcceptor searchMode)
        {
            Ms2ScanWithSpecificMass[] arrayOfSortedMs2Scans = GetScansInWindow(apexRT);
            double[] myScanPrecursorMasses = arrayOfSortedMs2Scans.Select(p => p.PrecursorMass).ToArray();
            foreach (AllowedIntervalWithNotch allowedIntervalWithNotch in searchMode.GetAllowedPrecursorMassIntervalsFromTheoreticalMass(peptideMonoisotopicMass).ToList())
            {
                DoubleRange allowedInterval = allowedIntervalWithNotch.AllowedInterval;
                int scanIndex = GetFirstScanWithMassOverOrEqual(allowedInterval.Minimum, myScanPrecursorMasses);
                if (scanIndex < arrayOfSortedMs2Scans.Length)
                {
                    var scanMass = myScanPrecursorMasses[scanIndex];
                    while (scanMass <= allowedInterval.Maximum)
                    {
                        var scan = arrayOfSortedMs2Scans[scanIndex];
                        yield return new ScanWithIndexAndNotchInfo(scan, allowedIntervalWithNotch.Notch, scanIndex);
                        scanIndex++;
                        if (scanIndex == arrayOfSortedMs2Scans.Length)
                        {
                            break;
                        }

                        scanMass = myScanPrecursorMasses[scanIndex];
                    }
                }
            }
        }
    }

    public class RecoveredMs2ScanWithSpecificMass : Ms2ScanWithSpecificMass
    {
        public MzPeak PeakClosestToDonor;
        public IsotopicEnvelope EnvelopeClosestToDonor;
 
        public RecoveredMs2ScanWithSpecificMass(
            MsDataScan dataScan, 
            double precursorMonoisotopicPeakMz, 
            int precursorCharge, 
            string fullFilePath,
            CommonParameters commonParam, 
            MzPeak peakClosestToDonor,
            IsotopicEnvelope envelopeClosestToDonor,
            IsotopicEnvelope[] neutralExperimentalFragments = null) :
            base(dataScan, precursorMonoisotopicPeakMz, precursorCharge, fullFilePath, commonParam,
                neutralExperimentalFragments)
        {
            PeakClosestToDonor = peakClosestToDonor;
            EnvelopeClosestToDonor = envelopeClosestToDonor;
        }
    }

    internal class DataScanWithIndex
    {
        public MsDataScan Scan { get; }
        public int Index { get; }

        internal DataScanWithIndex(MsDataScan scan, int scanIndex)
        {
            Scan = scan;
            Index = scanIndex;
        }
    }
}