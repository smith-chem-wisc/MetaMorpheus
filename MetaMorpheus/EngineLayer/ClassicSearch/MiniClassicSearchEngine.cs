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
using Easy.Common.Extensions;

namespace EngineLayer.ClassicSearch
{
    public class MiniClassicSearchEngine
    {
        private readonly SpectralLibrary SpectralLibrary;
        private MsDataFile _msDataFile;
        private readonly MsDataScan[] _ms2DataScansByRetentionTime;
        private object[] _ms2DataScanLocks;
        private readonly double[] _arrayOfMs2RTs;
        private CommonParameters FileSpecificParameters;
        private readonly string _fullFilePath;
        /// <summary>
        /// PeakRtApex +/- RetentionTimeWindowHalfWidth sets the retention time range for recovered spectra.
        /// </summary>
        public double RetentionTimeWindowHalfWidth { get; }

        // Each instance of MCSE will be specific to one file.
        public MiniClassicSearchEngine(
            MsDataFile dataFile,
            SpectralLibrary spectralLibrary,
            CommonParameters commonParameters,
            CommonParameters fileSpecificParameters,
            string fullFilePath,
            double retentionTimeWindowHalfWidth = 1.0)
        {
            SpectralLibrary = spectralLibrary;
            if (SpectralLibrary == null) throw new ArgumentNullException();
            FileSpecificParameters = fileSpecificParameters ?? commonParameters;
            RetentionTimeWindowHalfWidth = retentionTimeWindowHalfWidth;
            _fullFilePath = fullFilePath;

            _msDataFile = dataFile;
            _ms2DataScansByRetentionTime = _msDataFile
                .GetAllScansList()
                .Where(s => s.MsnOrder == 2)
                .Where(s => s.IsolationRange != null)
                .OrderBy(s => s.RetentionTime)
                .ToArray();
            _arrayOfMs2RTs = _ms2DataScansByRetentionTime.Select(s => s.RetentionTime).ToArray();
            _ms2DataScanLocks = new object[_ms2DataScansByRetentionTime.Length];
            for (int i = 0; i < _ms2DataScanLocks.Length; i++)
            {
                _ms2DataScanLocks[i] = new object();
            }
        }

        /// <summary>
        /// Searches all Ms2 scans in a 2 minute window around a peak apex for possible PSMs.
        /// Psms are scored against a peptideWithSetModifications that acted as a donor in MBR.
        /// Calculates traditional and spectral contrast angle scores
        /// </summary>
        /// <param name="donorPwsm"> Ms2 scans in window are searched for matches to this donor peptide</param>
        /// <param name="peakApexRT"> The center of the 2 minute window where the search occurs</param>
        /// <returns></returns>
        public List<PeptideSpectralMatch> SearchAroundPeak(PeptideWithSetModifications donorPwsm,
            ChromatographicPeak peak = null, double peakRetentionTime = -1, int peakCharge = -1)
        {
            if (peak?.Apex != null)
            {
                peakRetentionTime = peak.Apex.IndexedPeak.RetentionTime;
                peakCharge = peak.Apex.ChargeState;
            }
            if (peakRetentionTime < 0 | peakCharge < 0)
            {
                throw new ArgumentException("Peak can not be null!");
            }
            if (!SpectralLibrary.TryGetSpectrum(donorPwsm.FullSequence, peak.Apex.ChargeState,
                    out LibrarySpectrum donorSpectrum))
            {
                throw new MetaMorpheusException("Donor spectrum not found");
            }

            Dictionary<DissociationType, List<Product>> targetFragmentsForEachDissociationType = CreateFragmentDictionary();
            List<RecoveredMs2ScanWithSpecificMass> acceptableScans = RecoverScans(donorPwsm.MonoisotopicMass, peakRetentionTime, peakCharge);
            List<PeptideSpectralMatch> acceptablePsms = new();
            
            foreach (RecoveredMs2ScanWithSpecificMass scan in acceptableScans)
            {
                lock (_ms2DataScanLocks[scan.ScanIndex])
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

                    PeptideSpectralMatch psm = new PeptideSpectralMatch(donorPwsm, notch: -1, thisScore, scanIndex: -1,
                        scan, FileSpecificParameters, matchedIons);

                    CalculateSpectralAngle(psm, donorSpectrum);
                    
                    acceptablePsms.Add(psm);
                }
            }

            foreach (PeptideSpectralMatch psm in acceptablePsms)
            {
                psm.ResolveAllAmbiguities();
            }

            return acceptablePsms;
        }

        private Dictionary<DissociationType, List<Product>> CreateFragmentDictionary()
        {
            Dictionary<DissociationType, List<Product>> targetFragmentsForEachDissociationType = new();

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

            return targetFragmentsForEachDissociationType;
        }

        private List<RecoveredMs2ScanWithSpecificMass> RecoverScans(double donorPeptideMonoisotopicMass, double peakRT, int peakCharge)
        {
            double donorPeptideMz = donorPeptideMonoisotopicMass.ToMz(peakCharge);
            (int startIndex, int endIndex) rtSortedIndices = GetDataScansInWindow(peakRT);
            int scanIndex = rtSortedIndices.startIndex;
            List<RecoveredMs2ScanWithSpecificMass> recoveredScans = new();
            
            while (scanIndex < rtSortedIndices.endIndex)
            {
                lock (_ms2DataScanLocks[scanIndex])
                {
                    MsDataScan scan = _ms2DataScansByRetentionTime[scanIndex];
                    // This could be improved by sorting scans by isolation range, but I'm not sure how that would affect parrallelization
                    if (scan.IsolationRange.Contains(donorPeptideMz))
                    {
                        RecoveredMs2ScanWithSpecificMass recoveredScan =
                            DeconvoluteAcceptorScan(scan, scanIndex, donorPeptideMz, peakCharge);
                        if (recoveredScan != null)
                        {
                            recoveredScans.Add(recoveredScan);
                        }
                    }
                }
                scanIndex++;
            }

            return recoveredScans;
        }

        private (int startIndex, int endIndex) GetDataScansInWindow(double apexRT)
        {
            int startIndex = Array.BinarySearch(_arrayOfMs2RTs, apexRT - RetentionTimeWindowHalfWidth);
            if (startIndex < 0)
                startIndex = ~startIndex;
            int endIndex = Array.BinarySearch(_arrayOfMs2RTs, apexRT + RetentionTimeWindowHalfWidth);
            if (endIndex < 0)
                endIndex = ~endIndex;

            return (startIndex, endIndex);
        }

        private RecoveredMs2ScanWithSpecificMass DeconvoluteAcceptorScan(MsDataScan ms2scan, int scanIndex, double donorPeptideMz, int peakCharge)
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
            catch (MzLibException)
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
                ms2scan, scanIndex, closestDeconvolutedPrecursorMz, peakCharge, _fullFilePath, FileSpecificParameters,
                closestPrecursorPeak, closestPrecursorEnvelope, neutralExperimentalFragments: null); // neutralFragments are generated in constructor
        }

        public void CalculateSpectralAngle(PeptideSpectralMatch psm, LibrarySpectrum donorSpectrum)
        {
            SpectralSimilarity s = new SpectralSimilarity(psm.MsDataScan.MassSpectrum, donorSpectrum.XArray, donorSpectrum.YArray,
                    SpectralSimilarity.SpectrumNormalizationScheme.squareRootSpectrumSum, FileSpecificParameters.ProductMassTolerance.Value, false);
            psm.SpectralAngle = s.SpectralContrastAngle() ?? -1;
        }

    }

    public class RecoveredMs2ScanWithSpecificMass : Ms2ScanWithSpecificMass
    {
        public MzPeak PeakClosestToDonor { get; }
        public IsotopicEnvelope EnvelopeClosestToDonor { get; }
        public int ScanIndex { get; }
 
        public RecoveredMs2ScanWithSpecificMass(
            MsDataScan dataScan, 
            int scanIndex,
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
            ScanIndex = scanIndex;
        }
    }

}