using MassSpectrometry;
using MassSpectrometry.MzSpectra;
using MzLibUtil;
using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Linq;
using Chemistry;
using FlashLFQ;
using IsotopicEnvelope = MassSpectrometry.IsotopicEnvelope;


namespace EngineLayer.ClassicSearch
{
    public class MiniClassicSearchEngine
    {
        public MsDataFile InstanceSpecificMsDataFile { get; }
        public SpectralLibrary SpectralLibrary { get;  }
        /// <summary>
        /// PeakRtApex +/- RetentionTimeWindowHalfWidth sets the retention time range for recovered spectra.
        /// </summary>
        public double RetentionTimeWindowHalfWidth { get; }

        // no locks are implemented, because it's assumed that the the MsDataScans will not be modified in any way
        // violate this assumption at your own peril!!!
        private readonly MsDataScan[] _ms2DataScansByRetentionTime;
        private readonly double[] _arrayOfMs2RTs; 
        private readonly CommonParameters _fileSpecificParameters;
        private readonly string _fullFilePath;


        /// <summary>
        /// Every instance of MCSE is specific to one file.
        /// </summary>
        /// <param name="dataFile">File to be searched</param>
        /// <param name="spectralLibrary">Spectral library containing donor spectra </param>
        /// <param name="commonParameters"></param>
        /// <param name="fileSpecificParameters"></param>
        /// <param name="fullFilePath"></param>
        /// <param name="retentionTimeWindowHalfWidth"> The width of the retention time window for spectral recovery.
        /// Should be close to 1/2 the expected chromatographic peak width </param>
        /// <exception cref="ArgumentNullException"></exception>
        public MiniClassicSearchEngine(
            MsDataFile dataFile,
            SpectralLibrary spectralLibrary,
            CommonParameters commonParameters,
            CommonParameters fileSpecificParameters,
            string fullFilePath,
            double retentionTimeWindowHalfWidth = 1.0)
        {
            InstanceSpecificMsDataFile = dataFile;
            SpectralLibrary = spectralLibrary;
            if (SpectralLibrary == null) throw new ArgumentNullException();
            _fileSpecificParameters = fileSpecificParameters ?? commonParameters;
            _fullFilePath = fullFilePath;
            RetentionTimeWindowHalfWidth = retentionTimeWindowHalfWidth;

            _ms2DataScansByRetentionTime = InstanceSpecificMsDataFile
                .GetAllScansList()
                .Where(s => s.MsnOrder == 2)
                .Where(s => s.IsolationRange != null)
                .OrderBy(s => s.RetentionTime)
                .ToArray();

            _arrayOfMs2RTs = _ms2DataScansByRetentionTime.Select(s => s.RetentionTime).ToArray();
        }

        /// <summary>
        /// Recovers missed MS2 scans associated with MBR acceptors and calculates the spectral angle
        /// between recovered spectra and the donor spectrum.
        /// </summary>
        /// <param name="donorPwsm"> Donor peptide</param>
        /// <param name="peak"> Acceptor peak </param>
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
            if (!SpectralLibrary.TryGetSpectrum(donorPwsm.FullSequence, peakCharge,
                    out LibrarySpectrum donorSpectrum))
            {
                throw new MetaMorpheusException("Donor spectrum not found");
            }

            Dictionary<DissociationType, List<Product>> targetFragmentsForEachDissociationType = CreateFragmentDictionary();
            List<RecoveredMs2ScanWithSpecificMass> acceptableScans = FindScansInWindow(donorPwsm.MonoisotopicMass, peakRetentionTime, peakCharge);
            List<PeptideSpectralMatch> acceptablePsms = new();
            
            foreach (RecoveredMs2ScanWithSpecificMass scan in acceptableScans)
            {
                var dissociationType = _fileSpecificParameters.DissociationType == DissociationType.Autodetect & scan.TheScan.DissociationType != null 
                    ? scan.TheScan.DissociationType.Value 
                    : _fileSpecificParameters.DissociationType;
                if (!targetFragmentsForEachDissociationType.TryGetValue(dissociationType, out var peptideTheorProducts))
                {
                    //TODO: print some kind of warning here. the scan header dissociation type was unknown
                    continue;
                }

                // check if we've already generated theoretical fragments for this peptide+dissociation type
                if (peptideTheorProducts.Count == 0)
                {
                    donorPwsm.Fragment(dissociationType, _fileSpecificParameters.DigestionParams.FragmentationTerminus, peptideTheorProducts);
                    targetFragmentsForEachDissociationType[dissociationType] = peptideTheorProducts;
                }

                // match theoretical target ions to spectrum
                List<MatchedFragmentIon> matchedIons = MetaMorpheusEngine.MatchFragmentIons(scan, peptideTheorProducts, _fileSpecificParameters);

                // calculate the peptide's score
                double thisScore = MetaMorpheusEngine.CalculatePeptideScore(scan.TheScan, matchedIons);

                PeptideSpectralMatch psm = new PeptideSpectralMatch(donorPwsm, notch: -1, thisScore, scanIndex: -1,
                    scan, _fileSpecificParameters, matchedIons);

                CalculateSpectralAngle(psm, donorSpectrum);
                
                acceptablePsms.Add(psm);
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
            if (_fileSpecificParameters.DissociationType == DissociationType.Autodetect)
            {
                foreach (var item in GlobalVariables.AllSupportedDissociationTypes.Where(p => p.Value != DissociationType.Autodetect))
                {
                    targetFragmentsForEachDissociationType.Add(item.Value, new List<Product>());
                }
            }
            else
            {
                targetFragmentsForEachDissociationType.Add(_fileSpecificParameters.DissociationType, new List<Product>());
            }

            return targetFragmentsForEachDissociationType;
        }

        private List<RecoveredMs2ScanWithSpecificMass> FindScansInWindow(double donorPeptideMonoisotopicMass, double peakRT, int peakCharge)
        {
            double donorPeptideMz = donorPeptideMonoisotopicMass.ToMz(peakCharge);
            (int startIndex, int endIndex) rtSortedIndices = GetIndicesOfScansInRtWindow(peakRT);
            int scanIndex = rtSortedIndices.startIndex;
            List<RecoveredMs2ScanWithSpecificMass> scansInWindow = new();
            
            while (scanIndex < rtSortedIndices.endIndex)
            {
                MsDataScan scan = _ms2DataScansByRetentionTime[scanIndex];
                if (scan.IsolationRange.Contains(donorPeptideMz))
                {
                    RecoveredMs2ScanWithSpecificMass recoveredScan =
                        DeconvoluteAcceptorScan(scan, scanIndex, donorPeptideMz, peakCharge);
                    if (recoveredScan != null)
                    {
                        scansInWindow.Add(recoveredScan);
                    }
                }
                scanIndex++;
            }
            return scansInWindow;
        }

        private (int startIndex, int endIndex) GetIndicesOfScansInRtWindow(double apexRT)
        {
            int startIndex = Array.BinarySearch(_arrayOfMs2RTs, apexRT - RetentionTimeWindowHalfWidth);
            if (startIndex < 0)
                startIndex = ~startIndex;
            int endIndex = Array.BinarySearch(_arrayOfMs2RTs, apexRT + RetentionTimeWindowHalfWidth);
            if (endIndex < 0)
                endIndex = ~endIndex;

            return (startIndex, endIndex);
        }

        /// <summary>
        /// Returns one RecoveredMs2ScanWithSpecificMass. The "specific mass" is the m/z of the deconvoluted monoisotopic peak
        /// that is closest to the theoretical m/z of the donor.
        /// </summary>
        /// <param name="ms2scan"></param>
        /// <param name="scanIndex"></param>
        /// <param name="donorPeptideMz"></param>
        /// <param name="peakCharge"></param>
        /// <returns></returns>
        private RecoveredMs2ScanWithSpecificMass DeconvoluteAcceptorScan(MsDataScan ms2scan, int scanIndex, double donorPeptideMz, int peakCharge)
        {
            if (ms2scan.OneBasedPrecursorScanNumber == null)
            {
                return null;
            }

            MsDataScan precursorSpectrum = InstanceSpecificMsDataFile.GetOneBasedScan(ms2scan.OneBasedPrecursorScanNumber.Value);
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
                         _fileSpecificParameters.DeconvolutionMassTolerance.Value,
                         _fileSpecificParameters.DeconvolutionIntensityRatio))
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
                ms2scan, scanIndex, closestDeconvolutedPrecursorMz, peakCharge, _fullFilePath, _fileSpecificParameters,
                closestPrecursorPeak, closestPrecursorEnvelope, neutralExperimentalFragments: null); // neutralFragments are generated in constructor
        }

        public void CalculateSpectralAngle(PeptideSpectralMatch psm, LibrarySpectrum donorSpectrum)
        {
            SpectralSimilarity s = new SpectralSimilarity(psm.MsDataScan.MassSpectrum, donorSpectrum.XArray, donorSpectrum.YArray,
                    SpectralSimilarity.SpectrumNormalizationScheme.squareRootSpectrumSum, _fileSpecificParameters.ProductMassTolerance.Value, false);
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