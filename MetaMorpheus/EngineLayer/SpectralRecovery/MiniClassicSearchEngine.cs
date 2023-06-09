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
using Easy.Common.Extensions;
using ThermoFisher.CommonCore.Data;


namespace EngineLayer.SpectralRecovery
{
    public class MiniClassicSearchEngine
    {
        /// <summary>
        /// PeakRtApex +/- RetentionTimeWindowHalfWidth sets the retention time range for recovered spectra.
        /// </summary>
        public double RetentionTimeWindowHalfWidth { get; }
        public MsDataFile InstanceSpecificMsDataFile { get; }
        public SpectralLibrary SpectralLibrary { get;  }

        // no locks are implemented, because it's assumed that the the MsDataScans will not be modified in any way
        // violate this assumption at your own peril!!!
        private readonly MsDataScan[] _ms2DataScansByRetentionTime;
        private readonly double[] _arrayOfMs2RTs; 
        private readonly CommonParameters _fileSpecificParameters;
        private readonly string _fullFilePath;
        private readonly double[] _arrayOfAllScanRTs;
        private readonly bool _maxQuantAnalysis;

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
            double retentionTimeWindowHalfWidth = 1.0,
            bool maxQuantAnalysis = false)
        {
            InstanceSpecificMsDataFile = dataFile;
            SpectralLibrary = spectralLibrary;
            if (SpectralLibrary == null) throw new ArgumentNullException();
            _fileSpecificParameters = fileSpecificParameters ?? commonParameters;
            _fullFilePath = fullFilePath;
            RetentionTimeWindowHalfWidth = retentionTimeWindowHalfWidth;
            _maxQuantAnalysis = maxQuantAnalysis;

            _ms2DataScansByRetentionTime = InstanceSpecificMsDataFile
                .GetAllScansList()
                .Where(s => s.MsnOrder == 2)
                .Where(s => s.IsolationRange != null)
                .OrderBy(s => s.RetentionTime)
                .ToArray();

            _arrayOfMs2RTs = _ms2DataScansByRetentionTime.Select(s => s.RetentionTime).ToArray();
            _arrayOfAllScanRTs = dataFile.GetAllScansList().Select(s => s.RetentionTime).ToArray();
        }

        /// <summary>
        /// Recovers missed MS2 scans associated with MBR acceptors and calculates the spectral angle
        /// between recovered spectra and the donor spectrum.
        /// </summary>
        /// <param name="donorPwsm"> Donor peptide</param>
        /// <param name="peak"> Acceptor peak </param>
        /// <returns></returns>
        public List<SpectralRecoveryPSM> SearchAroundPeak(PeptideWithSetModifications donorPwsm,
            ChromatographicPeak peak = null, double peakRetentionTime = -1, int peakCharge = -1)
        {
            if (peak?.Apex != null)
            {
                peakRetentionTime = peak.Apex.IndexedPeak.RetentionTime;
                peakCharge = peak.Apex.ChargeState;
            } 
            // This is used for peaks read from MaxQuant
            else if(peak != null && peak.IsotopicEnvelopes.IsNotNullOrEmpty())
            {
                peakRetentionTime = peak.IsotopicEnvelopes[0].IndexedPeak.RetentionTime;
                peakCharge = peak.IsotopicEnvelopes[0].ChargeState;
            }

            if (peakRetentionTime < 0 | peakCharge < 0)
            {
                throw new ArgumentException("Peak can not be null!");
            }

            if (!SpectralLibrary.TryGetSpectrum(donorPwsm.FullSequence, peakCharge,
                    out LibrarySpectrum donorSpectrum))
            {
                return null;
            }

            Dictionary<DissociationType, List<Product>> targetFragmentsForEachDissociationType = CreateFragmentDictionary();
            List<RecoveredMs2ScanWithSpecificMass> acceptableScans = FindScansInWindow(donorPwsm.MonoisotopicMass, peakRetentionTime, peakCharge);
            List<SpectralRecoveryPSM> acceptablePsms = new();
            
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

                // For MaxQuant analysis, the peak doesn't contain the isotopic envelope. So, the precursor spectrum is scored
                // TODO: Refactor so this isn't spaghetti
                MzSpectrum precursorSpectrum = null;
                if (_maxQuantAnalysis)
                {
                    double retentionTime = peak.IsotopicEnvelopes.First().IndexedPeak.RetentionTime;
                    int oneBasedScanNumber = FindOneBasedScanNumber(retentionTime) + 1;
                    while (oneBasedScanNumber > 0)
                    {
                        MsDataScan precursorScan = InstanceSpecificMsDataFile.GetOneBasedScan(oneBasedScanNumber);
                        if (precursorScan.MsnOrder == 1)
                        {
                            precursorSpectrum = precursorScan.MassSpectrum;
                            break;
                        }

                        oneBasedScanNumber--;
                    }
                }

                SpectralRecoveryPSM psm = new SpectralRecoveryPSM(peak, donorPwsm, thisScore, scan,
                    _fileSpecificParameters, matchedIons, precursorSpectrum);

                CalculateSpectralAngle(psm, donorSpectrum);
                
                acceptablePsms.Add(psm);
            }

            if (!_maxQuantAnalysis)
            {
                foreach (var psm in acceptablePsms)
                {
                    psm.ResolveAllAmbiguities();
                }
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

        private int FindOneBasedScanNumber(double retentionTime)
        {
            int index = Array.BinarySearch(_arrayOfAllScanRTs, retentionTime);
            return index < 0 ? ~index : index;
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
                }
            }

            // TODO: Write more efficient method for getting precursors
            // Also, should probably check the isotopic envelope at the apex peak.
            return new RecoveredMs2ScanWithSpecificMass(
                ms2scan, scanIndex, closestDeconvolutedPrecursorMz, peakCharge, _fullFilePath, _fileSpecificParameters,
                closestPrecursorPeak, neutralExperimentalFragments: null); // neutralFragments are generated in constructor
        }

        public void CalculateSpectralAngle(PeptideSpectralMatch psm, LibrarySpectrum donorSpectrum)
        {
            if (psm.MsDataScan.MassSpectrum.Size == 0)
            {
                psm.SpectralAngle = -1;
                return; 
            }
            SpectralSimilarity s = new SpectralSimilarity(psm.MsDataScan.MassSpectrum, donorSpectrum.XArray, donorSpectrum.YArray,
                    SpectralSimilarity.SpectrumNormalizationScheme.squareRootSpectrumSum, _fileSpecificParameters.ProductMassTolerance.Value, false);
            psm.SpectralAngle = s.SpectralContrastAngle() ?? -1;
        }

    }

    public class RecoveredMs2ScanWithSpecificMass : Ms2ScanWithSpecificMass
    {
        public MzPeak PeakClosestToDonor { get; }

        public RecoveredMs2ScanWithSpecificMass(
            MsDataScan dataScan,
            int scanIndex,
            double precursorMonoisotopicPeakMz,
            int precursorCharge,
            string fullFilePath,
            CommonParameters commonParam,
            MzPeak peakClosestToDonor,
            IsotopicEnvelope[] neutralExperimentalFragments = null) :
            base(dataScan, precursorMonoisotopicPeakMz, precursorCharge, fullFilePath, commonParam,
                neutralExperimentalFragments)
        {
            PeakClosestToDonor = peakClosestToDonor;
        }
    }

}