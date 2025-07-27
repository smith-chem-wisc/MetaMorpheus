using MassSpectrometry;
using MzLibUtil;
using Proteomics;
using Omics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using Omics.Modifications;
using System.Collections.Concurrent;
using EngineLayer.Util;
using Omics;
using ThermoFisher.CommonCore.Data.Business;
using FlashLFQ;

namespace EngineLayer.DIA
{
    public class DIAEngine : MetaMorpheusEngine
    {
        public DIAparameters DIAparameters { get; set; } 
        public List<Ms2ScanWithSpecificMass> PseudoMs2Scans { get; set; } 
        protected override MetaMorpheusEngineResults RunSpecific()
        {

            return new MetaMorpheusEngineResults(this);
        }

        public DIAEngine(DIAparameters DIAparameters, CommonParameters commonParameters, List<(string FileName, CommonParameters Parameters)> fileSpecificParameters, List<string> nestedIds) :base(commonParameters, fileSpecificParameters, nestedIds)
        {
            DIAparameters = DIAparameters;
        }

        public static List<Ms2ScanWithSpecificMass> GetPseudoMs2Scans(MsDataFile dataFile, CommonParameters commonParameters, DIAparameters diaParam)
        {
            //read in scans
            var ms1Scans = dataFile.GetMS1Scans().ToArray();
            var ms2Scans = dataFile.GetAllScansList().Where(s => s.MsnOrder == 2).ToArray();
            var DIAScanWindowMap = ConstructMs2Groups(ms2Scans);

            //Get all MS1 XICs
            var Ms1PeakIndexingEngine = new Dictionary<(double min, double max), IndexingEngine<IIndexedPeak>>();

            var pseudoMs2Scans = new List<Ms2ScanWithSpecificMass>();
            return pseudoMs2Scans;
        }

        public static Dictionary<(double min, double max), List<MsDataScan>> ConstructMs2Groups(MsDataScan[] ms2Scans)
        {
            var DIAScanWindowMap = new Dictionary<(double min, double max), List<MsDataScan>>();
            foreach (var ms2 in ms2Scans)
            {
                (double min, double max) range = new(Math.Round(ms2.IsolationRange.Minimum, 2), Math.Round(ms2.IsolationRange.Maximum, 2));
                if (!DIAScanWindowMap.ContainsKey(range))
                {
                    DIAScanWindowMap[range] = new List<MsDataScan> { ms2 };
                }
                else
                {
                    DIAScanWindowMap[range].Add(ms2);
                }
            }
            return DIAScanWindowMap;
        }

        public static List<ExtractedIonChromatogram> GetAllXics(MsDataScan[] scans, IndexingEngineType indexingEngineType, Tolerance peakFindingTolerance, int maxMissedScanAllowed, double maxRTRange, int numPeakThreshold, DeconvolutionParameters deconParam = null)
        {
            switch (indexingEngineType)
            {
                case IndexingEngineType.MzPeak:
                    var mzPeakIndexingEngine = PeakIndexingEngine.InitializeIndexingEngine(scans);
                    return mzPeakIndexingEngine.GetAllXics(peakFindingTolerance, maxMissedScanAllowed, maxRTRange, numPeakThreshold);
                case IndexingEngineType.NeutralMass:
                    var neutralMassIndexingEngine = MassIndexingEngine.InitializeMassIndexingEngine(scans, deconParam);
                    return neutralMassIndexingEngine.GetAllXics(peakFindingTolerance, maxMissedScanAllowed, maxRTRange, numPeakThreshold);
                default:
                    throw new ArgumentException("Invalid indexing engine type specified.");
            }
        }

        public static List<PrecursorFragmentsGroup> GetAllPFGroupsInWindow(List<ExtractedIonChromatogram> precursorXics, List<ExtractedIonChromatogram> fragmentXics, DIAparameters diaParam)
        {
            return null;
        }
    }
}
