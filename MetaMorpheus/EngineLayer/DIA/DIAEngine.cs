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
using FlashLFQ;

namespace EngineLayer.DIA
{
    public class DIAEngine : MetaMorpheusEngine
    {
        public MsDataFile DataFile { get; set; }
        public DIAparameters DIAparams { get; set; } 
        public List<Ms2ScanWithSpecificMass> PseudoMs2Scans { get; set; } 
        protected override MetaMorpheusEngineResults RunSpecific()
        {
            //read in scans
            var ms1Scans = DataFile.GetMS1Scans().ToArray();
            var ms2Scans = DataFile.GetAllScansList().Where(s => s.MsnOrder == 2).ToArray();
            var DIAScanWindowMap = ConstructMs2Groups(ms2Scans);

            //Get all MS1 and MS2 XICs
            var allMs1Xics = new Dictionary<(double min, double max), List<ExtractedIonChromatogram>>();
            var allMs2Xics = new Dictionary<(double min, double max), List<ExtractedIonChromatogram>>();
            foreach (var ms2Group in DIAScanWindowMap)
            {
                allMs1Xics[ms2Group.Key] = DIAparams.Ms1XicConstructor.GetAllXics(ms1Scans, new MzRange(ms2Group.Key.min, ms2Group.Key.max));
                allMs2Xics[ms2Group.Key] = DIAparams.Ms2XicConstructor.GetAllXics(ms2Group.Value.ToArray());
            }

            //Precursor-fragment Grouping
            var allPfGroups = new List<PrecursorFragmentsGroup>();
            foreach (var ms2Group in DIAScanWindowMap.Keys)
            {
                var pfGroups = DIAparams.PfGroupingEngine.PrecursorFragmentGrouping(allMs1Xics[ms2Group], allMs2Xics[ms2Group]);
                allPfGroups.AddRange(pfGroups);
            }

            //Convert pfGroups to pseudo MS2 scans
            PseudoMs2Scans = new List<Ms2ScanWithSpecificMass>();
            int pfGroupIndex = 1;
            foreach (var pfGroup in allPfGroups)
            {
                pfGroup.PFgroupIndex = pfGroupIndex;
                var pseudoScan = PrecursorFragmentsGroup.GetPseudoMs2ScanFromPfGroup(pfGroup, DIAparams.PseudoMs2ConstructionType, CommonParameters, DataFile.FilePath);
                PseudoMs2Scans.Add(pseudoScan);
            }
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

            //Get all MS1 and MS2 XICs
            var ms1Xics = new List<ExtractedIonChromatogram>();
            var allMs2Xics = new Dictionary<(double min, double max), List<ExtractedIonChromatogram>>();
            foreach(var ms2Group in DIAScanWindowMap)
            {
                allMs2Xics[ms2Group.Key] = diaParam.Ms2XicConstructor.GetAllXics(ms2Group.Value.ToArray());
                diaParam.Ms2XicConstructor.XicSplineForAllXics(allMs2Xics[ms2Group.Key], commonParameters.MaxThreadsToUsePerFile);
            }

            //Precursor-fragment Grouping
            var allPfGroups = new List<PrecursorFragmentsGroup>();
            foreach(var ms2Group in DIAScanWindowMap)
            {
                var pfGroups = diaParam.PfGroupingEngine.PrecursorFragmentGrouping(ms1Xics, allMs2Xics[ms2Group.Key]);
                allPfGroups.AddRange(pfGroups);
            }

            //Convert pfGroups to pseudo MS2 scans
            var pseudoMs2Scans = new List<Ms2ScanWithSpecificMass>();
            int pfGroupIndex = 1;
            foreach (var pfGroup in allPfGroups)
            {
                pfGroup.PFgroupIndex = pfGroupIndex;
                var pseudoScan = PrecursorFragmentsGroup.GetPseudoMs2ScanFromPfGroup(pfGroup, diaParam.PseudoMs2ConstructionType, commonParameters, dataFile.FilePath);
                pseudoMs2Scans.Add(pseudoScan);
            }

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
    }
}
