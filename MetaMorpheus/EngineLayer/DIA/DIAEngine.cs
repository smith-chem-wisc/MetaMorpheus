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
    /// <summary>
    /// DIAEgine defines a workflow of generating DDA-like pseudo MS2 scans for DIA data analysis. It includes the processes of extracting precursor and 
    /// fragment XICs, grouping them into PrecursorFragmentsGroup objects, and constructing pseudo MS2 scans.
    /// <summary>
    public class DIAEngine : MetaMorpheusEngine
    {
        private readonly MsDataFile DataFile;
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

        public DIAEngine(DIAparameters DIAparameters, MsDataFile dataFile, CommonParameters commonParameters, List<(string FileName, CommonParameters Parameters)> fileSpecificParameters, List<string> nestedIds) :base(commonParameters, fileSpecificParameters, nestedIds)
        {
            DIAparams = DIAparameters;
            DataFile = dataFile;
        }

        /// <summary>
        /// Divide all MS2 scans into groups based on their isolation windows. 
        /// </summary>
        /// <returns> A dictionary where the key is (minumum of isolation window, maximum of isolation window), 
        /// and the values are MS2 scans that belongs to this isolation window. </returns>
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
