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
using ThermoFisher.CommonCore.Data.Business;
using Readers;

namespace EngineLayer.DIA
{
    /// <summary>
    /// DIAEngine defines a workflow of generating DDA-like pseudo MS2 scans for DIA data analysis. It includes the processes of extracting precursor and 
    /// fragment XICs, grouping them into PrecursorFragmentsGroup objects, and constructing pseudo MS2 scans.
    /// </summary>
    public class DIAEngine
    {
        private readonly MsDataFile DataFile;
        public CommonParameters CommonParams { get; set; }
        public DIAparameters DIAparams => CommonParams.DIAparameters;
        public int OneBasedScanNumber { get; set; } 
        public IEnumerable<Ms2ScanWithSpecificMass> GetPseudoMs2Scans()
        {
            //read in scans
            var ms1Scans = DataFile.GetMS1Scans().ToArray();
            var ms2Scans = DataFile.GetAllScansList().Where(s => s.MsnOrder == 2).ToArray();
            var DIAScanWindowMap = ConstructMs2Groups(ms2Scans);

            foreach (var kvp in DIAScanWindowMap)
            {
                var ms1Xics = DIAparams.Ms1XicConstructor.GetAllXics(ms1Scans, new MzRange(kvp.Key.min, kvp.Key.max));
                var ms2Xics = DIAparams.Ms2XicConstructor.GetAllXics(kvp.Value.ToArray());
                var pfGroups = DIAparams.PfGroupingEngine.PrecursorFragmentGrouping(ms1Xics, ms2Xics);

                foreach (var pfGroup in pfGroups)
                {
                    OneBasedScanNumber++;
                    pfGroup.PFgroupIndex = OneBasedScanNumber;
                    var pseudoScan = PrecursorFragmentsGroup.GetPseudoMs2ScanFromPfGroup(pfGroup, DIAparams.PseudoMs2ConstructionType, CommonParams, DataFile.FilePath);
                    yield return pseudoScan;
                }
            }
        }

        public DIAEngine(MsDataFile dataFile, CommonParameters commonParameters)
        {
            CommonParams = commonParameters;
            DataFile = dataFile;
        }

        /// <summary>
        /// Divide all MS2 scans into groups based on their isolation windows. 
        /// </summary>
        /// <returns> A dictionary where the key is (minumum of isolation window, maximum of isolation window), 
        /// and the values are MS2 scans that belongs to this isolation window. </returns>
        public static Dictionary<(double min, double max), List<MsDataScan>> ConstructMs2Groups(MsDataScan[] ms2Scans, int binSize = 100, double tolerance = 5)
        {
            var absoluteTolerance = new AbsoluteTolerance(tolerance);
            var diaWindowIDMap = new Dictionary<(double min, double max), int>();
            var diaIDScanMap = new SortedList<int, List<MsDataScan>> ();
            foreach (var ms2 in ms2Scans)
            {
                int roundedCenterMz = (int)Math.Round((ms2.IsolationRange.Maximum + ms2.IsolationRange.Minimum) * binSize / 2, 0); // Use the center of the isolation window as the identifier
                (double min, double max) range = new (Math.Round(ms2.IsolationRange.Minimum, 2), Math.Round(ms2.IsolationRange.Maximum, 2));
                int closestWindowIndex = GetClosestWindowIndex(roundedCenterMz, diaIDScanMap, absoluteTolerance);
                if (closestWindowIndex >= 0)
                {
                    int key = diaIDScanMap.Keys[closestWindowIndex];
                    diaIDScanMap[key].Add(ms2);
                }
                else
                {
                    diaIDScanMap[roundedCenterMz] = new List<MsDataScan> { ms2 };
                    diaWindowIDMap[range] = roundedCenterMz;
                }
            }
            var joinedMap = diaWindowIDMap.ToDictionary(kvp => kvp.Key, kvp => diaIDScanMap[kvp.Value]);
            return joinedMap;
        }

        private static int GetClosestWindowIndex(int mz, SortedList<int, List<MsDataScan>> map, Tolerance tolerance)
        {
            int closestIndex = -1;
            if (map.Count == 0) return closestIndex;

            //Find the closest windowID using binary search
            int index = Array.BinarySearch(map.Keys.ToArray(), mz);
            if (index >= 0)
            {
                closestIndex = index;
            }
            else
            {
                index = ~index; 
                if (index == 0)
                {
                    closestIndex = 0;
                }
                else if (index >= map.Count)
                {
                    closestIndex = map.Count - 1;
                }
                else
                {
                    if (Math.Abs(map.Keys[index] - mz) < Math.Abs(map.Keys[index - 1] - mz))
                    {
                        closestIndex = index;
                    }
                    else
                    {
                        closestIndex = index - 1;
                    }
                }
            }

            //check if the closest windowID is within tolerance
            if (tolerance.Within(mz, map.Keys[closestIndex]))
            {
                return closestIndex;
            }
            else
            {
                return -1;
            }
        }
    }
}
