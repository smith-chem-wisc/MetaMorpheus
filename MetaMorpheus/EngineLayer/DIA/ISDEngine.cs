using Easy.Common.Extensions;
using MassSpectrometry;
using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Text.RegularExpressions;

namespace EngineLayer.DIA
{
    /// <summary>
    /// ISDEngine inherits from DIAEngine and employs a similar workflow of generating XICs, grouping XICs and make pseudo Ms2 scans for ISD data analysis. 
    /// </summary>
    public class ISDEngine : DIAEngine
    {
        private readonly MsDataFile DataFile;

        public ISDEngine(MsDataFile dataFile, CommonParameters commonParameters) : base(dataFile, commonParameters)
        {
            DataFile = dataFile;
        }
       
        public override IEnumerable<Ms2ScanWithSpecificMass> GetPseudoMs2Scans()
        {
            //read in scans and isd scan pre-process
            var allScans = DataFile.GetAllScansList().ToArray();
            var isdVoltageMap = ConstructIsdGroups(allScans, out MsDataScan[] ms1Scans);
            ReLabelIsdScans(isdVoltageMap, allScans);

            //Get all MS1 XICs
            var allMs1Xics = DIAparams.Ms1XicConstructor.GetAllXicsWithXicSpline(ms1Scans);

            //Precursor-fragment Grouping
            var allPfGroups = new List<PrecursorFragmentsGroup>();
            if (DIAparams.CombineFragments)
            {
                var allMs2Xics = new Dictionary<int, List<ExtractedIonChromatogram>>();
                foreach (var kvp in isdVoltageMap)
                {
                    var ms2Xics = DIAparams.Ms2XicConstructor.GetAllXics(kvp.Value.ToArray());
                    allMs2Xics.Add(kvp.Key, ms2Xics);
                }
                allPfGroups = DIAparams.PfGroupingEngine.PrecursorFragmentGrouping(allMs1Xics, allMs2Xics.Values.SelectMany(p => p).ToList()).ToList();
            }
            else
            {
                foreach(var kvp in isdVoltageMap)
                {
                    var ms2Xics = DIAparams.Ms2XicConstructor.GetAllXics(kvp.Value.ToArray());
                    var pfGroups = DIAparams.PfGroupingEngine.PrecursorFragmentGrouping(allMs1Xics, ms2Xics);
                    allPfGroups.AddRange(pfGroups);
                }
            }

            foreach (var pfGroup in allPfGroups)
            {
                OneBasedScanNumber++;
                pfGroup.PFgroupIndex = OneBasedScanNumber;
                var pseudoScan = pfGroup.GetPseudoMs2ScanFromPfGroup(DIAparams.PseudoMs2ConstructionType, CommonParams, DataFile.FilePath);
                yield return pseudoScan;
            }
        }

        /// <summary>
        /// Relabel ISD scans to have correct Ms level, precursor scan number and isolation m/z/range.
        /// </summary>
        public static void ReLabelIsdScans(Dictionary<int, List<MsDataScan>> isdVoltageScanMap, MsDataScan[] ms1Scans)
        {
            for (int i = 0; i < isdVoltageScanMap.Count; i++)
            {
                foreach (var scan in isdVoltageScanMap.Values.ElementAt(i))
                {
                    scan.SetMsnOrder(2);
                    int oneBasedPrecursorScanNumber = scan.OneBasedScanNumber - i - 1;
                    var ms1Scan = ms1Scans.Where(s => s.OneBasedScanNumber == oneBasedPrecursorScanNumber).First();
                    scan.SetOneBasedPrecursorScanNumber(oneBasedPrecursorScanNumber);
                    scan.SetIsolationRange(ms1Scan.ScanWindowRange.Minimum, ms1Scan.ScanWindowRange.Maximum);
                    var scanWindowWidth = ms1Scan.ScanWindowRange.Maximum - ms1Scan.ScanWindowRange.Minimum;
                    scan.SetIsolationMz(ms1Scan.ScanWindowRange.Minimum + scanWindowWidth / 2);
                }
            }
        }

        /// <summary>
        /// Divide all ISD scans into groups based on their source energy level and extract real MS1 scans. 
        /// </summary>
        /// <returns> A dictionary where the key is ISD voltage and the values are MS2 scans. </returns>
        public static Dictionary<int, List<MsDataScan>> ConstructIsdGroups(MsDataScan[] scans, out MsDataScan[] ms1Scans)
        {
            var isdVoltageScanMap = new Dictionary<int, List<MsDataScan>>();
            //ISD voltage is written in the scan filter like "sid=15" in the raw file
            string pattern = $@"sid=(\d+)";
            foreach (var scan in scans)
            {
                int voltage = 0;
                var match = Regex.Match(scan.ScanFilter, pattern);
                if (match.Success) voltage = int.Parse(match.Groups[1].Value);
                if (!isdVoltageScanMap.ContainsKey(voltage))
                {
                    isdVoltageScanMap[voltage] = new List<MsDataScan> { scan };
                }
                else
                {
                    isdVoltageScanMap[voltage].Add(scan);
                }
            }
            var sortedMap = isdVoltageScanMap.OrderBy(kvp => kvp.Key).ToDictionary(kvp => kvp.Key, kvp => kvp.Value);
            //ISD runs start with MS1 scans so the first kvp is removed from the dictionary
            ms1Scans = sortedMap.First().Value.ToArray();
            sortedMap.Remove(sortedMap.First().Key);
            return sortedMap;
        }
    }
}
