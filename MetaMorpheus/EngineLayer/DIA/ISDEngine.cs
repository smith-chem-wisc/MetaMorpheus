using Easy.Common.Extensions;
using MassSpectrometry;
using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Text.RegularExpressions;
using System.Threading.Tasks;
using static Proteomics.RetentionTimePrediction.SSRCalc3;

namespace EngineLayer.DIA
{
    public class ISDEngine : DIAEngine
    {
        private readonly MsDataFile DataFile;
        public List<double> sourceVoltages { get; set; } 
        public ISDEngine(List<double> sourceVoltages, DIAparameters DIAparameters, MsDataFile dataFile, CommonParameters commonParameters, List<(string FileName, CommonParameters Parameters)> fileSpecificParameters, List<string> nestedIds) : base(DIAparameters, dataFile, commonParameters, fileSpecificParameters, nestedIds)
        {
        }

        protected override MetaMorpheusEngineResults RunSpecific()
        {
            //read in scans and isd scan pre-process
            var allScans = DataFile.GetAllScansList().ToArray();
            var isdVoltageMap = ConstructIsdGroups(allScans, out MsDataScan[] ms1Scans);
            ReLabelIsdScans(isdVoltageMap, allScans);

            //Get all MS1 and MS2 XICs
            var allMs1Xics = DIAparams.Ms1XicConstructor.GetAllXics(ms1Scans);
            var allMs2Xics = new Dictionary<double, List<ExtractedIonChromatogram>>();
            foreach (var ms2Group in isdVoltageMap)
            {
                allMs2Xics[ms2Group.Key] = DIAparams.Ms2XicConstructor.GetAllXics(ms2Group.Value.ToArray());
            }

            //Precursor-fragment Grouping
            var allPfGroups = new List<PrecursorFragmentsGroup>();
            if (DIAparams.CombineFragments)
            {
                allPfGroups = DIAparams.PfGroupingEngine.PrecursorFragmentGrouping(allMs1Xics, allMs2Xics.Values.SelectMany(p => p));
            }
            else
            {
                foreach (var ms2Group in isdVoltageMap.Keys)
                {
                    var pfGroups = DIAparams.PfGroupingEngine.PrecursorFragmentGrouping(allMs1Xics, allMs2Xics[ms2Group]);
                    allPfGroups.AddRange(pfGroups);
                }
            }
            
            //Convert pfGroups to pseudo MS2 scans
            PseudoMs2Scans = new List<Ms2ScanWithSpecificMass>();
            int pfGroupIndex = 1;
            foreach (var pfGroup in allPfGroups)
            {
                pfGroup.PFgroupIndex = pfGroupIndex;
                var pseudoScan = PrecursorFragmentsGroup.GetPseudoMs2ScanFromPfGroup(pfGroup, DIAparams.PseudoMs2ConstructionType, CommonParameters, DataFile.FilePath);
                PseudoMs2Scans.Add(pseudoScan);
                pfGroupIndex++;
            }

            return new MetaMorpheusEngineResults(this);
        }

        public static void ReLabelIsdScans(Dictionary<double, List<MsDataScan>> isdVoltageScanMap, MsDataScan[] ms1Scans)
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

        public static Dictionary<double, List<MsDataScan>> ConstructIsdGroups(MsDataScan[] scans, out MsDataScan[] ms1Scans)
        {
            var isdVoltageScanMap = new Dictionary<double, List<MsDataScan>>();
            string pattern = $@"sid=(\d+)";
            foreach (var scan in scans)
            {
                double voltage = 0;
                var match = Regex.Match(scan.ScanFilter, pattern);
                if (match.Success) voltage = double.Parse(match.Groups[1].Value);
                if (!isdVoltageScanMap.ContainsKey(voltage))
                {
                    isdVoltageScanMap[voltage] = new List<MsDataScan> { scan };
                }
                else
                {
                    isdVoltageScanMap[voltage].Add(scan);
                }
            }
            var soretedMap = isdVoltageScanMap.OrderBy(kvp => kvp.Key).ToDictionary(kvp => kvp.Key, kvp => kvp.Value);
            ms1Scans = soretedMap.First().Value.ToArray();
            soretedMap.Remove(soretedMap.First().Key);
            return soretedMap;
        }
    }
}
