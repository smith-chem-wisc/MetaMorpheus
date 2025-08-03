using MassSpectrometry;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Text.RegularExpressions;
using System.Threading.Tasks;

namespace EngineLayer.DIA
{
    public class ISDEngine : DIAEngine
    {
        public ISDEngine(DIAparameters DIAparameters, MsDataFile dataFile, CommonParameters commonParameters, List<(string FileName, CommonParameters Parameters)> fileSpecificParameters, List<string> nestedIds) : base(DIAparameters, dataFile, commonParameters, fileSpecificParameters, nestedIds)
        {
        }

        protected override MetaMorpheusEngineResults RunSpecific()
        {
            return new MetaMorpheusEngineResults(this);
        }

        public static void ReLabelIsdScans(Dictionary<double, List<MsDataScan>> isdScans, List<double> isdEnergies, MsDataScan[] allRawScans)
        {
            foreach(var voltage in isdEnergies)
            {
                int voltageIndex = 1;
                foreach(var scan in isdScans[voltage])
                {
                    scan.SetMsnOrder(2);
                    var ms1Scan = allRawScans[scan.OneBasedScanNumber - voltageIndex - 1];
                    scan.SetOneBasedPrecursorScanNumber(scan.OneBasedScanNumber - voltageIndex);
                    scan.SetIsolationRange(ms1Scan.ScanWindowRange.Minimum, ms1Scan.ScanWindowRange.Maximum);
                    var scanWindowWidth = ms1Scan.ScanWindowRange.Maximum - ms1Scan.ScanWindowRange.Minimum;
                    scan.SetIsolationMz(ms1Scan.ScanWindowRange.Minimum + scanWindowWidth / 2);
                }
                voltageIndex++;
            }
        }

        public static Dictionary<double, List<MsDataScan>> ConstructIsdGroups(MsDataScan[] scans)
        {
            var isdScanVoltageMap = new Dictionary<double, List<MsDataScan>>();
            string pattern = $@"sid=(\d+)";
            foreach (var scan in scans)
            {
                double voltage = 0;
                var match = Regex.Match(scan.ScanFilter, pattern);
                if (match.Success) voltage = double.Parse(match.Groups[1].Value);
                if (!isdScanVoltageMap.ContainsKey(voltage))
                {
                    isdScanVoltageMap[voltage] = new List<MsDataScan> { scan };
                }
                else
                {
                    isdScanVoltageMap[voltage].Add(scan);
                }
            }
            return isdScanVoltageMap;
        }
    }
}
