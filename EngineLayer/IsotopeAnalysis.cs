using Chemistry;
using IO.MzML;
using MassSpectrometry;
using MzLibUtil;
using NetSerializer;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using IO.ThermoRawFileReader;

namespace EngineLayer
{
    public class IsotopeAnalysis
    {
        public Dictionary<string, MsDataScan[]> _ms1Scans { get; private set; }

        public List<PeptideSpectralMatch> AllPsms { get; set; }

        public void ScoreEnvelopes()
        {
            var psmsByFile = AllPsms.OrderBy(p => p.ScanNumber).GroupBy(p => p.FullFilePath);
            foreach (var psmGroup in psmsByFile)
            {
                IEnumerable<int> precursorScanIndices = psmGroup.Where(p => p.PrecursorScanNumber != null).
                    Select(p => (int)p.PrecursorScanNumber).Distinct();
                Queue<int> precursorScanQueue = new Queue<int>(precursorScanIndices);
                Queue<MsDataScan> ms1Scans = null;
                for (int i = 0; i < _ms1Scans[psmGroup.Key].Length; i++ )
                {
                    if (_ms1Scans[psmGroup.Key][i] != null && _ms1Scans[psmGroup.Key][i].OneBasedScanNumber == precursorScanQueue.Peek())
                    {
                        ms1Scans.Enqueue(_ms1Scans[psmGroup.Key][i]);
                        precursorScanQueue.Dequeue();
                    }
                }

                foreach (PeptideSpectralMatch psm in psmGroup)
                {
                    
                }
                
            }
        }

        public bool ReadMS1Scans(string filePath, bool silent)
        {
            if (!silent)
            {
                Console.WriteLine("Reading spectra file");
            }

            MsDataScan[] msDataScans = null;

            // read spectra file
            string ext = Path.GetExtension(filePath).ToUpperInvariant();
            if (ext.Equals(".MZML"))
            {
                try
                {
                    msDataScans = Mzml.LoadAllStaticData(filePath).GetAllScansList()
                        .OrderBy(p => p.OneBasedScanNumber).ToArray();
                }
                catch (FileNotFoundException)
                {
                    if (!silent)
                    {
                        Console.WriteLine("\nCan't find .mzML file" + filePath + "\n");
                    }

                    return false;
                }
                catch (Exception e)
                {
                    if (!silent)
                    {
                        Console.WriteLine("Problem opening .mzML file " + filePath + "; " +
                                          e.Message);
                    }

                    return false;
                }

                for (int i = 0; i < msDataScans.Length; i++)
                {
                    if (msDataScans[i].MsnOrder > 1)
                    {
                        msDataScans[i] = null;
                    }
                }
            }
            else if (ext.Equals(".RAW"))
            {
                var tempList = new List<MsDataScan>();
                ThermoDynamicData dynamicConnection = null;

                try
                {
                    dynamicConnection = new ThermoDynamicData(filePath);

                    // use thermo dynamic connection to get the ms1 scans and then dispose of the connection
                    for (int i = 0; i < dynamicConnection.MsOrdersByScan.Length; i++)
                    {
                        if (dynamicConnection.MsOrdersByScan[i] == 1)
                        {
                            tempList.Add(dynamicConnection.GetOneBasedScanFromDynamicConnection(i + 1));
                        }
                        else
                        {
                            tempList.Add(null);
                        }
                    }

                    dynamicConnection.CloseDynamicConnection();
                }
                catch (FileNotFoundException)
                {
                    if (dynamicConnection != null)
                    {
                        dynamicConnection.CloseDynamicConnection();
                    }

                    if (!silent)
                    {
                        Console.WriteLine("\nCan't find .raw file" + filePath + "\n");
                    }

                    return false;
                }
                catch (Exception e)
                {
                    if (dynamicConnection != null)
                    {
                        dynamicConnection.CloseDynamicConnection();
                    }

                    if (!silent)
                    {
                        throw new MzLibException("FlashLFQ Error: Problem opening .raw file " + filePath + "; " + e.Message);
                    }
                }

                msDataScans = tempList.ToArray();
            }
            else
            {
                if (!silent)
                {
                    Console.WriteLine("Unsupported file type " + ext);
                    return false;
                }
            }

            _ms1Scans.Add(filePath, msDataScans);
            return true;
        }
    }
}
