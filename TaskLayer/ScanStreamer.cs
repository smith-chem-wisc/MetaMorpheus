using EngineLayer;
using MassSpectrometry;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace TaskLayer
{
    public class ScanStreamer
    {
        private DynamicDataConnection DataStream;
        private string FilePath;
        private CommonParameters commonParameters;
        private Queue<Ms2ScanWithSpecificMass> ForwardScanBuffer;
        private const int ForwardScanBufferSize = 1;

        private Dictionary<int, MsDataScan> ReverseScanBuffer = new Dictionary<int, MsDataScan>();
        private int ReverseScanBufferSize = 50;
        List<(double, int)> precursorBuffer = new List<(double, int)>();

        public bool IsDoneStreaming { get; private set; }
        public int CurrentScan;
        public int NumScans;
        public int calls = 0;

        public ScanStreamer(string filePath, CommonParameters commonParameters)
        {
            MyFileManager fileManager = new MyFileManager(false);
            DataStream = fileManager.OpenDynamicDataConnection(filePath);
            ForwardScanBuffer = new Queue<Ms2ScanWithSpecificMass>(ForwardScanBufferSize);
            CurrentScan = 1;

            this.commonParameters = commonParameters;
            FilePath = filePath;
            IsDoneStreaming = false;
        }

        public Ms2ScanWithSpecificMass NextScan()
        {
            calls++;
            if (IsDoneStreaming)
            {
                return null;
            }

            lock (ForwardScanBuffer)
            {
                if (ForwardScanBuffer.Count > 0)
                {
                    return ForwardScanBuffer.Dequeue();
                }
                else
                {
                    FillForwardBuffer();
                    return NextScan();
                }
            }
        }

        private void FillForwardBuffer()
        {
            while (ForwardScanBuffer.Count < ForwardScanBufferSize)
            {
                // get the scan from the file
                MsDataScan scan = DataStream.GetOneBasedScanFromDynamicConnection(CurrentScan);

                if (scan == null)
                {
                    // done streaming data file
                    DataStream.CloseDynamicConnection();
                    IsDoneStreaming = true;
                    return;
                }

                // add the scan to the buffer
                ReverseScanBuffer.Add(scan.OneBasedScanNumber, scan);

                // if this is a fragmentation scan, get + return precursor info
                if (scan.MsnOrder > 1 && scan.OneBasedPrecursorScanNumber.HasValue)
                {
                    // get the precursor scan
                    int precursorScanNum = scan.OneBasedPrecursorScanNumber.Value;
                    MsDataScan precursorScan;

                    if (!ReverseScanBuffer.TryGetValue(precursorScanNum, out precursorScan))
                    {
                        precursorScan = DataStream.GetOneBasedScanFromDynamicConnection(precursorScanNum);
                    }

                    // deconvolute isolation window + fragmentation scan

                    foreach (Ms2ScanWithSpecificMass scanWithPrecursor in MetaMorpheusTask.DeconvolutePrecursors(scan, precursorScan, commonParameters, FilePath, precursorBuffer))
                    {
                        ForwardScanBuffer.Enqueue(scanWithPrecursor);
                    }
                }

                // remove old scans from buffer that are unlikely to be used later
                while (ReverseScanBuffer.Count > ReverseScanBufferSize)
                {
                    int minScanNumberInBuffer = ReverseScanBuffer.Min(p => p.Key);
                    ReverseScanBuffer.Remove(minScanNumberInBuffer);
                }

                CurrentScan++;
            }
        }
    }
}
