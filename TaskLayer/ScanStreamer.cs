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
        private int CurrentScan;
        private readonly string FilePath;
        private readonly CommonParameters commonParameters;
        private readonly FilteringParams filteringParams;

        // buffer stuff
        private readonly int ForwardScanBufferSize;
        private readonly int BackScanBufferSize;
        private readonly Queue<Ms2ScanWithSpecificMass> ForwardScanBuffer;
        private readonly Dictionary<int, MsDataScan> BackScanBuffer;
        private readonly List<(double, int)> precursorBuffer;

        /// <summary>
        /// Returns true if there are more scans in the data stream; otherwise, returns false.
        /// </summary>
        public bool IsDoneStreaming { get; private set; }

        /// <summary>
        /// Allows thread-safe streaming of an MS data file to Ms2ScanWithSpecificMass objects. 
        /// Useful for large data files that can't be loaded into memory all at once (e.g., timsTOF).
        /// </summary>
        public ScanStreamer(string filePath, CommonParameters commonParameters)
        {
            DataStream = MyFileManager.OpenDynamicDataConnection(filePath);
            filteringParams = MyFileManager.GetFilterParamsFromCommonParams(commonParameters);
            CurrentScan = 1;

            ForwardScanBufferSize = 1;
            BackScanBufferSize = 100;
            ForwardScanBuffer = new Queue<Ms2ScanWithSpecificMass>(ForwardScanBufferSize * 50);
            BackScanBuffer = new Dictionary<int, MsDataScan>(BackScanBufferSize * 2);
            precursorBuffer = new List<(double, int)>(ForwardScanBufferSize * 50);

            this.commonParameters = commonParameters;
            FilePath = filePath;
            IsDoneStreaming = false;
        }

        /// <summary>
        /// Retrieves the next Ms2ScanWithSpecificMass from the data stream. Returns null at the end of the stream.
        /// </summary>
        public Ms2ScanWithSpecificMass NextScan()
        {
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
                    FillBuffers();
                    return NextScan();
                }
            }
        }

        /// <summary>
        /// Gets the next few of Ms2ScanWithSpecificMass from the MS data file and adds them to the forward buffer.
        /// Also fills the back scan buffer w/ MsDataScans for retrieval for precursor info later.
        /// </summary>
        private void FillBuffers()
        {
            while (ForwardScanBuffer.Count < ForwardScanBufferSize)
            {
                // get the scan from the file
                MsDataScan scan = DataStream.GetOneBasedScanFromDynamicConnection(CurrentScan, filteringParams);

                if (scan == null)
                {
                    // done streaming data file
                    DataStream.CloseDynamicConnection();
                    IsDoneStreaming = true;
                    return;
                }

                // add the scan to the backscan buffer. this makes getting precursor scans easier/faster
                BackScanBuffer.Add(scan.OneBasedScanNumber, scan);

                // if this is a fragmentation scan, get + return precursor info
                if (scan.MsnOrder > 1 && scan.OneBasedPrecursorScanNumber.HasValue)
                {
                    // get the precursor scan
                    int precursorScanNum = scan.OneBasedPrecursorScanNumber.Value;
                    MsDataScan precursorScan;

                    if (!BackScanBuffer.TryGetValue(precursorScanNum, out precursorScan))
                    {
                        precursorScan = DataStream.GetOneBasedScanFromDynamicConnection(precursorScanNum, filteringParams);
                    }

                    // deconvolute precursor isolation window + this scan's fragments
                    foreach (Ms2ScanWithSpecificMass scanWithPrecursor in MetaMorpheusTask.DeconvolutePrecursors(scan, precursorScan, commonParameters, FilePath, precursorBuffer))
                    {
                        ForwardScanBuffer.Enqueue(scanWithPrecursor);
                    }
                }

                // remove old scans from buffer that are unlikely to be used later as precursor scans
                while (BackScanBuffer.Count > BackScanBufferSize)
                {
                    int minScanNumberInBuffer = BackScanBuffer.Min(p => p.Key);
                    BackScanBuffer.Remove(minScanNumberInBuffer);
                }

                CurrentScan++;
            }
        }
    }
}
