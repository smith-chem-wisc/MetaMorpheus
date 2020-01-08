using EngineLayer;
using IO.Mgf;
using IO.MzML;
using MassSpectrometry;
using System;
using System.Collections.Generic;
using System.IO;
using System.IO.Compression;
using ThermoRawFileReader;

namespace TaskLayer
{
    public class MyFileManager
    {
        private readonly bool DisposeOfFileWhenDone;
        private readonly Dictionary<string, MsDataFile> MyMsDataFiles = new Dictionary<string, MsDataFile>();
        private readonly object FileLoadingLock = new object();

        public MyFileManager(bool disposeOfFileWhenDone)
        {
            DisposeOfFileWhenDone = disposeOfFileWhenDone;
        }

        public static event EventHandler<StringEventArgs> WarnHandler;

        public bool SeeIfOpen(string path)
        {
            return (MyMsDataFiles.ContainsKey(path) && MyMsDataFiles[path] != null);
        }

        public MsDataFile LoadFile(string origDataFile, CommonParameters commonParameters)
        {
            FilteringParams filter = new FilteringParams(commonParameters.NumberOfPeaksToKeepPerWindow, commonParameters.MinimumAllowedIntensityRatioToBasePeak, commonParameters.WindowWidthThomsons, commonParameters.NumberOfWindows, commonParameters.NormalizePeaksAccrossAllWindows, commonParameters.TrimMs1Peaks, commonParameters.TrimMsMsPeaks);

            if (commonParameters.DissociationType == DissociationType.LowCID || commonParameters.ChildScanDissociationType == DissociationType.LowCID)
            {
                filter = null;
            }

            if (MyMsDataFiles.TryGetValue(origDataFile, out MsDataFile value) && value != null)
            {
                return value;
            }

            // By now know that need to load this file!!!
            lock (FileLoadingLock) // Lock because reading is sequential
            {
                if (Path.GetExtension(origDataFile).Equals(".mzML", StringComparison.OrdinalIgnoreCase))
                {
                    MyMsDataFiles[origDataFile] = Mzml.LoadAllStaticData(origDataFile, filter, commonParameters.MaxThreadsToUsePerFile);
                }
                else if (Path.GetExtension(origDataFile).Equals(".mgf", StringComparison.OrdinalIgnoreCase))
                {
                    MyMsDataFiles[origDataFile] = Mgf.LoadAllStaticData(origDataFile, filter);
                }
                else
                {
                    MyMsDataFiles[origDataFile] = ThermoRawFileReaderData.LoadAllStaticData(origDataFile, filter, maxThreads: 1);
                }

                return MyMsDataFiles[origDataFile];
            }
        }

        internal void DoneWithFile(string origDataFile)
        {
            if (DisposeOfFileWhenDone)
            {
                MyMsDataFiles[origDataFile] = null;
            }
        }

        private void Warn(string v)
        {
            WarnHandler?.Invoke(this, new StringEventArgs(v, null));
        }

        public static void CompressDirectory(DirectoryInfo directorySelected)
        {
            foreach (FileInfo fileToCompress in directorySelected.GetFiles())
            {
                CompressFile(fileToCompress);
            }
        }

        public static void CompressFile(FileInfo fileToCompress)
        {
            using (FileStream originalFileStream = fileToCompress.OpenRead())
            {
                if ((File.GetAttributes(fileToCompress.FullName) &
                   FileAttributes.Hidden) != FileAttributes.Hidden & fileToCompress.Extension != ".gz")
                {
                    using (FileStream compressedFileStream = File.Create(fileToCompress.FullName + ".gz"))
                    {
                        using (GZipStream compressionStream = new GZipStream(compressedFileStream,
                           CompressionMode.Compress))
                        {
                            originalFileStream.CopyTo(compressionStream);
                        }
                    }
                }
            }

            if (File.Exists(fileToCompress.FullName))
            {
                File.Delete(fileToCompress.FullName);
            }
        }

        public static void DecompressDirectory(DirectoryInfo directorySelected)
        {
            foreach (FileInfo fileToDecompress in directorySelected.GetFiles())
            {
                DecompressFile(fileToDecompress);
            }
        }

        public static void DecompressFile(FileInfo fileToDecompress)
        {
            using (FileStream originalFileStream = fileToDecompress.OpenRead())
            {
                string currentFileName = fileToDecompress.FullName;
                string newFileName = currentFileName.Remove(currentFileName.Length - fileToDecompress.Extension.Length);

                using (FileStream decompressedFileStream = File.Create(newFileName))
                {
                    using (GZipStream decompressionStream = new GZipStream(originalFileStream, CompressionMode.Decompress))
                    {
                        decompressionStream.CopyTo(decompressedFileStream);
                        Console.WriteLine("Decompressed: {0}", fileToDecompress.Name);
                    }
                }
            }

            if (File.Exists(fileToDecompress.FullName))
            {
                File.Delete(fileToDecompress.FullName);
            }
        }
    }
}