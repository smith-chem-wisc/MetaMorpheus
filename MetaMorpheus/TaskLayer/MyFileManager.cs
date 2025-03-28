using EngineLayer;
using MassSpectrometry;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.IO;
using System.IO.Compression;
using Readers;
using System.Reflection.Metadata.Ecma335;

namespace TaskLayer
{
    public class MyFileManager
    {
        private readonly bool DisposeOfFileWhenDone;
        private readonly ConcurrentDictionary<string, MsDataFile> MyMsDataFiles = new ConcurrentDictionary<string, MsDataFile>();

        public MyFileManager(bool disposeOfFileWhenDone)
        {
            DisposeOfFileWhenDone = disposeOfFileWhenDone;
        }

        public static event EventHandler<StringEventArgs> WarnHandler;

        public bool SeeIfOpen(string path)
        {
            return MyMsDataFiles.TryGetValue(path, out var file);
        }

        public MsDataFile LoadFile(string origDataFile, CommonParameters commonParameters)
        {
            FilteringParams filter = new FilteringParams(
                commonParameters.NumberOfPeaksToKeepPerWindow,
                commonParameters.MinimumAllowedIntensityRatioToBasePeak,
                commonParameters.WindowWidthThomsons,
                commonParameters.NumberOfWindows,
                commonParameters.NormalizePeaksAccrossAllWindows,
                commonParameters.TrimMs1Peaks,
                commonParameters.TrimMsMsPeaks);

            if (commonParameters.DissociationType == DissociationType.LowCID || commonParameters.MS2ChildScanDissociationType == DissociationType.LowCID || commonParameters.MS3ChildScanDissociationType == DissociationType.LowCID)
            {
                filter = null;
            }

            if (MyMsDataFiles.TryGetValue(origDataFile, out MsDataFile value))
                return value;

            return MyMsDataFiles.GetOrAdd(origDataFile,
                MsDataFileReader.GetDataFile(origDataFile).LoadAllStaticData(filter, commonParameters.MaxThreadsToUsePerFile));
        }

        internal void DoneWithFile(string origDataFile)
        {
            if (DisposeOfFileWhenDone)
            {
                // This would only return false if the file was not in the dictionary
                MyMsDataFiles.TryRemove(origDataFile, out var file);
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