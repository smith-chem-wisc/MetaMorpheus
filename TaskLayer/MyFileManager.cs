using EngineLayer;
using IO.MzML;
using IO.Mgf;

#if NETFRAMEWORK

using IO.Thermo;

#else
#endif

using MassSpectrometry;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.IO.Compression;

namespace TaskLayer
{
    public class MyFileManager
    {
        public enum ThermoMsFileReaderVersionCheck { DllsNotFound, IncorrectVersion, CorrectVersion, SomeDllsMissing };

        private readonly bool DisposeOfFileWhenDone;
        private readonly Dictionary<string, MsDataFile> MyMsDataFiles = new Dictionary<string, MsDataFile>();
        private readonly object FileLoadingLock = new object();
        private const string AssumedThermoMsFileReaderDllPath = @"C:\Program Files\Thermo\MSFileReader";
        private const string DesiredFileIoVersion = "3.0";
        private const string DesiredFregistryVersion = "3.0";
        private const string DesiredXRawFileVersion = "3.0.29.0";

        public MyFileManager(bool disposeOfFileWhenDone)
        {
            DisposeOfFileWhenDone = disposeOfFileWhenDone;
        }

        public static event EventHandler<StringEventArgs> WarnHandler;

        public bool SeeIfOpen(string path)
        {
            return (MyMsDataFiles.ContainsKey(path) && MyMsDataFiles[path] != null);
        }

        public static ThermoMsFileReaderVersionCheck ValidateThermoMsFileReaderVersion()
        {
            string fileIoAssumedPath = Path.Combine(AssumedThermoMsFileReaderDllPath, "Fileio_x64.dll");
            string fregistryAssumedPath = Path.Combine(AssumedThermoMsFileReaderDllPath, "fregistry_x64.dll");
            string xRawFileAssumedPath = Path.Combine(AssumedThermoMsFileReaderDllPath, "XRawfile2_x64.dll");

            if (File.Exists(fileIoAssumedPath) && File.Exists(fregistryAssumedPath) && File.Exists(xRawFileAssumedPath))
            {
                string fileIoVersion = FileVersionInfo.GetVersionInfo(fileIoAssumedPath).FileVersion;
                string fregistryVersion = FileVersionInfo.GetVersionInfo(fregistryAssumedPath).FileVersion;
                string xRawFileVersion = FileVersionInfo.GetVersionInfo(xRawFileAssumedPath).FileVersion;

                if (fileIoVersion.Equals(DesiredFileIoVersion) && fregistryVersion.Equals(DesiredFregistryVersion) && xRawFileVersion.Equals(DesiredXRawFileVersion))
                {
                    return ThermoMsFileReaderVersionCheck.CorrectVersion;
                }
                else
                {
                    return ThermoMsFileReaderVersionCheck.IncorrectVersion;
                }
            }
            else if (File.Exists(fileIoAssumedPath) || File.Exists(fregistryAssumedPath) || File.Exists(xRawFileAssumedPath))
            {
                return ThermoMsFileReaderVersionCheck.SomeDllsMissing;
            }

            return ThermoMsFileReaderVersionCheck.DllsNotFound;
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
#if NETFRAMEWORK
                    MyMsDataFiles[origDataFile] = ThermoStaticData.LoadAllStaticData(origDataFile, filter);
#else
                    Warn("No capability for reading " + origDataFile);
#endif
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