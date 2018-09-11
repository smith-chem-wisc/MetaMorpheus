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
            else if(File.Exists(fileIoAssumedPath) || File.Exists(fregistryAssumedPath) || File.Exists(xRawFileAssumedPath))
            {
                return ThermoMsFileReaderVersionCheck.SomeDllsMissing;
            }

            return ThermoMsFileReaderVersionCheck.DllsNotFound;
        }

        public MsDataFile LoadFile(string origDataFile, int? topNpeaks, double? minRatio, bool trimMs1Peaks, bool trimMsMsPeaks, CommonParameters commonParameters)
        {
            FilteringParams filter = new FilteringParams(topNpeaks, minRatio, 1, trimMs1Peaks, trimMsMsPeaks);
            if (MyMsDataFiles.TryGetValue(origDataFile, out MsDataFile value) && value != null)
                return value;

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
                MyMsDataFiles[origDataFile] = null;
        }

        private void Warn(string v)
        {
            WarnHandler?.Invoke(this, new StringEventArgs(v, null));
        }
    }
}