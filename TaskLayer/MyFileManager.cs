using EngineLayer;
using IO.MzML;

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
        #region Public Fields

        public enum ThermoMsFileReaderVersionCheck { DllsNotFound, IncorrectVersion, CorrectVersion };

        #endregion Public Fields

        #region Private Fields

        private readonly bool disposeOfFileWhenDone;
        private readonly Dictionary<string, IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>>> myMsDataFiles = new Dictionary<string, IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>>>();
        private readonly object fileLoadingLock = new object();
        private const string AssumedThermoMsFileReaderDllPath = @"C:\Program Files\Thermo\MSFileReader";
        private const string DesiredFileIoVersion = "3.0";
        private const string DesiredFregistryVersion = "3.0";
        private const string DesiredXRawFileVersion = "3.0.29.0";

        #endregion Private Fields

        #region Public Constructors

        public MyFileManager(bool disposeOfFileWhenDone)
        {
            this.disposeOfFileWhenDone = disposeOfFileWhenDone;
        }

        #endregion Public Constructors

        #region Public Events

        public static event EventHandler<StringEventArgs> WarnHandler;

        #endregion Public Events

        #region Public Methods

        public bool SeeIfOpen(string path)
        {
            return (myMsDataFiles.ContainsKey(path) && myMsDataFiles[path] != null);
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

            return ThermoMsFileReaderVersionCheck.DllsNotFound;
        }

        #endregion Public Methods

        #region Internal Methods

        internal IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> LoadFile(string origDataFile, int? topNpeaks, double? minRatio, bool trimMs1Peaks, bool trimMsMsPeaks)
        {
            FilteringParams filter = new FilteringParams(topNpeaks, minRatio, null, trimMs1Peaks, trimMsMsPeaks);
            if (myMsDataFiles.TryGetValue(origDataFile, out IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> value) && value != null)
                return value;

            // By now know that need to load this file!!!
            lock (fileLoadingLock) // Lock because reading is sequential
                if (Path.GetExtension(origDataFile).Equals(".mzML", StringComparison.OrdinalIgnoreCase))
                {
                    myMsDataFiles[origDataFile] = Mzml.LoadAllStaticData(origDataFile, filter);
                }
                else
                {
#if NETFRAMEWORK
                    myMsDataFiles[origDataFile] = ThermoStaticData.LoadAllStaticData(origDataFile, filter);
#else
                    Warn("No capability for reading " + origDataFile);
#endif
                }
            return myMsDataFiles[origDataFile];
        }

        internal void DoneWithFile(string origDataFile)
        {
            if (disposeOfFileWhenDone)
                myMsDataFiles[origDataFile] = null;
        }

        #endregion Internal Methods

        #region Private Methods

        private void Warn(string v)
        {
            WarnHandler?.Invoke(this, new StringEventArgs(v, null));
        }

        #endregion Private Methods
    }
}