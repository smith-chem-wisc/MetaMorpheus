using EngineLayer;
using IO.MzML;

#if NET461

using IO.Thermo;

#else
#endif

using MassSpectrometry;
using System;
using System.Collections.Generic;
using System.IO;

namespace TaskLayer
{
    public class MyFileManager
    {
        #region Private Fields

        private readonly bool disposeOfFileWhenDone;
        private readonly Dictionary<string, IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>>> myMsDataFiles = new Dictionary<string, IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>>>();
        private readonly object fileLoadingLock = new object();

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

        #region Internal Methods

        internal IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> LoadFile(string origDataFile, int? topNpeaks, double? minRatio, bool trimMs1Peaks, bool trimMsMsPeaks)
        {
            if (myMsDataFiles.TryGetValue(origDataFile, out IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> value) && value != null)
                return value;

            // By now know that need to load this file!!!
            lock (fileLoadingLock) // Lock because reading is sequential
                if (Path.GetExtension(origDataFile).Equals(".mzML", StringComparison.OrdinalIgnoreCase))
                    myMsDataFiles[origDataFile] = Mzml.LoadAllStaticData(origDataFile, topNpeaks, minRatio, trimMs1Peaks, trimMsMsPeaks);
                else

#if NET461
                    myMsDataFiles[origDataFile] = ThermoStaticData.LoadAllStaticData(origDataFile, topNpeaks, minRatio, trimMs1Peaks, trimMsMsPeaks);
#else
                    Warn("No capability for reading " + origDataFile);
#endif

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