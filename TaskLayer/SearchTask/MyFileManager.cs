using IO.MzML;
using IO.Thermo;
using MassSpectrometry;
using System;
using System.Collections.Generic;
using System.IO;

namespace TaskLayer
{
    internal class MyFileManager
    {
        #region Private Fields

        private readonly bool disposeOfFileWhenDone;
        private Dictionary<string, IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>>> myMsDataFiles = new Dictionary<string, IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>>>();

        private object fileLoadingLock = new object();

        #endregion Private Fields

        #region Public Constructors

        public MyFileManager(bool disposeOfFileWhenDone)
        {
            this.disposeOfFileWhenDone = disposeOfFileWhenDone;
        }

        #endregion Public Constructors

        #region Internal Methods

        internal IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> LoadFile(string origDataFile)
        {
            if (myMsDataFiles.TryGetValue(origDataFile, out IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> value) && value != null)
                return value;

            // By now know that need to load this file!!!
            lock (fileLoadingLock) // Lock because reading is sequential
                if (Path.GetExtension(origDataFile).Equals(".mzML", StringComparison.InvariantCultureIgnoreCase))
                    myMsDataFiles[origDataFile] = Mzml.LoadAllStaticData(origDataFile);
                else
                    myMsDataFiles[origDataFile] = ThermoStaticData.LoadAllStaticData(origDataFile);

            return myMsDataFiles[origDataFile];
        }

        internal void DoneWithFile(string origDataFile)
        {
            if (disposeOfFileWhenDone)
                myMsDataFiles[origDataFile] = null;
        }

        #endregion Internal Methods
    }
}