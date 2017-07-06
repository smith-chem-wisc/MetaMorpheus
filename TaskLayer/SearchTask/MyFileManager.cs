using IO.MzML;
using IO.Thermo;
using MassSpectrometry;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace TaskLayer
{
    internal class MyFileManager
    {

        #region Private Fields

        private Dictionary<string, IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>>> myMsDataFiles = new Dictionary<string, IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>>>();
        private Dictionary<string, bool> inUse = new Dictionary<string, bool>();

        private object fileLoadingLock = new object();

        #endregion Private Fields

        #region Public Constructors

        public MyFileManager()
        {
        }

        #endregion Public Constructors

        #region Internal Methods

        internal IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> LoadFile(string origDataFile)
        {
            IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> value;
            if (myMsDataFiles.TryGetValue(origDataFile, out value) && value != null)
                return value;

            // By now know that need to load this file!!!

            var success = false;
            lock (fileLoadingLock) // Lock because reading is sequential
                while (!success)
                {
                    try
                    {
                        if (Path.GetExtension(origDataFile).Equals(".mzML"))
                            myMsDataFiles[origDataFile] = Mzml.LoadAllStaticData(origDataFile);
                        else
                            myMsDataFiles[origDataFile] = ThermoStaticData.LoadAllStaticData(origDataFile);
                        success = true;
                    }
                    catch (OutOfMemoryException)
                    {
                        var notInUse = inUse.First(b => !b.Value);
                        myMsDataFiles[notInUse.Key] = null;
                    }
                }

            inUse[origDataFile] = true;
            return myMsDataFiles[origDataFile];
        }

        internal void DoneWithFile(string origDataFile)
        {
            inUse[origDataFile] = false;
        }

        #endregion Internal Methods

    }
}