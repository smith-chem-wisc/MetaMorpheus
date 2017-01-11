using InternalLogic;
using InternalLogicCalibration;
using IO.MzML;
using IO.Thermo;
using MassSpectrometry;
using Spectra;
using System;
using System.Collections.Generic;
using System.IO;
using System.Text.RegularExpressions;

namespace InternalLogicWithFileIO
{
    public static class mzCalIO
    {
        public static string elementsLocation = @"elements.dat";

        private static int GetLastNumberFromString(string s)
        {
            return Convert.ToInt32(Regex.Match(s, @"\d+$").Value);
        }

        public static void Load()
        {
            UsefulProteomicsDatabases.Loaders.LoadElements(elementsLocation);
        }

        //public static void MzmlOutput(SoftwareLockMassParams p)
        //{
        //    p.po.status("Creating _indexedmzMLConnection, and putting data in it");
        //    var path = Path.Combine(Path.GetDirectoryName(p.myMsDataFile.FilePath), Path.GetFileNameWithoutExtension(p.myMsDataFile.FilePath) + p.paramString + "-Calibrated.mzML");
        //    MzmlMethods.CreateAndWriteMyIndexedMZmlwithCalibratedSpectra(p.myMsDataFile, path);
        //    //p.po.SucessfullyFinishedWritingFile(new SingleFileEventArgs(path));
        //    //p.myTaskResults.newSpectra.Add(path);
        //}
    }
}