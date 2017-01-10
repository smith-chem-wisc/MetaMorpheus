using IO.MzML;
using IO.Thermo;
using MassSpectrometry;
using Spectra;
using System;
using System.Collections.Generic;
using System.IO;
using System.Text.RegularExpressions;

namespace MetaMorpheusLogic
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

        public static void MzmlOutput(SoftwareLockMassParams p)
        {
            p.po.status("Creating _indexedmzMLConnection, and putting data in it");
            var path = Path.Combine(Path.GetDirectoryName(p.myMsDataFile.FilePath), Path.GetFileNameWithoutExtension(p.myMsDataFile.FilePath) + p.paramString + "-Calibrated.mzML");
            MzmlMethods.CreateAndWriteMyIndexedMZmlwithCalibratedSpectra(p.myMsDataFile, path);
            p.po.SucessfullyFinishedWritingFile(new SingleFileEventArgs(path));
            p.myTaskResults.newSpectra.Add(path);
        }

        public static SoftwareLockMassParams GetReady(string origDataFile, List<NewPsmWithFDR> psms, Tolerance searchfragmentTolerance, AllTasksParams po, MyTaskResults myTaskResults)
        {
            IMsDataFile<IMzSpectrum<MzPeak>> myMsDataFile;
            if (Path.GetExtension(origDataFile).Equals(".mzML"))
            {
                myMsDataFile = new Mzml(origDataFile);
                myMsDataFile.Open();
            }
            else
            {
                myMsDataFile = new ThermoRawFile(origDataFile);
                myMsDataFile.Open();
            }
            int randomSeed = 1;
            // TODO: fix the tolerance calculation below
            var a = new SoftwareLockMassParams(myMsDataFile, randomSeed, searchfragmentTolerance.Value * 2, myTaskResults);

            a.postProcessing = MzmlOutput;
            a.identifications = psms;
            a.mzRange = new DoubleRange(0, 0);

            //a.MS1spectraToWatch.Add(22557);

            //a.MS2spectraToWatch.Add(22564);

            a.matchesToExclude = new HashSet<int>();

            a.po = po;
            return a;
        }
    }
}