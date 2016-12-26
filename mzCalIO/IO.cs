using IO.MzML;
using IO.Thermo;
using MassSpectrometry;
using mzCal;
using MzIdentML;
using Proteomics;
using Spectra;
using System;
using System.Collections.Generic;
using System.IO;
using System.Text.RegularExpressions;

namespace mzCalIO
{
    public static class mzCalIO
    {
        public static UsefulProteomicsDatabases.Generated.unimod unimodDeserialized;
        public static UsefulProteomicsDatabases.Generated.obo psimodDeserialized;
        public static Dictionary<int, ChemicalFormulaModification> uniprotDeseralized;

        public static string unimodLocation = @"unimod_tables.xml";
        public static string psimodLocation = @"PSI-MOD.obo.xml";
        public static string elementsLocation = @"elements.dat";
        public static string uniprotLocation = @"ptmlist.txt";

        private static int GetLastNumberFromString(string s)
        {
            return Convert.ToInt32(Regex.Match(s, @"\d+$").Value);
        }

        public static void Load()
        {
            UsefulProteomicsDatabases.Loaders.LoadElements(elementsLocation);
            unimodDeserialized = UsefulProteomicsDatabases.Loaders.LoadUnimod(unimodLocation);
            psimodDeserialized = UsefulProteomicsDatabases.Loaders.LoadPsiMod(psimodLocation);
            uniprotDeseralized = UsefulProteomicsDatabases.Loaders.LoadUniprot(uniprotLocation);
        }

        public static string getFormulaFromDictionary(string dictionary, string acession)
        {
            if (dictionary == "UNIMOD")
            {
                string unimodAcession = acession;
                var indexToLookFor = GetLastNumberFromString(unimodAcession) - 1;
                while (unimodDeserialized.modifications[indexToLookFor].record_id != GetLastNumberFromString(unimodAcession))
                    indexToLookFor--;
                return Regex.Replace(unimodDeserialized.modifications[indexToLookFor].composition, @"[\s()]", ""); ;
            }
            else if (dictionary == "PSI-MOD")
            {
                string psimodAcession = acession;
                UsefulProteomicsDatabases.Generated.oboTerm ksadklfj = (UsefulProteomicsDatabases.Generated.oboTerm)psimodDeserialized.Items[GetLastNumberFromString(psimodAcession) + 2];

                if (GetLastNumberFromString(psimodAcession) != GetLastNumberFromString(ksadklfj.id))
                    throw new Exception("Error in reading psi-mod file, acession mismatch!");
                else
                {
                    foreach (var a in ksadklfj.xref_analog)
                    {
                        if (a.dbname == "DiffFormula")
                        {
                            return Regex.Replace(a.name, @"[\s()]", "");
                        }
                    }
                    return uniprotDeseralized[GetLastNumberFromString(psimodAcession)].ThisChemicalFormula.Formula;
                }
            }
            else if (dictionary == "v.tsv" && acession == "v.tsv:35")
                return "O";
            else if (dictionary == "f.tsv" && acession == "f.tsv:4")
                return "CH3NO";
            else
                //Not familiar with modification dictionary
                return null;
        }

        public static void MzmlOutput(SoftwareLockMassParams p)
        {
            p.OnOutput("Creating _indexedmzMLConnection, and putting data in it");
            var path = Path.Combine(Path.GetDirectoryName(p.myMsDataFile.FilePath), Path.GetFileNameWithoutExtension(p.myMsDataFile.FilePath) + p.paramString + "-Calibrated.mzML");
            MzmlMethods.CreateAndWriteMyIndexedMZmlwithCalibratedSpectra(p.myMsDataFile, path);
            p.OnFileFinished(path);
        }

        public static SoftwareLockMassParams GetReady(string origDataFile, string mzidName)
        {
            IMsDataFile<IMzSpectrum<MzPeak>> myMsDataFile;
            bool deconvolute = true;
            if (Path.GetExtension(origDataFile).Equals(".mzML"))
            {
                myMsDataFile = new Mzml(origDataFile);
                myMsDataFile.Open();
            }
            else
            {
                myMsDataFile = new ThermoRawFile(origDataFile);
                myMsDataFile.Open();
                if (((ThermoRawFile)myMsDataFile).monoisotopicPrecursorSelectionEnabled)
                    deconvolute = false;
            }
            int randomSeed = 1;
            var identifications = new MzidIdentifications(mzidName);
            var a = new SoftwareLockMassParams(myMsDataFile, randomSeed, deconvolute, identifications.fragmentTolerance.Value * 2);

            a.postProcessing = MzmlOutput;
            a.getFormulaFromDictionary = getFormulaFromDictionary;
            a.identifications = identifications;
            a.mzRange = new DoubleRange(0, 0);

            //a.MS1spectraToWatch.Add(22557);

            //a.MS2spectraToWatch.Add(22564);

            a.matchesToExclude = new HashSet<int>();
            return a;
        }
    }
}