using EngineLayer;
using MassSpectrometry;
using NUnit.Framework;
using Proteomics;
using Proteomics.AminoAcidPolymer;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using TaskLayer;
using UsefulProteomicsDatabases;

namespace Test
{
    [TestFixture]
    public static class SilacTest
    {
        [Test]
        public static void TestSilacNoLightProtein()
        {
            //The concern with multiple mods per label is the conversions back and forth between "light" and "heavy" labels
            Residue heavyArginine = new Residue("c", 'c', "c", Chemistry.ChemicalFormula.ParseFormula("C6H12N{15}4O"), ModificationSites.All); //+4 arginine
            Residue heavierArginine = new Residue("d", 'd', "d", Chemistry.ChemicalFormula.ParseFormula("C{13}6H12N{15}4O"), ModificationSites.All); //+10 arginine
            Residue.AddNewResiduesToDictionary(new List<Residue> { heavyArginine }); //These should be added in the  search task, but we need to add this one earlier so that we can create a heavy pwsm

            Residue lightArginine = Residue.GetResidue('R');

            SilacLabel heavyLabel = new SilacLabel(lightArginine.Letter, heavyArginine.Letter, heavyArginine.ThisChemicalFormula.Formula, heavyArginine.MonoisotopicMass - lightArginine.MonoisotopicMass);
            SilacLabel heavierLabel = new SilacLabel(lightArginine.Letter, heavierArginine.Letter, heavierArginine.ThisChemicalFormula.Formula, heavierArginine.MonoisotopicMass - lightArginine.MonoisotopicMass);

            SearchTask task = new SearchTask
            {
                SearchParameters = new SearchParameters
                {
                    SilacLabels = new List<SilacLabel> { heavyLabel, heavierLabel }
                },
                CommonParameters = new CommonParameters(digestionParams: new DigestionParams(generateUnlabeledProteinsForSilac: false)) //this is the important part of the unit test
            };

            PeptideWithSetModifications heavyPeptide = new PeptideWithSetModifications("PEPTIDEc", new Dictionary<string, Modification>());

            List<double> massDifferences = new List<double> { heavierArginine.MonoisotopicMass - heavyArginine.MonoisotopicMass };
            MsDataFile myMsDataFile1 = new TestDataFile(heavyPeptide, massDifferences);
            string mzmlName = @"silac.mzML";
            IO.MzML.MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(myMsDataFile1, mzmlName, false);

            string xmlName = "SilacDb.xml";
            Protein theProtein = new Protein("PEPTIDER", "accession1");
            ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), new List<Protein> { theProtein }, xmlName);

            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestSilac");
            Directory.CreateDirectory(outputFolder);
            var theStringResult = task.RunTask(outputFolder, new List<DbForTask> { new DbForTask(xmlName, false) }, new List<string> { mzmlName }, "taskId1").ToString();

            //test proteins
            string[] output = File.ReadAllLines(TestContext.CurrentContext.TestDirectory + @"/TestSilac/AllProteinGroups.tsv");
            Assert.AreEqual(output.Length, 2);
            Assert.IsTrue(output[0].Contains("Modification Info List\tIntensity_silac(R+3.988)\tIntensity_silac(R+10.008)")); //test that two files were made and no light file
            Assert.IsTrue(output[1].Contains("875000\t437500")); //test the heavier intensity is half that of the heavy (per the raw file)

            //test peptides
            output = File.ReadAllLines(TestContext.CurrentContext.TestDirectory + @"/TestSilac/AllQuantifiedPeptides.tsv");
            Assert.AreEqual(output.Length, 2);
            Assert.IsTrue(output[0].Contains("Organism\tIntensity_silac(R+3.988)\tIntensity_silac(R+10.008)")); //test the two files were made and no light file
            Assert.IsTrue(output[1].Contains("875000\t437500")); //test intensity

            //test peaks
            output = File.ReadAllLines(TestContext.CurrentContext.TestDirectory + @"/TestSilac/AllQuantifiedPeaks.tsv");
            Assert.AreEqual(output.Length, 3);
            Assert.IsTrue(output[1].Contains("silac\t")); //test the filename was NOT modified (it was for proteins, but we don't want it for peptides)
            Assert.IsTrue(output[2].Contains("silac\t"));//test the filename was NOT modified (it was for proteins, but we don't want it for peptides)
            Assert.IsTrue(output[1].Contains("PEPTIDER(+3.988)\t")); //test light sequence was not modified
            Assert.IsTrue(output[2].Contains("PEPTIDER(+10.008)\t")); //test heavy sequence was output correctly (do NOT want "PEPTIDEa")
            Assert.IsTrue(output[1].Contains("959.44")); //test light mass
            Assert.IsTrue(output[2].Contains("965.46")); //test heavy mass

            //test PSMs
            output = File.ReadAllLines(TestContext.CurrentContext.TestDirectory + @"/TestSilac/AllPSMs.psmtsv");
            Assert.IsTrue(output[1].Contains("959.44")); //test the correct monoisotopic mass
            Assert.IsTrue(output[1].Contains("PEPTIDER(+3.988)")); //test the correct psm
            Assert.IsTrue(output[1].Contains("silac\t")); //test the filename was NOT modified (it was for proteins, but we don't want it for peptides)

            //Clear the old files
            Directory.Delete(outputFolder, true);
            File.Delete(xmlName);
            File.Delete(mzmlName);
        }

        [Test]
        public static void TestSilacMultipleModsPerCondition()
        {
            //The concern with multiple mods per label is the conversions back and forth between "light" and "heavy" labels
            Residue heavyLysine = new Residue("a", 'a', "a", Chemistry.ChemicalFormula.ParseFormula("C{13}6H12N{15}2O"), ModificationSites.All); //+8 lysine
            Residue heavyArginine = new Residue("b", 'b', "b", Chemistry.ChemicalFormula.ParseFormula("C{13}6H12N4O"), ModificationSites.All); //+6 arginine
            Residue lightLysine = Residue.GetResidue('K');
            Residue lightArginine = Residue.GetResidue('R');

            SilacLabel krLabel = new SilacLabel(lightLysine.Letter, heavyLysine.Letter, heavyLysine.ThisChemicalFormula.Formula, heavyLysine.MonoisotopicMass - lightLysine.MonoisotopicMass);
            krLabel.AddAdditionalSilacLabel(new SilacLabel(lightArginine.Letter, heavyArginine.Letter, heavyArginine.ThisChemicalFormula.Formula, heavyArginine.MonoisotopicMass - lightArginine.MonoisotopicMass));

            SearchTask task = new SearchTask
            {
                SearchParameters = new SearchParameters
                {
                    SilacLabels = new List<SilacLabel> { krLabel }
                },
                CommonParameters = new CommonParameters(digestionParams: new DigestionParams(minPeptideLength: 2))
            };

            PeptideWithSetModifications lightPeptide = new PeptideWithSetModifications("SEQENEWITHAKANDANR", new Dictionary<string, Modification>());

            List<double> massDifferences = new List<double> { (heavyLysine.MonoisotopicMass + heavyArginine.MonoisotopicMass) - (lightLysine.MonoisotopicMass + lightArginine.MonoisotopicMass) };
            MsDataFile myMsDataFile1 = new TestDataFile(lightPeptide, massDifferences);
            string mzmlName = @"silac.mzML";
            IO.MzML.MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(myMsDataFile1, mzmlName, false);

            string xmlName = "SilacDb.xml";
            Protein theProtein = new Protein("MPRTEINRSEQENEWITHAKANDANRANDSMSTFF", "accession1");
            ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), new List<Protein> { theProtein }, xmlName);

            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestSilac");
            Directory.CreateDirectory(outputFolder);
            task.RunTask(outputFolder, new List<DbForTask> { new DbForTask(xmlName, false) }, new List<string> { mzmlName }, "taskId1");

            //test proteins
            string[] output = File.ReadAllLines(TestContext.CurrentContext.TestDirectory + @"/TestSilac/AllProteinGroups.tsv");
            Assert.AreEqual(output.Length, 2);
            Assert.IsTrue(output[0].Contains("Intensity_silac\tIntensity_silac(K+8.014 & R+6.020)")); //test that two files were made
            Assert.IsTrue(output[1].Contains("875000\t437500")); //test the heavy intensity is half that of the light (per the raw file)

            //test peptides
            output = File.ReadAllLines(TestContext.CurrentContext.TestDirectory + @"/TestSilac/AllQuantifiedPeptides.tsv");
            Assert.AreEqual(output.Length, 2);
            Assert.IsTrue(output[1].Contains("SEQENEWITHAKANDANR\taccession1\t"));//test the sequence and accession were not modified
            Assert.IsTrue(output[1].Contains("875000")); //test intensity
            Assert.IsFalse(output[1].Contains("SEQENEWITHAK(+8.014)ANDANR(+6.020)")); //test the sequence was not doubled modified
            Assert.IsTrue(output[1].Contains("437500")); //test intensity

            //test peaks
            output = File.ReadAllLines(TestContext.CurrentContext.TestDirectory + @"/TestSilac/AllQuantifiedPeaks.tsv");
            Assert.AreEqual(output.Length, 3);
            Assert.IsTrue(output[1].Contains("silac\t")); //test the filename was NOT modified (it was for proteins, but we don't want it for peptides)
            Assert.IsTrue(output[2].Contains("silac\t"));//test the filename was NOT modified (it was for proteins, but we don't want it for peptides)
            Assert.IsTrue(output[1].Contains("SEQENEWITHAKANDANR\t")); //test light sequence was not modified
            Assert.IsTrue(output[2].Contains("SEQENEWITHAK(+8.014)ANDANR(+6.020)\t")); //test heavy sequence was output correctly (do NOT want "PEPTIDEa")
            Assert.IsTrue(output[1].Contains("2111.96")); //test light mass
            Assert.IsTrue(output[2].Contains("2125.99")); //test heavy mass
            Assert.IsTrue(output[2].Contains("accession1")); //test heavy accesssion is light in output


            ///Test for when an additional label is the only label on a peptide
            ///Usually crashes in mzId
            //Delete old files
            File.Delete(mzmlName);
            Directory.Delete(outputFolder, true);

            lightPeptide = new PeptideWithSetModifications("ANDANR", new Dictionary<string, Modification>()); //has the additional, but not the original
            massDifferences = new List<double> { (heavyArginine.MonoisotopicMass) - (lightArginine.MonoisotopicMass) };
            myMsDataFile1 = new TestDataFile(lightPeptide, massDifferences, true);
            mzmlName = @"silac.mzML";
            IO.MzML.MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(myMsDataFile1, mzmlName, false);

            Directory.CreateDirectory(outputFolder);
            task.RunTask(outputFolder, new List<DbForTask> { new DbForTask(xmlName, false) }, new List<string> { mzmlName }, "taskId1");

            //Clear the old files
            Directory.Delete(outputFolder, true);
            File.Delete(xmlName);
            File.Delete(mzmlName);


        }

        [Test]
        public static void TestSilacQuantification()
        {
            //make heavy residue and add to search task
            Residue heavyLysine = new Residue("a", 'a', "a", Chemistry.ChemicalFormula.ParseFormula("C{13}6H12N{15}2O"), ModificationSites.All); //+8 lysine
            Residue lightLysine = Residue.GetResidue('K');

            SearchTask task = new SearchTask
            {
                SearchParameters = new SearchParameters
                {
                    SilacLabels = new List<SilacLabel> { new SilacLabel(lightLysine.Letter, heavyLysine.Letter, heavyLysine.ThisChemicalFormula.Formula, heavyLysine.MonoisotopicMass - lightLysine.MonoisotopicMass) }
                }
            };

            PeptideWithSetModifications lightPeptide = new PeptideWithSetModifications("PEPTIDEK", new Dictionary<string, Modification>());

            List<double> massDifferences = new List<double> { heavyLysine.MonoisotopicMass - lightLysine.MonoisotopicMass };
            MsDataFile myMsDataFile1 = new TestDataFile(lightPeptide, massDifferences);
            string mzmlName = @"silac.mzML";
            IO.MzML.MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(myMsDataFile1, mzmlName, false);

            string xmlName = "SilacDb.xml";
            Protein theProtein = new Protein("PEPTIDEK", "accession1");
            ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), new List<Protein> { theProtein }, xmlName);

            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestSilac");
            Directory.CreateDirectory(outputFolder);
            var theStringResult = task.RunTask(outputFolder, new List<DbForTask> { new DbForTask(xmlName, false) }, new List<string> { mzmlName }, "taskId1").ToString();

            Assert.IsTrue(theStringResult.Contains("All target PSMS within 1% FDR: 1")); //it's not a psm, it's a MBR feature

            ///Normal Peptide
            //test proteins
            string[] output = File.ReadAllLines(TestContext.CurrentContext.TestDirectory + @"/TestSilac/AllProteinGroups.tsv");
            Assert.AreEqual(output.Length, 2);
            Assert.IsTrue(output[0].Contains("Intensity_silac\tIntensity_silac(K+8.014)")); //test that two files were made
            Assert.IsTrue(output[1].Contains("875000\t437500")); //test the heavy intensity is half that of the light (per the raw file)

            //test peptides
            output = File.ReadAllLines(TestContext.CurrentContext.TestDirectory + @"/TestSilac/AllQuantifiedPeptides.tsv");
            Assert.AreEqual(output.Length, 2);
            Assert.IsTrue(output[1].Contains("PEPTIDEK\taccession1\t"));//test the sequence and accession were not modified
            Assert.IsTrue(output[1].Contains("875000")); //test intensity
            Assert.IsFalse(output[1].Contains("PEPTIDEK(+8.014)")); //test the sequence was not doubled modified
            Assert.IsTrue(output[1].Contains("437500")); //test intensity

            //test peaks
            output = File.ReadAllLines(TestContext.CurrentContext.TestDirectory + @"/TestSilac/AllQuantifiedPeaks.tsv");
            Assert.AreEqual(output.Length, 3);
            Assert.IsTrue(output[1].Contains("silac\t")); //test the filename was NOT modified (it was for proteins, but we don't want it for peptides)
            Assert.IsTrue(output[2].Contains("silac\t"));//test the filename was NOT modified (it was for proteins, but we don't want it for peptides)
            Assert.IsTrue(output[1].Contains("PEPTIDEK\t")); //test light sequence was not modified
            Assert.IsTrue(output[2].Contains("PEPTIDEK(+8.014)\t")); //test heavy sequence was output correctly (do NOT want "PEPTIDEa")
            Assert.IsTrue(output[1].Contains("927.45")); //test light mass
            Assert.IsTrue(output[2].Contains("935.46")); //test heavy mass

            ///Ambiguous base sequence peptide
            //Clear the old files
            Directory.Delete(outputFolder, true);
            File.Delete(xmlName);
            File.Delete(mzmlName);

            //make a heavy peptide
            massDifferences = new List<double> { heavyLysine.MonoisotopicMass - lightLysine.MonoisotopicMass };
            myMsDataFile1 = new TestDataFile(lightPeptide, massDifferences, true);
            IO.MzML.MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(myMsDataFile1, mzmlName, false);

            //make an ambiguous database
            Protein theProtein2 = new Protein("PEPTLDEKPEPTIDEK", "accession2");
            ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), new List<Protein> { theProtein, theProtein2 }, xmlName);

            Directory.CreateDirectory(outputFolder);
            theStringResult = task.RunTask(outputFolder, new List<DbForTask> { new DbForTask(xmlName, false) }, new List<string> { mzmlName }, "taskId1").ToString();

            output = File.ReadAllLines(TestContext.CurrentContext.TestDirectory + @"/TestSilac/AllPSMs.psmtsv");
            Assert.IsTrue(output[1].Contains("silac\t")); //test the filename was NOT modified (it was for proteins, but we don't want it for peptides)
            Assert.IsTrue(output[1].Contains("PEPTIDEK(+8.014)|PEPTLDEK(+8.014)|PEPTIDEK(+8.014)")
                || output[1].Contains("PEPTIDEK(+8.014)|PEPTIDEK(+8.014)|PEPTLDEK(+8.014)")
                || output[1].Contains("PEPTLDEK(+8.014)|PEPTIDEK(+8.014)|PEPTIDEK(+8.014)")); //test the heavy ambiguous peptides were all found
            //Need the options, because output isn't consistent as of 3/26/19

            ///Ambiguous proteinGroup
            //Clear the old files
            Directory.Delete(outputFolder, true);
            File.Delete(xmlName);

            //make an ambiguous database
            theProtein2 = new Protein("PEPTIDEK", "accession2");
            ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), new List<Protein> { theProtein, theProtein2 }, xmlName);

            Directory.CreateDirectory(outputFolder);
            theStringResult = task.RunTask(outputFolder, new List<DbForTask> { new DbForTask(xmlName, false) }, new List<string> { mzmlName }, "taskId1").ToString();

            output = File.ReadAllLines(TestContext.CurrentContext.TestDirectory + @"/TestSilac/AllPSMs.psmtsv");
            Assert.IsTrue(output[1].Contains("accession1|accession2")
                || output[1].Contains("accession2|accession1")); //test the heavy ambiguous peptides were all found
            //Need the options, because output isn't consistent as of 3/26/19
            Assert.IsTrue(output[1].Contains("\tPEPTIDEK(+8.014)\t")); //test the heavy ambiguous peptides were all found

            //delete files
            Directory.Delete(outputFolder, true);
            File.Delete(xmlName);
            File.Delete(mzmlName);
        }

        [Test]
        public static void TestSilacWhenProteinIsMissing()
        {
            //make heavy residue and add to search task
            Residue heavyLysine = new Residue("a", 'a', "a", Chemistry.ChemicalFormula.ParseFormula("C{13}6H12N{15}2O"), ModificationSites.All); //+8 lysine
            Residue lightLysine = Residue.GetResidue('K');

            SearchTask task = new SearchTask
            {
                SearchParameters = new SearchParameters
                {
                    SilacLabels = new List<SilacLabel> { new SilacLabel(lightLysine.Letter, heavyLysine.Letter, heavyLysine.ThisChemicalFormula.Formula, heavyLysine.MonoisotopicMass - lightLysine.MonoisotopicMass) },
                    NoOneHitWonders = true
                    //The NoOneHitWonders=true doesn't really seem like a SILAC test, but we're testing that there's no crash if a quantified peptide's proteinGroup isn't quantified
                    //This happens if somebody messed with parsimony (picked TDS) or from requiring two peptides per protein (and we're only finding one). We're testing the second case here.
                }
            };

            PeptideWithSetModifications lightPeptide = new PeptideWithSetModifications("PEPTIDEK", new Dictionary<string, Modification>());

            List<double> massDifferences = new List<double> { heavyLysine.MonoisotopicMass - lightLysine.MonoisotopicMass };
            MsDataFile myMsDataFile1 = new TestDataFile(lightPeptide, massDifferences);
            string mzmlName = @"silac.mzML";
            IO.MzML.MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(myMsDataFile1, mzmlName, false);

            string xmlName = "SilacDb.xml";
            Protein theProtein = new Protein("PEPTIDEK", "accession1");
            ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), new List<Protein> { theProtein }, xmlName);

            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestSilac");
            Directory.CreateDirectory(outputFolder);
            var theStringResult = task.RunTask(outputFolder, new List<DbForTask> { new DbForTask(xmlName, false) }, new List<string> { mzmlName }, "taskId1").ToString();
        }
    }
}