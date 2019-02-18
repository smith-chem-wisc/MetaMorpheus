using EngineLayer;
using MassSpectrometry;
using NUnit.Framework;
using Proteomics;
using Proteomics.AminoAcidPolymer;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.IO;
using System.Text;
using TaskLayer;
using UsefulProteomicsDatabases;

namespace Test
{
    [TestFixture]
    public static class SilacTest
    {
        [Test]
        public static void TestSilacQuantification()
        {
            Residue heavyLysine = new Residue("a", 'a', "a", Chemistry.ChemicalFormula.ParseFormula("C{13}6H12N{15}2O"), ModificationSites.All); //+8 lysine
            Residue lightLysine = Residue.GetResidue('K');
            List<SilacLabel> silacLabels = new List<Proteomics.SilacLabel>
            {
                new SilacLabel(lightLysine.Letter, lightLysine.Letter, heavyLysine.ThisChemicalFormula.Formula, heavyLysine.MonoisotopicMass - lightLysine.MonoisotopicMass)
            };
            SearchTask task = new SearchTask
            {
                CommonParameters = new CommonParameters(),
                SearchParameters = new SearchParameters
                {
                    SilacLabels = silacLabels
                }
            };

            List<(string, MetaMorpheusTask)> taskList = new List<(string, MetaMorpheusTask)> { ("1", task) };

            PeptideWithSetModifications lightPeptide = new PeptideWithSetModifications("PEPTIDEK", new Dictionary<string, Proteomics.Modification>());

            List<double> massDifferences = new List<double>{ heavyLysine.MonoisotopicMass - lightLysine.MonoisotopicMass };
            MsDataFile myMsDataFile1 = new TestDataFile(lightPeptide, massDifferences);
            string mzmlName = @"silac.mzML";
            IO.MzML.MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(myMsDataFile1, mzmlName, false);

            string xmlName = "SilacDb.xml";
            {
                Protein theProtein = new Protein("PEPTIDEK", "accession1");
                ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), new List<Protein> { theProtein }, xmlName);
            }

            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestSilac");
            Directory.CreateDirectory(outputFolder);
            var theStringResult = task.RunTask(outputFolder, new List<DbForTask> { new DbForTask(xmlName, false) }, new List<string> { mzmlName }, "taskId1").ToString();

            Assert.IsTrue(theStringResult.Contains("All target PSMS within 1% FDR: 1"));
            Directory.Delete(outputFolder, true);
            File.Delete(xmlName);
            File.Delete(mzmlName);
        }
    }
}
