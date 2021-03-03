using Chemistry;
using EngineLayer;
using EngineLayer.ClassicSearch;
using MzLibUtil;
using NUnit.Framework;
using Proteomics;
using Proteomics.Fragmentation;
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
    internal static class AmbiguityTest
    {
        [Test]
        public static void TestResolveAmbiguities()
        {
            Protease protease = new Protease("Custom Protease4", CleavageSpecificity.Full, null, null, new List<DigestionMotif> { new DigestionMotif("K", null, 1, "") });
            ProteaseDictionary.Dictionary.Add(protease.Name, protease);
            CommonParameters CommonParameters_t = new CommonParameters(
                dissociationType: MassSpectrometry.DissociationType.HCD,
                digestionParams: new DigestionParams(
                    protease: protease.Name,
                    minPeptideLength: 1),
                scoreCutoff: 1,
                reportAllAmbiguity: true);
            var fsp_t = new List<(string fileName, CommonParameters fileSpecificParameters)>();
            fsp_t.Add(("", CommonParameters_t));

            CommonParameters CommonParameters_f = new CommonParameters(
                dissociationType: MassSpectrometry.DissociationType.HCD,
                digestionParams: new DigestionParams(
                    protease: protease.Name,
                    minPeptideLength: 1),
                scoreCutoff: 1,
                reportAllAmbiguity: false);
            var fsp_f = new List<(string fileName, CommonParameters fileSpecificParameters)>();
            fsp_f.Add(("", CommonParameters_t));

            var myMsDataFile = new TestDataFile();
            var variableModifications = new List<Modification>();
            var fixedModifications = new List<Modification>();
            var proteinList = new List<Protein> { new Protein("MNNKNKNKQQQ", "Prot1"), new Protein("MNNNKQQQ", "Prot2") };
            var searchModes = new SinglePpmAroundZeroSearchMode(5);

            Tolerance DeconvolutionMassTolerance = new PpmTolerance(5);

            var listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(myMsDataFile, null, new CommonParameters()).OrderBy(b => b.PrecursorMass).ToArray();
            
            PeptideSpectralMatch[] allPsmsArray_withAmbiguity = new PeptideSpectralMatch[listOfSortedms2Scans.Length];
            
            PeptideSpectralMatch[] allPsmsArray_withOutAmbiguity = new PeptideSpectralMatch[listOfSortedms2Scans.Length];

            new ClassicSearchEngine(allPsmsArray_withAmbiguity, listOfSortedms2Scans, variableModifications, fixedModifications, null, null, null, 
                proteinList, searchModes, CommonParameters_t, fsp_t, null, new List<string>()).Run(); //report all ambiguity TRUE
            new ClassicSearchEngine(allPsmsArray_withOutAmbiguity, listOfSortedms2Scans, variableModifications, fixedModifications, null, null, null,
                proteinList, searchModes, CommonParameters_f, fsp_f, null, new List<string>()).Run(); //report all ambiguity FALSE

            Assert.AreEqual("QQQ", allPsmsArray_withAmbiguity[0].BaseSequence);
            Assert.AreEqual("QQQ", allPsmsArray_withOutAmbiguity[0].BaseSequence);
            Assert.IsTrue(allPsmsArray_withAmbiguity[0].ProteinLength == null);
            Assert.IsTrue(allPsmsArray_withOutAmbiguity[0].ProteinLength != null);
            Assert.IsTrue(allPsmsArray_withAmbiguity[0].OneBasedStartResidueInProtein == null);
            Assert.IsTrue(allPsmsArray_withOutAmbiguity[0].OneBasedStartResidueInProtein != null);
        }

        [Test]
        public static void TestContaminantAmbiguity()
        {
            //create an ms file and a database for the peptide
            Protein targetProtein = new Protein("PEPTIDE", "target");
            string xmlName = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\PEPTIDE.xml");
            ProteinDbWriter.WriteXmlDatabase(null, new List<Protein> { targetProtein }, xmlName);
            PeptideWithSetModifications pepWithSetMods = targetProtein.Digest(new DigestionParams(), null, null).First();
            TestDataFile msFile = new TestDataFile(pepWithSetMods);
            string mzmlName = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\PEPTIDE.mzML");
            IO.MzML.MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(msFile, mzmlName, false);
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestContaminantAmbiguityOutput");

            //run a full modern search using two databases (the same database) but one is called a target and the other is called a contaminant
            //KEEP BOTH TARGET AND CONTAMINANT
            SearchParameters modernSearchParams = new SearchParameters();
            modernSearchParams.SearchType = SearchType.Modern;
            modernSearchParams.TCAmbiguity = TargetContaminantAmbiguity.RenameProtein;
            SearchTask modernTask = new SearchTask();
            modernTask.SearchParameters = modernSearchParams;

            EverythingRunnerEngine engine = new EverythingRunnerEngine(new List<(string, MetaMorpheusTask)> { ("task1", modernTask) }, new List<string> { mzmlName }, new List<DbForTask> { new DbForTask(xmlName, false), new DbForTask(xmlName, true) }, outputFolder);
            engine.Run();
            //run the modern search again now that it's reading the index instead of writing it.
            engine.Run();

            //check that the psm file shows it's both a target and a contaminant
            string psmLine = File.ReadAllLines(Path.Combine(outputFolder, "task1", "AllPSMs.psmtsv"))[1];
            string[] splitLine = psmLine.Split('\t');
            Assert.IsTrue(splitLine[30].Equals("N|Y")); //column "Contaminant"
            Assert.IsTrue(splitLine[37].Equals("T|C")); //column "Decoy/Contaminant/Target"


            //KEEP ONLY TARGET
            modernSearchParams = new SearchParameters();
            modernSearchParams.SearchType = SearchType.Modern;
            modernSearchParams.TCAmbiguity = TargetContaminantAmbiguity.RemoveContaminant;
            modernTask = new SearchTask();
            modernTask.SearchParameters = modernSearchParams;

            engine = new EverythingRunnerEngine(new List<(string, MetaMorpheusTask)> { ("task1", modernTask) }, new List<string> { mzmlName }, new List<DbForTask> { new DbForTask(xmlName, false), new DbForTask(xmlName, true) }, outputFolder);
            engine.Run();
            //run the modern search again now that it's reading the index instead of writing it.
            engine.Run();

            //check that the psm file shows it's both a target and a contaminant
            psmLine = File.ReadAllLines(Path.Combine(outputFolder, "task1", "AllPSMs.psmtsv"))[1];
            splitLine = psmLine.Split('\t');
            Assert.IsTrue(splitLine[30].Equals("N")); //column "Contaminant"
            Assert.IsTrue(splitLine[37].Equals("T")); //column "Decoy/Contaminant/Target"


            //KEEP ONLY CONTAMINANT
            modernSearchParams = new SearchParameters();
            modernSearchParams.SearchType = SearchType.Modern;
            modernSearchParams.TCAmbiguity = TargetContaminantAmbiguity.RemoveTarget;
            modernTask = new SearchTask();
            modernTask.SearchParameters = modernSearchParams;

            engine = new EverythingRunnerEngine(new List<(string, MetaMorpheusTask)> { ("task1", modernTask) }, new List<string> { mzmlName }, new List<DbForTask> { new DbForTask(xmlName, false), new DbForTask(xmlName, true) }, outputFolder);
            engine.Run();
            //run the modern search again now that it's reading the index instead of writing it.
            engine.Run();

            //check that the psm file shows it's both a target and a contaminant
            psmLine = File.ReadAllLines(Path.Combine(outputFolder, "task1", "AllPSMs.psmtsv"))[1];
            splitLine = psmLine.Split('\t');
            Assert.IsTrue(splitLine[30].Equals("Y")); //column "Contaminant"
            Assert.IsTrue(splitLine[37].Equals("C")); //column "Decoy/Contaminant/Target"


            Directory.Delete(outputFolder, true);
        }

        public class SanitizeTestClass : SearchTask
        {
            [Test]
            public static void TestDatabaseSanitizationPrioritizesXml()
            {
                //if a fasta and a xml database are both used and called targets, prioritize the xml
                Protein fastaProtein = new Protein("APEPTIDE", "Test");
                ModificationMotif.TryGetMotif("A", out ModificationMotif motif);
                Modification mod = new Modification(_originalId: "acetylation", _modificationType: "testModType", _target: motif, _chemicalFormula: ChemicalFormula.ParseFormula("C2H2O1"), _locationRestriction: "Anywhere.");
                Protein xmlProtein = new Protein("APEPTIDE", "Test", oneBasedModifications: new Dictionary<int, List<Modification>> { { 1, new List<Modification> { mod } } });
                List<Protein> proteins = new List<Protein> { fastaProtein, xmlProtein };
                SanitizeProteinDatabase(proteins, TargetContaminantAmbiguity.RemoveTarget);
                Assert.IsTrue(proteins.Count == 1);
                Assert.IsTrue(proteins.First().OneBasedPossibleLocalizedModifications.Count != 0);

                //reverse order and try again
                proteins = new List<Protein> { xmlProtein, fastaProtein };
                SanitizeProteinDatabase(proteins, TargetContaminantAmbiguity.RemoveTarget);
                Assert.IsTrue(proteins.Count == 1);
                Assert.IsTrue(proteins.First().OneBasedPossibleLocalizedModifications.Count != 0);

                //same with no mods
                xmlProtein = new Protein("APEPTIDE", "Test", proteolysisProducts: new List<ProteolysisProduct> { new ProteolysisProduct(1, 3, "zrWuzHere") });
                proteins = new List<Protein> { xmlProtein, fastaProtein };
                SanitizeProteinDatabase(proteins, TargetContaminantAmbiguity.RemoveTarget);
                Assert.IsTrue(proteins.Count == 1);
                Assert.IsTrue(proteins.First().ProteolysisProducts.Count() != 0);
            }

            [Test]
            public static void TestDatabaseSanitizationRemovesDecoys()
            {
                ModificationMotif.TryGetMotif("A", out ModificationMotif motif);
                Modification mod = new Modification(_originalId: "acetylation", _modificationType: "testModType", _target: motif, _chemicalFormula: ChemicalFormula.ParseFormula("C2H2O1"), _locationRestriction: "Anywhere.");
                Protein targetProtein = new Protein("APEPTIDE", "Test", oneBasedModifications: new Dictionary<int, List<Modification>> { { 1, new List<Modification> { mod } } });
                Protein contaminantProtein = new Protein("APEPTIDE", "Test", isContaminant: true);
                Protein tDecoyProtein = new Protein("AEDITPEP", "DECOY_Test", oneBasedModifications: new Dictionary<int, List<Modification>> { { 1, new List<Modification> { mod } } });
                Protein cDecoyProtein = new Protein("AEDITPEP", "DECOY_Test");
                List<Protein> proteins = new List<Protein> { targetProtein, contaminantProtein, tDecoyProtein, cDecoyProtein }; //two decoys, one for target, one for contaminant
                SanitizeProteinDatabase(proteins, TargetContaminantAmbiguity.RemoveContaminant);
                Assert.IsTrue(proteins.Count == 2);
                Assert.IsTrue(!proteins[0].IsContaminant); //forward is target
                Assert.IsTrue(proteins[0].OneBasedPossibleLocalizedModifications.Count == 1);
                Assert.IsTrue(proteins[1].OneBasedPossibleLocalizedModifications.Count == 1);

                proteins = new List<Protein> { targetProtein, contaminantProtein, tDecoyProtein, cDecoyProtein }; //two decoys, one for target, one for contaminant
                SanitizeProteinDatabase(proteins, TargetContaminantAmbiguity.RemoveTarget);
                Assert.IsTrue(proteins.Count == 2);
                Assert.IsTrue(proteins[0].IsContaminant); //forward is contaminant
                Assert.IsTrue(proteins[0].OneBasedPossibleLocalizedModifications.Count == 0); //contaminant doesn't have mods                
                Assert.IsTrue(proteins[1].OneBasedPossibleLocalizedModifications.Count == 1); //would be cool if this was 0 and we could keep the correct decoy, but kind of an edge case
            }
        }
    }
}