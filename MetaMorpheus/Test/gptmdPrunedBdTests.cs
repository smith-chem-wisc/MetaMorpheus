using EngineLayer;
using MassSpectrometry;
using NUnit.Framework;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using Omics.BioPolymer;
using Omics.Modifications;
using TaskLayer;
using UsefulProteomicsDatabases;
using Omics;

namespace Test
{
    [TestFixture]
    public static class GptmdPrunedDbTests
    {
        // want a psm whose base sequence is not ambigous but full sequence is (ptm is not localized): make sure this does not make it in DB

        [Test]
        public static void TestPrunedGeneration()
        {
            //Create GPTMD Task
            //Create Search Task
            GptmdTask task1 = new GptmdTask
            {
                CommonParameters = new CommonParameters(),
                GptmdParameters = new GptmdParameters
                {
                    ListOfModsGptmd = GlobalVariables.AllModsKnown.Where(b =>
                        b.ModificationType.Equals("Common Artifact")
                        || b.ModificationType.Equals("Common Biological")
                        || b.ModificationType.Equals("Metal")
                        || b.ModificationType.Equals("Less Common")
                        ).Select(b => (b.ModificationType, b.IdWithMotif)).ToList()
                }
            };

            SearchTask task2 = new SearchTask
            {
                CommonParameters = new CommonParameters(),

                SearchParameters = new SearchParameters
                {
                    DoParsimony = true,
                    SearchTarget = true,
                    WritePrunedDatabase = true,
                    SearchType = SearchType.Classic
                }
            };
            List<(string, MetaMorpheusTask)> taskList = new List<(string, MetaMorpheusTask)> { ("task1", task1), ("task2", task2) };
            string mzmlName = @"TestData\PrunedDbSpectra.mzml";
            string fastaName = @"TestData\DbForPrunedDb.fasta";
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestPrunedGeneration");
            var engine = new EverythingRunnerEngine(taskList, new List<string> { mzmlName }, new List<DbForTask> { new DbForTask(fastaName, false) }, outputFolder);
            engine.Run();
            string final = Path.Combine(MySetUpClass.outputFolder, "task2", "DbForPrunedDbGPTMDproteinPruned.xml");
            List<Protein> proteins = ProteinDbLoader.LoadProteinXML(final, true, DecoyType.Reverse, new List<Modification>(), false, new List<string>(), out var ok);
            // ensures that protein out put contains the correct number of proteins to match the following conditions.
            // all proteins in DB have baseSequence!=null (not ambiguous)
            // all proteins that belong to a protein group are written to DB
            Assert.That(proteins.Count, Is.EqualTo(18));
            int totalNumberOfMods = proteins.Sum(p => p.OneBasedPossibleLocalizedModifications.Count + p.SequenceVariations.Sum(sv => sv.OneBasedModifications.Count));


            //tests that modifications are being done correctly
            Assert.That(totalNumberOfMods, Is.EqualTo(4));
            Directory.Delete(outputFolder, true);
        }

        //test if prunedDatabase matches expected output
        [Test]
        public static void TestPrunedDatabase()
        {
            //Create Search Task
            SearchTask task1 = new SearchTask
            {
                SearchParameters = new SearchParameters
                {
                    WritePrunedDatabase = true,
                    SearchTarget = true,
                    MassDiffAcceptorType = MassDiffAcceptorType.Exact,
                    ModsToWriteSelection = new Dictionary<string, int>
                    {
                        {"ConnorModType", 1}
                    }
                },
                CommonParameters = new CommonParameters(digestionParams: new DigestionParams(minPeptideLength: 5))
            };

            //add task to task list
            List<(string, MetaMorpheusTask)> taskList = new List<(string, MetaMorpheusTask)>
            {
               ("task1", task1)
            };

            ModificationMotif.TryGetMotif("P", out ModificationMotif motif);

            var connorMod = new Modification(_originalId: "ConnorMod on P", _modificationType: "ConnorModType", _target: motif, _locationRestriction: "Anywhere.", _monoisotopicMass: 10);

            GlobalVariables.AddMods(new List<Modification>
            {
                connorMod
            }, false);

            //create modification lists
            List<Modification> variableModifications = GlobalVariables.AllModsKnown.OfType<Modification>()
                .Where(b => task1.CommonParameters.ListOfModsVariable.Contains((b.ModificationType, b.IdWithMotif))).ToList();

            //add modification to Protein object
            var dictHere = new Dictionary<int, List<Modification>>();
            Modification modToAdd = connorMod;
            Modification modToAdd2 = connorMod;
            dictHere.Add(1, new List<Modification> { modToAdd });
            dictHere.Add(3, new List<Modification> { modToAdd2 });

            //protein Creation (One with mod and one without)
            Protein TestProteinWithMod = new Protein("PEPTID", "accession1", "organism", new List<Tuple<string, string>>(), dictHere);

            //First Write XML Database
            string xmlName = "okkk.xml";

            //Add Mod to list and write XML input database
            var Hash = new HashSet<Tuple<int, Modification>>
            {
                new Tuple<int, Modification>(3, modToAdd)
            };

            ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), new List<Protein> { TestProteinWithMod }, xmlName);

            //now write MZML file
            var protein = ProteinDbLoader.LoadProteinXML(xmlName, true,
                DecoyType.Reverse, new List<Modification>(), false, new List<string>(), out Dictionary<string, Modification> ok);

            //Dictionary 'ok' contains unknown modifications. There are no unknown modifications in this test.
            Assert.That(ok.Count, Is.EqualTo(0));
            //One protein is read from the .xml database and one decoy is created. Therefore, the list of proteins contains 2 entries.
            Assert.That(protein.Count, Is.EqualTo(2));
            //The original database had two localized mods on the protein. Therefore. both protein and decoy should have two mods.
            Assert.That(protein[0].OneBasedPossibleLocalizedModifications.Count, Is.EqualTo(2));
            List<int> foundResidueIndicies = protein[0].OneBasedPossibleLocalizedModifications.Select(k => k.Key).ToList();
            List<int> expectedResidueIndices = new List<int>() { 1, 3 };
            Assert.That(foundResidueIndicies, Is.EquivalentTo(expectedResidueIndices));
            Assert.That(protein[1].OneBasedPossibleLocalizedModifications.Count, Is.EqualTo(2));
            foundResidueIndicies = protein[1].OneBasedPossibleLocalizedModifications.Select(k => k.Key).ToList();
            expectedResidueIndices = new List<int>() { 4, 6 }; //originally modified residues are now at the end in the decoy
            Assert.That(foundResidueIndicies, Is.EquivalentTo(expectedResidueIndices));

            var thisOk = ok;//for debugging
            var commonParamsAtThisPoint = task1.CommonParameters.DigestionParams; //for debugging

            var digestedList = protein[0].Digest(task1.CommonParameters.DigestionParams, new List<Modification> { },
                variableModifications).ToList();
            Assert.That(digestedList.Count, Is.EqualTo(4));

            //Set Peptide with 1 mod at position 3
            var pepWithSetMods1 = digestedList[1];

            //Finally Write MZML file
            Assert.That(pepWithSetMods1.FullSequence, Is.EqualTo("PEP[ConnorModType:ConnorMod on P]TID"));//this might be base sequence
            MsDataFile myMsDataFile = new TestDataFile(new List<IBioPolymerWithSetMods> { pepWithSetMods1 });
            string mzmlName = @"hello.mzML";
            Readers.MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(myMsDataFile, mzmlName, false);

            //run!
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestPrunedDatabase");
            var engine = new EverythingRunnerEngine(taskList, new List<string> { mzmlName },
                new List<DbForTask> { new DbForTask(xmlName, false) }, outputFolder);
            engine.Run();

            string final = Path.Combine(MySetUpClass.outputFolder, "task1", "okkkpruned.xml");

            var proteins = ProteinDbLoader.LoadProteinXML(final, true, DecoyType.Reverse, new List<Modification>(), false, new List<string>(), out ok);
            //check length
            Assert.That(proteins[0].OneBasedPossibleLocalizedModifications.Count, Is.EqualTo(1));
            //check location (key)
            Assert.That(proteins[0].OneBasedPossibleLocalizedModifications.ContainsKey(3), Is.EqualTo(true));
            List<Modification> listOfMods = proteins[0].OneBasedPossibleLocalizedModifications[3];
            //check Type, count, ID
            Assert.That(listOfMods[0].ModificationType, Is.EqualTo("ConnorModType"));
            Assert.That(listOfMods[0].IdWithMotif, Is.EqualTo("ConnorMod on P"));
            Assert.That(listOfMods.Count, Is.EqualTo(1));
            Directory.Delete(outputFolder, true);
            File.Delete(xmlName);
            File.Delete(mzmlName);
        }

        [Test]
        public static void TestUserModSelectionInPrunedDB()
        {
            List<(string, string)> listOfModsFixed = new List<(string, string)> { ("Common Fixed", "Carbamidomethyl on C"), ("Common Fixed", "Carbamidomethyl on U") };
            //Create Search Task
            SearchTask task5 = new SearchTask
            {
                SearchParameters = new SearchParameters
                {
                    WritePrunedDatabase = true,
                    SearchTarget = true,
                    MassDiffAcceptorType = MassDiffAcceptorType.Exact,
                },
                CommonParameters = new CommonParameters(listOfModsFixed: listOfModsFixed)
            };

            task5.SearchParameters.ModsToWriteSelection["Mod"] = 0;
            task5.SearchParameters.ModsToWriteSelection["Common Fixed"] = 1;
            task5.SearchParameters.ModsToWriteSelection["Glycan"] = 2;
            task5.SearchParameters.ModsToWriteSelection["missing"] = 3;

            //add task 1 to task list
            List<(string, MetaMorpheusTask)> taskList = new List<(string, MetaMorpheusTask)> { ("task5", task5) };
            ModificationMotif.TryGetMotif("P", out ModificationMotif motif);
            ModificationMotif.TryGetMotif("E", out ModificationMotif motif2);

            var connorMod = new Modification(_originalId: "ModToNotAppear", _modificationType: "Mod", _target: motif, _locationRestriction: "Anywhere.", _monoisotopicMass: 10);
            var connorMod2 = new Modification(_originalId: "Default(Mod in DB and Observed)", _modificationType: "Common Fixed", _target: motif, _locationRestriction: "Anywhere.", _monoisotopicMass: 10);
            var connorMod3 = new Modification(_originalId: "ModToAlwaysAppear", _modificationType: "Glycan", _target: motif, _locationRestriction: "Anywhere.", _monoisotopicMass: 10);
            var connorMod4 = new Modification(_originalId: "ModObservedNotinDB", _modificationType: "missing", _target: motif2, _locationRestriction: "Anywhere.", _monoisotopicMass: 5);

            GlobalVariables.AddMods(new List<Modification>
            {
                connorMod,
                connorMod2,
                connorMod3,
                connorMod4
            }, false);

            //create modification lists
            List<Modification> variableModifications = GlobalVariables.AllModsKnown.OfType<Modification>().Where(b => task5.CommonParameters.ListOfModsVariable.Contains
            ((b.ModificationType, b.IdWithMotif))).ToList();
            List<Modification> fixedModifications = GlobalVariables.AllModsKnown.OfType<Modification>().Where(b => task5.CommonParameters.ListOfModsFixed.Contains
            ((b.ModificationType, b.IdWithMotif))).ToList();

            //add modification to Protein object
            var dictHere = new Dictionary<int, List<Modification>>();
            Modification modToAdd = connorMod;
            Modification modToAdd2 = connorMod2;
            Modification modToAdd3 = connorMod3;
            Modification modToAdd4 = connorMod4;

            //add Fixed modifcation so can test if mod that is observed and not in DB
            fixedModifications.Add(connorMod4);
            listOfModsFixed.Add((connorMod4.ModificationType, connorMod4.IdWithMotif));

            dictHere.Add(1, new List<Modification> { modToAdd });
            dictHere.Add(2, new List<Modification> { modToAdd2 }); //default
            dictHere.Add(3, new List<Modification> { modToAdd3 }); //Alway Appear

            var dictHere2 = new Dictionary<int, List<Modification>>
            {
                { 1, new List<Modification> { modToAdd } },
                { 2, new List<Modification> { modToAdd2 } }, //default
                { 3, new List<Modification> { modToAdd3 } }, //Alway Appear
                { 4, new List<Modification> { modToAdd4 } }//observed
            };

            //protein Creation (One with mod and one without)
            Protein TestProteinWithModForDB = new Protein("PPPPPPPPPPE", "accession1", "organism", new List<Tuple<string, string>>(), dictHere);
            Protein TestProteinWithModObsevred = new Protein("PPPPPPPPPPE", "accession1", "organism", new List<Tuple<string, string>>(), dictHere2);

            //First Write XML Database
            string xmlName = "selectedMods.xml";
            string xmlName2 = "selectedModsObvs.xml";

            //Add Mod to list and write XML input database
            Dictionary<string, HashSet<Tuple<int, Modification>>> modList = new Dictionary<string, HashSet<Tuple<int, Modification>>>();
            var Hash = new HashSet<Tuple<int, Modification>>
            {
                new Tuple<int, Modification>(1, modToAdd),
                new Tuple<int, Modification>(2, modToAdd2),
                new Tuple<int, Modification>(3, modToAdd3),
                new Tuple<int, Modification>(4, modToAdd4), //Observed Only
            };

            modList.Add("test", Hash);
            ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), new List<Protein> { TestProteinWithModForDB }, xmlName);

            //Add Observed Only
            modList.Add("test2", Hash);
            ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), new List<Protein> { TestProteinWithModObsevred }, xmlName2);

            //now create MZML data
            var protein = ProteinDbLoader.LoadProteinXML(xmlName2, true, DecoyType.Reverse, new List<Modification>(), false, new List<string>(), out Dictionary<string, Modification> ok);
            var digestedList = protein[0].Digest(task5.CommonParameters.DigestionParams, fixedModifications, variableModifications).ToList();

            //Set Peptide with 1 mod at position 3
            var pepWithSetMods1 = digestedList[0];
            var pepWithSetMods2 = digestedList[1];
            var pepWithSetMods3 = digestedList[2];
            var pepWithSetMods4 = digestedList[3];
            var pepWithSetMods5 = digestedList[4];

            //CUSTOM PEP
            MsDataFile myMsDataFile = new TestDataFile(new List<IBioPolymerWithSetMods>
            { pepWithSetMods1, pepWithSetMods2, pepWithSetMods3, pepWithSetMods4, pepWithSetMods5 });
            string mzmlName = @"newMzml.mzML";
            Readers.MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(myMsDataFile, mzmlName, false);

            //make sure this runs correctly
            //run!
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestUserModSelectionInPrunedDB");
            var engine = new EverythingRunnerEngine(taskList, new List<string> { mzmlName }, new List<DbForTask> { new DbForTask(xmlName, false) }, outputFolder);
            engine.Run();
            string final = Path.Combine(MySetUpClass.outputFolder, "task5", "selectedModspruned.xml");
            var proteins = ProteinDbLoader.LoadProteinXML(final, true, DecoyType.Reverse, new List<Modification>(), false, new List<string>(), out ok);
            var Dlist = proteins[0].GetVariantBioPolymers().SelectMany(vp => vp.Digest(task5.CommonParameters.DigestionParams, fixedModifications, variableModifications)).ToList();
            Assert.That(Dlist[0].NumFixedMods, Is.EqualTo(1));

            //check length
            Assert.That(proteins[0].OneBasedPossibleLocalizedModifications.Count, Is.EqualTo(3));
            List<Modification> listOfLocalMods = new List<Modification>();
            listOfLocalMods.AddRange(proteins[0].OneBasedPossibleLocalizedModifications[2]);
            listOfLocalMods.AddRange(proteins[0].OneBasedPossibleLocalizedModifications[3]);
            listOfLocalMods.AddRange(proteins[0].OneBasedPossibleLocalizedModifications[11]);

            //check Type, count, ID
            Assert.That(listOfLocalMods[0].ModificationType, Is.EqualTo("Common Fixed"));
            Assert.That(listOfLocalMods[2].ModificationType, Is.EqualTo("missing"));
            Assert.That(!listOfLocalMods.Contains(connorMod)); //make sure that mod set not to show up is not in mod list

            Assert.That(listOfLocalMods[0].IdWithMotif, Is.EqualTo("Default(Mod in DB and Observed) on P"));
            Assert.That(listOfLocalMods[1].IdWithMotif, Is.EqualTo("ModToAlwaysAppear on P"));
            //Makes sure Mod that was not in the DB but was observed is in pruned DB
            Assert.That(listOfLocalMods[2].IdWithMotif, Is.EqualTo("ModObservedNotinDB on E"));
            Assert.That(listOfLocalMods.Count, Is.EqualTo(3));
            Directory.Delete(outputFolder, true);
            File.Delete(mzmlName);
            File.Delete(xmlName);
            File.Delete(xmlName2);
        }

        [Test]
        public static void TestProteinPrunedWithModSelection()
        {
            var modToWrite = GlobalVariables.AllModsKnown.Where(p => p.ModificationType == "UniProt" && p.Target.ToString() == "T").First();
            var modToNotWrite = GlobalVariables.AllModsKnown.Where(p => p.ModificationType == "Common Artifact" && p.Target.ToString() == "X").First();

            var protein1 = new Protein("PEPTIDEKPEPT", "1", oneBasedModifications: new Dictionary<int, List<Modification>> { { 1, new List<Modification> { modToNotWrite } }, { 12, new List<Modification> { modToWrite } } });
            var protein2 = new Protein("PEPTIDPEPT", "2", oneBasedModifications: new Dictionary<int, List<Modification>> { { 1, new List<Modification> { modToNotWrite } }, { 10, new List<Modification> { modToWrite } } });

            string path = @"temp";

            ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), new List<Protein> { protein1, protein2 }, path);

            Directory.CreateDirectory(Path.Combine(TestContext.CurrentContext.TestDirectory, @"PrunedDbTest"));

            Dictionary<string, HashSet<Tuple<int, Modification>>> modList = new Dictionary<string, HashSet<Tuple<int, Modification>>>();
            var Hash = new HashSet<Tuple<int, Modification>>
            {
                new Tuple<int, Modification>(1, modToWrite),
                new Tuple<int, Modification>(2, modToNotWrite),

            };

            var db = ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), new List<Protein> { protein1, protein2 }, Path.Combine(TestContext.CurrentContext.TestDirectory, @"PrunedDbTest/fakeDb.xml"));

            var peptideObserved = protein1.Digest(new DigestionParams(minPeptideLength: 1), new List<Modification>(), new List<Modification>())
            .Where(p => p.BaseSequence == "PEPT" && p.AllModsOneIsNterminus.Count > 0).First();
            PostSearchAnalysisParameters testPostTaskParameters = new PostSearchAnalysisParameters();
            CommonParameters commonParam = new CommonParameters();
            double[,] noiseData = new double[10000, 10000];
            noiseData[0,0] = 1.0; 
            List<Omics.Fragmentation.MatchedFragmentIon> matchedFragmentIons = new List<Omics.Fragmentation.MatchedFragmentIon>() { };
            MzSpectrum spectrum = new MzSpectrum(noiseData);
            MsDataScan scan = new MsDataScan(spectrum , 1, 1, true, Polarity.Unknown, 2, new MzLibUtil.MzRange(10, 1000), "", MZAnalyzerType.Orbitrap, 10000, null, noiseData, "");
            testPostTaskParameters.BioPolymerList = new List<IBioPolymer>() { protein1, protein2 };
            testPostTaskParameters.AllSpectralMatches = new List<SpectralMatch> { new PeptideSpectralMatch(peptideObserved, 0, 20, 1, new Ms2ScanWithSpecificMass(scan, 100, 1, @"", commonParam), commonParam, matchedFragmentIons) };
            testPostTaskParameters.SearchParameters = new SearchParameters();
            testPostTaskParameters.SearchParameters.WritePrunedDatabase = true;
            testPostTaskParameters.SearchParameters.DoLabelFreeQuantification = false;
            testPostTaskParameters.SearchParameters.WriteMzId = false;
            testPostTaskParameters.DatabaseFilenameList = new List<DbForTask>() { new DbForTask(Path.Combine(TestContext.CurrentContext.TestDirectory, @"PrunedDbTest/fakeDb.xml"), false) };
            testPostTaskParameters.OutputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"PrunedDbTest");
            Directory.CreateDirectory(Path.Combine(TestContext.CurrentContext.TestDirectory, @"PrunedDbTest/individual"));
            testPostTaskParameters.IndividualResultsOutputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"PrunedDbTest/individual");
            int[] stuffForSpectraFile = new int[2];
            stuffForSpectraFile[0] = 10;
            stuffForSpectraFile[1] = 10;
            Dictionary<string, int[]> numSpectraPerFile = new Dictionary<string, int[]>();
            numSpectraPerFile.Add("", stuffForSpectraFile);
            testPostTaskParameters.NumMs2SpectraPerFile = numSpectraPerFile;

            MsDataFile myMsDataFile = new TestDataFile(new List<IBioPolymerWithSetMods>
            { peptideObserved});
            string mzmlName = @"newMzml.mzML";
            Readers.MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(myMsDataFile, mzmlName, false);
                        
            modList.Add("test", Hash);
            
            testPostTaskParameters.CurrentRawFileList = new List<string>() { mzmlName };

            SearchTask task5 = new SearchTask
            {
                SearchParameters = new SearchParameters
                {
                    WritePrunedDatabase = true,
                    SearchTarget = true,
                    MassDiffAcceptorType = MassDiffAcceptorType.Exact,
                },
                CommonParameters = new CommonParameters()
            };

            var test = task5.RunTask(Path.Combine(TestContext.CurrentContext.TestDirectory, @"PrunedDbTest"), new List<DbForTask>() { new DbForTask(Path.Combine(TestContext.CurrentContext.TestDirectory, @"PrunedDbTest/fakeDb.xml"),false)}, new List<string>() { mzmlName }, "name");
            testPostTaskParameters.SearchTaskResults = test;
           
            PostSearchAnalysisTask testPostTask = new PostSearchAnalysisTask();
            testPostTask.Parameters = testPostTaskParameters;
            testPostTask.CommonParameters = commonParam;
            testPostTask.FileSpecificParameters = new List<(string FileName, CommonParameters Parameters)> { ("newMzMl.mzml", commonParam) };
            testPostTask.Run();

            var proteinsLoaded = ProteinDbLoader.LoadProteinXML(path, true, DecoyType.None, GlobalVariables.AllModsKnown, false, new List<string>(), out var unknownMods);

            // assert that mods on proteins are the same before/after task is run
            Assert.That(protein1.Equals(proteinsLoaded.First(p => p.Accession == "1")));
            Assert.That(protein2.Equals(proteinsLoaded.First(p => p.Accession == "2")));

            // assert that protein pruned DB has correct proteins mods
            var proteinPruned = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, @"PrunedDbTest/fakeDbproteinPruned.xml"), true, DecoyType.None, GlobalVariables.AllModsKnown, false, new List<string>(), out var unknownMods1);
            Assert.That(proteinPruned.Count().Equals(1));
            Assert.That(proteinPruned.FirstOrDefault().OneBasedPossibleLocalizedModifications.Count().Equals(1));
            // assert that mod-pruned DB has correct proteins and mods
            var modPruned = ProteinDbLoader.LoadProteinXML(Path.Combine(TestContext.CurrentContext.TestDirectory, @"PrunedDbTest/fakeDbpruned.xml"), true, DecoyType.None, GlobalVariables.AllModsKnown, false, new List<string>(), out var unknownMods2);
            Assert.That(modPruned.Count().Equals(2));
            Assert.That(modPruned.ElementAt(0).OneBasedPossibleLocalizedModifications.Count().Equals(1));
            Assert.That(modPruned.ElementAt(1).OneBasedPossibleLocalizedModifications.Count().Equals(1));
        }
        [Test]
        public static void TestProteinPrunedWithModSelectionAndVariants()
        {
            // PURPOSE
            // This test verifies that:
            // 1) We can construct a small input DB containing target proteins plus variant-applied isoforms.
            // 2) A synthetic PSM derived from a variant isoform leads PostSearchAnalysisTask to write:
            //    - a mod-pruned DB (...pruned.xml) containing all proteins with only confidently localized (and configured) mods
            //    - a protein-pruned DB (...proteinPruned.xml) containing only the proteins/consensus variants actually supported by confident PSMs
            // 3) Baseline properties of proteins and modifications remain stable when read back from the written DBs.

            // ------------------------------------------------------------
            // Arrange
            // ------------------------------------------------------------

            // Pick two existing modifications from the global variable set:
            // - modToWrite (UniProt on T) is intended to survive pruning in realistic scenarios.
            // - modToNotWrite (Common Artifact on X) is a control and is generally filtered.
            var modToWrite = GlobalVariables.AllModsKnown
                .First(p => p.ModificationType == "UniProt" && p.Target.ToString() == "T");
            var modToNotWrite = GlobalVariables.AllModsKnown
                .First(p => p.ModificationType == "Common Artifact" && p.Target.ToString() == "X");

            // A variant-localized modification map to attach to the sequence variation
            var variantMods = new Dictionary<int, List<Modification>>
            {
                { 1, new List<Modification> { modToNotWrite } }
            };

            // Define a single missense variant for protein1: V(4) -> T(4).
            var variants = new List<SequenceVariation>
            {
                new SequenceVariation(
                    4, 4,
                    "V", "T", "",
                    // VCF-like annotation string
                    "20\t41168825\t.\tT\tC\t14290.77\t.\tANN=C|missense_variant|MODERATE|PLCG1|ENSG00000124181|transcript|ENST00000244007.7|protein_coding|22/33|c.2438T>C|p.Ile813Thr|2635/5285|2438/3876|813/1291||\tGT:AD:DP:GQ:PL\t1/1:1,392:393:99:14319,1142,0",
                    variantMods)
            };

            // Construct two proteins:
            // - protein1 has two base-level modifications and one sequence variation (V->T at 4).
            // - protein2 is a control protein with its own base-level mods, no variants.
            var protein1 = new Protein(
                "PEPVIDEKPEPT",
                "1",
                oneBasedModifications: new Dictionary<int, List<Modification>>
                {
                    { 1,  new List<Modification> { modToNotWrite } },
                    { 12, new List<Modification> { modToWrite } }
                },
                sequenceVariations: variants);

            var protein2 = new Protein(
                "PEPIDPEPT",
                "2",
                oneBasedModifications: new Dictionary<int, List<Modification>>
                {
                    { 1, new List<Modification> { modToNotWrite } },
                    { 9, new List<Modification> { modToWrite } }
                });

            // Generate variant-applied isoforms for protein1 (variants only; exclude consensus)
            var protein1Variants = protein1
                .GetVariantBioPolymers(1, 0, maxSequenceVariantIsoforms: 2)
                .Where(v => v.AppliedSequenceVariations.Count > 0)
                .ToList();

            Assert.That(protein1Variants.Count == 1, "Test setup failed: expected exactly one variant isoform");
            Assert.That(protein1Variants[0].BaseSequence == "PEPTIDEKPEPT", "Test setup failed: unexpected variant isoform sequence");

            // Build a working set: targets + applied-variant isoforms.
            var workingProteins = new List<Protein> { protein1, protein2 };
            workingProteins.AddRange(protein1Variants);

            // I/O locations (unique for this test to avoid cross-test contamination)
            string inputDbFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, "PrunedDbTestVariant");
            string inputDbPath = Path.Combine(inputDbFolder, "fakeDb.xml");
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, "PrunedDbTestVariant_Output");
            string rawPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "newMzml.mzML");

            Directory.CreateDirectory(inputDbFolder);
            Directory.CreateDirectory(outputFolder);

            // Write the full input DB (targets + variant isoforms)
            ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), workingProteins, inputDbPath);

            // Create an observed peptide from a variant-applied isoform to simulate evidence coming from the variant.
            var observedVariantIsoform = protein1Variants.First();
            var observedPeptide = observedVariantIsoform
                .Digest(new DigestionParams(minPeptideLength: 1), new List<Modification>(), new List<Modification>())
                .First(p => p.BaseSequence == "PEPT");

            // Build a minimal mzML with a single synthetic MS2 spectrum for the observed peptide
            MsDataFile myMsDataFile = new TestDataFile(new List<IBioPolymerWithSetMods> { observedPeptide });
            Readers.MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(myMsDataFile, rawPath, false);

            // Prepare a synthetic PSM used by PostSearchAnalysis
            var commonParams = new CommonParameters();
            double[,] noiseData = new double[100, 100];
            noiseData[0, 0] = 1.0;
            var spectrum = new MzSpectrum(noiseData);
            var scan = new MsDataScan(
                spectrum, 1, 1, true, Polarity.Unknown, 2,
                new MzLibUtil.MzRange(10, 1000), "", MZAnalyzerType.Orbitrap, 10000, null, noiseData, "");
            var matchedIons = new List<Omics.Fragmentation.MatchedFragmentIon>();
            var psm = new PeptideSpectralMatch(
                observedPeptide, 0, 20, 1,
                new Ms2ScanWithSpecificMass(scan, 100, 1, @"", commonParams),
                commonParams, matchedIons);

            // Configure a SearchTask (needed to produce SearchTaskResults for PostSearchAnalysis)
            var searchTask = new SearchTask
            {
                SearchParameters = new SearchParameters
                {
                    WritePrunedDatabase = true,
                    SearchTarget = true,
                    MassDiffAcceptorType = MassDiffAcceptorType.Exact
                },
                CommonParameters = new CommonParameters()
            };

            // ------------------------------------------------------------
            // Act
            // ------------------------------------------------------------

            // Run the search end-to-end to produce results for PostSearchAnalysis
            var searchResults = searchTask.RunTask(
                outputFolder,
                new List<DbForTask> { new DbForTask(inputDbPath, false) },
                new List<string> { rawPath },
                "name");

            // Prepare and run PostSearchAnalysis (this writes the pruned databases)
            var postParams = new PostSearchAnalysisParameters
            {
                BioPolymerList = workingProteins.Cast<IBioPolymer>().ToList(),
                AllSpectralMatches = new List<SpectralMatch> { psm },
                SearchParameters = new SearchParameters
                {
                    WritePrunedDatabase = true,
                    DoLabelFreeQuantification = false,
                    WriteMzId = false,
                },
                DatabaseFilenameList = new List<DbForTask> { new DbForTask(inputDbPath, false) },
                OutputFolder = outputFolder,
                IndividualResultsOutputFolder = Path.Combine(outputFolder, "individual"),
                NumMs2SpectraPerFile = new Dictionary<string, int[]> { { "", new[] { 10, 10 } } },
                CurrentRawFileList = new List<string> { rawPath },
                SearchTaskResults = searchResults
            };

            Directory.CreateDirectory(postParams.IndividualResultsOutputFolder);

            var postTask = new PostSearchAnalysisTask
            {
                Parameters = postParams,
                CommonParameters = commonParams,
                FileSpecificParameters = new List<(string FileName, CommonParameters Parameters)>
        {
            ("newMzMl.mzml", commonParams)
        }
            };

            postTask.Run();

            // ------------------------------------------------------------
            // Assert
            // ------------------------------------------------------------

            // Sanity: load back the original DB (not pruned) and ensure baseline content matches what we wrote.
            var proteinsLoaded = ProteinDbLoader.LoadProteinXML(
                inputDbPath, true, DecoyType.None, GlobalVariables.AllModsKnown, false, new List<string>(), out var unknownMods);

            // Paths for pruned outputs (names come from WritePrunedDatabase)
            string modPrunedPath = Path.Combine(outputFolder, "fakeDbpruned.xml");
            string proteinPrunedPath = Path.Combine(outputFolder, "fakeDbproteinPruned.xml");

            Assert.Multiple(() =>
            {
                // Input DB expectations
                Assert.That(File.Exists(inputDbPath), "Input DB not written");
                Assert.That(proteinsLoaded, Is.Not.Null);
                Assert.That(proteinsLoaded.Count, Is.GreaterThanOrEqualTo(2), "Input DB must contain at least protein1, protein2, plus any variant isoforms");
                Assert.That(unknownMods.Count, Is.EqualTo(0), "Unexpected unknown modifications while loading input DB");

                // Pruned outputs exist
                Assert.That(File.Exists(modPrunedPath), "Mod-pruned DB was not written");
                Assert.That(File.Exists(proteinPrunedPath), "Protein-pruned DB was not written");

                // Basic integrity checks on pruned outputs (use current LoadProteinXML overload)
                var proteinPruned = ProteinDbLoader.LoadProteinXML(
                    proteinPrunedPath, true, DecoyType.None, GlobalVariables.AllModsKnown, false, new List<string>(), out var _);
                var modPruned = ProteinDbLoader.LoadProteinXML(
                    modPrunedPath, true, DecoyType.None, GlobalVariables.AllModsKnown, false, new List<string>(), out var _);

                // Build accession sets for clear diagnostics
                var proteinPrunedAccessions = proteinPruned.Select(p => p.Accession).ToHashSet();
                var modPrunedAccessions = modPruned.Select(p => p.Accession).ToHashSet();

                // Protein-pruned: should contain the supported consensus protein (accession "1")
                Assert.That(proteinPruned.Any(), "Protein-pruned DB is empty");
                Assert.That(proteinPrunedAccessions.Contains("1"),
                    $"Protein-pruned DB missing supported accession '1'. Actual: [{string.Join(", ", proteinPrunedAccessions)}]");

                // Also verify pruning behavior on that protein: exactly 1 kept site
                var prunedProt1 = proteinPruned.First(p => p.Accession == "1");
                Assert.That(prunedProt1.OneBasedPossibleLocalizedModifications.Count, Is.EqualTo(1),
                    $"Protein-pruned '1' expected 1 kept mod site; actual: {prunedProt1.OneBasedPossibleLocalizedModifications.Count}");

                // Mod-pruned: should include both consensus proteins "1" and "2"
                Assert.That(modPruned.Any(), "Mod-pruned DB is empty");
                Assert.That(modPrunedAccessions.Contains("1"),
                    $"Mod-pruned DB missing accession '1'. Actual: [{string.Join(", ", modPrunedAccessions)}]");
                Assert.That(modPrunedAccessions.Contains("2"),
                    $"Mod-pruned DB missing accession '2'. Actual: [{string.Join(", ", modPrunedAccessions)}]");

                // Optional: ensure each has at least 1 kept site (do not require exact count to avoid reference/collapse nuances)
                Assert.That(modPruned.First(p => p.Accession == "1").OneBasedPossibleLocalizedModifications.Count, Is.GreaterThanOrEqualTo(1),
                    "Mod-pruned '1' expected >=1 kept mod site");
                Assert.That(modPruned.First(p => p.Accession == "2").OneBasedPossibleLocalizedModifications.Count, Is.GreaterThanOrEqualTo(1),
                    "Mod-pruned '2' expected >=1 kept mod site");

                // Extra validation: the pruned files are non-empty and well-formed
                Assert.That(new FileInfo(modPrunedPath).Length, Is.GreaterThan(0), "mod-pruned DB file size is zero");
                Assert.That(new FileInfo(proteinPrunedPath).Length, Is.GreaterThan(0), "protein-pruned DB file size is zero");
            });

            // Cleanup
            try
            {
                if (File.Exists(rawPath)) File.Delete(rawPath);
                if (File.Exists(inputDbPath)) File.Delete(inputDbPath);
                if (Directory.Exists(inputDbFolder)) Directory.Delete(inputDbFolder, true);
                if (Directory.Exists(outputFolder)) Directory.Delete(outputFolder, true);
            }
            catch
            {
                // ignore cleanup errors (e.g., locked by IDE)
            }
        }
        [Test]
        public static void TestContaminantGPTMD()
        {
            GptmdTask task1 = new GptmdTask
            {
                CommonParameters = new CommonParameters(),
                GptmdParameters = new GptmdParameters
                {
                    ListOfModsGptmd = GlobalVariables.AllModsKnown.Where(b =>
                        b.ModificationType.Equals("Common Artifact")
                        || b.ModificationType.Equals("Common Biological")
                        || b.ModificationType.Equals("Metal")
                        || b.ModificationType.Equals("Less Common")
                        ).Select(b => (b.ModificationType, b.IdWithMotif)).ToList()
                }
            };

            List<(string, MetaMorpheusTask)> taskList = new List<(string, MetaMorpheusTask)> { ("task1", task1) };
            string mzmlName = @"TestData\PrunedDbSpectra.mzml";
            string fastaName = @"TestData\DbForPrunedDb.fasta";
            string contaminantName = @"DatabaseTests\ProteaseModTest.fasta";
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestPrunedGeneration");
            var engine = new EverythingRunnerEngine(taskList, new List<string> { mzmlName },
                new List<DbForTask> { new DbForTask(fastaName, false), new DbForTask(contaminantName, true) },
                outputFolder);
            engine.Run();
            string final = Path.Combine(MySetUpClass.outputFolder, "task1", "DbForPrunedDbGPTMD.xml");
            string contaminantFinal = Path.Combine(MySetUpClass.outputFolder, "task1", "ProteaseModTestGPTMD.xml");
            Assert.That(File.Exists(final));
            Assert.That(File.Exists(contaminantFinal));
            Directory.Delete(outputFolder, true);
        }
    }
}