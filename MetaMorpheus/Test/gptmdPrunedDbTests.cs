using EngineLayer;
using MassSpectrometry;
using NUnit.Framework;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using Omics.Modifications;
using TaskLayer;
using UsefulProteomicsDatabases;
using Omics;
using System.Threading.Tasks;
using Transcriptomics;
using UsefulProteomicsDatabases.Transcriptomics;
using System.Windows.Automation.Peers;

namespace Test
{
    [TestFixture]
    public static class GptmdPrunedDbTests
    {

        private static string _testingDirectory;
        [OneTimeSetUp]
        public static void OneTimeSetUp()
        {
            _testingDirectory = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestPrunedGenerationWriter");
            if (Directory.Exists(_testingDirectory))
                Directory.Delete(_testingDirectory, true);
            Directory.CreateDirectory(_testingDirectory);
        }

        [OneTimeTearDown]
        public static void OneTimeTearDown()
        {
            Directory.Delete(_testingDirectory, true);
        }


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
            Assert.That(totalNumberOfMods, Is.EqualTo(8));
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
            Dictionary<string, HashSet<Tuple<int, Modification>>> modList = new Dictionary<string, HashSet<Tuple<int, Modification>>>();
            var Hash = new HashSet<Tuple<int, Modification>>
            {
                new Tuple<int, Modification>(3, modToAdd)
            };
            modList.Add("test", Hash);
            ProteinDbWriter.WriteXmlDatabase(modList, new List<Protein> { TestProteinWithMod }, xmlName);

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
            PeptideWithSetModifications pepWithSetMods1 = digestedList[1];

            //Finally Write MZML file
            Assert.That(pepWithSetMods1.FullSequence, Is.EqualTo("PEP[ConnorModType:ConnorMod on P]TID"));//this might be base sequence
            MsDataFile myMsDataFile = new TestDataFile(new List<PeptideWithSetModifications> { pepWithSetMods1 });
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
            ProteinDbWriter.WriteXmlDatabase(modList, new List<Protein> { TestProteinWithModForDB }, xmlName);

            //Add Observed Only
            modList.Add("test2", Hash);
            ProteinDbWriter.WriteXmlDatabase(modList, new List<Protein> { TestProteinWithModObsevred }, xmlName2);

            //now create MZML data
            var protein = ProteinDbLoader.LoadProteinXML(xmlName2, true, DecoyType.Reverse, new List<Modification>(), false, new List<string>(), out Dictionary<string, Modification> ok);
            var digestedList = protein[0].Digest(task5.CommonParameters.DigestionParams, fixedModifications, variableModifications).ToList();

            //Set Peptide with 1 mod at position 3
            PeptideWithSetModifications pepWithSetMods1 = digestedList[0];
            PeptideWithSetModifications pepWithSetMods2 = digestedList[1];
            PeptideWithSetModifications pepWithSetMods3 = digestedList[2];
            PeptideWithSetModifications pepWithSetMods4 = digestedList[3];
            PeptideWithSetModifications pepWithSetMods5 = digestedList[4];

            //CUSTOM PEP
            MsDataFile myMsDataFile = new TestDataFile(new List<PeptideWithSetModifications>
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
            var Dlist = proteins[0].GetVariantProteins().SelectMany(vp => vp.Digest(task5.CommonParameters.DigestionParams, fixedModifications, variableModifications)).ToList();
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

            var db = ProteinDbWriter.WriteXmlDatabase(modList, new List<Protein> { protein1, protein2 }, Path.Combine(TestContext.CurrentContext.TestDirectory, @"PrunedDbTest/fakeDb.xml"));

            var peptideObserved = protein1.Digest(new DigestionParams(minPeptideLength: 1), new List<Modification>(), new List<Modification>())
            .Where(p => p.BaseSequence == "PEPT" && p.AllModsOneIsNterminus.Count > 0).First();
            PostSearchAnalysisParameters testPostTaskParameters = new PostSearchAnalysisParameters();
            CommonParameters commonParam = new CommonParameters();
            double[,] noiseData = new double[10000, 10000];
            noiseData[0,0] = 1.0; 
            List<Omics.Fragmentation.MatchedFragmentIon> matchedFragmentIons = new List<Omics.Fragmentation.MatchedFragmentIon>() { };
            MzSpectrum spectrum = new MzSpectrum(noiseData);
            MsDataScan scan = new MsDataScan(spectrum , 1, 1, true, Polarity.Unknown, 2, new MzLibUtil.MzRange(10, 1000), "", MZAnalyzerType.Orbitrap, 10000, null, noiseData, "");
            testPostTaskParameters.ProteinList = new List<Protein>() { protein1, protein2 };
            testPostTaskParameters.AllPsms = new List<SpectralMatch> { new PeptideSpectralMatch(peptideObserved, 0, 20, 1, new Ms2ScanWithSpecificMass(scan, 100, 1, @"", commonParam), commonParam, matchedFragmentIons) };
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

            MsDataFile myMsDataFile = new TestDataFile(new List<PeptideWithSetModifications>
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
            var modToWrite = GlobalVariables.AllModsKnown.Where(p => p.ModificationType == "UniProt" && p.Target.ToString() == "T").First();
            var modToNotWrite = GlobalVariables.AllModsKnown.Where(p => p.ModificationType == "Common Artifact" && p.Target.ToString() == "X").First();
            Dictionary<int, List<Modification>> variantMods = new Dictionary<int, List<Modification>>();
            variantMods.Add(1, new List<Modification>(){ modToNotWrite});

            List<SequenceVariation> variants = new List<SequenceVariation> { new SequenceVariation(4, 4, "V", "T", @"20\t41168825\t.\tT\tC\t14290.77\t.\tANN=C|missense_variant|MODERATE|PLCG1|ENSG00000124181|transcript|ENST00000244007.7|protein_coding|22/33|c.2438T>C|p.Ile813Thr|2635/5285|2438/3876|813/1291||\tGT:AD:DP:GQ:PL\t1/1:1,392:393:99:14319,1142,0", variantMods) };

            var protein1 = new Protein("PEPVIDEKPEPT", "1", oneBasedModifications: new Dictionary<int, List<Modification>> { { 1, new List<Modification> { modToNotWrite } }, { 12, new List<Modification> { modToWrite } } }, sequenceVariations: variants);
            var protein2 = new Protein("PEPIDPEPT", "2", oneBasedModifications: new Dictionary<int, List<Modification>> { { 1, new List<Modification> { modToNotWrite } }, { 9, new List<Modification> { modToWrite } } });
            var protein1Variants = protein1.GetVariantProteins(1,0);           
            
            string path = @"temp";

            var proteinList = new List<Protein> { protein1, protein2};
            proteinList.AddRange(protein1Variants);
          

            ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), proteinList, path);

            Directory.CreateDirectory(Path.Combine(TestContext.CurrentContext.TestDirectory, @"PrunedDbTestVariant"));

            Dictionary<string, HashSet<Tuple<int, Modification>>> modList = new Dictionary<string, HashSet<Tuple<int, Modification>>>();
            var Hash = new HashSet<Tuple<int, Modification>>
            {
                new Tuple<int, Modification>(1, modToWrite),
                new Tuple<int, Modification>(2, modToNotWrite),

            };

            var db = ProteinDbWriter.WriteXmlDatabase(modList,  proteinList , Path.Combine(TestContext.CurrentContext.TestDirectory, @"PrunedDbTestVariant/fakeDb.xml"));

            var peptideObserved = protein1Variants.First().Digest(new DigestionParams(minPeptideLength: 1), new List<Modification>(), new List<Modification>())
            .Where(p => p.BaseSequence == "PEPT").First();
            PostSearchAnalysisParameters testPostTaskParameters = new PostSearchAnalysisParameters();
            CommonParameters commonParam = new CommonParameters();
            double[,] noiseData = new double[10000, 10000];
            noiseData[0, 0] = 1.0;
            List<Omics.Fragmentation.MatchedFragmentIon> matchedFragmentIons = new List<Omics.Fragmentation.MatchedFragmentIon>() { };
            MzSpectrum spectrum = new MzSpectrum(noiseData);
            MsDataScan scan = new MsDataScan(spectrum, 1, 1, true, Polarity.Unknown, 2, new MzLibUtil.MzRange(10, 1000), "", MZAnalyzerType.Orbitrap, 10000, null, noiseData, "");
            testPostTaskParameters.ProteinList = proteinList;
            testPostTaskParameters.AllPsms = new List<SpectralMatch> { new PeptideSpectralMatch(peptideObserved, 0, 20, 1, new Ms2ScanWithSpecificMass(scan, 100, 1, @"", commonParam), commonParam, matchedFragmentIons) };
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

            MsDataFile myMsDataFile = new TestDataFile(new List<PeptideWithSetModifications>
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

            var test = task5.RunTask(Path.Combine(TestContext.CurrentContext.TestDirectory, @"PrunedDbTest"), new List<DbForTask>() { new DbForTask(Path.Combine(TestContext.CurrentContext.TestDirectory, @"PrunedDbTest/fakeDb.xml"), false) }, new List<string>() { mzmlName }, "name");
            testPostTaskParameters.SearchTaskResults = test;

            PostSearchAnalysisTask testPostTask = new PostSearchAnalysisTask();
            testPostTask.Parameters = testPostTaskParameters;
            testPostTask.CommonParameters = commonParam;
            testPostTask.FileSpecificParameters = new List<(string FileName, CommonParameters Parameters)> { ("newMzMl.mzml", commonParam) };
            testPostTask.Run();

            var proteinsLoaded = ProteinDbLoader.LoadProteinXML(path, true, DecoyType.None, GlobalVariables.AllModsKnown, false, new List<string>(), out var unknownMods);

            // assert that mods on proteins are the same before/after task is run            
            Assert.That(protein1Variants.First().Accession, Is.EqualTo(proteinsLoaded.First().Accession));
            Assert.That(protein1Variants.First().OneBasedPossibleLocalizedModifications.Count(), Is.EqualTo(proteinsLoaded.First().OneBasedPossibleLocalizedModifications.Count()));
            Assert.That(protein2.OneBasedPossibleLocalizedModifications.Count(), Is.EqualTo(proteinsLoaded.ElementAt(1).OneBasedPossibleLocalizedModifications.Count()));

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
        [NonParallelizable]
        public static async Task TestWriteDataAsync_Protein()
        {
            string fastaPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "DbForPrunedDb.fasta");
            var proteins = ProteinDbLoader.LoadProteinFasta(fastaPath, true, DecoyType.None,  false,  out _);
            string outputPath = Path.Combine(_testingDirectory, "test_output_protein.xml");

            bool eventTriggered = false;
            var eventListener = new EventHandler<SingleFileEventArgs>((sender, e) => 
            {
                eventTriggered = true;
                Assert.That(e, Is.TypeOf<SingleFileEventArgs>());
                Assert.That(e.WrittenFile, Is.EqualTo(outputPath));
            });

            PrunedDatabaseWriter.FinishedWritingFileHandler += eventListener;

            await PrunedDatabaseWriter.WriteDataAsync(outputPath, proteins);
            Assert.That(eventTriggered, Is.True);

            PrunedDatabaseWriter.FinishedWritingFileHandler -= eventListener;
        }

        [Test]
        [NonParallelizable]
        public static async Task TestWriteDataAsync_RNA()
        {
            string fastaPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "Transcriptomics", "TestData", "ModomicsUnmodifiedTrimmed.fasta");
            var rnas = RnaDbLoader.LoadRnaFasta(fastaPath, true, DecoyType.None, false, out _);
            string outputPath = Path.Combine(_testingDirectory, "test_output_rna.xml");

            bool eventTriggered = false;
            var eventListener = new EventHandler<SingleFileEventArgs>((sender, e) =>
            {
                eventTriggered = true;
                Assert.That(e, Is.TypeOf<SingleFileEventArgs>());
                Assert.That(e.WrittenFile, Is.EqualTo(outputPath));
            });

            PrunedDatabaseWriter.FinishedWritingFileHandler += eventListener;

            await PrunedDatabaseWriter.WriteDataAsync(outputPath, rnas);
            Assert.That(eventTriggered, Is.True);

            PrunedDatabaseWriter.FinishedWritingFileHandler -= eventListener;
        }

        [Test]
        public static void TestWriteData_InvalidType()
        {
            var bioPolymers = new List<IBioPolymer> { new Protein("MNNNKQQQ", null), new RNA("AUGCUA", null) };
            string outputPath = "test_output_invalid.xml";

            Assert.Throws<ArgumentException>(() => PrunedDatabaseWriter.WriteData(outputPath, bioPolymers));
        }

        [Test]
        public static void TestGetModificationsToWrite()
        {
            var modsToWrite = new Dictionary<string, int>
            {
                { "Less Common", 1},
                { "UniProt", 2 },
                { "Common Biological", 3 }
            };

            var (modificationsToWriteIfBoth, modificationsToWriteIfInDatabase, modificationsToWriteIfObserved) 
                = PrunedDatabaseWriter.GetModificationsToWrite(modsToWrite);

            Assert.That(modificationsToWriteIfBoth, Is.Not.Null);
            Assert.That(modificationsToWriteIfInDatabase, Is.Not.Null);
            Assert.That(modificationsToWriteIfObserved, Is.Not.Null);

            // Assuming GlobalVariables.AllModsKnown contains modifications of types "modType1", "modType2", and "modType3"
            Assert.That(modificationsToWriteIfBoth.Any(mod => mod.ModificationType == "Less Common"), Is.True);
            Assert.That(modificationsToWriteIfInDatabase.Any(mod => mod.ModificationType == "UniProt"), Is.True);
            Assert.That(modificationsToWriteIfObserved.Any(mod => mod.ModificationType == "Common Biological"), Is.True);
        }


        private static Protein _testProtein = new Protein("PEPTIDEKPEPTK", "1");

        private class TestSpectralMatch : SpectralMatch
        {
            // Construct the bare minimum of our objects needed to run this section of PostSearchAnalysisTask
            public TestSpectralMatch(string fullSequence, int startResidue = 0, int endResidue = 13,
                string filePath = "default", string protease = "trypsin")
                : base(
                    new PeptideWithSetModifications(fullSequence, GlobalVariables.AllModsKnownDictionary,
                        p: _testProtein, oneBasedStartResidueInProtein: startResidue,
                        oneBasedEndResidueInProtein: endResidue, digestionParams: new DigestionParams(protease)),
                    0, 10, 0,
                    new Ms2ScanWithSpecificMass(
                        new MsDataScan(
                            new MzSpectrum([], [], false), 0, 0,
                            true, Polarity.Positive, 0, default, "", 0,
                            0, 0,
                            new double[0, 0], ""), 0, 0, filePath, new CommonParameters(), []),
                    new CommonParameters(), [])
            {
                FdrInfo = new EngineLayer.FdrAnalysis.FdrInfo
                {
                    PEP = 0,
                    PEP_QValue = 0,
                    QValue = 0,
                    QValueNotch = 0
                };
                ResolveAllAmbiguities();

                // adjust modification index 
                foreach (var bioPolymer in BestMatchingBioPolymersWithSetMods)
                {
                    var allMods = bioPolymer.Peptide.AllModsOneIsNterminus.ToArray();
                    bioPolymer.Peptide.AllModsOneIsNterminus.Clear();
                    foreach (var modification in allMods)
                    {
                        bioPolymer.Peptide.AllModsOneIsNterminus.Add(modification.Key + startResidue-1, modification.Value);
                    }
                }
            }
        };

        [Test]
        public static void GetProteinToConfidentModifiedSequences_SingleEvidence_OneMod()
        {
            // Arrange
            var psms = new List<SpectralMatch>
            {
                new TestSpectralMatch("PEPTIDE[Common Biological:Carboxylation on E]K", 0 , 8)
            };
            var postSearchAnalysisTask = new PostSearchAnalysisTask()
            {
                CommonParameters = new(),
                FileSpecificParameters = new()
            };

            // Act
            var result = postSearchAnalysisTask.GetProteinToConfidentModifiedSequences(psms);

            // Assert
            Assert.That(result, Is.Not.Null);
            Assert.That(result.Count, Is.EqualTo(1));
            Assert.That(result[_testProtein].Count, Is.EqualTo(1));
        }

        [Test]
        public static void GetProteinToConfidentModifiedSequences_SingleEvidence_TwoMod()
        {
            // Arrange
            var psms = new List<SpectralMatch>
            {
                new TestSpectralMatch("PEPTIDE[Common Biological:Carboxylation on E]K", 0 , 8),
                new TestSpectralMatch("PE[Common Biological:Carboxylation on E]PTIDEK", 0 , 8)
            };
            var postSearchAnalysisTask = new PostSearchAnalysisTask()
            {
                CommonParameters = new(),
                FileSpecificParameters = new()
            };

            // Act
            var result = postSearchAnalysisTask.GetProteinToConfidentModifiedSequences(psms);

            // Assert
            Assert.That(result, Is.Not.Null);
            Assert.That(result.Count, Is.EqualTo(1));
            Assert.That(result[_testProtein].Count, Is.EqualTo(2));
        }

        [Test]
        public static void GetProteinToConfidentModifiedSequences_MultipleEvidence_AllFail()
        {
            // Arrange
            var psms = new List<SpectralMatch>
            {
                new TestSpectralMatch("PEPTIDE[Common Biological:Carboxylation on E]K", 0 , 8),
                new TestSpectralMatch("PE[Common Biological:Carboxylation on E]PTIDEK", 0 , 8),
                new TestSpectralMatch("PEPTIDEKPE[Common Biological:Carboxylation on E]PTK", 0 , 13),
                new TestSpectralMatch("PEPTIDEKPE[Common Biological:Carboxylation on E]PTK", 0 , 13)
            };
            var postSearchAnalysisTask = new PostSearchAnalysisTask()
            {
                CommonParameters = new(),
                FileSpecificParameters = new()
            };

            // Act
            var result = postSearchAnalysisTask.GetProteinToConfidentModifiedSequences(psms, 2);

            // Assert
            Assert.That(result, Is.Not.Null);
            Assert.That(result.Count, Is.EqualTo(1));
            Assert.That(result[_testProtein], Is.Empty);
        }

        [Test]
        public static void GetProteinToConfidentModifiedSequences_MultipleEvidence_MissedCleavage_OneSharedModSimple()
        {
            // Arrange
            var psms = new List<SpectralMatch>
            {
                new TestSpectralMatch("PEPTIDE[Common Biological:Carboxylation on E]K", 0 , 8),
                new TestSpectralMatch("PE[Common Biological:Carboxylation on E]PTIDEK", 0 , 8),
                new TestSpectralMatch("PEPTIDEKPE[Common Biological:Carboxylation on E]PTK", 0 , 13),
                new TestSpectralMatch("PE[Common Biological:Carboxylation on E]PTIDEKPEPTK", 0 , 13)
            };
            var postSearchAnalysisTask = new PostSearchAnalysisTask()
            {
                CommonParameters = new(),
                FileSpecificParameters = new()
            };

            // Act
            var result = postSearchAnalysisTask.GetProteinToConfidentModifiedSequences(psms, 2);

            // Assert
            Assert.That(result, Is.Not.Null);
            Assert.That(result.Count, Is.EqualTo(1));
            Assert.That(result[_testProtein].Count, Is.EqualTo(1));

            string[] fullSequences = result.SelectMany(p => p.Value.Select(p => p.FullSequence)).ToArray();
            Assert.That(fullSequences, Does.Contain("PE[Common Biological:Carboxylation on E]PTIDEK"));
        }

        [Test]
        public static void GetProteinToConfidentModifiedSequences_MultipleEvidence_MissedCleavage_TwoSharedMod_AlternatingTerm()
        {
            // Arrange
            var psms = new List<SpectralMatch>
            {
                new TestSpectralMatch("PEPT[Common Biological:Phosphorylation on T]K", 8, 13),
                new TestSpectralMatch("PE[Common Biological:Carboxylation on E]PTIDEK", 0 , 8),
                new TestSpectralMatch("PEPTIDEKPE[Common Biological:Carboxylation on E]PTK", 0 , 13),

                // both of the below satisfy the criteria of covering all modifications, but only one should be selected
                new TestSpectralMatch("PE[Common Biological:Carboxylation on E]PTIDEKPEPT[Common Biological:Phosphorylation on T]K", 0 , 13),
                new TestSpectralMatch("PE[Common Biological:Carboxylation on E]PTIDEKPE[Common Biological:Carboxylation on E]PT[Common Biological:Phosphorylation on T]K", 0 , 13)
            };
            var postSearchAnalysisTask = new PostSearchAnalysisTask()
            {
                CommonParameters = new(),
                FileSpecificParameters = new()
            };

            // Act
            var result = postSearchAnalysisTask.GetProteinToConfidentModifiedSequences(psms, 2);

            // Assert
            Assert.That(result, Is.Not.Null);
            Assert.That(result.Count, Is.EqualTo(1));
            Assert.That(result[_testProtein].Count, Is.EqualTo(1));

            string[] fullSequences = result.SelectMany(p => p.Value.Select(p => p.FullSequence)).ToArray();
            Assert.That(fullSequences, Does.Contain("PE[Common Biological:Carboxylation on E]PTIDEKPEPT[Common Biological:Phosphorylation on T]K"));
        }

        [Test]
        public static void GetProteinToConfidentModifiedSequences_MultipleEvidence_DifferentDissociationTypes()
        {
            // Arrange
            var psms = new List<SpectralMatch>
            { 
                // two mods found in the same condition, one mod of those found in a different contdition
                new TestSpectralMatch("PE[Common Biological:Carboxylation on E]PT[Common Biological:Phosphorylation on T]IDEK", 0, 8, "hcd"),
                new TestSpectralMatch("PE[Common Biological:Carboxylation on E]PT[Common Biological:Phosphorylation on T]IDEK", 0, 8, "hcd"),
                new TestSpectralMatch("PE[Common Biological:Carboxylation on E]PTIDEK", 0, 8, "hcd"),
                new TestSpectralMatch("PE[Common Biological:Carboxylation on E]PTIDEK", 0 , 8, "etd"),
                new TestSpectralMatch("PEPTIDEK", 0 , 8),
            };
            var postSearchAnalysisTask = new PostSearchAnalysisTask()
            {
                CommonParameters = new(),
                FileSpecificParameters = new()
                {
                    ("hcd", new CommonParameters(dissociationType: DissociationType.HCD)),
                    ("etd", new CommonParameters(dissociationType: DissociationType.ETD))
                },
            };

            // Act
            var result = postSearchAnalysisTask.GetProteinToConfidentModifiedSequences(psms, 2);

            // Assert
            Assert.That(result, Is.Not.Null);
            Assert.That(result.Count, Is.EqualTo(1));
            Assert.That(result[_testProtein].Count, Is.EqualTo(1));

            string[] fullSequences = result.SelectMany(p => p.Value.Select(p => p.FullSequence)).ToArray();
            Assert.That(fullSequences, Does.Contain("PE[Common Biological:Carboxylation on E]PTIDEK"));
        }

        [Test]
        public static void GetProteinToConfidentModifiedSequences_MultipleEvidence_DifferentDissociationTypes_MoreMods()
        {
            // Arrange
            var psms = new List<SpectralMatch>
            { 
                // three mods found in the same condition, one mod of those found in a different condition, two found in the same condition in a missed cleavage product
                new TestSpectralMatch("PE[Common Biological:Carboxylation on E]PT[Common Biological:Phosphorylation on T]IDE[Common Biological:Carboxylation on E]K", 0, 8, "hcd"),
                new TestSpectralMatch("PE[Common Biological:Carboxylation on E]PT[Common Biological:Phosphorylation on T]IDE[Common Biological:Carboxylation on E]KPEPTK", 0, 13, "hcd"),
                new TestSpectralMatch("PE[Common Biological:Carboxylation on E]PTIDEK", 0, 8, "hcd"),
                new TestSpectralMatch("PE[Common Biological:Carboxylation on E]PTIDEK", 0 , 8, "etd"),
                new TestSpectralMatch("PEPTIDEK", 0 , 8),
            };
            var postSearchAnalysisTask = new PostSearchAnalysisTask()
            {
                CommonParameters = new(),
                FileSpecificParameters = new()
                {
                    ("hcd", new CommonParameters(dissociationType: DissociationType.HCD)),
                    ("etd", new CommonParameters(dissociationType: DissociationType.ETD))
                },
            };

            // Act
            var result = postSearchAnalysisTask.GetProteinToConfidentModifiedSequences(psms, 2);

            // Assert
            Assert.That(result, Is.Not.Null);
            Assert.That(result.Count, Is.EqualTo(1));
            Assert.That(result[_testProtein].Count, Is.EqualTo(1));

            string[] fullSequences = result.SelectMany(p => p.Value.Select(p => p.FullSequence)).ToArray();
            Assert.That(fullSequences, Does.Contain("PE[Common Biological:Carboxylation on E]PT[Common Biological:Phosphorylation on T]IDE[Common Biological:Carboxylation on E]K"));
        }

        //[Test]
        //public static void GetProteinToConfidentModifiedSequences_MultipleEvidence_DifferentDigestionAgents()
        //{
        //    // Arrange
        //    var psms = new List<SpectralMatch>
        //    {
        //        // Add test SpectralMatch objects with different digestion agents here
        //    };
        //    var postSearchAnalysisTask = new PostSearchAnalysisTask();

        //    // Act
        //    var result = postSearchAnalysisTask.GetProteinToConfidentModifiedSequences(psms, 1);

        //    // Assert
        //    Assert.That(result, Is.Not.Null);
        //    // Add more assertions to verify the correctness of the result
        //}
    }
}