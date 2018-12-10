using EngineLayer;
using MassSpectrometry;
using MzLibUtil;
using Nett;
using NUnit.Framework;
using Proteomics;
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
    public static class MyTaskTest
    {
        [Test]
        public static void TestEverythingRunner()
        {
            foreach (var modFile in Directory.GetFiles(@"Mods"))
                GlobalVariables.AddMods(PtmListLoader.ReadModsFromFile(modFile, out var fmww), false);

            CalibrationTask task1 = new CalibrationTask
            {
                CommonParameters = new CommonParameters(digestionParams: new DigestionParams(maxMissedCleavages: 0, minPeptideLength: 1, initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain)),

                CalibrationParameters = new CalibrationParameters
                {
                    WriteIntermediateFiles = true,
                    NumFragmentsNeededForEveryIdentification = 6,
                }
            };

            GptmdTask task2 = new GptmdTask
            {
                CommonParameters = new CommonParameters()
            };

            SearchTask task3 = new SearchTask
            {
                CommonParameters = new CommonParameters(),

                SearchParameters = new SearchParameters
                {
                    DoParsimony = true,
                    SearchTarget = true,
                    SearchType = SearchType.Modern
                }
            };

            SearchTask task4 = new SearchTask
            {
                CommonParameters = new CommonParameters(),

                SearchParameters = new SearchParameters
                {
                    SearchType = SearchType.Modern,
                }
            };

            List<(string, MetaMorpheusTask)> taskList = new List<(string, MetaMorpheusTask)> {
                ("task1", task1),
                ("task2", task2),
                ("task3", task3),
                ("task4", task4),};

            List<Modification> variableModifications = GlobalVariables.AllModsKnown.OfType<Modification>().Where(b => task1.CommonParameters.ListOfModsVariable.Contains((b.ModificationType, b.IdWithMotif))).ToList();
            List<Modification> fixedModifications = GlobalVariables.AllModsKnown.OfType<Modification>().Where(b => task1.CommonParameters.ListOfModsFixed.Contains((b.ModificationType, b.IdWithMotif))).ToList();

            // Generate data for files
            Protein ParentProtein = new Protein("MPEPTIDEKANTHE", "accession1");

            var digestedList = ParentProtein.Digest(task1.CommonParameters.DigestionParams, fixedModifications, variableModifications).ToList();

            Assert.AreEqual(3, digestedList.Count);

            PeptideWithSetModifications pepWithSetMods1 = digestedList[0];

            PeptideWithSetModifications pepWithSetMods2 = digestedList[2];

            var dictHere = new Dictionary<int, List<Modification>>();
            ModificationMotif.TryGetMotif("E", out ModificationMotif motif);
            dictHere.Add(3, new List<Modification> { new Modification(_originalId: "21", _target: motif, _locationRestriction: "Anywhere.", _monoisotopicMass: 21.981943) });
            Protein ParentProteinToNotInclude = new Protein("MPEPTIDEK", "accession2", "organism", new List<Tuple<string, string>>(), dictHere);
            digestedList = ParentProteinToNotInclude.Digest(task1.CommonParameters.DigestionParams, fixedModifications, variableModifications).ToList();

            MsDataFile myMsDataFile = new TestDataFile(new List<PeptideWithSetModifications> { pepWithSetMods1, pepWithSetMods2, digestedList[1] });

            Protein proteinWithChain = new Protein("MAACNNNCAA", "accession3", "organism", new List<Tuple<string, string>>(), new Dictionary<int, List<Modification>>(), new List<ProteolysisProduct> { new ProteolysisProduct(4, 8, "chain") }, "name2", "fullname2");

            string mzmlName = @"ok.mzML";
            IO.MzML.MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(myMsDataFile, mzmlName, false);
            string xmlName = "okk.xml";
            ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), new List<Protein> { ParentProtein, proteinWithChain }, xmlName);

            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestEverythingRunner");
            // RUN!
            var engine = new EverythingRunnerEngine(taskList, new List<string> { mzmlName }, new List<DbForTask> { new DbForTask(xmlName, false) }, outputFolder);
            engine.Run();
            File.Delete(Path.Combine(TestContext.CurrentContext.TestDirectory, mzmlName));
            File.Delete(Path.Combine(TestContext.CurrentContext.TestDirectory, xmlName));
            Directory.Delete(outputFolder, true);
        }

        [Test]
        public static void TestMultipleFilesRunner()
        {
            foreach (var modFile in Directory.GetFiles(@"Mods"))
                GlobalVariables.AddMods(PtmListLoader.ReadModsFromFile(modFile, out var fmww), false);

            CalibrationTask task1 = new CalibrationTask
            {
                CommonParameters = new CommonParameters
                (
                    digestionParams: new DigestionParams(maxMissedCleavages: 0, minPeptideLength: 1, initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain),
                    listOfModsVariable: new List<(string, string)> { ("Common Variable", "Oxidation on M") },
                    listOfModsFixed: new List<(string, string)> { ("Common Fixed", "Carbamidomethyl on C") },
                    productMassTolerance: new AbsoluteTolerance(0.01)
                ),
                CalibrationParameters = new CalibrationParameters
                {
                    NumFragmentsNeededForEveryIdentification = 6,
                }
            };
            GptmdTask task2 = new GptmdTask
            {
                CommonParameters = new CommonParameters
                (
                    digestionParams: new DigestionParams(),
                    productMassTolerance: new AbsoluteTolerance(0.01)
                ),
            };

            SearchTask task3 = new SearchTask
            {
                CommonParameters = new CommonParameters(),

                SearchParameters = new SearchParameters
                {
                    DoParsimony = true,
                    SearchTarget = true,
                    SearchType = SearchType.Modern,
                }
            };
            SearchTask task4 = new SearchTask
            {
                CommonParameters = new CommonParameters(),

                SearchParameters = new SearchParameters
                {
                    SearchType = SearchType.Modern,
                }
            };
            List<(string, MetaMorpheusTask)> taskList = new List<(string, MetaMorpheusTask)> {
                ("task1", task1),
                ("task2", task2),
                ("task3", task3),
                ("task4", task4),};

            List<Modification> variableModifications = GlobalVariables.AllModsKnown.OfType<Modification>().Where(b => task1.CommonParameters.ListOfModsVariable.Contains((b.ModificationType, b.IdWithMotif))).ToList();
            List<Modification> fixedModifications = GlobalVariables.AllModsKnown.OfType<Modification>().Where(b => task1.CommonParameters.ListOfModsFixed.Contains((b.ModificationType, b.IdWithMotif))).ToList();

            // Generate data for files
            Protein ParentProtein = new Protein("MPEPTIDEKANTHE", "accession1");

            var digestedList = ParentProtein.Digest(task1.CommonParameters.DigestionParams, fixedModifications, variableModifications).ToList();

            Assert.AreEqual(3, digestedList.Count);

            PeptideWithSetModifications pepWithSetMods1 = digestedList[0];

            PeptideWithSetModifications pepWithSetMods2 = digestedList[2];

            var dictHere = new Dictionary<int, List<Modification>>();
            ModificationMotif.TryGetMotif("E", out ModificationMotif motif);
            dictHere.Add(3, new List<Modification> { new Modification(_originalId: "21", _modificationType: "myModType", _target: motif, _locationRestriction: "Anywhere.", _monoisotopicMass: 21.981943) });
            Protein ParentProteinToNotInclude = new Protein("MPEPTIDEK", "accession2", "organism", new List<Tuple<string, string>>(), dictHere);
            digestedList = ParentProteinToNotInclude.Digest(task1.CommonParameters.DigestionParams, fixedModifications, variableModifications).ToList();
            Assert.AreEqual(4, digestedList.Count);

            MsDataFile myMsDataFile1 = new TestDataFile(new List<PeptideWithSetModifications> { pepWithSetMods1, pepWithSetMods2, digestedList[1] });

            string mzmlName1 = @"ok1.mzML";
            IO.MzML.MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(myMsDataFile1, mzmlName1, false);

            MsDataFile myMsDataFile2 = new TestDataFile(new List<PeptideWithSetModifications> { pepWithSetMods1, pepWithSetMods2, digestedList[1] });

            string mzmlName2 = @"ok2.mzML";
            IO.MzML.MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(myMsDataFile2, mzmlName2, false);

            Protein proteinWithChain1 = new Protein("MAACNNNCAA", "accession3", "organism", new List<Tuple<string, string>>(), new Dictionary<int, List<Modification>>(), new List<ProteolysisProduct> { new ProteolysisProduct(4, 8, "chain") }, "name2", "fullname2", false, false, new List<DatabaseReference>(), new List<SequenceVariation>(), null);
            Protein proteinWithChain2 = new Protein("MAACNNNCAA", "accession3", "organism", new List<Tuple<string, string>>(), new Dictionary<int, List<Modification>>(), new List<ProteolysisProduct> { new ProteolysisProduct(4, 8, "chain") }, "name2", "fullname2", false, false, new List<DatabaseReference>(), new List<SequenceVariation>(), null);

            string xmlName = "okk.xml";
            ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), new List<Protein> { ParentProtein, proteinWithChain1, proteinWithChain2 }, xmlName);

            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestMultipleFilesRunner");
            // RUN!
            var engine = new EverythingRunnerEngine(taskList, new List<string> { mzmlName1, mzmlName2 }, new List<DbForTask> { new DbForTask(xmlName, false) }, outputFolder);
            engine.Run();
            Directory.Delete(outputFolder, true);
            File.Delete(xmlName);
            File.Delete(mzmlName1);
            File.Delete(mzmlName2);
        }

        [Test]
        public static void MakeSureFdrDoesntSkip()
        {
            MetaMorpheusTask task = new SearchTask
            {
                CommonParameters = new CommonParameters
                (
                    digestionParams: new DigestionParams(minPeptideLength: 2),

                    scoreCutoff: 1,
                    deconvolutionIntensityRatio: 999,
                    deconvolutionMassTolerance: new PpmTolerance(50)
                ),
                SearchParameters = new SearchParameters
                {
                    DecoyType = DecoyType.None,
                    MassDiffAcceptorType = MassDiffAcceptorType.Open,
                }
            };

            string xmlName = "MakeSureFdrDoesntSkip.xml";

            {
                Protein theProtein = new Protein("MG", "accession1");
                ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), new List<Protein> { theProtein }, xmlName);
            }

            string mzmlName = @"MakeSureFdrDoesntSkip.mzML";


            var theProteins = ProteinDbLoader.LoadProteinXML(xmlName, true, DecoyType.Reverse, new List<Modification>(), false, new List<string>(), out Dictionary<string, Modification> ok);

            List<Modification> fixedModifications = new List<Modification>();

            var targetDigested = theProteins[0].Digest(task.CommonParameters.DigestionParams, fixedModifications, GlobalVariables.AllModsKnown.OfType<Modification>().ToList()).ToList();

            PeptideWithSetModifications targetGood = targetDigested.First();

            TestDataFile myMsDataFile = new TestDataFile(new List<PeptideWithSetModifications> { targetGood }, true);

            var ii = myMsDataFile.GetOneBasedScan(1).MassSpectrum.YArray.ToList();

            ii.Add(1);
            ii.Add(1);
            ii.Add(1);
            ii.Add(1);

            var intensities = ii.ToArray();

            var mm = myMsDataFile.GetOneBasedScan(1).MassSpectrum.XArray.ToList();

            var hah = 104.35352;
            mm.Add(hah);
            mm.Add(hah + 1);
            mm.Add(hah + 2);

            var mz = mm.ToArray();

            Array.Sort(mz, intensities);

            myMsDataFile.ReplaceFirstScanArrays(mz, intensities);

            IO.MzML.MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(myMsDataFile, mzmlName, false);
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestMakeSureFdrDoesntSkip");
            Directory.CreateDirectory(outputFolder);

            // RUN!
            var theStringResult = task.RunTask(outputFolder, new List<DbForTask> { new DbForTask(xmlName, false) }, new List<string> { mzmlName }, "taskId1").ToString();
            Assert.IsTrue(theStringResult.Contains("All target PSMS within 1% FDR: 1"));
            Directory.Delete(outputFolder, true);
            File.Delete(xmlName);
            File.Delete(mzmlName);
        }

        [Test]
        public static void MakeSureGptmdTaskMatchesExactMatches()
        {
            MetaMorpheusTask task1;

            {
                ModificationMotif.TryGetMotif("T", out ModificationMotif motif);
                Modification myNewMod = new Modification(_originalId: "ok", _modificationType: "okType", _target: motif, _locationRestriction: "Anywhere.", _monoisotopicMass: 229);

                GlobalVariables.AddMods(new List<Modification> { myNewMod }, false);
                task1 = new GptmdTask
                {
                    CommonParameters = new CommonParameters
                    (
                        digestionParams: new DigestionParams(initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain),
                        listOfModsVariable: new List<(string, string)>(),
                        listOfModsFixed: new List<(string, string)>(),
                        scoreCutoff: 1,
                        precursorMassTolerance: new AbsoluteTolerance(1)
                    ),

                    GptmdParameters = new GptmdParameters
                    {
                        ListOfModsGptmd = new List<(string, string)> { ("okType", "ok on T") },
                    }
                };
            }

            string xmlName = "sweetness.xml";

            {
                Protein theProtein = new Protein("MPEPTIDEKANTHE", "accession1");
                ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), new List<Protein> { theProtein }, xmlName);
            }

            string mzmlName = @"ok.mzML";

            {
                var theProteins = ProteinDbLoader.LoadProteinXML(xmlName, true, DecoyType.Reverse, new List<Modification>(), false, new List<string>(), out Dictionary<string, Modification> ok);

                List<Modification> fixedModifications = new List<Modification>();

                var targetDigested = theProteins[0].Digest(task1.CommonParameters.DigestionParams, fixedModifications, GlobalVariables.AllModsKnown.OfType<Modification>().Where(b => b.OriginalId.Equals("ok")).ToList()).ToList();

                ModificationMotif.TryGetMotif("T", out ModificationMotif motif);
                PeptideWithSetModifications targetGood = targetDigested[0];

                PeptideWithSetModifications targetWithUnknownMod = targetDigested[1];
                MsDataFile myMsDataFile = new TestDataFile(new List<PeptideWithSetModifications> { targetGood, targetWithUnknownMod }, true);

                IO.MzML.MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(myMsDataFile, mzmlName, false);
            }
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestMakeSureGptmdTaskMatchesExactMatchesTest");
            Directory.CreateDirectory(outputFolder);

            // RUN!
            var theStringResult = task1.RunTask(outputFolder, new List<DbForTask> { new DbForTask(xmlName, false) }, new List<string> { mzmlName }, "taskId1").ToString();
            Assert.IsTrue(theStringResult.Contains("Modifications added: 1"));
            Directory.Delete(outputFolder, true);
            File.Delete(xmlName);
            File.Delete(mzmlName);
            Directory.Delete(Path.Combine(TestContext.CurrentContext.TestDirectory, @"Task Settings"), true);
        }

        [Test]
        public static void TestPeptideCount()
        {
            SearchTask testPeptides = new SearchTask
            {
                CommonParameters = new CommonParameters
                (

                    digestionParams: new DigestionParams(minPeptideLength: 5)
                ),
                SearchParameters = new SearchParameters
                {
                    WritePrunedDatabase = true,
                    SearchTarget = true,
                    MassDiffAcceptorType = MassDiffAcceptorType.Exact
                }
            };

            List<(string, MetaMorpheusTask)> taskList = new List<(string, MetaMorpheusTask)>
            {
               ("TestPeptides", testPeptides)
            };

            ModificationMotif.TryGetMotif("P", out ModificationMotif motif);

            var testUniqeMod = new Modification(_originalId: "testPeptideMod", _modificationType: "mt", _target: motif, _locationRestriction: "Anywhere.", _monoisotopicMass: 10);
            GlobalVariables.AddMods(new List<Modification>
            {
                testUniqeMod
            }, false);

            //create modification lists

            List<Modification> variableModifications = GlobalVariables.AllModsKnown.OfType<Modification>().Where
                (b => testPeptides.CommonParameters.ListOfModsVariable.Contains((b.ModificationType, b.IdWithMotif))).ToList();

            //add modification to Protein object
            var modDictionary = new Dictionary<int, List<Modification>>();
            Modification modToAdd = testUniqeMod;
            modDictionary.Add(1, new List<Modification> { modToAdd });
            modDictionary.Add(3, new List<Modification> { modToAdd });

            //protein Creation (One with mod and one without)
            Protein TestProtein = new Protein("PEPTID", "accession1", "organism", new List<Tuple<string, string>>(), modDictionary);

            //First Write XML Database

            string xmlName = "singleProteinWithTwoMods.xml";

            //Add Mod to list and write XML input database
            Dictionary<string, HashSet<Tuple<int, Modification>>> modList = new Dictionary<string, HashSet<Tuple<int, Modification>>>();
            var Hash = new HashSet<Tuple<int, Modification>>
            {
                new Tuple<int, Modification>(3, modToAdd)
            };
            modList.Add("test", Hash);
            ProteinDbWriter.WriteXmlDatabase(modList, new List<Protein> { TestProtein }, xmlName);

            //now write MZML file
            var protein = ProteinDbLoader.LoadProteinXML(xmlName, true, DecoyType.Reverse, new List<Modification>(), false, new List<string>(), out Dictionary<string, Modification> ok);
            var setList1 = protein[0].Digest(testPeptides.CommonParameters.DigestionParams, new List<Modification> { }, variableModifications).ToList();
            Assert.AreEqual(4, setList1.Count);

            //Finally Write MZML file
            MsDataFile myMsDataFile = new TestDataFile(new List<PeptideWithSetModifications> { setList1[0], setList1[1], setList1[2], setList1[3], setList1[0], setList1[1] });
            string mzmlName = @"singleProteinWithRepeatedMods.mzML";
            IO.MzML.MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(myMsDataFile, mzmlName, false);

            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestMultipleFilesRunner");
            var engine = new EverythingRunnerEngine(taskList, new List<string> { mzmlName }, new List<DbForTask> { new DbForTask(xmlName, false) }, outputFolder);
            engine.Run();

            string line;

            bool foundD = false;
            using (StreamReader file = new StreamReader(Path.Combine(MySetUpClass.outputFolder, "TestPeptides", "results.txt")))
            {
                while ((line = file.ReadLine()) != null)
                {
                    if (line.Contains("All target peptides within 1% FDR: 4"))
                    {
                        foundD = true;
                    }
                }
            }
            Assert.IsTrue(foundD);
            Directory.Delete(outputFolder, true);
            File.Delete(mzmlName);
            File.Delete(xmlName);
        }

        [Test]
        public static void TestFileOutput()
        {
            string thisTaskOutputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestFileOutput");

            SearchTask task = Toml.ReadFile<SearchTask>(Path.Combine(TestContext.CurrentContext.TestDirectory, @"SlicedSearchTaskConfig.toml"), MetaMorpheusTask.tomlConfig);
            task.SearchParameters.DecoyType = DecoyType.None;

            DbForTask db = new DbForTask(Path.Combine(TestContext.CurrentContext.TestDirectory, @"sliced-db.fasta"), false);
            DbForTask db2 = new DbForTask(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", @"DbForPrunedDb.fasta"), false);
            string raw = Path.Combine(TestContext.CurrentContext.TestDirectory, @"sliced-raw.mzML");
            string raw2 = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", @"PrunedDbSpectra.mzml");
            EverythingRunnerEngine singleMassSpectraFile = new EverythingRunnerEngine(new List<(string, MetaMorpheusTask)> { ("SingleMassSpectraFileOutput", task) }, new List<string> { raw }, new List<DbForTask> { db }, thisTaskOutputFolder);
            EverythingRunnerEngine multipleMassSpectraFiles = new EverythingRunnerEngine(new List<(string, MetaMorpheusTask)> { ("MultipleMassSpectraFileOutput", task) }, new List<string> { raw, raw2 }, new List<DbForTask> { db, db2 }, thisTaskOutputFolder);

            singleMassSpectraFile.Run();
            multipleMassSpectraFiles.Run();

            // test single file output
            HashSet<string> expectedFiles = new HashSet<string> {
                "AllPeptides.psmtsv", "AllProteinGroups.tsv", "AllPSMs.psmtsv", "AllPSMs_FormattedForPercolator.tsv", "AllQuantifiedPeaks.tsv",
                "AllQuantifiedPeptides.tsv", "prose.txt", "results.txt" };

            HashSet<string> files = new HashSet<string>(Directory.GetFiles(Path.Combine(thisTaskOutputFolder, "SingleMassSpectraFileOutput")).Select(v => Path.GetFileName(v)));

            // these 2 lines are for debug purposes, so you can see which files you're missing (if any)
            var missingFiles = expectedFiles.Except(files);
            var extraFiles = files.Except(expectedFiles);

            // test that output is what's expected
            Assert.That(files.SetEquals(expectedFiles));

            // test multi file output
            files = new HashSet<string>(Directory.GetFiles(Path.Combine(thisTaskOutputFolder, "MultipleMassSpectraFileOutput")).Select(v => Path.GetFileName(v)));
            missingFiles = expectedFiles.Except(files);
            extraFiles = files.Except(expectedFiles);

            Assert.That(files.SetEquals(expectedFiles));

            expectedFiles = new HashSet<string> {
                "PrunedDbSpectra.mzID", "PrunedDbSpectra_PSMs.psmtsv", "PrunedDbSpectra_PSMsFormattedForPercolator.tsv", "PrunedDbSpectra_Peptides.psmtsv", "PrunedDbSpectra_ProteinGroups.tsv", "PrunedDbSpectra_QuantifiedPeaks.tsv",
                "sliced-raw.mzID", "sliced-raw_PSMs.psmtsv", "sliced-raw_PSMsFormattedForPercolator.tsv", "sliced-raw_Peptides.psmtsv", "sliced-raw_ProteinGroups.tsv", "sliced-raw_QuantifiedPeaks.tsv" };

            string individualFilePath = Path.Combine(thisTaskOutputFolder, "MultipleMassSpectraFileOutput", "Individual File Results");
            Assert.That(Directory.Exists(individualFilePath));

            files = new HashSet<string>(Directory.GetFiles(individualFilePath).Select(v => Path.GetFileName(v)));
            missingFiles = expectedFiles.Except(files);
            extraFiles = files.Except(expectedFiles);

            Assert.That(files.SetEquals(expectedFiles));

            files = new HashSet<string>(Directory.GetFiles(Path.Combine(thisTaskOutputFolder, "Task Settings")).Select(v => Path.GetFileName(v)));
            expectedFiles = new HashSet<string> {
                "MultipleMassSpectraFileOutputconfig.toml", "SingleMassSpectraFileOutputconfig.toml" };
            Assert.That(files.SetEquals(expectedFiles));
            Directory.Delete(thisTaskOutputFolder, true);
        }

        /// <summary>
        /// This tests for a bug in annotating mods in the search task. The situation is that if you search with a fasta database (no mods annotated),
        /// and then do GPTMD, then search with the GPTMD database, the resulting PSM will have a UniProt mod annotated on it.
        /// Also, if GPTMD has a mod with the same name as a UniProt mod, the annotated PSM will be ambiguous between
        /// the UniProt and the MetaMorpheus modification.
        /// </summary>
        [Test]
        public static void TestUniprotNamingConflicts()
        {
            // write the mod
            var outputDir = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestUniprotNamingConflicts");
            Directory.CreateDirectory(outputDir);
            string modToWrite = "Custom List\nID   Hydroxyproline\nTG   P\nPP   Anywhere.\nMT   Biological\nCF   O1\n" + @"//";
            var filePath = Path.Combine(GlobalVariables.DataDir, @"Mods", @"hydroxyproline.txt");
            File.WriteAllLines(filePath, new string[] { modToWrite });

            // read the mod
            GlobalVariables.AddMods(PtmListLoader.ReadModsFromFile(filePath, out var fmww), false);
            Assert.That(GlobalVariables.AllModsKnown.Where(v => v.IdWithMotif == "Hydroxyproline on P").Count() == 1);

            // should have an error message...
            Assert.That(GlobalVariables.ErrorsReadingMods.Where(v => v.Contains("Hydroxyproline")).Count() > 0);
            Directory.Delete(outputDir, true);
        }

        /// <summary>
        /// Tests that pepXML is written
        /// 
        /// TODO: Assert pepXML properties
        /// </summary>
        [Test]
        public static void TestPepXmlOutput()
        {
            SearchTask search = new SearchTask
            {
                SearchParameters = new SearchParameters
                {
                    WritePepXml = true
                }
            };

            List<(string, MetaMorpheusTask)> taskList = new List<(string, MetaMorpheusTask)> { ("TestPepXmlOutput", search) };

            string mzmlName = @"TestData\PrunedDbSpectra.mzml";
            string fastaName = @"TestData\DbForPrunedDb.fasta";
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestPepXmlOutput");

            var engine = new EverythingRunnerEngine(taskList, new List<string> { mzmlName }, new List<DbForTask> { new DbForTask(fastaName, false) }, outputFolder);
            engine.Run();

            string outputPepXmlPath = Path.Combine(outputFolder, @"TestPepXmlOutput\Individual File Results\PrunedDbSpectra.pep.XML");
            Assert.That(File.Exists(outputPepXmlPath));
            Directory.Delete(outputFolder, true);
        }

        [Test]
        public static void TestModernAndClassicSearch()
        {
            SearchTask classicSearch = new SearchTask();

            SearchTask modernSearch = new SearchTask
            {
                SearchParameters = new SearchParameters
                {
                    SearchType = SearchType.Modern
                }
            };
            List<int> counts = new List<int>();

            List<(string, MetaMorpheusTask)> taskList = new List<(string, MetaMorpheusTask)> { ("ClassicSearch", classicSearch), ("ModernSearch", modernSearch) };

            string mzmlName = @"TestData\PrunedDbSpectra.mzml";
            string fastaName = @"TestData\DbForPrunedDb.fasta";
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestPepXmlOutput");

            var engine = new EverythingRunnerEngine(taskList, new List<string> { mzmlName }, new List<DbForTask> { new DbForTask(fastaName, false) }, outputFolder);
            engine.Run();

            string classicPath = Path.Combine(outputFolder, @"ClassicSearch\AllPSMs.psmtsv");
            var classicPsms = File.ReadAllLines(classicPath).ToList();

            string modernPath = Path.Combine(outputFolder, @"ModernSearch\AllPSMs.psmtsv");
            var modernPsms = File.ReadAllLines(modernPath).ToList();
            counts.Add(modernPsms.Count);

            Assert.That(modernPsms.SequenceEqual(classicPsms));
            Directory.Delete(outputFolder, true);
        }
    }
}
