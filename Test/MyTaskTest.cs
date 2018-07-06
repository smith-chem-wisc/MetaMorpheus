using EngineLayer;
using MassSpectrometry;
using MzLibUtil;
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
            // Setup tasks

            foreach (var modFile in Directory.GetFiles(@"Mods"))
            {
                GlobalVariables.AddMods(PtmListLoader.ReadModsFromFile(modFile));
            }

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
            List<(string, MetaMorpheusTask)> taskList = new List<(string, MetaMorpheusTask)>
            {
                ("task1", task1),
                ("task2", task2),
                ("task3", task3),
                ("task4", task4),
            };

            List<ModificationWithMass> variableModifications = GlobalVariables.AllModsKnown.OfType<ModificationWithMass>().Where(b => task1.CommonParameters.ListOfModsVariable.Contains((b.modificationType, b.id))).ToList();
            List<ModificationWithMass> fixedModifications = GlobalVariables.AllModsKnown.OfType<ModificationWithMass>().Where(b => task1.CommonParameters.ListOfModsFixed.Contains((b.modificationType, b.id))).ToList();

            // Generate data for files
            Protein ParentProtein = new Protein("MPEPTIDEKANTHE", "accession1");
            var digestedList = ParentProtein.Digest(task1.CommonParameters.DigestionParams, fixedModifications, variableModifications).ToList();
            Assert.AreEqual(3, digestedList.Count);

            PeptideWithSetModifications pepWithSetMods1 = digestedList[0];
            PeptideWithSetModifications pepWithSetMods2 = digestedList[2];

            var dictHere = new Dictionary<int, List<Modification>>();
            ModificationMotif.TryGetMotif("E", out ModificationMotif motif);
            dictHere.Add(3, new List<Modification> { new ModificationWithMass("21", null, motif, TerminusLocalization.Any, 21.981943) });
            Protein ParentProteinToNotInclude = new Protein("MPEPTIDEK", "accession2", "organism", new List<Tuple<string, string>>(), dictHere);
            digestedList = ParentProteinToNotInclude.Digest(task1.CommonParameters.DigestionParams, fixedModifications, variableModifications).ToList();

            MsDataFile myMsDataFile = new TestDataFile(new List<PeptideWithSetModifications> { pepWithSetMods1, pepWithSetMods2, digestedList[1] });

            Protein proteinWithChain = new Protein("MAACNNNCAA", "accession3", "organism", new List<Tuple<string, string>>(), new Dictionary<int, List<Modification>>(), new List<ProteolysisProduct> { new ProteolysisProduct(4, 8, "chain") }, "name2", "fullname2");

            // Write the files

            string mzmlName = @"ok.mzML";
            IO.MzML.MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(myMsDataFile, mzmlName, false);
            string xmlName = "okk.xml";
            ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), new List<Protein> { ParentProtein, proteinWithChain }, xmlName);

            // RUN!

            var engine = new EverythingRunnerEngine(taskList, new List<string> { mzmlName }, new List<DbForTask> { new DbForTask(xmlName, false) }, Environment.CurrentDirectory);
            engine.Run();
        }

        [Test]
        public static void TestMultipleFilesRunner()
        {
            // Setup tasks

            foreach (var modFile in Directory.GetFiles(@"Mods"))
            {
                GlobalVariables.AddMods(PtmListLoader.ReadModsFromFile(modFile));
            }

            CalibrationTask task1 = new CalibrationTask
            {
                CommonParameters = new CommonParameters
                (
                    digestionParams: new DigestionParams(maxMissedCleavages: 0, minPeptideLength: 1, initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain),
                    listOfModsVariable: new List<(string, string)> { ("Common Variable", "Oxidation of M") },
                    listOfModsFixed: new List<(string, string)> { ("Common Fixed", "Carbamidomethyl of C") },
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
            List<(string, MetaMorpheusTask)> taskList = new List<(string, MetaMorpheusTask)>
            {
                ("task1", task1),
                ("task2", task2),
                ("task3", task3),
                ("task4", task4),
            };
            List<ModificationWithMass> variableModifications = GlobalVariables.AllModsKnown.OfType<ModificationWithMass>().Where(b => task1.CommonParameters.ListOfModsVariable.Contains((b.modificationType, b.id))).ToList();
            List<ModificationWithMass> fixedModifications = GlobalVariables.AllModsKnown.OfType<ModificationWithMass>().Where(b => task1.CommonParameters.ListOfModsFixed.Contains((b.modificationType, b.id))).ToList();

            // Generate data for files
            Protein ParentProtein = new Protein("MPEPTIDEKANTHE", "accession1");

            var digestedList = ParentProtein.Digest(task1.CommonParameters.DigestionParams, fixedModifications, variableModifications).ToList();

            Assert.AreEqual(3, digestedList.Count);

            PeptideWithSetModifications pepWithSetMods1 = digestedList[0];

            PeptideWithSetModifications pepWithSetMods2 = digestedList[2];

            var dictHere = new Dictionary<int, List<Modification>>();
            ModificationMotif.TryGetMotif("E", out ModificationMotif motif);
            dictHere.Add(3, new List<Modification> { new ModificationWithMass("21", null, motif, TerminusLocalization.Any, 21.981943) });
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

            // RUN!

            var engine = new EverythingRunnerEngine(taskList, new List<string> { mzmlName1, mzmlName2 }, new List<DbForTask> { new DbForTask(xmlName, false) }, Environment.CurrentDirectory);
            engine.Run();
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

            // Generate protein and write to file

            Protein theProtein = new Protein("MG", "accession1");
            ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), new List<Protein> { theProtein }, xmlName);

            // Generate and write the mzml

            string mzmlName = @"MakeSureFdrDoesntSkip.mzML";
            var theProteins = ProteinDbLoader.LoadProteinXML(xmlName, true, DecoyType.Reverse, new List<Modification>(), false, new List<string>(), out Dictionary<string, Modification> ok);
            var targetDigested = theProteins[0].Digest(task.CommonParameters.DigestionParams, new List<ModificationWithMass>(), GlobalVariables.AllModsKnown.OfType<ModificationWithMass>().ToList()).ToList();
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

            // RUN!

            var theStringResult = task.RunTask(TestContext.CurrentContext.TestDirectory, new List<DbForTask> { new DbForTask(xmlName, false) }, new List<string> { mzmlName }, "taskId1").ToString();
            Assert.IsTrue(theStringResult.Contains("All target PSMS within 1% FDR: 1"));
        }

        [Test]
        public static void MakeSureGptmdTaskMatchesExactMatches()
        {
            MetaMorpheusTask task1;

            // Setup tasks

            ModificationMotif.TryGetMotif("T", out ModificationMotif motif);
            GlobalVariables.AddMods(new List<ModificationWithMass> { new ModificationWithMass("ok", "okType", motif, TerminusLocalization.Any, 229) });
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
                    ListOfModsGptmd = new List<(string, string)> { ("okType", "ok") },
                }
            };

            // Generate protein and write to file
            string xmlName = "sweetness.xml";
            Protein theProtein = new Protein("MPEPTIDEKANTHE", "accession1");
            ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), new List<Protein> { theProtein }, xmlName);

            // Generate protein and write to file
            string mzmlName = @"ok.mzML";
            var theProteins = ProteinDbLoader.LoadProteinXML(xmlName, true, DecoyType.Reverse, new List<Modification>(), false, new List<string>(), out var ok);
            var targetDigested = theProteins[0].Digest(task1.CommonParameters.DigestionParams, new List<ModificationWithMass>(),
                GlobalVariables.AllModsKnown.OfType<ModificationWithMass>().Where(b => b.id.Equals("ok")).ToList()).ToList();

            PeptideWithSetModifications targetGood = targetDigested[0];
            PeptideWithSetModifications targetWithUnknownMod = targetDigested[1];
            MsDataFile myMsDataFile = new TestDataFile(new List<PeptideWithSetModifications> { targetGood, targetWithUnknownMod }, true);

            IO.MzML.MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(myMsDataFile, mzmlName, false);

            // RUN!
            var theStringResult = task1.RunTask(TestContext.CurrentContext.TestDirectory, new List<DbForTask> { new DbForTask(xmlName, false) }, new List<string> { mzmlName }, "taskId1").ToString();
            Assert.IsTrue(theStringResult.Contains("Modifications added: 1"));
        }

        [Test]
        public static void TestPeptideCount()
        {
            // Setup

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

            var testUniqeMod = new ModificationWithMass("testPeptideMod", "mt", motif, TerminusLocalization.Any, 10);
            GlobalVariables.AddMods(new List<ModificationWithLocation>
            {
                testUniqeMod
            });

            // Mod setup and protein creation: protein Creation (One with mod and one without)

            List<ModificationWithMass> variableModifications = GlobalVariables.AllModsKnown.OfType<ModificationWithMass>().Where
                (b => testPeptides.CommonParameters.ListOfModsVariable.Contains((b.modificationType, b.id))).ToList();
            var modDictionary = new Dictionary<int, List<Modification>>();
            ModificationWithMass modToAdd = testUniqeMod;
            modDictionary.Add(1, new List<Modification> { modToAdd });
            modDictionary.Add(3, new List<Modification> { modToAdd });
            Protein TestProtein = new Protein("PEPTID", "accession1", "organism", new List<Tuple<string, string>>(), modDictionary);

            //First Write XML Database
            string xmlName = "singleProteinWithTwoMods.xml";
            Dictionary<string, HashSet<Tuple<int, Modification>>> modList = new Dictionary<string, HashSet<Tuple<int, Modification>>>();
            var Hash = new HashSet<Tuple<int, Modification>> { new Tuple<int, Modification>(3, modToAdd) };
            modList.Add("test", Hash);
            ProteinDbWriter.WriteXmlDatabase(modList, new List<Protein> { TestProtein }, xmlName);

            //now write MZML file
            var protein = ProteinDbLoader.LoadProteinXML(xmlName, true, DecoyType.Reverse, new List<Modification>(), false, new List<string>(), out Dictionary<string, Modification> ok);
            var setList1 = protein[0].Digest(testPeptides.CommonParameters.DigestionParams, new List<ModificationWithMass> { }, variableModifications).ToList();
            Assert.AreEqual(4, setList1.Count);

            //Finally Write MZML file
            MsDataFile myMsDataFile = new TestDataFile(new List<PeptideWithSetModifications> { setList1[0], setList1[1], setList1[2], setList1[3], setList1[0], setList1[1] });
            string mzmlName = @"singleProteinWithRepeatedMods.mzML";
            IO.MzML.MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(myMsDataFile, mzmlName, false);

            // RUN!

            string outputFolderInThisTest = MySetUpClass.outputFolder;
            var engine = new EverythingRunnerEngine(taskList, new List<string> { mzmlName }, new List<DbForTask> { new DbForTask(xmlName, false) }, Environment.CurrentDirectory);
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
        }
    }
}