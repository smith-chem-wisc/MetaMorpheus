using EngineLayer;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Proteomics;
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
        #region Public Fields

        public static bool hasPrunedRun;

        #endregion Public Fields

        #region Public Methods

        [Test]
        public static void TestEverythingRunner()
        {
            #region Setup tasks

            foreach (var modFile in Directory.GetFiles(@"Mods"))
                GlobalEngineLevelSettings.AddMods(PtmListLoader.ReadModsFromFile(modFile));

            CalibrationTask task1 = new CalibrationTask
            {
                CommonParameters = new CommonParameters
                {
                    ConserveMemory = false,
                    DigestionParams = new DigestionParams
                    {
                        MaxMissedCleavages = 0,
                        MinPeptideLength = null,
                        InitiatorMethionineBehavior = InitiatorMethionineBehavior.Retain
                    },
                },
                CalibrationParameters = new CalibrationParameters
                {
                    WriteIntermediateFiles = true,
                    NumFragmentsNeededForEveryIdentification = 6,
                }
            };
            GptmdTask task2 = new GptmdTask
            {
                CommonParameters = new CommonParameters
                {
                    ConserveMemory = false
                },
            };

            SearchTask task3 = new SearchTask
            {
                CommonParameters = new CommonParameters
                {
                    ConserveMemory = false
                },
                SearchParameters = new SearchParameters
                {
                    DoParsimony = true,
                    SearchType = SearchType.Modern
                }
            };

            SearchTask task4 = new SearchTask
            {
                CommonParameters = new CommonParameters
                {
                    ConserveMemory = false
                },
                SearchParameters = new SearchParameters
                {
                    SearchType = SearchType.Modern,
                }
            };
            List<Tuple<string, MetaMorpheusTask>> taskList = new List<Tuple<string, MetaMorpheusTask>> {
                new Tuple<string, MetaMorpheusTask>("task1", task1),
                new Tuple<string, MetaMorpheusTask>("task2", task2),
                new Tuple<string, MetaMorpheusTask>("task3", task3),
                new Tuple<string, MetaMorpheusTask>("task4", task4),};

            #endregion Setup tasks

            List<ModificationWithMass> variableModifications = GlobalEngineLevelSettings.AllModsKnown.OfType<ModificationWithMass>().Where(b => task1.CommonParameters.ListOfModsVariable.Contains(new Tuple<string, string>(b.modificationType, b.id))).ToList();
            List<ModificationWithMass> fixedModifications = GlobalEngineLevelSettings.AllModsKnown.OfType<ModificationWithMass>().Where(b => task1.CommonParameters.ListOfModsFixed.Contains(new Tuple<string, string>(b.modificationType, b.id))).ToList();

            // Generate data for files
            Protein ParentProtein = new Protein("MPEPTIDEKANTHE", "accession1");

            var digestedList = ParentProtein.Digest(task1.CommonParameters.DigestionParams, fixedModifications, variableModifications).ToList();

            Assert.AreEqual(3, digestedList.Count);

            PeptideWithSetModifications pepWithSetMods1 = digestedList[0];

            PeptideWithSetModifications pepWithSetMods2 = digestedList[2];

            var dictHere = new Dictionary<int, List<Modification>>();
            ModificationMotif.TryGetMotif("E", out ModificationMotif motif);
            dictHere.Add(3, new List<Modification> { new ModificationWithMass("21", null, motif, TerminusLocalization.Any, 21.981943) });
            Protein ParentProteinToNotInclude = new Protein("MPEPTIDEK", "accession2", new List<Tuple<string, string>>(), dictHere);
            digestedList = ParentProteinToNotInclude.Digest(task1.CommonParameters.DigestionParams, fixedModifications, variableModifications).ToList();

            IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> myMsDataFile = new TestDataFile(new List<PeptideWithSetModifications> { pepWithSetMods1, pepWithSetMods2, digestedList[1] });

            Protein proteinWithChain = new Protein("MAACNNNCAA", "accession3", new List<Tuple<string, string>>(), new Dictionary<int, List<Modification>>(), new List<ProteolysisProduct> { new ProteolysisProduct(4, 8, "chain") }, "name2", "fullname2");

            #region Write the files

            string mzmlName = @"ok.mzML";
            IO.MzML.MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(myMsDataFile, mzmlName, false);
            string xmlName = "okk.xml";
            ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), new List<Protein> { ParentProtein, proteinWithChain }, xmlName);

            #endregion Write the files

            // RUN!
            var engine = new EverythingRunnerEngine(taskList, new List<string> { mzmlName }, new List<DbForTask> { new DbForTask(xmlName, false) }, Environment.CurrentDirectory);
            engine.Run();
        }

        [Test]
        public static void TestMultipleFilesRunner()
        {
            #region Setup tasks

            foreach (var modFile in Directory.GetFiles(@"Mods"))
                GlobalEngineLevelSettings.AddMods(PtmListLoader.ReadModsFromFile(modFile));

            CalibrationTask task1 = new CalibrationTask
            {
                CommonParameters = new CommonParameters
                {
                    DigestionParams = new DigestionParams
                    {
                        MaxMissedCleavages = 0,
                        MinPeptideLength = null,
                        InitiatorMethionineBehavior = InitiatorMethionineBehavior.Retain
                    },
                    ListOfModsVariable = new List<Tuple<string, string>> { new Tuple<string, string>("Common Variable", "Oxidation of M") },
                    ListOfModsFixed = new List<Tuple<string, string>> { new Tuple<string, string>("Common Fixed", "Carbamidomethyl of C") },
                    ListOfModsLocalize = GlobalEngineLevelSettings.AllModsKnown.Select(b => new Tuple<string, string>(b.modificationType, b.id)).ToList(),
                    ProductMassTolerance = new AbsoluteTolerance(0.01)
                },
                CalibrationParameters = new CalibrationParameters
                {
                    NumFragmentsNeededForEveryIdentification = 6,
                }
            };
            GptmdTask task2 = new GptmdTask
            {
                CommonParameters = new CommonParameters
                {
                    DigestionParams = new DigestionParams
                    {
                        Protease = GlobalEngineLevelSettings.ProteaseDictionary["trypsin"],
                    },
                    ProductMassTolerance = new AbsoluteTolerance(0.01)
                },
            };

            SearchTask task3 = new SearchTask
            {
                CommonParameters = new CommonParameters
                {
                    ConserveMemory = false
                },
                SearchParameters = new SearchParameters
                {
                    DoParsimony = true,
                    SearchType = SearchType.Modern,
                }
            };
            SearchTask task4 = new SearchTask
            {
                CommonParameters = new CommonParameters
                {
                    ConserveMemory = false
                },
                SearchParameters = new SearchParameters
                {
                    SearchType = SearchType.Modern,
                }
            };
            List<Tuple<string, MetaMorpheusTask>> taskList = new List<Tuple<string, MetaMorpheusTask>> {
                new Tuple<string, MetaMorpheusTask>("task1", task1),
                new Tuple<string, MetaMorpheusTask>("task2", task2),
                new Tuple<string, MetaMorpheusTask>("task3", task3),
                new Tuple<string, MetaMorpheusTask>("task4", task4),};

            #endregion Setup tasks

            List<ModificationWithMass> variableModifications = GlobalEngineLevelSettings.AllModsKnown.OfType<ModificationWithMass>().Where(b => task1.CommonParameters.ListOfModsVariable.Contains(new Tuple<string, string>(b.modificationType, b.id))).ToList();
            List<ModificationWithMass> fixedModifications = GlobalEngineLevelSettings.AllModsKnown.OfType<ModificationWithMass>().Where(b => task1.CommonParameters.ListOfModsFixed.Contains(new Tuple<string, string>(b.modificationType, b.id))).ToList();

            // Generate data for files
            Protein ParentProtein = new Protein("MPEPTIDEKANTHE", "accession1");

            var digestedList = ParentProtein.Digest(task1.CommonParameters.DigestionParams, fixedModifications, variableModifications).ToList();

            Assert.AreEqual(3, digestedList.Count);

            PeptideWithSetModifications pepWithSetMods1 = digestedList[0];

            PeptideWithSetModifications pepWithSetMods2 = digestedList[2];

            var dictHere = new Dictionary<int, List<Modification>>();
            ModificationMotif.TryGetMotif("E", out ModificationMotif motif);
            dictHere.Add(3, new List<Modification> { new ModificationWithMass("21", null, motif, TerminusLocalization.Any, 21.981943) });
            Protein ParentProteinToNotInclude = new Protein("MPEPTIDEK", "accession2", new List<Tuple<string, string>>(), dictHere);
            digestedList = ParentProteinToNotInclude.Digest(task1.CommonParameters.DigestionParams, fixedModifications, variableModifications).ToList();
            Assert.AreEqual(4, digestedList.Count);

            IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> myMsDataFile1 = new TestDataFile(new List<PeptideWithSetModifications> { pepWithSetMods1, pepWithSetMods2, digestedList[1] });

            string mzmlName1 = @"ok1.mzML";
            IO.MzML.MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(myMsDataFile1, mzmlName1, false);

            IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> myMsDataFile2 = new TestDataFile(new List<PeptideWithSetModifications> { pepWithSetMods1, pepWithSetMods2, digestedList[1] });

            string mzmlName2 = @"ok2.mzML";
            IO.MzML.MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(myMsDataFile2, mzmlName2, false);

            Protein proteinWithChain1 = new Protein("MAACNNNCAA", "accession3", new List<Tuple<string, string>>(), new Dictionary<int, List<Modification>>(), new List<ProteolysisProduct> { new ProteolysisProduct(4, 8, "chain") }, "name2", "fullname2", false, false, new List<DatabaseReference>(), new List<SequenceVariation>(), null);
            Protein proteinWithChain2 = new Protein("MAACNNNCAA", "accession3", new List<Tuple<string, string>>(), new Dictionary<int, List<Modification>>(), new List<ProteolysisProduct> { new ProteolysisProduct(4, 8, "chain") }, "name2", "fullname2", false, false, new List<DatabaseReference>(), new List<SequenceVariation>(), null);

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
                {
                    DigestionParams = new DigestionParams
                    {
                        MinPeptideLength = 2,
                    },
                    ScoreCutoff = 1,
                    DeconvolutionIntensityRatio = 999,
                    DeconvolutionMassTolerance = new PpmTolerance(50),
                },
                SearchParameters = new SearchParameters
                {
                    DecoyType = DecoyType.None,
                    MassDiffAcceptorType = MassDiffAcceptorType.Open,
                }
            };

            string xmlName = "MakeSureFdrDoesntSkip.xml";

            #region Generate protein and write to file

            {
                Protein theProtein = new Protein("MG", "accession1");
                ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), new List<Protein> { theProtein }, xmlName);
            }

            #endregion Generate protein and write to file

            string mzmlName = @"MakeSureFdrDoesntSkip.mzML";

            #region Generate and write the mzml

            {
                var theProteins = ProteinDbLoader.LoadProteinXML(xmlName, true, DecoyType.Reverse, new List<Modification>(), false, new List<string>(), out Dictionary<string, Modification> ok);

                List<ModificationWithMass> fixedModifications = new List<ModificationWithMass>();

                var targetDigested = theProteins[0].Digest(task.CommonParameters.DigestionParams, fixedModifications, GlobalEngineLevelSettings.AllModsKnown.OfType<ModificationWithMass>().ToList()).ToList();

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
            }

            #endregion Generate and write the mzml

            // RUN!
            var theStringResult = task.RunTask(TestContext.CurrentContext.TestDirectory, new List<DbForTask> { new DbForTask(xmlName, false) }, new List<string> { mzmlName }, "taskId1").ToString();
            Assert.IsTrue(theStringResult.Contains("All target PSMS within 1% FDR: 1"));
        }

        [Test]
        public static void MakeSureGptmdTaskMatchesExactMatches()
        {
            MetaMorpheusTask task1;

            #region Setup tasks

            {
                ModificationMotif.TryGetMotif("T", out ModificationMotif motif);
                GlobalEngineLevelSettings.AddMods(new List<ModificationWithMass> { new ModificationWithMass("ok", "okType", motif, TerminusLocalization.Any, 229) });
                task1 = new GptmdTask
                {
                    CommonParameters = new CommonParameters
                    {
                        ConserveMemory = false,
                        DigestionParams = new DigestionParams
                        {
                            InitiatorMethionineBehavior = InitiatorMethionineBehavior.Retain,
                        },
                        ListOfModsVariable = new List<Tuple<string, string>>(),
                        ListOfModsFixed = new List<Tuple<string, string>>(),
                        ScoreCutoff = 1,
                        PrecursorMassTolerance = new AbsoluteTolerance(1)
                    },

                    GptmdParameters = new GptmdParameters
                    {
                        ListOfModsGptmd = new List<Tuple<string, string>> { new Tuple<string, string>("okType", "ok") },
                    }
                };
            }

            #endregion Setup tasks

            string xmlName = "sweetness.xml";

            #region Generate protein and write to file

            {
                Protein theProtein = new Protein("MPEPTIDEKANTHE", "accession1");
                ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), new List<Protein> { theProtein }, xmlName);
            }

            #endregion Generate protein and write to file

            string mzmlName = @"ok.mzML";

            #region Generate and write the mzml

            {
                var theProteins = ProteinDbLoader.LoadProteinXML(xmlName, true, DecoyType.Reverse, new List<Modification>(), false, new List<string>(), out Dictionary<string, Modification> ok);

                List<ModificationWithMass> fixedModifications = new List<ModificationWithMass>();

                var targetDigested = theProteins[0].Digest(task1.CommonParameters.DigestionParams, fixedModifications, GlobalEngineLevelSettings.AllModsKnown.OfType<ModificationWithMass>().ToList()).ToList();

                ModificationMotif.TryGetMotif("T", out ModificationMotif motif);
                PeptideWithSetModifications targetGood = targetDigested[0];

                PeptideWithSetModifications targetWithUnknownMod = targetDigested[1];
                IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> myMsDataFile = new TestDataFile(new List<PeptideWithSetModifications> { targetGood, targetWithUnknownMod }, true);

                IO.MzML.MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(myMsDataFile, mzmlName, false);
            }

            #endregion Generate and write the mzml

            // RUN!
            var theStringResult = task1.RunTask(TestContext.CurrentContext.TestDirectory, new List<DbForTask> { new DbForTask(xmlName, false) }, new List<string> { mzmlName }, "taskId1").ToString();
            Assert.IsTrue(theStringResult.Contains("Modifications added: 1"));
        }

        //test if prunedDatabase matches expected output
        [Test]
        public static void TestPrunedDatabase()
        {
            hasPrunedRun = true;

            #region setup

            //Create Search Task
            SearchTask task1 = new SearchTask
            {
                SearchParameters = new SearchParameters
                {
                    WritePrunedDatabase = true,
                    MassDiffAcceptorType = MassDiffAcceptorType.Exact
                }
            };

            //add task 1 to task list
            List<Tuple<string, MetaMorpheusTask>> taskList = new List<Tuple<string, MetaMorpheusTask>> {
               new Tuple<string, MetaMorpheusTask>("task1", task1)};

            ModificationMotif.TryGetMotif("P", out ModificationMotif motif);

            var connorMod = new ModificationWithMass("ConnorMod", "ConnorModType", motif, TerminusLocalization.Any, 10);

            GlobalEngineLevelSettings.AddMods(new List<ModificationWithLocation>
            {
                connorMod
            });

            #endregion setup

            #region Protein and Mod Creation

            //create modification lists
            List<ModificationWithMass> variableModifications = GlobalEngineLevelSettings.AllModsKnown.OfType<ModificationWithMass>().Where(b => task1.CommonParameters.ListOfModsVariable.Contains(new Tuple<string, string>(b.modificationType, b.id))).ToList();

            //add modification to Protein object
            var dictHere = new Dictionary<int, List<Modification>>();
            ModificationWithMass modToAdd = connorMod;
            ModificationWithMass modToAdd2 = connorMod;
            dictHere.Add(1, new List<Modification> { modToAdd });
            dictHere.Add(3, new List<Modification> { modToAdd2 });

            //protein Creation (One with mod and one without)
            Protein TestProtein = new Protein("PEPTID", "accession1");
            Protein TestProteinWithMod = new Protein("PEPTID", "accession1", new List<Tuple<string, string>>(), dictHere);

            #endregion Protein and Mod Creation

            #region XML File

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

            #endregion XML File

            #region MZML File

            //now write MZML file
            var protein = ProteinDbLoader.LoadProteinXML(xmlName, true, DecoyType.Reverse, new List<Modification>(), false, new List<string>(), out Dictionary<string, Modification> ok);
            var digestedList = protein[0].Digest(task1.CommonParameters.DigestionParams, new List<ModificationWithMass> { }, variableModifications).ToList();
            Assert.AreEqual(4, digestedList.Count);

            //Set Peptide with 1 mod at position 3
            PeptideWithSetModifications pepWithSetMods1 = digestedList[1];

            //Finally Write MZML file
            Assert.AreEqual("PEP[ConnorModType:ConnorMod]TID", pepWithSetMods1.Sequence);
            IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> myMsDataFile = new TestDataFile(new List<PeptideWithSetModifications> { pepWithSetMods1 });
            string mzmlName = @"hello.mzML";
            IO.MzML.MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(myMsDataFile, mzmlName, false);

            #endregion MZML File

            //run!
            var engine = new EverythingRunnerEngine(taskList, new List<string> { mzmlName }, new List<DbForTask> { new DbForTask(xmlName, false) }, Environment.CurrentDirectory);
            engine.Run();

            string outputFolderInThisTest = MySetUpClass.outputFolder;
            string final = Path.Combine(MySetUpClass.outputFolder, "task1", "okkkpruned.xml");
            //string[] files = Directory.GetFiles(fileAtPath);
            //string file = fileAtPath;
            var proteins = ProteinDbLoader.LoadProteinXML(final, true, DecoyType.Reverse, new List<Modification>(), false, new List<string>(), out ok);
            //check length
            Assert.AreEqual(proteins[0].OneBasedPossibleLocalizedModifications.Count, 1);
            //check location (key)
            Assert.AreEqual(proteins[0].OneBasedPossibleLocalizedModifications.ContainsKey(3), true);
            List<Modification> listOfMods = new List<Modification>();
            listOfMods = proteins[0].OneBasedPossibleLocalizedModifications[3];
            //check Type, count, ID
            Assert.AreEqual(listOfMods[0].modificationType, "ConnorModType");
            Assert.AreEqual(listOfMods[0].id, "ConnorMod");
            Assert.AreEqual(listOfMods.Count, 1);
        }

        [Test]
        public static void TestUniquePeptideCount()
        {
            #region setup

            SearchTask testUnique = new SearchTask
            {
                CommonParameters = new CommonParameters
                {
                    ListOfModsLocalize = new List<Tuple<string, string>> { new Tuple<string, string>("ConnorModType", "ConnorMod") },
                },
                SearchParameters = new SearchParameters
                {
                    WritePrunedDatabase = true,
                    MassDiffAcceptorType = MassDiffAcceptorType.Exact
                }
            };

            List<Tuple<string, MetaMorpheusTask>> taskList = new List<Tuple<string, MetaMorpheusTask>> {
               new Tuple<string, MetaMorpheusTask>("TestUnique", testUnique)};

            ModificationMotif.TryGetMotif("P", out ModificationMotif motif);

            var testUniqeMod = new ModificationWithMass("testUniqeMod", "mt", motif, TerminusLocalization.Any, 10);
            GlobalEngineLevelSettings.AddMods(new List<ModificationWithLocation>
            {
                testUniqeMod
            });

            #endregion setup

            #region mod setup and protein creation

            //create modification lists

            List<ModificationWithMass> variableModifications = GlobalEngineLevelSettings.AllModsKnown.OfType<ModificationWithMass>().Where(b => testUnique.CommonParameters.ListOfModsVariable.Contains(new Tuple<string, string>(b.modificationType, b.id))).ToList();

            //add modification to Protein object
            var modDictionary = new Dictionary<int, List<Modification>>();
            ModificationWithMass modToAdd = testUniqeMod;
            modDictionary.Add(1, new List<Modification> { modToAdd });
            modDictionary.Add(3, new List<Modification> { modToAdd });

            //protein Creation (One with mod and one without)
            Protein TestProtein = new Protein("PEPTID", "accession1", new List<Tuple<string, string>>(), modDictionary);

            #endregion mod setup and protein creation

            #region XML setup

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

            #endregion XML setup

            #region MZML setup

            //now write MZML file
            var protein = ProteinDbLoader.LoadProteinXML(xmlName, true, DecoyType.Reverse, new List<Modification>(), false, new List<string>(), out Dictionary<string, Modification> ok);
            var setList1 = protein[0].Digest(testUnique.CommonParameters.DigestionParams, new List<ModificationWithMass> { }, variableModifications).ToList();
            Assert.AreEqual(4, setList1.Count);

            //Finally Write MZML file
            IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> myMsDataFile = new TestDataFile(new List<PeptideWithSetModifications> { setList1[0], setList1[1], setList1[2], setList1[3], setList1[0], setList1[1] });
            string mzmlName = @"singleProteinWithRepeatedMods.mzML";
            IO.MzML.MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(myMsDataFile, mzmlName, false);

            #endregion MZML setup

            #region run

            string outputFolderInThisTest = MySetUpClass.outputFolder;
            var engine = new EverythingRunnerEngine(taskList, new List<string> { mzmlName }, new List<DbForTask> { new DbForTask(xmlName, false) }, Environment.CurrentDirectory);
            engine.Run();

            string line;

            bool foundD = false;
            using (StreamReader file = new StreamReader(Path.Combine(MySetUpClass.outputFolder, "TestUnique", "results.txt")))
            {
                while ((line = file.ReadLine()) != null)
                {
                    if (line.Contains("Unique target peptides within 1% FDR: 4"))
                    {
                        foundD = true;
                    }
                }
            }
            Assert.IsTrue(foundD);

            #endregion run
        }

        #endregion Public Methods
    }
}