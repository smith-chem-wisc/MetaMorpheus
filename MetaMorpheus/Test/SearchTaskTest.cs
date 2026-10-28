using EngineLayer;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework; 
using Proteomics;
using Omics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using EngineLayer.DatabaseLoading;
using Omics;
using Omics.Digestion;
using Omics.Modifications;
using Readers;
using TaskLayer;
using UsefulProteomicsDatabases;

namespace Test
{
    [TestFixture]
    public static class SearchTaskTest
    {
        /// <summary>
        /// Tests each type of mass difference acceptor type to make sure values are assigned properly
        /// </summary>
        [Test]
        public static void MassDiffAceptorTest()
        {
            SearchTask searchTask = new SearchTask();
            var result = SearchTask.GetMassDiffAcceptor(searchTask.CommonParameters.PrecursorMassTolerance, searchTask.SearchParameters.MassDiffAcceptorType, searchTask.SearchParameters.CustomMdac);
            Assert.That(result.FileNameAddition.Equals("1mm"));

            result = SearchTask.GetMassDiffAcceptor(searchTask.CommonParameters.PrecursorMassTolerance, MassDiffAcceptorType.TwoMM, searchTask.SearchParameters.CustomMdac);
            Assert.That(result.FileNameAddition.Equals("2mm"));

            result = SearchTask.GetMassDiffAcceptor(searchTask.CommonParameters.PrecursorMassTolerance, MassDiffAcceptorType.ThreeMM, searchTask.SearchParameters.CustomMdac);
            Assert.That(result.FileNameAddition.Equals("3mm"));

            result = SearchTask.GetMassDiffAcceptor(searchTask.CommonParameters.PrecursorMassTolerance, MassDiffAcceptorType.ModOpen, searchTask.SearchParameters.CustomMdac);
            Assert.That(result.FileNameAddition.Equals("-187andUp"));

            result = SearchTask.GetMassDiffAcceptor(searchTask.CommonParameters.PrecursorMassTolerance, MassDiffAcceptorType.Open, searchTask.SearchParameters.CustomMdac);
            Assert.That(result.FileNameAddition.Equals("OpenSearch"));

            result = SearchTask.GetMassDiffAcceptor(searchTask.CommonParameters.PrecursorMassTolerance, MassDiffAcceptorType.Custom, "custom ppmAroundZero 4");
            Assert.That(result.FileNameAddition.Equals("4ppmAroundZero"));

            result = SearchTask.GetMassDiffAcceptor(searchTask.CommonParameters.PrecursorMassTolerance, MassDiffAcceptorType.Exact, searchTask.SearchParameters.CustomMdac);
            Assert.That(result.FileNameAddition.Equals("5ppmAroundZero"));

            result = SearchTask.GetMassDiffAcceptor(searchTask.CommonParameters.PrecursorMassTolerance, MassDiffAcceptorType.PlusOrMinusThreeMM, searchTask.SearchParameters.CustomMdac);
            Assert.That(result.FileNameAddition.Equals("PlusOrMinus3Da"));

            result = SearchTask.GetMassDiffAcceptor(searchTask.CommonParameters.PrecursorMassTolerance, MassDiffAcceptorType.MostAbundant_Exact, searchTask.SearchParameters.CustomMdac);
            Assert.That(result.NumNotches == 1);

            result = SearchTask.GetMassDiffAcceptor(searchTask.CommonParameters.PrecursorMassTolerance, MassDiffAcceptorType.MostAbundant_PlusMinusOne, searchTask.SearchParameters.CustomMdac);
            Assert.That(result.NumNotches == 3);

            result = SearchTask.GetMassDiffAcceptor(searchTask.CommonParameters.PrecursorMassTolerance, MassDiffAcceptorType.MostAbundant_PlusMinusTwo, searchTask.SearchParameters.CustomMdac);
            Assert.That(result.NumNotches == 5);
        }

        /// <summary>
        /// Tests to make sure custom mass difference acceptor inputs are parsed properly
        /// </summary>
        [Test]
        public static void ParseSearchModeTest()
        {
            SearchTask searchTask = new SearchTask();
            var result = SearchTask.GetMassDiffAcceptor(searchTask.CommonParameters.PrecursorMassTolerance, MassDiffAcceptorType.Custom, "TestCustom dot 5 ppm 0,1.0029,2.0052");
            Assert.That(result.NumNotches == 3);

            result = SearchTask.GetMassDiffAcceptor(searchTask.CommonParameters.PrecursorMassTolerance, MassDiffAcceptorType.Custom, "TestCustom dot 5 da 0,1.0029,2.0052");
            Assert.That(result.NumNotches == 3);

            result = SearchTask.GetMassDiffAcceptor(searchTask.CommonParameters.PrecursorMassTolerance, MassDiffAcceptorType.Custom, "TestCustom interval [0,5];[0,5]");
            Assert.That(result.NumNotches == 1);

            result = SearchTask.GetMassDiffAcceptor(searchTask.CommonParameters.PrecursorMassTolerance, MassDiffAcceptorType.Custom, "TestCustom OpenSearch 5");
            Assert.That(result.FileNameAddition.Equals("OpenSearch"));

            result = SearchTask.GetMassDiffAcceptor(searchTask.CommonParameters.PrecursorMassTolerance, MassDiffAcceptorType.Custom, "TestCustom daltonsAroundZero 5");
            Assert.That(result.FileNameAddition.Equals("5daltonsAroundZero"));

            result = SearchTask.GetMassDiffAcceptor(searchTask.CommonParameters.PrecursorMassTolerance, MassDiffAcceptorType.Custom, "TestCustom ppmAroundZero 5");
            Assert.That(result.FileNameAddition.Equals("5ppmAroundZero"));

            Assert.That(() => SearchTask.GetMassDiffAcceptor(searchTask.CommonParameters.PrecursorMassTolerance, MassDiffAcceptorType.Custom, "TestCustom Test 5"),
                Throws.TypeOf<MetaMorpheusException>());
        }

        /// <summary>
        /// Ensures that the minimum peptide length is observed (KLEDHPK)
        /// Ensures semispecific search finds peptides that were cleaved correctly during the first digestion (precursor index is made and searched correctly) (KEDEEDKFDAMGNK)
        /// </summary>
        [Test]
        public static void SemiSpecificFullAndSmallMatches()
        {
            
            SearchTask searchTask = new SearchTask()
            {

                SearchParameters = new SearchParameters
                {
                    WriteMzId = true,
                    SearchType = SearchType.NonSpecific,
                    LocalFdrCategories = new List<FdrCategory>
                        {
                            FdrCategory.FullySpecific,
                            FdrCategory.SemiSpecific
                        }
                },
                CommonParameters = new CommonParameters(addCompIons: true, scoreCutoff: 11,
                    digestionParams: new DigestionParams(minPeptideLength: 7, searchModeType: CleavageSpecificity.Semi, fragmentationTerminus: FragmentationTerminus.N))
            };

            string myFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\tinySemi.mgf");
            string myDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\semiTest.fasta");
            DbForTask db = new DbForTask(myDatabase, false);

            List<(string, MetaMorpheusTask)> taskList = new List<(string, MetaMorpheusTask)> { ("TestSemiSpecificSmall", searchTask) };


            var engine = new EverythingRunnerEngine(taskList, new List<string> { myFile }, new List<DbForTask> { new DbForTask(myDatabase, false) }, Environment.CurrentDirectory);
            engine.Run();

            string outputPath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestSemiSpecificSmall\AllPSMs.psmtsv");
            var output = File.ReadAllLines(outputPath);
            Assert.That(output.Length == 3);

            var mzId = File.ReadAllLines(Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestSemiSpecificSmall\tinySemi.mzID"));
            Assert.That(mzId[115].Equals("          <cvParam name=\"mzML format\" cvRef=\"PSI-MS\" accession=\"MS:1000584\" />"));
            Assert.That(mzId[118].Equals("          <cvParam name=\"mzML unique identifier\" cvRef=\"PSI-MS\" accession=\"MS:1001530\" />"));
            Assert.That(mzId[97].Equals("        <cvParam name=\"pep:FDR threshold\" value=\"0.01\" cvRef=\"PSI-MS\" accession=\"MS:1001448\" />"));

        }

        /// <summary>
        /// Ensures semispecific search runs and outputs properly
        /// </summary>
        [Test]
        public static void SemiSpecificTest()
        {
            List<FragmentationTerminus> terminiToTest = new List<FragmentationTerminus>
            {
                FragmentationTerminus.N,
                FragmentationTerminus.C
            };
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestSemiSpecific");
            string myFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\PrunedDbSpectra.mzml");
            string myDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\DbForPrunedDb.fasta");
            foreach (FragmentationTerminus fragTerm in terminiToTest)
            {
                SearchTask searchTask = new SearchTask()
                {
                    SearchParameters = new SearchParameters
                    {
                        SearchType = SearchType.NonSpecific,
                        LocalFdrCategories = new List<FdrCategory>
                        {
                            FdrCategory.FullySpecific,
                            FdrCategory.SemiSpecific
                        }
                    },
                    CommonParameters = new CommonParameters(scoreCutoff: 4, addCompIons: true,
                    digestionParams: new DigestionParams(searchModeType: CleavageSpecificity.Semi, fragmentationTerminus: fragTerm))
                };

                List<(string, MetaMorpheusTask)> taskList = new List<(string, MetaMorpheusTask)> { ("TestSemiSpecific", searchTask) };

                var engine = new EverythingRunnerEngine(taskList, new List<string> { myFile }, new List<DbForTask> { new DbForTask(myDatabase, false) }, outputFolder);
                engine.Run();

                string outputPath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestSemiSpecific\TestSemiSpecific\AllPSMs.psmtsv");
                var output = File.ReadAllLines(outputPath);
                Assert.That(output.Length == 12); //if N is only producing 11 lines, then the c is not being searched with it. //Possibly 13 lines if decon changes because of missed mono
            }
            Directory.Delete(outputFolder, true);
        }


        /// <summary>
        /// Ensure internal fragment ions are being matched correctly and can disambiguate ambiguous proteoforms
        /// </summary>
        [Test]
        public static void InternalFragmentIonTest()
        {
            SearchTask searchTask = new SearchTask()
            {

                SearchParameters = new SearchParameters
                {
                    MinAllowedInternalFragmentLength = 1
                },
                CommonParameters = new CommonParameters(
                   digestionParams: new DigestionParams("top-down"),
                   listOfModsVariable: new List<(string, string)> {
                       ("Common Variable", "Oxidation on M"),
                       ("Common Biological", "Acetylation on K"),
                       ("Common Biological", "Acetylation on X")  })
            };

            string myFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\InternalTest.mgf");
            string myDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\InternalTest.fasta");
            DbForTask db = new DbForTask(myDatabase, false);

            List<(string, MetaMorpheusTask)> taskList = new List<(string, MetaMorpheusTask)> { ("TestInternal", searchTask) };


            var engine = new EverythingRunnerEngine(taskList, new List<string> { myFile }, new List<DbForTask> { new DbForTask(myDatabase, false) }, Environment.CurrentDirectory);
            engine.Run();

            string outputPath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestInternal\AllPSMs.psmtsv");
            //var output = File.ReadAllLines(outputPath);
            //read the psmtsv
            List<PsmFromTsv> psms = SpectrumMatchTsvReader.ReadPsmTsv(outputPath, out var warning);
            Assert.That(psms.Count == 1);
            //check that it's been disambiguated
            Assert.That(!(psms[0].FullSequence.Contains("|")));
            Assert.That(psms[0].MatchedIons.First().Intensity, Is.EqualTo(161210));
            Assert.That(psms[0].MatchedIons.First().Mz, Is.EqualTo(585.25292));
            Assert.That(psms[0].MatchedIons.First().Charge, Is.EqualTo(1));
            Assert.That(psms[0].MatchedIons[4].Intensity, Is.EqualTo(131546));
            Assert.That(psms[0].MatchedIons[4].Mz, Is.EqualTo(782.84816));
            Assert.That(psms[0].MatchedIons[4].Charge, Is.EqualTo(2));
            int numTotalFragments = psms[0].MatchedIons.Count;


            //test again but no variable acetyl on K. Make sure that internal fragments are still searched even without ambiguity
            searchTask = new SearchTask()
            {

                SearchParameters = new SearchParameters
                {
                    MinAllowedInternalFragmentLength = 1
                },
                CommonParameters = new CommonParameters(
                   digestionParams: new DigestionParams("top-down"),
                   listOfModsVariable: new List<(string, string)> {
                       ("Common Variable", "Oxidation on M"),
                       ("Common Biological", "Acetylation on X")  })
            };
            taskList = new List<(string, MetaMorpheusTask)> { ("TestInternal", searchTask) };
            engine = new EverythingRunnerEngine(taskList, new List<string> { myFile }, new List<DbForTask> { new DbForTask(myDatabase, false) }, Environment.CurrentDirectory);
            engine.Run();
            psms = SpectrumMatchTsvReader.ReadPsmTsv(outputPath, out warning);
            Assert.That(psms.Count == 1);
            Assert.That(psms[0].MatchedIons.Count == numTotalFragments);

            Directory.Delete(Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestInternal"), true);
        }

        /// <summary>
        /// Tests that normalization in a search task works properly with an Experimental Design file read in,
        /// and skips quantification when that file is absent
        /// </summary>
        [Test]
        public static void PostSearchNormalizeTest()
        {
            SearchTask searchTask = new SearchTask()
            {
                SearchParameters = new SearchParameters
                {
                    Normalize = true
                },
                CommonParameters = new( precursorDeconParams:new IsoDecDeconvolutionParameters())
                {
                }
            };

            string myFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\PrunedDbSpectra.mzml");
            string myDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\DbForPrunedDb.fasta");
            string folderPath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestNormalizationExperDesign");
            string experimentalDesignFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\ExperimentalDesign.tsv");
            DbForTask db = new DbForTask(myDatabase, false);

            try
            {
                using (StreamWriter output = new StreamWriter(experimentalDesignFile))
                {
                    output.WriteLine("FileName\tCondition\tBiorep\tFraction\tTechrep");
                    output.WriteLine("PrunedDbSpectra.mzml" + "\t" + "condition" + "\t" + "1" + "\t" + "1" + "\t" + "1");
                }

                // run the task
                Directory.CreateDirectory(folderPath);
                searchTask.RunTask(folderPath, new List<DbForTask> { db }, new List<string> { myFile }, "normal");

                Directory.Delete(folderPath, true);

                // delete the exper design and try again. this should skip quantification
                File.Delete(experimentalDesignFile);

                // run the task
                Directory.CreateDirectory(folderPath);
                searchTask.RunTask(folderPath, new List<DbForTask> { db }, new List<string> { myFile }, "normal");

                // PSMs should be present but no quant output
                Assert.That(!File.Exists(Path.Combine(folderPath, "AllQuantifiedPeptides.tsv")));
                Assert.That(File.Exists(Path.Combine(folderPath, "AllPSMs.psmtsv")));

                Directory.Delete(folderPath, true);
            }
            finally
            {
                // Always remove the ExperimentalDesign.tsv we wrote into shared TestData\, even
                // on assertion/exception, so subsequent tests in the same suite run don't pick
                // up a stale file (which would silently break quantification for them).
                if (File.Exists(experimentalDesignFile)) File.Delete(experimentalDesignFile);
                if (Directory.Exists(folderPath)) Directory.Delete(folderPath, true);
            }
        }

        // Malformed exp design (Normalize=false) -> quant skipped, protein-groups TSV has no dynamic quant columns.
        // Filename remains AllQuantifiedProteinGroups.tsv because it's driven by DoLabelFreeQuantification, not by quant success.
        [Test]
        public static void PostSearchMalformedExperimentalDesignSkipsQuant()
        {
            SearchTask searchTask = new SearchTask()
            {
                SearchParameters = new SearchParameters
                {
                    Normalize = false,
                    DoParsimony = true
                },
                CommonParameters = new(precursorDeconParams: new IsoDecDeconvolutionParameters())
            };

            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestMalformedExpDesign");
            string inputFolder = Path.Combine(outputFolder, "inputs");
            Directory.CreateDirectory(inputFolder);
            string mzmlPath = Path.Combine(inputFolder, "PrunedDbSpectra.mzml");
            string fastaPath = Path.Combine(inputFolder, "DbForPrunedDb.fasta");
            string expDesignPath = Path.Combine(inputFolder, "ExperimentalDesign.tsv");
            File.Copy(Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\PrunedDbSpectra.mzml"), mzmlPath, true);
            File.Copy(Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\DbForPrunedDb.fasta"), fastaPath, true);

            // Lists an unrelated file so ReadExperimentalDesign emits errors.
            using (StreamWriter w = new StreamWriter(expDesignPath))
            {
                w.WriteLine("FileName\tCondition\tBiorep\tFraction\tTechrep");
                w.WriteLine("UnrelatedFile.mzml\tcondition\t1\t1\t1");
            }

            try
            {
                searchTask.RunTask(outputFolder, new List<DbForTask> { new DbForTask(fastaPath, false) }, new List<string> { mzmlPath }, "normal");

                Assert.That(File.Exists(Path.Combine(outputFolder, "AllQuantifiedPeptides.tsv")), Is.False);

                string pgPath = Path.Combine(outputFolder, "AllQuantifiedProteinGroups.tsv");
                Assert.That(File.Exists(pgPath), Is.True);
                var lines = File.ReadAllLines(pgPath);
                Assert.That(lines.Length, Is.GreaterThan(1));
                var header = lines[0];
                // Quantification is skipped, so nothing intensity-derived is reported. Spectral counts
                // and count-based occupancy need no quantification and are reported per spectra file.
                Assert.That(header.Contains("SpectralCount_"), Is.True);
                Assert.That(header.Contains("CountOccupancy_"), Is.True);
                Assert.That(header.Contains("Intensity_"), Is.False);
                Assert.That(header.Contains("IntensityOccupancy_"), Is.False);
                Assert.That(lines.Select(l => l.Split('\t').Length).AllSame(), Is.True);
            }
            finally
            {
                if (Directory.Exists(outputFolder)) Directory.Delete(outputFolder, true);
            }
        }

        // Count-only occupancy across more than one spectra file. Each protein group used to derive
        // its columns from the files its own PSMs came from, so a group seen in only one of two files
        // wrote a short row while the writer emitted the first group's header - sliding Number of PSMs
        // and every column after it two positions left. The single-file guard above cannot see this.
        [Test]
        public static void PostSearchCountOnlyMultiFileWritesOneSchema()
        {
            var uniqueToFirst = new PeptideWithSetModifications("PEPTIDEK", new Dictionary<string, Modification>());
            var uniqueToSecond = new PeptideWithSetModifications("MYRIADPEPTIDER", new Dictionary<string, Modification>());
            var inBothFiles = new PeptideWithSetModifications("ANOTHERLONGPEPTIDEK", new Dictionary<string, Modification>());

            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestCountOnlyMultiFile");
            string inputFolder = Path.Combine(outputFolder, "inputs");
            Directory.CreateDirectory(inputFolder);

            string mzmlPathA = Path.Combine(inputFolder, "countOnlyA.mzML");
            string mzmlPathB = Path.Combine(inputFolder, "countOnlyB.mzML");
            MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(
                new TestDataFile(new List<IBioPolymerWithSetMods> { uniqueToFirst, inBothFiles }), mzmlPathA, false);
            MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(
                new TestDataFile(new List<IBioPolymerWithSetMods> { uniqueToSecond, inBothFiles }), mzmlPathB, false);

            string xmlPath = Path.Combine(inputFolder, "countOnlyDb.xml");
            ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(),
                new List<Protein>
                {
                    new Protein(uniqueToFirst.BaseSequence, "accA"),
                    new Protein(uniqueToSecond.BaseSequence, "accB"),
                    new Protein(inBothFiles.BaseSequence, "accShared")
                }, xmlPath);

            SearchTask searchTask = new SearchTask()
            {
                SearchParameters = new SearchParameters
                {
                    DoLabelFreeQuantification = false,
                    DoParsimony = true
                }
            };

            try
            {
                searchTask.RunTask(outputFolder, new List<DbForTask> { new DbForTask(xmlPath, false) },
                    new List<string> { mzmlPathA, mzmlPathB }, "normal");

                string pgPath = Path.Combine(outputFolder, "AllProteinGroups.tsv");
                Assert.That(File.Exists(pgPath), Is.True);

                var lines = File.ReadAllLines(pgPath);
                Assert.That(lines.Length, Is.GreaterThan(2));

                var headerFields = lines[0].Split('\t');
                Assert.That(headerFields, Does.Contain("SpectralCount_countOnlyA"));
                Assert.That(headerFields, Does.Contain("SpectralCount_countOnlyB"));
                Assert.That(headerFields, Does.Contain("CountOccupancy_countOnlyA"));
                Assert.That(headerFields, Does.Contain("CountOccupancy_countOnlyB"));
                Assert.That(headerFields.Any(f => f.StartsWith("Intensity_")), Is.False);

                Assert.That(lines.Select(l => l.Split('\t').Length).AllSame(), Is.True);

                // Guards the guard: if every group happened to be seen in both files the row widths
                // would agree even with per-group schemas, and the assertion above would prove nothing.
                int countA = Array.IndexOf(headerFields, "SpectralCount_countOnlyA");
                int countB = Array.IndexOf(headerFields, "SpectralCount_countOnlyB");
                var rows = lines.Skip(1).Select(l => l.Split('\t')).ToList();
                Assert.That(rows.Any(r => (r[countA] == "0") != (r[countB] == "0")), Is.True,
                    "every protein group was observed in both files, so the shared schema went untested");
            }
            finally
            {
                if (Directory.Exists(outputFolder)) Directory.Delete(outputFolder, true);
            }
        }

        // No exp design + Normalize=false -> defaults built, LFQ runs, dynamic columns appear.
        // Complements PostSearchNormalizeTest (no exp design + Normalize=true -> skip).
        [Test]
        public static void PostSearchNoExpDesignNoNormalizeRunsQuant()
        {
            SearchTask searchTask = new SearchTask()
            {
                SearchParameters = new SearchParameters
                {
                    Normalize = false,
                    DoParsimony = true
                },
                CommonParameters = new(precursorDeconParams: new IsoDecDeconvolutionParameters())
            };

            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestNoExpDesignNoNormalize");
            string inputFolder = Path.Combine(outputFolder, "inputs");
            Directory.CreateDirectory(inputFolder);
            string mzmlPath = Path.Combine(inputFolder, "PrunedDbSpectra.mzml");
            string fastaPath = Path.Combine(inputFolder, "DbForPrunedDb.fasta");
            File.Copy(Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\PrunedDbSpectra.mzml"), mzmlPath, true);
            File.Copy(Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\DbForPrunedDb.fasta"), fastaPath, true);

            try
            {
                searchTask.RunTask(outputFolder, new List<DbForTask> { new DbForTask(fastaPath, false) }, new List<string> { mzmlPath }, "normal");

                string pgPath = Path.Combine(outputFolder, "AllQuantifiedProteinGroups.tsv");
                Assert.That(File.Exists(pgPath), Is.True);
                var lines = File.ReadAllLines(pgPath);
                Assert.That(lines.Length, Is.GreaterThan(1));
                var header = lines[0];

                Assert.That(header.Contains("SpectralCount_"), Is.True);
                Assert.That(header.Contains("CountOccupancy_"), Is.True);
                Assert.That(lines.Select(l => l.Split('\t').Length).AllSame(), Is.True);
            }
            finally
            {
                if (Directory.Exists(outputFolder)) Directory.Delete(outputFolder, true);
            }
        }

        // Per-file (individual) protein-group TSVs must carry the same quant/occupancy columns as the
        // combined AllQuantifiedProteinGroups.tsv, computed from each file's own intensities. Guards the
        // regression where ConstructSubsetProteinGroup left the subset's SampleGroupResults unpopulated,
        // silently dropping every quant column from the individual-file output.
        [Test]
        public static void PostSearchIndividualFileProteinGroupsHaveQuantColumns()
        {
            SearchTask searchTask = new SearchTask()
            {
                SearchParameters = new SearchParameters
                {
                    Normalize = false,
                    DoParsimony = true,
                    WriteIndividualFiles = true
                },
                CommonParameters = new(precursorDeconParams: new IsoDecDeconvolutionParameters())
            };

            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestIndividualPgQuantColumns");
            string inputFolder = Path.Combine(outputFolder, "inputs");
            Directory.CreateDirectory(inputFolder);
            string fastaPath = Path.Combine(inputFolder, "DbForPrunedDb.fasta");
            File.Copy(Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\DbForPrunedDb.fasta"), fastaPath, true);

            // Two spectra files so individual-file results are written (requires more than one file).
            string mzml1 = Path.Combine(inputFolder, "PrunedDbSpectra1.mzml");
            string mzml2 = Path.Combine(inputFolder, "PrunedDbSpectra2.mzml");
            File.Copy(Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\PrunedDbSpectra.mzml"), mzml1, true);
            File.Copy(Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\PrunedDbSpectra.mzml"), mzml2, true);

            try
            {
                searchTask.RunTask(outputFolder, new List<DbForTask> { new DbForTask(fastaPath, false) },
                    new List<string> { mzml1, mzml2 }, "normal");

                // Combined file: quant ran, so dynamic columns are present (sanity check).
                string combinedPath = Path.Combine(outputFolder, "AllQuantifiedProteinGroups.tsv");
                Assert.That(File.Exists(combinedPath), Is.True);
                Assert.That(File.ReadAllLines(combinedPath)[0].Contains("Intensity_"), Is.True);

                // Individual-file protein-group TSVs must also carry the quant/occupancy columns,
                // each computed from that file's own intensity.
                string individualFolder = Path.Combine(outputFolder, "Individual File Results");
                var individualPgFiles = Directory.GetFiles(individualFolder, "*_ProteinGroups.tsv");
                Assert.That(individualPgFiles.Length, Is.EqualTo(2));

                foreach (var pgFile in individualPgFiles)
                {
                    var lines = File.ReadAllLines(pgFile);
                    Assert.That(lines.Length, Is.GreaterThan(1), $"{pgFile} has no data rows");
                    string header = lines[0];

                    Assert.That(header.Contains("SpectralCount_"), Is.True, $"{pgFile} missing SpectralCount_ column");
                    Assert.That(header.Contains("Intensity_"), Is.True, $"{pgFile} missing Intensity_ column");
                    Assert.That(header.Contains("CountOccupancy_"), Is.True, $"{pgFile} missing CountOccupancy_ column");
                    Assert.That(header.Contains("IntensityOccupancy_"), Is.True, $"{pgFile} missing IntensityOccupancy_ column");

                    // header and every data row must agree on the column count
                    Assert.That(lines.Select(l => l.Split('\t').Length).AllSame(), Is.True, $"{pgFile} has inconsistent column counts");
                }
            }
            finally
            {
                if (Directory.Exists(outputFolder)) Directory.Delete(outputFolder, true);
            }
        }

        /// <summary>
        /// The per-file protein tables list each group's unique and shared peptides, as the combined table does. They
        /// left both cells empty: the per-file subsets were never given their peptide text, while "Number of Unique
        /// Peptides" (counted from the same peptides) was filled. So a row could say 5 unique peptides and list none.
        /// </summary>
        [Test]
        public static void PostSearchIndividualFileProteinGroupsListTheirPeptides()
        {
            SearchTask searchTask = new SearchTask()
            {
                SearchParameters = new SearchParameters
                {
                    Normalize = false,
                    DoParsimony = true,
                    WriteIndividualFiles = true
                },
                CommonParameters = new(precursorDeconParams: new IsoDecDeconvolutionParameters())
            };

            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestIndividualPgPeptideCells");
            string inputFolder = Path.Combine(outputFolder, "inputs");
            Directory.CreateDirectory(inputFolder);
            string fastaPath = Path.Combine(inputFolder, "DbForPrunedDb.fasta");
            File.Copy(Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\DbForPrunedDb.fasta"), fastaPath, true);
            string mzml1 = Path.Combine(inputFolder, "PrunedDbSpectra1.mzml");
            string mzml2 = Path.Combine(inputFolder, "PrunedDbSpectra2.mzml");
            File.Copy(Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\PrunedDbSpectra.mzml"), mzml1, true);
            File.Copy(Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\PrunedDbSpectra.mzml"), mzml2, true);

            try
            {
                searchTask.RunTask(outputFolder, new List<DbForTask> { new DbForTask(fastaPath, false) },
                    new List<string> { mzml1, mzml2 }, "normal");

                var individualPgFiles = Directory.GetFiles(Path.Combine(outputFolder, "Individual File Results"), "*_ProteinGroups.tsv");
                Assert.That(individualPgFiles.Length, Is.EqualTo(2));

                foreach (var pgFile in individualPgFiles)
                {
                    var lines = File.ReadAllLines(pgFile);
                    var header = lines[0].Split('\t').ToList();
                    int unique = header.IndexOf("Unique Peptides");
                    int shared = header.IndexOf("Shared Peptides");
                    int numberUnique = header.IndexOf("Number of Unique Peptides");
                    var rows = lines.Skip(1).Select(l => l.Split('\t')).ToList();

                    Assert.That(rows.Count(r => r[unique].Length > 0), Is.GreaterThan(0), $"{pgFile} lists no unique peptide at all");
                    foreach (var row in rows)
                    {
                        int listedUnique = row[unique].Split('|', StringSplitOptions.RemoveEmptyEntries).Length;
                        Assert.That(listedUnique, Is.EqualTo(int.Parse(row[numberUnique])),
                            $"{pgFile} {row[0]}: the Unique Peptides cell and Number of Unique Peptides disagree");
                        Assert.That(row[unique].Length + row[shared].Length, Is.GreaterThan(0),
                            $"{pgFile} {row[0]}: a group with peptides lists none");
                    }
                }
            }
            finally
            {
                if (Directory.Exists(outputFolder)) Directory.Delete(outputFolder, true);
            }
        }

        /// <summary>
        /// Each per-file protein table lists the peptides identified in that file, not the peptides of the whole search.
        /// The test above searches two copies of one file, so its per-file lists equal the combined list and it cannot
        /// tell them apart. Here the two files identify different peptides of the same protein: protB's YSLTPEGK is in
        /// file 1 only and its GFHWSEAR in file 2 only, while LVNELTEFAK (also in protC) is in both.
        /// </summary>
        [Test]
        public static void PostSearchIndividualFileProteinGroupsListOnlyTheirOwnFilesPeptides()
        {
            SearchTask searchTask = new SearchTask()
            {
                CommonParameters = new CommonParameters(
                    scoreCutoff: 1,
                    digestionParams: new DigestionParams(
                        maxMissedCleavages: 0,
                        minPeptideLength: 5,
                        initiatorMethionineBehavior: InitiatorMethionineBehavior.Retain)),
                SearchParameters = new SearchParameters
                {
                    DoParsimony = true,
                    DoLabelFreeQuantification = false,
                    WriteIndividualFiles = true
                }
            };

            Protein protA = new Protein("AEDLNTKVQGLSASK", "protA");
            Protein protB = new Protein("YSLTPEGKLVNELTEFAKGFHWSEAR", "protB");
            Protein protC = new Protein("LVNELTEFAKDNWQSQR", "protC");
            var peptides = new[] { protA, protB, protC }
                .SelectMany(p => p.Digest(searchTask.CommonParameters.DigestionParams, new List<Modification>(), new List<Modification>()))
                .GroupBy(p => p.BaseSequence)
                .ToDictionary(g => g.Key, g => (IBioPolymerWithSetMods)g.First());

            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestIndividualPgOwnPeptides");
            string inputFolder = Path.Combine(outputFolder, "inputs");
            string dbPath = Path.Combine(inputFolder, "ownPeptides.xml");
            string mzml1 = Path.Combine(inputFolder, "ownPeptides1.mzML");
            string mzml2 = Path.Combine(inputFolder, "ownPeptides2.mzML");

            // accession -> (unique peptides, shared peptides), as one protein table lists them
            static Dictionary<string, (HashSet<string> Unique, HashSet<string> Shared)> ReadPeptideCells(string path)
            {
                var lines = File.ReadAllLines(path);
                var header = lines[0].Split('\t').ToList();
                int unique = header.IndexOf("Unique Peptides");
                int shared = header.IndexOf("Shared Peptides");
                return lines.Skip(1).Select(l => l.Split('\t')).ToDictionary(r => r[0], r => (
                    r[unique].Split('|', StringSplitOptions.RemoveEmptyEntries).ToHashSet(),
                    r[shared].Split('|', StringSplitOptions.RemoveEmptyEntries).ToHashSet()));
            }

            try
            {
                Directory.CreateDirectory(inputFolder);
                ProteinDbWriter.WriteXmlDatabase(new Dictionary<string, HashSet<Tuple<int, Modification>>>(), new List<Protein> { protA, protB, protC }, dbPath);
                Readers.MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(new TestDataFile(
                    new[] { "AEDLNTK", "VQGLSASK", "YSLTPEGK", "LVNELTEFAK" }.Select(s => peptides[s]).ToList()), mzml1, false);
                Readers.MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(new TestDataFile(
                    new[] { "GFHWSEAR", "LVNELTEFAK", "DNWQSQR" }.Select(s => peptides[s]).ToList()), mzml2, false);

                searchTask.RunTask(outputFolder, new List<DbForTask> { new DbForTask(dbPath, false) },
                    new List<string> { mzml1, mzml2 }, "normal");

                var combined = ReadPeptideCells(Path.Combine(outputFolder, "AllProteinGroups.tsv"));
                string individualFolder = Path.Combine(outputFolder, "Individual File Results");
                var file1 = ReadPeptideCells(Path.Combine(individualFolder, "ownPeptides1_ProteinGroups.tsv"));
                var file2 = ReadPeptideCells(Path.Combine(individualFolder, "ownPeptides2_ProteinGroups.tsv"));

                // the combined table lists protB's peptides from both files
                Assert.That(combined["protB"].Unique, Is.EquivalentTo(new[] { "YSLTPEGK", "GFHWSEAR" }));
                Assert.That(combined["protB"].Shared, Is.EquivalentTo(new[] { "LVNELTEFAK" }));

                // each per-file table lists only what its own file identified
                Assert.That(file1["protA"].Unique, Is.EquivalentTo(new[] { "AEDLNTK", "VQGLSASK" }));
                Assert.That(file1["protB"].Unique, Is.EquivalentTo(new[] { "YSLTPEGK" }));
                Assert.That(file1["protB"].Shared, Is.EquivalentTo(new[] { "LVNELTEFAK" }));
                Assert.That(file2["protB"].Unique, Is.EquivalentTo(new[] { "GFHWSEAR" }));
                Assert.That(file2["protB"].Shared, Is.EquivalentTo(new[] { "LVNELTEFAK" }));
                Assert.That(file2["protC"].Unique, Is.EquivalentTo(new[] { "DNWQSQR" }));

                // and no per-file row lists a peptide the other file alone identified
                Assert.That(file1.Values.SelectMany(v => v.Unique.Concat(v.Shared)).Intersect(new[] { "GFHWSEAR", "DNWQSQR" }), Is.Empty);
                Assert.That(file2.Values.SelectMany(v => v.Unique.Concat(v.Shared)).Intersect(new[] { "AEDLNTK", "VQGLSASK", "YSLTPEGK" }), Is.Empty);
            }
            finally
            {
                if (Directory.Exists(outputFolder)) Directory.Delete(outputFolder, true);
            }
        }

        /// <summary>
        /// Test that we don't get a crash if protein groups are not constructed
        /// </summary>
        [Test]
        public static void ProteinGroupsNoParsimonyTest()
        {
            SearchTask searchTask = new SearchTask()
            {
                SearchParameters = new SearchParameters
                {
                    DoParsimony = false
                },
            };

            string myFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\PrunedDbSpectra.mzml");
            string myDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\DbForPrunedDb.fasta");
            string folderPath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestProteinGroupsNoParsimony");

            DbForTask db = new DbForTask(myDatabase, false);
            Directory.CreateDirectory(folderPath);

            searchTask.RunTask(folderPath, new List<DbForTask> { db }, new List<string> { myFile }, "normal");
            Directory.Delete(folderPath, true);
        }

        /// <summary>
        /// Test ensures pruned databases are written when contaminant DB is searched
        /// </summary>
        [Test]
        public static void PrunedDbWithContaminantsTest()
        {
            SearchTask searchTask = new SearchTask()
            {
                SearchParameters = new SearchParameters
                {
                    WritePrunedDatabase = true
                },
            };

            string myFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\PrunedDbSpectra.mzml");
            string myDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\DbForPrunedDb.fasta");
            string folderPath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestNormalization");
            string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\ExperimentalDesign.tsv");

            // contaminant DB
            DbForTask db = new DbForTask(myDatabase, true);
            Directory.CreateDirectory(folderPath);

            searchTask.RunTask(folderPath, new List<DbForTask> { db }, new List<string> { myFile }, "normal");

            Assert.That(File.ReadAllLines(Path.Combine(folderPath, @"DbForPrunedDbproteinPruned.xml")).Length > 0);
            Assert.That(File.ReadAllLines(Path.Combine(folderPath, @"DbForPrunedDbPruned.xml")).Length > 0);
            Directory.Delete(folderPath, true);
        }

        /// <summary>
        /// Test ensures peptide FDR is calculated and that it doesn't output PSM FDR results
        /// </summary>
        [Test]
        public static void PeptideFDRTest()
        {
            string myFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\PrunedDbSpectra.mzml");
            string myFile2 = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\PrunedDbSpectra2.mzml");
            if (!File.Exists(myFile2)) { File.Copy(myFile, myFile2); }
            string myDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\DbForPrunedDb.fasta");
            string folderPath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestPeptideFDR");
            DbForTask db = new DbForTask(myDatabase, true);
            Directory.CreateDirectory(folderPath);

            // search something with multiple hits of the same peptide to see if peptide FDR is calculated at the end
            new SearchTask().RunTask(folderPath, new List<DbForTask> { db }, new List<string> { myFile, myFile2 }, "normal");
            List<string> columns = null;
            int cumDecoys = 0;
            int cumTargets = 0;
            double finalQValue = 0;
            foreach (string line in File.ReadAllLines(Path.Combine(folderPath, @"AllPeptides.psmtsv")))
            {
                string[] lineline = line.Split('\t');
                if (line.StartsWith("File Name")) // header
                {
                    columns = lineline.ToList();
                }
                else
                {
                    // since each PSM has a duplicate, these counts will be 1,3,5,7, etc. if peptide FDR isn't calculated
                    // if peptide FDR is calculated, they will be 1,2,3,4, etc. as expected
                    if (lineline[columns.IndexOf("Decoy/Contaminant/Target")] == "D")
                    {
                        Assert.That(++cumDecoys, Is.EqualTo(int.Parse(lineline[columns.IndexOf("Cumulative Decoy")])));
                    }
                    else
                    {
                        Assert.That(++cumTargets, Is.EqualTo(int.Parse(lineline[columns.IndexOf("Cumulative Target")])));
                    }

                    finalQValue = Math.Max(finalQValue, (double)cumDecoys / (double)cumTargets);
                }
                
            }

            // test that the final q-value follows the (target / decoy) formula
            // intermediate q-values no longer always follow this formula, so I'm not testing them here
            Assert.That(0.5, Is.EqualTo(finalQValue).Within(.0005));
            Directory.Delete(folderPath, true);
        }

        /// <summary>
        /// Tests interval mass difference acceptor type to make sure values are assigned properly
        /// </summary>
        [Test]
        public static void IntervalMassDiffAceptorTest()
        {
            var result = new IntervalMassDiffAcceptor("-187andUp", new List<DoubleRange> { new DoubleRange(-187, double.PositiveInfinity) });
            Assert.That(result.Accepts(2.0, 2.0) == 0);
            Assert.That(result.Accepts(double.PositiveInfinity, 2.0) == 0);
            result.GetAllowedPrecursorMassIntervalsFromTheoreticalMass(2.0);
            result.GetAllowedPrecursorMassIntervalsFromObservedMass(2.0);
            result.ToString();
            result.ToProseString();
            result = new IntervalMassDiffAcceptor("-187andUp", new List<DoubleRange> { new DoubleRange(-5, 0) });
            Assert.That(result.Accepts(2.0, 4.0) == 0);
            Assert.That(result.Accepts(2.0, 10.0) == -1);
        }

        /// <summary>
        /// Tests interval mass difference acceptor type to make sure values are assigned properly
        /// </summary>
        [Test]
        public static void SingleAbsoluteAroundZeroSearchMode()
        {
            var result = new SingleAbsoluteAroundZeroSearchMode(2.0);
            Assert.That(result.Accepts(2.0, 2.0) == 0);
            result.GetAllowedPrecursorMassIntervalsFromTheoreticalMass(2.0);
            result.GetAllowedPrecursorMassIntervalsFromObservedMass(2.0);
            result.ToString();
            result.ToProseString();
        }

        [Test]
        public static void TestMzIdentMlWriterWithUniprotResId()
        {
            Protein protein = new Protein("PEPTIDE", "", databaseFilePath: "temp");

            // The mod has to have a RESID reference and NO Unimod one, because GetUnimodCvParam prefers
            // Unimod and returns before it reaches the RESID branch this test exists to cover.
            // "FMN phosphoryl threonine on T" used to qualify. It stopped when the modification data moved
            // to the mzLib package, whose ptmlist.txt is a newer UniProt release that added
            // "DR Unimod; 442." to that entry -- so it now writes UNIMOD:442 and the RESID branch went
            // untested. This one still has RESID and PSI-MOD only, and targets the same residue.
            Modification uniProtMod = GlobalVariables.AllModsKnown.First(p =>
                p.IdWithMotif == "(E)-2,3-didehydrobutyrine on T"
                && p.ModificationType == "UniProt"
                && p.Target.ToString() == "T"
                && p.DatabaseReference.ContainsKey("RESID")
                && p.LocationRestriction == "Anywhere.");

            string resIdAccession = uniProtMod.DatabaseReference["RESID"].First();
            var peptide = protein.Digest(new DigestionParams(), new List<Modification> { uniProtMod }, new List<Modification>()).First();

            MsDataScan dfb = new MsDataScan(new MzSpectrum(new double[] { 1 }, new double[] { 1 }, false), 0, 1, true, Polarity.Positive, double.NaN, null,
                null, MZAnalyzerType.Orbitrap, double.NaN, null, null, "scan=1", double.NaN, null, null, double.NaN, null, DissociationType.AnyActivationType, 0, null);
            Ms2ScanWithSpecificMass scan = new Ms2ScanWithSpecificMass(dfb, 2, 0, "File", new CommonParameters());

            var psm = new PeptideSpectralMatch(peptide, 0, 1, 0, scan, new CommonParameters(), new List<MatchedFragmentIon>());
            psm.ResolveAllAmbiguities();
            psm.SetFdrValues(0, 0, 0, 0, 0, 0, 0, 0);

            string path = Path.Combine(TestContext.CurrentContext.TestDirectory, "ResIdOutput.mzID");

            MzIdentMLWriter.WriteMzIdentMl(new List<SpectralMatch> { psm }, new List<ProteinGroup>(), new List<Modification>(),
                new List<Modification>(), new List<SilacLabel>(), new List<DigestionAgent>(), new PpmTolerance(20), new PpmTolerance(20),
                0, path, true);

            var file = File.ReadAllLines(path);
            bool found = false;
            foreach (var line in file)
            {
                if (line.Contains("(E)-2,3-didehydrobutyrine on T") && line.Contains("RESID:" + resIdAccession))
                {
                    found = true;
                }
            }
            Assert.That(found);

            File.Delete(path);
        }

        [Test]
        public static void TestMzIdentMlWriterWithUniprotPsiMod()
        {
            Protein protein = new Protein("PEPTIDE", "", databaseFilePath: "temp");

            ModificationMotif.TryGetMotif("T", out var motif);

            Modification fakeMod = new Modification(_originalId: "FAKE", _accession: "FAKE_MOD_ACCESSION", _modificationType: "fake",
                _target: motif, _locationRestriction: "Anywhere.", _monoisotopicMass: 0,
                _databaseReference: new Dictionary<string, IList<string>> { { "PSI-MOD", new List<string> { "FAKE_MOD_ACCESSION" } } });

            string resIdAccession = fakeMod.DatabaseReference["PSI-MOD"].First();
            var peptide = protein.Digest(new DigestionParams(), new List<Modification> { fakeMod }, new List<Modification>()).First();

            MsDataScan dfb = new MsDataScan(new MzSpectrum(new double[] { 1 }, new double[] { 1 }, false), 0, 1, true, Polarity.Positive, double.NaN, null,
                null, MZAnalyzerType.Orbitrap, double.NaN, null, null, "scan=1", double.NaN, null, null, double.NaN, null, DissociationType.AnyActivationType, 0, null);
            Ms2ScanWithSpecificMass scan = new Ms2ScanWithSpecificMass(dfb, 2, 0, "File", new CommonParameters());

            var psm = new PeptideSpectralMatch(peptide, 0, 1, 0, scan, new CommonParameters(), new List<MatchedFragmentIon>());
            psm.ResolveAllAmbiguities();
            psm.SetFdrValues(0, 0, 0, 0, 0, 0, 0, 0);

            string path = Path.Combine(TestContext.CurrentContext.TestDirectory, "ResIdOutput.mzID");

            MzIdentMLWriter.WriteMzIdentMl(new List<SpectralMatch> { psm }, new List<ProteinGroup>(), new List<Modification>(),
                new List<Modification>(), new List<SilacLabel>(), new List<DigestionAgent>(), new PpmTolerance(20), new PpmTolerance(20),
                0, path, true);

            var file = File.ReadAllLines(path);
            bool found = false;
            foreach (var line in file)
            {
                if (line.Contains("FAKE on T") && line.Contains("PSI-MOD:" + resIdAccession))
                {
                    found = true;
                }
            }
            Assert.That(found);

            // test again w/ NOT appending motifs onto mod names
            MzIdentMLWriter.WriteMzIdentMl(new List<SpectralMatch> { psm }, new List<ProteinGroup>(), new List<Modification>(),
                new List<Modification>(), new List<SilacLabel>(), new List<DigestionAgent>(), new PpmTolerance(20), new PpmTolerance(20),
                0, path, false);

            file = File.ReadAllLines(path);
            found = false;
            foreach (var line in file)
            {
                if (line.Contains("FAKE") && !line.Contains(" on ") && line.Contains("PSI-MOD:" + resIdAccession))
                {
                    found = true;
                }
            }
            Assert.That(found);

            File.Delete(path);
        }

        /// <summary>
        /// 
        /// </summary>
        [Test]
        public static void TestAutodetectDissocationTypeFromScanHeader()
        {
            SearchTask searchTask = new SearchTask()
            {
                SearchParameters = new SearchParameters
                {
                    DoLabelFreeQuantification = false // quant disabled just to save some time
                },

                // use DissociationType.Autodetect as the dissociation type. this signals to the search that the dissociation type
                // should be taken from the scan header on a scan-specific basis
                CommonParameters = new CommonParameters(dissociationType: DissociationType.Autodetect, pepQValueThreshold: 0.01, qValueThreshold: 1.0)
            };

            string myFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\SmallCalibratible_Yeast.mzML");
            string myDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\smalldb.fasta");
            string folderPath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestAutodetectDissocationTypeFromScanHeader");

            DbForTask db = new DbForTask(myDatabase, false);

            // run the task
            var autoTaskFolder = Path.Combine(folderPath, @"Autodetect");
            Directory.CreateDirectory(autoTaskFolder);
            searchTask.RunTask(autoTaskFolder, new List<DbForTask> { db }, new List<string> { myFile }, "");

            // run identical task but select the HCD dissociation type instead of autodetect
            var hcdTaskFolder = Path.Combine(folderPath, @"HCD");
            Directory.CreateDirectory(hcdTaskFolder);
            searchTask = new SearchTask()
            {
                SearchParameters = new SearchParameters
                {
                    DoLabelFreeQuantification = false
                },

                CommonParameters = new CommonParameters(dissociationType: DissociationType.HCD)
            };

            searchTask.RunTask(hcdTaskFolder, new List<DbForTask> { db }, new List<string> { myFile }, "");

            // check search results
            var psmFileAutodetect = File.ReadAllLines(Path.Combine(autoTaskFolder, "AllPSMs.psmtsv"));
            var psmFileHcd = File.ReadAllLines(Path.Combine(hcdTaskFolder, "AllPSMs.psmtsv"));

            Assert.That(psmFileAutodetect.Length == psmFileHcd.Length);

            for (int i = 0; i < psmFileAutodetect.Length; i++)
            {
                Assert.That(psmFileAutodetect[i].Equals(psmFileHcd[i]));
            }

            // clean up
            Directory.Delete(folderPath, true);
        }

        [Test]
        public static void TestPepFilteringFewerThan100Psms()
        {
            SearchTask searchTask = new SearchTask()
            {
                SearchParameters = new SearchParameters
                {
                    DoLabelFreeQuantification = false // quant disabled just to save some time
                },

                // Use pepQvalue filtering with fewer than 100 psms
                CommonParameters = new CommonParameters(dissociationType: DissociationType.HCD,
                    pepQValueThreshold: 0.02, qValueThreshold: 1.0)
            };

            string myFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\SmallCalibratible_Yeast.mzML");
            string myDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\smalldb.fasta");
            string folderPath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestPepFiltering");

            DbForTask db = new DbForTask(myDatabase, false);

            // run the task
            var pepTaskFolder = Path.Combine(folderPath, @"pepTest");
            Directory.CreateDirectory(pepTaskFolder);
            searchTask.RunTask(pepTaskFolder, new List<DbForTask> { db },
                new List<string> { myFile }, "");

            string resultsFile = Path.Combine(pepTaskFolder, "results.txt");
            string[] results = File.ReadAllLines(resultsFile);
            // Each sentence sits just above the count it qualifies; the old text ("Results were filtered by q-value.")
            // said a filter ran.
            const string notFiltered = "PEP filtering was requested, but there were fewer than 100 PSMs (85, counting decoys), so results were not filtered.";
            Assert.That(results[5], Is.EqualTo(notFiltered));
            Assert.That(results[6], Is.EqualTo("All target PSMs with q-value <= 1: 84"));
            Assert.That(results[7], Is.EqualTo(notFiltered));
            Assert.That(results[8], Does.StartWith("All target peptides with q-value <= 1: "));

            // clean up
            Directory.Delete(folderPath, true);
        }
        /// <summary>
        /// Tests the contents of the results.txt file generated by a SearchTask to ensure that
        /// PSM, peptide, and protein group counts are correctly reported with the appropriate
        /// filter type and threshold values.
        /// 
        /// This test specifically validates the behavior when PEP Q-value filtering is requested
        /// (pepQValueThreshold: 0.02) but cannot be performed due to insufficient PSMs (fewer than 100).
        /// In this scenario, the system should fall back to standard Q-value filtering with a 
        /// threshold of 1.0 (as specified by qValueThreshold).
        /// 
        /// The test verifies:
        /// 1. Combined results across all files show correct counts with "q-value <= 1" format
        /// 2. Individual file results (mouseOne, mouseTwo) show correct per-file counts
        /// 3. Protein groups are reported using the same filter type and threshold as PSMs and peptides
        ///    (both combined and per-file results should use consistent formatting)
        /// 4. Precursor and MS2 scan counts are accurately reported for each file and combined
        /// 
        /// Expected output format in results.txt:
        /// - "All target PSMs with q-value <= 1: [count]"
        /// - "All target peptides with q-value <= 1: [count]"  
        /// - "All target protein groups with q-value <= 1: [count]"
        /// - "[fileName] - Target protein groups with q-value <= 1: [count]"
        /// </summary>
        [Test]
        public static void TestSearchTaskResultsTextContents()
        {
            SearchTask searchTask = new SearchTask()
            {
                SearchParameters = new SearchParameters
                {
                    DoLabelFreeQuantification = false, // quant disabled just to save some time
                    WriteIndividualFiles = true // the per-file counts are checked below
                },

                // Use pepQvalue filtering with fewer than 100 psms
                CommonParameters = new CommonParameters(dissociationType: DissociationType.HCD,
        pepQValueThreshold: 0.02, qValueThreshold: 1.0)
            };

            string myFileOne = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\mouseOne.mzML");
            string myFileTwo = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\mouseTwo.mzML");
            string myDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\mouseOne.xml");
            string folderPath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestPepResultsOutput");

            DbForTask db = new DbForTask(myDatabase, false);

            // run the task
            var pepTaskFolder = Path.Combine(folderPath, @"pepTest");
            Directory.CreateDirectory(pepTaskFolder);
            searchTask.RunTask(pepTaskFolder, new List<DbForTask> { db },
                new List<string> { myFileOne, myFileTwo }, "");

            string resultsFile = Path.Combine(pepTaskFolder, "results.txt");
            string[] results = File.ReadAllLines(resultsFile);

            // Lines 5 and 7 say PEP filtering was not done, each just above the count it qualifies.
            Assert.That(results[5], Does.StartWith("PEP filtering was requested, but there were fewer than 100 PSMs ("));
            Assert.That(results[6], Is.EqualTo("All target PSMs with q-value <= 1: 46"));
            Assert.That(results[7], Is.EqualTo(results[5]));
            Assert.That(results[8], Is.EqualTo("All target peptides with q-value <= 1: 14"));
            Assert.That(results[9], Is.EqualTo("All target protein groups with q-value <= 1: 14"));
            Assert.That(results[10], Is.EqualTo("All Precursors: 123"));
            Assert.That(results[11], Is.EqualTo("All MS2 Scans: 42"));

            Assert.That(results[13], Is.EqualTo("mouseOne - Target PSMs with q-value <= 1: 22"));
            Assert.That(results[14], Is.EqualTo("mouseOne - Target peptides with q-value <= 1: 14"));
            Assert.That(results[15], Is.EqualTo("mouseOne - Target protein groups with q-value <= 1: 14"));
            Assert.That(results[16], Is.EqualTo("mouseOne - Precursors: 57"));
            Assert.That(results[17], Is.EqualTo("mouseOne - MS2 Scans: 20"));

            Assert.That(results[19], Is.EqualTo("mouseTwo - Target PSMs with q-value <= 1: 24"));
            Assert.That(results[20], Is.EqualTo("mouseTwo - Target peptides with q-value <= 1: 14"));
            Assert.That(results[21], Is.EqualTo("mouseTwo - Target protein groups with q-value <= 1: 14"));
            Assert.That(results[22], Is.EqualTo("mouseTwo - Precursors: 66"));
            Assert.That(results[23], Is.EqualTo("mouseTwo - MS2 Scans: 22"));

            // clean up
            Directory.Delete(folderPath, true);
        }

        /// <summary>
        /// A PSM that is ambiguous when parsimony starts is excluded from it, but parsimony's own
        /// trimming can leave it unambiguous. Quantification then accepts it, so it must still be
        /// matched to the protein group its proteins landed in rather than reported as UNDEFINED.
        /// </summary>
        [Test]
        public static void QuantifiedPeptideResolvedByParsimonyGetsItsProteinGroup()
        {
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestIssue2218");
            string spectraFile = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "TaGe_SA_HeLa_04_subset_longestSeq.mzML");
            string database = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "hela_snip_for_unitTest.fasta");

            var searchTask = new SearchTask();

            new EverythingRunnerEngine(
                new List<(string, MetaMorpheusTask)> { ("Task", searchTask) },
                new List<string> { spectraFile },
                new List<DbForTask> { new DbForTask(database, false) },
                outputFolder).Run();

            string[] lines = File.ReadAllLines(Path.Combine(outputFolder, "Task", "AllQuantifiedPeptides.tsv"));
            int proteinGroupColumn = Array.IndexOf(lines[0].Split('\t'), "Protein Groups");
            Assert.That(proteinGroupColumn, Is.GreaterThanOrEqualTo(0));

            var undefined = lines.Skip(1)
                .Select(line => line.Split('\t'))
                .Where(cells => cells[proteinGroupColumn] == "UNDEFINED")
                .ToList();
            Assert.That(undefined, Is.Empty, "quantified peptides with no protein group: "
                + string.Join(", ", undefined.Select(cells => cells[0])));

            // this is the peptide parsimony disambiguates; it belongs to the beta-tubulin group
            var tubulin = lines.Skip(1)
                .Select(line => line.Split('\t'))
                .Single(cells => cells[1] == "EAESCDCLQGFQLTHSLGGGTGSGMGTLLISK");
            Assert.That(tubulin[proteinGroupColumn], Does.Contain("P04350"));

            Directory.Delete(outputFolder, true);
        }
    }
}
