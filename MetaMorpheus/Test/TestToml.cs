using Chemistry;
using EngineLayer;
using EngineLayer.ClassicSearch;
using MassSpectrometry;
using MzLibUtil;
using Nett;
using NUnit.Framework; using Assert = NUnit.Framework.Legacy.ClassicAssert;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using Omics.Digestion;
using Omics.Fragmentation;
using SpectralAveraging;
using TaskLayer;
using UsefulProteomicsDatabases;

namespace Test
{
    [TestFixture]
    public static class TestToml
    {
        [Test]
        public static void TestTomlFunction()
        {
            SearchTask searchTask = new()
            {
                CommonParameters = new CommonParameters(
                    productMassTolerance: new PpmTolerance(666),
                    listOfModsFixed: new List<(string, string)> { ("Common Fixed", "Carbamidomethyl on C"), ("Common Fixed", "Carbamidomethyl on U") },
                    listOfModsVariable: new List<(string, string)> { ("Common Variable", "Oxidation on M") }),
            };
            Toml.WriteFile(searchTask, "SearchTask.toml", MetaMorpheusTask.tomlConfig);
            var searchTaskLoaded = Toml.ReadFile<SearchTask>("SearchTask.toml", MetaMorpheusTask.tomlConfig);

            Assert.AreEqual(searchTask.CommonParameters.DeconvolutionMassTolerance.ToString(), searchTaskLoaded.CommonParameters.DeconvolutionMassTolerance.ToString());
            Assert.AreEqual(searchTask.CommonParameters.ProductMassTolerance.ToString(), searchTaskLoaded.CommonParameters.ProductMassTolerance.ToString());
            Assert.AreEqual(searchTask.CommonParameters.PrecursorMassTolerance.ToString(), searchTaskLoaded.CommonParameters.PrecursorMassTolerance.ToString());

            Assert.AreEqual(searchTask.CommonParameters.ListOfModsFixed.Count(), searchTaskLoaded.CommonParameters.ListOfModsFixed.Count());
            Assert.AreEqual(searchTask.CommonParameters.ListOfModsFixed.First().Item1, searchTaskLoaded.CommonParameters.ListOfModsFixed.First().Item1);
            Assert.AreEqual(searchTask.CommonParameters.ListOfModsFixed.First().Item2, searchTaskLoaded.CommonParameters.ListOfModsFixed.First().Item2);

            Assert.AreEqual(searchTask.CommonParameters.ListOfModsVariable.Count(), searchTaskLoaded.CommonParameters.ListOfModsVariable.Count());

            Assert.AreEqual(searchTask.SearchParameters.MassDiffAcceptorType, searchTaskLoaded.SearchParameters.MassDiffAcceptorType);
            Assert.AreEqual(searchTask.SearchParameters.CustomMdac, searchTaskLoaded.SearchParameters.CustomMdac);

            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestConsistency");
            string myFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\PrunedDbSpectra.mzml");
            string myDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\DbForPrunedDb.fasta");

            var engine = new EverythingRunnerEngine(new List<(string, MetaMorpheusTask)> { ("Search", searchTask) }, new List<string> { myFile }, new List<DbForTask> { new DbForTask(myDatabase, false) }, outputFolder);
            engine.Run();
            var engineToml = new EverythingRunnerEngine(new List<(string, MetaMorpheusTask)> { ("SearchTOML", searchTaskLoaded) }, new List<string> { myFile }, new List<DbForTask> { new DbForTask(myDatabase, false) }, outputFolder);
            engineToml.Run();

            var results = File.ReadAllLines(Path.Combine(outputFolder, @"Search\AllPSMs.psmtsv"));
            var resultsToml = File.ReadAllLines(Path.Combine(outputFolder, @"SearchTOML\AllPSMs.psmtsv"));
            Assert.That(results.SequenceEqual(resultsToml));

            CalibrationTask calibrationTask = new();
            Toml.WriteFile(calibrationTask, "CalibrationTask.toml", MetaMorpheusTask.tomlConfig);
            _ = Toml.ReadFile<CalibrationTask>("CalibrationTask.toml", MetaMorpheusTask.tomlConfig);

            GptmdTask gptmdTask = new GptmdTask();
            Toml.WriteFile(gptmdTask, "GptmdTask.toml", MetaMorpheusTask.tomlConfig);
            var gptmdTaskLoaded = Toml.ReadFile<GptmdTask>("GptmdTask.toml", MetaMorpheusTask.tomlConfig);

            var gptmdEngine = new EverythingRunnerEngine(new List<(string, MetaMorpheusTask)> { ("GPTMD", gptmdTask) }, new List<string> { myFile }, new List<DbForTask> { new DbForTask(myDatabase, false) }, outputFolder);
            gptmdEngine.Run();
            var gptmdEngineToml = new EverythingRunnerEngine(new List<(string, MetaMorpheusTask)> { ("GPTMDTOML", gptmdTaskLoaded) }, new List<string> { myFile }, new List<DbForTask> { new DbForTask(myDatabase, false) }, outputFolder);
            gptmdEngineToml.Run();

            XLSearchTask xLSearchTask = new();
            Toml.WriteFile(xLSearchTask, "XLSearchTask.toml", MetaMorpheusTask.tomlConfig);
            var xLSearchTaskLoaded = Toml.ReadFile<XLSearchTask>("XLSearchTask.toml", MetaMorpheusTask.tomlConfig);

            string myFileXl = Path.Combine(TestContext.CurrentContext.TestDirectory, @"XlTestData\BSA_DSSO_ETchD6010.mgf");
            string myDatabaseXl = Path.Combine(TestContext.CurrentContext.TestDirectory, @"XlTestData\BSA.fasta");

            var xlEngine = new EverythingRunnerEngine(new List<(string, MetaMorpheusTask)> { ("XLSearch", xLSearchTask) }, new List<string> { myFileXl }, new List<DbForTask> { new DbForTask(myDatabaseXl, false) }, outputFolder);
            xlEngine.Run();
            var xlEngineToml = new EverythingRunnerEngine(new List<(string, MetaMorpheusTask)> { ("XLSearchTOML", xLSearchTaskLoaded) }, new List<string> { myFileXl }, new List<DbForTask> { new DbForTask(myDatabaseXl, false) }, outputFolder);
            xlEngineToml.Run();

            var xlResults = File.ReadAllLines(Path.Combine(outputFolder, @"XLSearch\XL_Intralinks.tsv"));
            var xlResultsToml = File.ReadAllLines(Path.Combine(outputFolder, @"XLSearchTOML\XL_Intralinks.tsv"));

            SpectralAveragingTask averagingTask = new SpectralAveragingTask(new SpectralAveragingParameters()
                { NumberOfScansToAverage = 15, BinSize = 2});
            Toml.WriteFile(averagingTask, "averagingTask.toml", MetaMorpheusTask.tomlConfig);
            var readInAveragingTask = Toml.ReadFile<SpectralAveragingTask>("averagingTask.toml", MetaMorpheusTask.tomlConfig);
            Assert.That(averagingTask.Parameters.NumberOfScansToAverage, Is.EqualTo(readInAveragingTask.Parameters.NumberOfScansToAverage).And.EqualTo(15));
            Assert.That(averagingTask.Parameters.BinSize, Is.EqualTo(readInAveragingTask.Parameters.BinSize).And.EqualTo(2));


            Assert.That(xlResults.SequenceEqual(xlResultsToml));
            Directory.Delete(outputFolder, true);
            File.Delete(Path.Combine(TestContext.CurrentContext.TestDirectory, @"GptmdTask.toml"));
            File.Delete(Path.Combine(TestContext.CurrentContext.TestDirectory, @"XLSearchTask.toml"));
            File.Delete(Path.Combine(TestContext.CurrentContext.TestDirectory, @"SearchTask.toml"));
            File.Delete(Path.Combine(TestContext.CurrentContext.TestDirectory, @"CalibrationTask.toml"));
            File.Delete(Path.Combine(TestContext.CurrentContext.TestDirectory, @"averagingTask.toml"));
        }

        [Test]
        public static void TestTomlForSpecficFiles()
        {
            var fileSpecificToml = Toml.ReadFile(Path.Combine(TestContext.CurrentContext.TestDirectory, "testFileSpecfic.toml"), MetaMorpheusTask.tomlConfig);
            var tomlSettingsList = fileSpecificToml.ToDictionary(p => p.Key);
            Assert.AreEqual(tomlSettingsList["Protease"].Value.Get<string>(), "Asp-N");
            Assert.AreEqual(tomlSettingsList["DissociationType"].Value.Get<string>(), "ETD");
            Assert.IsFalse(tomlSettingsList.ContainsKey("maxMissedCleavages"));
            Assert.IsFalse(tomlSettingsList.ContainsKey("InitiatorMethionineBehavior"));

            FileSpecificParameters f = new(fileSpecificToml);

            Assert.AreEqual("Asp-N", f.Protease.Name);
            Assert.AreEqual(DissociationType.ETD, f.DissociationType);
            Assert.IsNull(f.MaxMissedCleavages);

            CommonParameters c = MetaMorpheusTask.SetAllFileSpecificCommonParams(new CommonParameters(), f);

            Assert.AreEqual("Asp-N", c.DigestionParams.Protease.Name);
            Assert.AreEqual(DissociationType.ETD, c.DissociationType);
            Assert.AreEqual(2, c.DigestionParams.MaxMissedCleavages);
        }

        [Test]
        public static void TestBadFileSpecificProtease()
        {
            //this test checks for a catch statement (or some other handling) for file-specific toml loading
            //create a toml with a protease that doesn't exist in the protease.tsv dictionary
            string proteaseNotInDictionary = "aaa"; //arbitrary. If somebody adds a protease with this name, use a different name
            string proteaseInDictionary = "trypsin"; //just make sure we are doing this right
            Assert.IsTrue(!ProteaseDictionary.Dictionary.ContainsKey(proteaseNotInDictionary));
            Assert.IsTrue(ProteaseDictionary.Dictionary.ContainsKey(proteaseInDictionary));

            //write the toml
            //let's use the datafile ok.mgf (arbitrary)
            File.WriteAllLines(Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\ok.toml"), new string[] { "Protease = " + '"' + proteaseNotInDictionary + '"' });

            //create a task with this, we want the run to work and just ignore the bad toml
            SearchTask task = new SearchTask();
            //just test it doesn't crash (i.e. the crash is handled)
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"BadProteaseTest");
            Directory.CreateDirectory(outputFolder);
            task.RunTask(outputFolder, 
                new List<DbForTask> { new DbForTask(Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\okk.xml"), false) },
                new List<string>{Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\ok.mgf") },
                outputFolder);

            //Clear result files
            Directory.Delete(outputFolder, true);
        }

        [Test]
        public static void TestOldModsInToml()
        {
            //These mods are not formatted correctly. This should throw a warning in MetaMorpheusTask
            List<(string, string)> variableMods = new() {("","") };
            List<(string, string)> fixedMods = new() { ("", "") };

            CommonParameters commonParameters = new(listOfModsFixed: fixedMods, listOfModsVariable: variableMods);

            SearchTask searchTask = new()
            {
                SearchParameters = new SearchParameters
                {
                    Normalize = true
                },
                CommonParameters = commonParameters
            };

            string myFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\PrunedDbSpectra.mzml");
            string myDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\DbForPrunedDb.fasta");
            string folderPath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestNormalizationExperDesign");
            string experimentalDesignFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\ExperimentalDesign.tsv");
            using (StreamWriter output = new StreamWriter(experimentalDesignFile))
            {
                output.WriteLine("FileName\tCondition\tBiorep\tFraction\tTechrep");
                output.WriteLine("PrunedDbSpectra.mzml" + "\t" + "condition" + "\t" + "1" + "\t" + "1" + "\t" + "1");
            }
            DbForTask db = new(myDatabase, false);

            // run the task
            Directory.CreateDirectory(folderPath);
            

            bool wasCalled = false;
            MetaMorpheusTask.WarnHandler += (o, e) => wasCalled = true;

            //Executing this run task will load modifications in the MetaMorpheusTask. Mods are written incorrectly, which will throw a warning.
            searchTask.RunTask(folderPath, new List<DbForTask> { db }, new List<string> { myFile }, "normal");

            Directory.Delete(folderPath, true);

            //
            Assert.IsTrue(wasCalled);
        }

        [Test]
        public static void FileSpecificParametersTest()
        {
            var filePath = Path.Combine(TestContext.CurrentContext.TestDirectory, "testFileParams.toml");

            var fileSpecificToml = Toml.ReadFile(filePath, MetaMorpheusTask.tomlConfig);

            FileSpecificParameters fsp = new(fileSpecificToml);
            Assert.AreEqual(DissociationType.CID, fsp.DissociationType);
            Assert.AreEqual(0, fsp.MaxMissedCleavages);
            Assert.AreEqual(0, fsp.MaxModsForPeptide);
            Assert.AreEqual(0, fsp.MaxPeptideLength);
            Assert.AreEqual(0, fsp.MinPeptideLength);
            Assert.AreEqual(5.0d, fsp.PrecursorMassTolerance.Value);
            Assert.AreEqual(5.0d, fsp.ProductMassTolerance.Value);
            Assert.AreEqual("Asp-N", fsp.Protease.Name);
            Assert.AreEqual("HPLC", fsp.SeparationType.ToString());

            filePath = Path.Combine(TestContext.CurrentContext.TestDirectory, "testFileParams_bad.toml");

            var fileSpecificTomlBad = Toml.ReadFile(filePath, MetaMorpheusTask.tomlConfig);

            Assert.Throws<MetaMorpheusException>(() => new FileSpecificParameters(fileSpecificTomlBad));
        }

        /// <summary>
        /// This test ensures that all properties within FileSpecificParams shares both a name and a type with
        /// a properties in CommonParameters or DigestionParameters and replaces the untestable method
        /// ValidateFileSpecificVariableNames from FileSpecificParameters
        /// </summary>
        [Test]
        public static void TestFileSpecificAndCommonParametersNameEquality()
        {
            var filePath = Path.Combine(TestContext.CurrentContext.TestDirectory, "testFileParams.toml");
            var fileSpecificToml = Toml.ReadFile(filePath, MetaMorpheusTask.tomlConfig);

            FileSpecificParameters fileSpecificParameters = new(fileSpecificToml);
            CommonParameters commonParameters = new();

            // foreach property in File Specific Parameters, ensure common parameters has a property with the same name
            foreach (var fileSpecificProperty in fileSpecificParameters.GetType().GetProperties())
            {
                string fileSpecificPropertyName = fileSpecificProperty.Name;
                var commonProperty = commonParameters.GetType().GetProperty(fileSpecificPropertyName);
                if (commonProperty is null)
                {
                    // if not found in common parameters, check digestion parameters
                    commonProperty = commonParameters.DigestionParams.GetType().GetProperty(fileSpecificPropertyName);
                    if (commonProperty is null)
                        Assert.Fail("Common Parameters does not have a property with the name " + fileSpecificProperty);
                }

                string commonPropertyName = commonProperty.Name;

                Assert.That(commonPropertyName, Is.EqualTo(fileSpecificPropertyName));
                Assert.That(commonProperty.GetType(), Is.EqualTo(fileSpecificProperty.GetType()));
            }
        }

        [Test]
        public static void TestDigestionParamsTomlReadingWriting()
        {
            var digestionParams = new DigestionParams("top-down", 4, 5, 12345, 2012,
                InitiatorMethionineBehavior.Undefined, 4, CleavageSpecificity.Semi, FragmentationTerminus.Both, false,
                true, true);
            var commonParams = new CommonParameters(digestionParams: digestionParams);
            var searchTask = new SearchTask();
            searchTask.CommonParameters = commonParams;

            // check that digestion params are correct in search task
            Assert.That(searchTask.CommonParameters.DigestionParams.FragmentationTerminus, Is.EqualTo(digestionParams.FragmentationTerminus));
            Assert.That(searchTask.CommonParameters.DigestionParams.SearchModeType, Is.EqualTo(digestionParams.SearchModeType));
            Assert.That(searchTask.CommonParameters.DigestionParams.InitiatorMethionineBehavior, Is.EqualTo(digestionParams.InitiatorMethionineBehavior));
            Assert.That(searchTask.CommonParameters.DigestionParams.MaxMissedCleavages, Is.EqualTo(digestionParams.MaxMissedCleavages));
            Assert.That(searchTask.CommonParameters.DigestionParams.MaxModificationIsoforms, Is.EqualTo(digestionParams.MaxModificationIsoforms));
            Assert.That(searchTask.CommonParameters.DigestionParams.MinPeptideLength, Is.EqualTo(digestionParams.MinPeptideLength));
            Assert.That(searchTask.CommonParameters.DigestionParams.MaxPeptideLength, Is.EqualTo(digestionParams.MaxPeptideLength));
            Assert.That(searchTask.CommonParameters.DigestionParams.Protease.Name, Is.EqualTo(digestionParams.Protease.Name));
            Assert.That(searchTask.CommonParameters.DigestionParams.GeneratehUnlabeledProteinsForSilac, Is.EqualTo(digestionParams.GeneratehUnlabeledProteinsForSilac));
            Assert.That(searchTask.CommonParameters.DigestionParams.KeepNGlycopeptide, Is.EqualTo(digestionParams.KeepNGlycopeptide));
            Assert.That(searchTask.CommonParameters.DigestionParams.KeepOGlycopeptide, Is.EqualTo(digestionParams.KeepOGlycopeptide));

            // write and read file 
            string filePath = Path.Combine(TestContext.CurrentContext.TestDirectory, "testDigestionParams.toml");
            Toml.WriteFile(searchTask, filePath, MetaMorpheusTask.tomlConfig);
            var searchTaskLoaded = Toml.ReadFile<SearchTask>(filePath, MetaMorpheusTask.tomlConfig);

            // check that digestion params are correct in search task
            Assert.That(searchTaskLoaded.CommonParameters.DigestionParams.FragmentationTerminus, Is.EqualTo(digestionParams.FragmentationTerminus));
            Assert.That(searchTaskLoaded.CommonParameters.DigestionParams.SearchModeType, Is.EqualTo(digestionParams.SearchModeType));
            Assert.That(searchTaskLoaded.CommonParameters.DigestionParams.InitiatorMethionineBehavior, Is.EqualTo(digestionParams.InitiatorMethionineBehavior));
            Assert.That(searchTaskLoaded.CommonParameters.DigestionParams.MaxMissedCleavages, Is.EqualTo(digestionParams.MaxMissedCleavages));
            Assert.That(searchTaskLoaded.CommonParameters.DigestionParams.MaxModificationIsoforms, Is.EqualTo(digestionParams.MaxModificationIsoforms));
            Assert.That(searchTaskLoaded.CommonParameters.DigestionParams.MinPeptideLength, Is.EqualTo(digestionParams.MinPeptideLength));
            Assert.That(searchTaskLoaded.CommonParameters.DigestionParams.MaxPeptideLength, Is.EqualTo(digestionParams.MaxPeptideLength));
            Assert.That(searchTaskLoaded.CommonParameters.DigestionParams.Protease.Name, Is.EqualTo(digestionParams.Protease.Name));
            Assert.That(searchTaskLoaded.CommonParameters.DigestionParams.GeneratehUnlabeledProteinsForSilac, Is.EqualTo(digestionParams.GeneratehUnlabeledProteinsForSilac));
            Assert.That(searchTaskLoaded.CommonParameters.DigestionParams.KeepNGlycopeptide, Is.EqualTo(digestionParams.KeepNGlycopeptide));
            Assert.That(searchTaskLoaded.CommonParameters.DigestionParams.KeepOGlycopeptide, Is.EqualTo(digestionParams.KeepOGlycopeptide));

            File.Delete(filePath);
        }
        
    }
}