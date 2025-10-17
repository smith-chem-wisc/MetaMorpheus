using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using Easy.Common.Extensions;
using EngineLayer;
using GuiFunctions;
using MassSpectrometry;
using MzLibUtil;
using Nett;
using NUnit.Framework;
using NUnit.Framework.Legacy;
using Readers;
using TaskLayer;
using Transcriptomics.Digestion;
using UsefulProteomicsDatabases;

namespace Test.MetaDraw
{
    [ExcludeFromCodeCoverage]
    internal class FragmentReanalysis
    {
        [Test]
        public static void TestFragmentationReanalysisViewModel_DefaultProperties()
        {
            var viewModel = new FragmentationReanalysisViewModel();

            // check default parameters on loading
            Assert.That(viewModel.SelectedDissociationType, Is.EqualTo(DissociationType.HCD));
            Assert.That(!viewModel.Persist);
            Assert.That(!viewModel.UseInternalIons);
            Assert.That(viewModel.MinInternalIonLength, Is.EqualTo(10));
            Assert.That(viewModel.DissociationTypes.Count(), Is.EqualTo(7));
            Assert.That(viewModel.ProductIonMassTolerance, Is.EqualTo(20));

            var productsToUse = viewModel.PossibleProducts.Where(p => p.Use).Select(p => p.ProductType).ToList();
            var hcdProducts = Omics.Fragmentation.Peptide.DissociationTypeCollection.ProductsFromDissociationType[DissociationType.HCD];
            Assert.That(productsToUse.Count, Is.EqualTo(hcdProducts.Count));
            CollectionAssert.AreEqual(productsToUse, hcdProducts);
        }

        [Test]
        public static void TestFragmentationReanalysisViewModel_DissociationTypeSelectionChanged()
        {
            var viewModel = new FragmentationReanalysisViewModel();
            Assert.That(viewModel.SelectedDissociationType, Is.EqualTo(DissociationType.HCD));


            foreach (var dissociationType in viewModel.DissociationTypes)
            {
                var products = Omics.Fragmentation.Peptide.DissociationTypeCollection.ProductsFromDissociationType[dissociationType];
                viewModel.SelectedDissociationType = dissociationType;
                var productsToUse = viewModel.PossibleProducts.Where(p => p.Use).Select(p => p.ProductType).ToList();
                CollectionAssert.AreEquivalent(products, productsToUse);

                foreach (var viewModelPossibleProduct in viewModel.PossibleProducts)
                {
                    Assert.That(viewModelPossibleProduct.TypeString, Is.EqualTo(viewModelPossibleProduct.ProductType.ToString()));
                }
            }
        }

        [Test]
        public static void TestFragmentationReanalysisViewModel_DissociationTypeError()
        {
            var viewModel = new FragmentationReanalysisViewModel();
            Assert.That(viewModel.SelectedDissociationType, Is.EqualTo(DissociationType.HCD));
            viewModel.SelectedDissociationType = DissociationType.BIRD;
            Assert.That(viewModel.SelectedDissociationType, Is.EqualTo(DissociationType.HCD));

            viewModel = new FragmentationReanalysisViewModel(false);
            Assert.That(viewModel.SelectedDissociationType, Is.EqualTo(DissociationType.CID));
            viewModel.SelectedDissociationType = DissociationType.BIRD;
            Assert.That(viewModel.SelectedDissociationType, Is.EqualTo(DissociationType.HCD));
        }

        [Test]
        public static void TestFragmentationReanalysisViewModel_RematchIons()
        {
            var viewModel = new FragmentationReanalysisViewModel();

            // run a quick search 
            var myTomlPath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\Task1-SearchTaskconfig.toml");
            var searchTaskLoaded = Toml.ReadFile<SearchTask>(myTomlPath, MetaMorpheusTask.tomlConfig);
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\TestConsistency");
            string myFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\TaGe_SA_A549_3_snip.mzML");
            string myDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\TaGe_SA_A549_3_snip.fasta");
            var engineToml = new EverythingRunnerEngine(new List<(string, MetaMorpheusTask)> { ("SearchTOML", searchTaskLoaded) }, new List<string> { myFile }, new List<DbForTask> { new DbForTask(myDatabase, false) }, outputFolder);
            engineToml.Run();
            string psmFile = Path.Combine(outputFolder, @"SearchTOML\AllPSMs.psmtsv");
            var dataFile = MsDataFileReader.GetDataFile(myFile);

            // parse out psm and its respective scan
            List<PsmFromTsv> parsedPsms = SpectrumMatchTsvReader.ReadPsmTsv(psmFile, out var warnings);
            var psmToResearch = parsedPsms.First();
            var scan = dataFile.GetOneBasedScan(psmToResearch.Ms2ScanNumber);

            // increase product ion tolerance and ensure more ions are found
            viewModel.ProductIonMassTolerance = 200;
            var newMatchedIons = viewModel.MatchIonsWithNewTypes(scan, psmToResearch);
            Assert.That(psmToResearch.MatchedIons.Count < newMatchedIons.Count);

            // decrease product ion tolerance and ensure fewer ions are found
            viewModel.ProductIonMassTolerance = 2;
            newMatchedIons = viewModel.MatchIonsWithNewTypes(scan, psmToResearch);
            Assert.That(newMatchedIons.Count < psmToResearch.MatchedIons.Count);
            viewModel.ProductIonMassTolerance = 150; // Set big to ensure intersection with XCorr results. 

            // ensure ambiguous psms get kicked out before reanalysis
            var ambiguousPsm = parsedPsms.First(p => p.FullSequence.Contains('|'));
            var newIons = viewModel.MatchIonsWithNewTypes(scan, ambiguousPsm);
            CollectionAssert.AreEqual(ambiguousPsm.MatchedIons, newIons);

            viewModel.PossibleProducts.ForEach(p => p.Use = true);
            newMatchedIons = viewModel.MatchIonsWithNewTypes(scan, psmToResearch);

            // searching with additional ions should yield more 
            Assert.That(psmToResearch.MatchedIons.Count, Is.LessThan(newMatchedIons.Count));
            // all original ions should be retained
            var intersect = newMatchedIons.Select(p => p.Annotation).Intersect(psmToResearch.MatchedIons.Select(p => p.Annotation));
            Assert.That(intersect.Count(), Is.EqualTo(psmToResearch.MatchedIons.Count));

            // perform the same operation but also with internal ions
            viewModel.UseInternalIons = true;
            viewModel.MinInternalIonLength = 3;

            // searching with even more ions should yield even more results but still contain all original results
            var internalIonNewIons = viewModel.MatchIonsWithNewTypes(scan, psmToResearch);
            Assert.That(psmToResearch.MatchedIons.Count, Is.LessThan(internalIonNewIons.Count));
            Assert.That(newMatchedIons.Count, Is.LessThan(internalIonNewIons.Count));

            // all original ions should be retained
            intersect = internalIonNewIons.Select(p => p.Annotation).Intersect(psmToResearch.MatchedIons.Select(p => p.Annotation));
            Assert.That(intersect.Count(), Is.EqualTo(psmToResearch.MatchedIons.Count));

            // clean up
            Directory.Delete(outputFolder, true);
        }

        [Test]
        [NonParallelizable]
        public static void TestFragmentationReanalysisViewModel_RematchIons_RNA()
        {
            GlobalVariables.AnalyteType = AnalyteType.Oligo;
            var viewModel = new FragmentationReanalysisViewModel(false);
            viewModel.SelectedDissociationType = DissociationType.CID;

            // run a quick search 
            var searchTaskLoaded = new SearchTask()
            {
                SearchParameters = new RnaSearchParameters
                {
                    DecoyType = DecoyType.Reverse,
                    MassDiffAcceptorType = MassDiffAcceptorType.Custom,
                    CustomMdac = "Custom interval [-5,5]",
                    DisposeOfFileWhenDone = true
                },
                CommonParameters = new CommonParameters
                (
                    dissociationType: DissociationType.CID,
                    deconvolutionMaxAssumedChargeState: -20,
                    deconvolutionIntensityRatio: 3,
                    deconvolutionMassTolerance: new PpmTolerance(4),
                    precursorMassTolerance: new PpmTolerance(5),
                    productMassTolerance: new PpmTolerance(20),
                    scoreCutoff: 5,
                    totalPartitions: 1,
                    maxThreadsToUsePerFile: 1,
                    doPrecursorDeconvolution: true,
                    useProvidedPrecursorInfo: false,
                    digestionParams: new RnaDigestionParams()
                ),
            };
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\TestRefragmentRNA");
            var myFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"Transcriptomics\TestData\GUACUG_NegativeMode_Sliced.mzML");
            string myDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, @"Transcriptomics\TestData\6mer.fasta");
            var engineToml = new EverythingRunnerEngine(new List<(string, MetaMorpheusTask)> { ("SearchTOML", searchTaskLoaded) }, new List<string> { myFile }, new List<DbForTask> { new DbForTask(myDatabase, false) }, outputFolder);
            engineToml.Run();
            string psmFile = Path.Combine(outputFolder, @"SearchTOML\AllOSMs.osmtsv");
            var dataFile = MsDataFileReader.GetDataFile(myFile);

            // parse out psm and its respective scan
            List<OsmFromTsv> parsedOsms = SpectrumMatchTsvReader.ReadOsmTsv(psmFile, out var warnings);
            var psmToResearch = parsedOsms.First();
            var scan = dataFile.GetOneBasedScan(psmToResearch.Ms2ScanNumber);

            // increase product ion tolerance and ensure more ions are found
            viewModel.ProductIonMassTolerance = 200;
            var newMatchedIons = viewModel.MatchIonsWithNewTypes(scan, psmToResearch);
            Assert.That(psmToResearch.MatchedIons.Count < newMatchedIons.Count);

            // decrease product ion tolerance and ensure fewer ions are found
            viewModel.ProductIonMassTolerance = 2;
            newMatchedIons = viewModel.MatchIonsWithNewTypes(scan, psmToResearch);
            Assert.That(newMatchedIons.Count < psmToResearch.MatchedIons.Count);
            viewModel.ProductIonMassTolerance = 20;

            viewModel.PossibleProducts.ForEach(p => p.Use = true);
            newMatchedIons = viewModel.MatchIonsWithNewTypes(scan, psmToResearch);

            // clean up
            Directory.Delete(outputFolder, true);
            GlobalVariables.AnalyteType = AnalyteType.Peptide;
        }
    }
}
