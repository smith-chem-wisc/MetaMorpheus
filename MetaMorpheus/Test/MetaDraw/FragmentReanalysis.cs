using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Easy.Common.Extensions;
using EngineLayer;
using GuiFunctions;
using MassSpectrometry;
using Nett;
using NUnit.Framework;
using NUnit.Framework.Legacy;
using Assert = NUnit.Framework.Legacy.ClassicAssert;
using Omics.Fragmentation;
using pepXML.Generated;
using Readers;
using TaskLayer;

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
            Assert.AreEqual(DissociationType.HCD, viewModel.SelectedDissociationType);
            Assert.False(viewModel.Persist);
            Assert.False(viewModel.UseInternalIons);
            Assert.AreEqual(10, viewModel.MinInternalIonLength);
            Assert.AreEqual(7, viewModel.DissociationTypes.Count());
            Assert.AreEqual(20, viewModel.ProductIonMassTolerance);

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
            Assert.That(viewModel.SelectedDissociationType, Is.EqualTo(DissociationType.LowCID));
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
            List<PsmFromTsv> parsedPsms = PsmTsvReader.ReadTsv(psmFile, out var warnings);
            var psmToResearch = parsedPsms.First();
            var scan = dataFile.GetOneBasedScan(psmToResearch.Ms2ScanNumber);

            // increase product ion tolerance and ensure more ions are found
            viewModel.ProductIonMassTolerance = 200;
            var newMatchedIons = viewModel.MatchIonsWithNewTypes(scan, psmToResearch);
            Assert.Less(psmToResearch.MatchedIons.Count, newMatchedIons.Count);

            // decrease product ion tolerance and ensure fewer ions are found
            viewModel.ProductIonMassTolerance = 2;
            newMatchedIons = viewModel.MatchIonsWithNewTypes(scan, psmToResearch);
            Assert.Less(newMatchedIons.Count, psmToResearch.MatchedIons.Count);
            viewModel.ProductIonMassTolerance = 20;

            // ensure ambiguous psms get kicked out before reanalysis
            var ambiguousPsm = parsedPsms.First(p => p.FullSequence.Contains('|'));
            var newIons = viewModel.MatchIonsWithNewTypes(scan, ambiguousPsm);
            CollectionAssert.AreEqual(ambiguousPsm.MatchedIons, newIons);

            viewModel.PossibleProducts.ForEach(p => p.Use = true);
            newMatchedIons = viewModel.MatchIonsWithNewTypes(scan, psmToResearch);

            // searching with additional ions should yield more 
            Assert.That(psmToResearch.MatchedIons.Count, Is.LessThan(newMatchedIons.Count));
            // all original ions should be retained
            var intersect = newMatchedIons.Intersect(psmToResearch.MatchedIons);
            Assert.That(intersect.Count(), Is.EqualTo(psmToResearch.MatchedIons.Count));

            // perform the same operation but also with internal ions
            viewModel.UseInternalIons = true;
            viewModel.MinInternalIonLength = 3;

            // searching with even more ions should yield even more results but still contain all original results
            var internalIonNewIons = viewModel.MatchIonsWithNewTypes(scan, psmToResearch);
            Assert.That(psmToResearch.MatchedIons.Count, Is.LessThan(internalIonNewIons.Count));
            Assert.That(newMatchedIons.Count, Is.LessThan(internalIonNewIons.Count));

            // all original ions should be retained
            intersect = internalIonNewIons.Intersect(psmToResearch.MatchedIons);
            Assert.That(intersect.Count(), Is.EqualTo(psmToResearch.MatchedIons.Count));

            // clean up
            Directory.Delete(outputFolder, true);
        }
    }
}
