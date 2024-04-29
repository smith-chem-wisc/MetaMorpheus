using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using GuiFunctions;
using MassSpectrometry;
using NUnit.Framework;
using Org.BouncyCastle.Bcpg;
using pepXML.Generated;

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
            Assert.False(viewModel.Persist);
            Assert.False(viewModel.UseInternalIons);
            Assert.That(viewModel.MinInternalIonLength, Is.EqualTo(10));
            Assert.That(viewModel.DissociationTypes.Count(), Is.EqualTo(7));

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
            }
        }
    }
}
