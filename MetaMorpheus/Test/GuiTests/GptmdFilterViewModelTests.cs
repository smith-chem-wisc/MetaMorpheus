using System;
using GuiFunctions;
using EngineLayer;
using NUnit.Framework;
using Omics;
using System.Collections.Generic;
using Omics.Fragmentation;
using Omics.Modifications;

namespace Test.GuiTests
{
    public class DummyGptmdFilter : IGptmdFilter
    {
        public static string GetFilterTypeName(IGptmdFilter filter) => "Dummy";
        public bool Passes(IBioPolymerWithSetMods candidatePeptide, SpectralMatch psm, double newScore, double originalScore, List<MatchedFragmentIon> matchedIons, int peptideOneBasedModSite, int peptideLength, Modification modAttemptingToAdd) => true;
    }

    [TestFixture]
    public class GptmdFilterViewModelTests
    {
        [Test]
        public void Constructor_SetsPropertiesCorrectly()
        {
            var filter = new DummyGptmdFilter();
            var name = "TestFilter";
            var summary = "A test summary";
            var viewModel = new GptmdFilterViewModel(filter, name, summary, isSelected: false);

            Assert.That(viewModel.Filter, Is.EqualTo(filter));
            Assert.That(viewModel.Name, Is.EqualTo(name));
            Assert.That(viewModel.Summary, Is.EqualTo(summary));
            Assert.That(viewModel.IsSelected, Is.False);
        }

        [Test]
        public void Constructor_DefaultIsSelectedTrue()
        {
            var filter = new DummyGptmdFilter();
            var viewModel = new GptmdFilterViewModel(filter, "n", "s");

            Assert.That(viewModel.IsSelected, Is.True);
        }

        [Test]
        public void IsSelected_Property_SetAndGet()
        {
            var filter = new DummyGptmdFilter();
            var viewModel = new GptmdFilterViewModel(filter, "n", "s", isSelected: false);

            Assert.That(viewModel.IsSelected, Is.False);

            viewModel.IsSelected = true;
            Assert.That(viewModel.IsSelected, Is.True);

            viewModel.IsSelected = false;
            Assert.That(viewModel.IsSelected, Is.False);
        }

        [Test]
        public void Name_And_Summary_AreReadOnly()
        {
            var filter = new DummyGptmdFilter();
            var name = "ReadonlyName";
            var summary = "ReadonlySummary";
            var viewModel = new GptmdFilterViewModel(filter, name, summary);

            Assert.That(viewModel.Name, Is.EqualTo(name));
            Assert.That(viewModel.Summary, Is.EqualTo(summary));
        }

        [Test]
        public void Filter_Property_IsReadOnly()
        {
            var filter = new DummyGptmdFilter();
            var viewModel = new GptmdFilterViewModel(filter, "n", "s");

            Assert.That(viewModel.Filter, Is.EqualTo(filter));
        }
    }
}
