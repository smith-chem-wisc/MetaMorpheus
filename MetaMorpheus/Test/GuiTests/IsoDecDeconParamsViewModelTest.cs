using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using GuiFunctions;
using MassSpectrometry;

namespace Test.GuiTests
{
    [TestFixture]
    internal class IsoDecDeconParamsViewModelTest
    {
        private IsoDecDeconParamsViewModel _viewModel;

        [SetUp]
        public void Setup()
        {
            _viewModel = new IsoDecDeconParamsViewModel(new IsoDecDeconvolutionParameters());
        }

        [Test]
        public void TestPhaseRes()
        {
            _viewModel.PhaseRes = 4;
            Assert.That(_viewModel.PhaseRes, Is.EqualTo(4));

            _viewModel.PhaseRes = 8;
            Assert.That(_viewModel.PhaseRes, Is.EqualTo(8));

            _viewModel.PhaseRes = 5;
            Assert.That(_viewModel.PhaseRes, Is.EqualTo(8)); // Should not change
        }

        [Test]
        public void TestCssThreshold()
        {
            _viewModel.CssThreshold = 0.5f;
            Assert.That(_viewModel.CssThreshold, Is.EqualTo(0.5f));

            _viewModel.CssThreshold = 1.1f;
            Assert.That(_viewModel.CssThreshold, Is.EqualTo(0.5f)); // Should not change

            _viewModel.CssThreshold = -0.1f;
            Assert.That(_viewModel.CssThreshold, Is.EqualTo(0.5f)); // Should not change
        }

        [Test]
        public void TestMatchTolerance()
        {
            _viewModel.MatchTolerance = 10f;
            Assert.That(_viewModel.MatchTolerance, Is.EqualTo(10f));

            _viewModel.MatchTolerance = -1f;
            Assert.That(_viewModel.MatchTolerance, Is.EqualTo(10f)); // Should not change
        }

        [Test]
        public void TestMaximumShift()
        {
            _viewModel.MaximumShift = 5;
            Assert.That(_viewModel.MaximumShift, Is.EqualTo(5));

            _viewModel.MaximumShift = -1;
            Assert.That(_viewModel.MaximumShift, Is.EqualTo(5)); // Should not change
        }

        [Test]
        public void TestMzWindowForIsotopeDistributionMinimum()
        {
            _viewModel.MzWindowForIsotopeDistributionMinimum = -0.5f;
            Assert.That(_viewModel.MzWindowForIsotopeDistributionMinimum, Is.EqualTo(-0.5f));

            _viewModel.MzWindowForIsotopeDistributionMinimum = 0.5f;
            Assert.That(_viewModel.MzWindowForIsotopeDistributionMinimum, Is.EqualTo(-0.5f)); // Should not change
        }

        [Test]
        public void TestMzWindowForIsotopeDistributionMaximum()
        {
            _viewModel.MzWindowForIsotopeDistributionMaximum = 0.5f;
            Assert.That(_viewModel.MzWindowForIsotopeDistributionMaximum, Is.EqualTo(0.5f));

            _viewModel.MzWindowForIsotopeDistributionMaximum = -0.5f;
            Assert.That(_viewModel.MzWindowForIsotopeDistributionMaximum, Is.EqualTo(0.5f)); // Should not change
        }

        [Test]
        public void TestKnockdownRounds()
        {
            _viewModel.KnockdownRounds = 3;
            Assert.That(_viewModel.KnockdownRounds, Is.EqualTo(3));

            _viewModel.KnockdownRounds = -1;
            Assert.That(_viewModel.KnockdownRounds, Is.EqualTo(3)); // Should not change
        }

        [Test]
        public void TestMinAreaCovered()
        {
            _viewModel.MinAreaCovered = 0.5f;
            Assert.That(_viewModel.MinAreaCovered, Is.EqualTo(0.5f));

            _viewModel.MinAreaCovered = 1.1f;
            Assert.That(_viewModel.MinAreaCovered, Is.EqualTo(0.5f)); // Should not change

            _viewModel.MinAreaCovered = -0.1f;
            Assert.That(_viewModel.MinAreaCovered, Is.EqualTo(0.5f)); // Should not change
        }

        [Test]
        public void TestDataThreshold()
        {
            _viewModel.DataThreshold = 0.5f;
            Assert.That(_viewModel.DataThreshold, Is.EqualTo(0.5f));

            _viewModel.DataThreshold = 1.1f;
            Assert.That(_viewModel.DataThreshold, Is.EqualTo(0.5f)); // Should not change

            _viewModel.DataThreshold = -0.1f;
            Assert.That(_viewModel.DataThreshold, Is.EqualTo(0.5f)); // Should not change
        }

        [Test]
        public void TestReportMultipleMonoisotopicMasses()
        {
            _viewModel.ReportMultipleMonoisotopicMasses = true;
            Assert.That(_viewModel.ReportMultipleMonoisotopicMasses, Is.True);

            _viewModel.ReportMultipleMonoisotopicMasses = false;
            Assert.That(_viewModel.ReportMultipleMonoisotopicMasses, Is.False);
        }
    }
}
