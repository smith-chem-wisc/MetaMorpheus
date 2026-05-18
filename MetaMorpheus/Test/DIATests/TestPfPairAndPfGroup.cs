using EngineLayer;
using EngineLayer.DIA;
using MassSpectrometry;
using MathNet.Numerics.Statistics;
using MzLibUtil;
using NUnit.Framework;
using NUnit.Framework.Interfaces;
using System;
using System.Collections.Generic;
using System.Data.Entity.Core;
using System.Linq;
using TaskLayer;
namespace Test.DIATests
{
    public class TestPfPairAndPfGroup
    {
        [Test]
        public static void TestCorrelationCalculation()
        {
            //Create two perfectly aligned XICs
            var peakList1 = new List<IIndexedPeak>();
            double[] intensityMultipliers = { 1, 2, 3, 2, 1 };
            for (int i = 0; i < intensityMultipliers.Length; i++)
            {
                peakList1.Add(new IndexedMassSpectralPeak(intensity: 1e5 * intensityMultipliers[i], retentionTime: 1 + i / 10, zeroBasedScanIndex: i, mz: 500.0));
            }
            var xic1 = new ExtractedIonChromatogram(peakList1);
            var peakList2 = new List<IIndexedPeak>();
            for (int i = 0; i < intensityMultipliers.Length; i++)
            {
                peakList2.Add(new IndexedMassSpectralPeak(intensity: 1e6 * intensityMultipliers[i], retentionTime: 1 + i / 10, zeroBasedScanIndex: i, mz: 501.0));
            }
            var xic2 = new ExtractedIonChromatogram(peakList2);
            double corr = PrecursorFragmentsGroup.CalculateXicCorrelation(xic1, xic2);
            //The Pearson's correlation should be 1
            Assert.That(corr, Is.EqualTo(1.0).Within(1e-6));

            //XICs with completely opposite trends
            var peakList3 = new List<IIndexedPeak>();
            double[] intensityMultipliers2 = { 3, 2, 1, 2, 3 };
            for (int i = 0; i < intensityMultipliers2.Length; i++)
            {
                peakList3.Add(new IndexedMassSpectralPeak(intensity: 1e6 * intensityMultipliers2[i], retentionTime: 1 + i / 10, zeroBasedScanIndex: i, mz: 501.0));
            }
            var xic3 = new ExtractedIonChromatogram(peakList3);
            corr = PrecursorFragmentsGroup.CalculateXicCorrelation(xic1, xic3);
            //The Pearson's correlation should be -1
            Assert.That(corr, Is.EqualTo(-1.0).Within(1e-6));

            //XICs with insufficient overlap points
            var peakList4 = new List<IIndexedPeak>();
            for (int i = 0; i < intensityMultipliers.Length; i++)
            {
                peakList4.Add(new IndexedMassSpectralPeak(intensity: 1e6 * intensityMultipliers[i], retentionTime: 1 + i / 10, zeroBasedScanIndex: i + 4, mz: 502.0));
            }
            var xic4 = new ExtractedIonChromatogram(peakList4);
            corr = PrecursorFragmentsGroup.CalculateXicCorrelation(xic1, xic4);
            Assert.That(corr, Is.NaN);

            //Test with missing points in XICs
            var peakList5 = new List<IIndexedPeak>();
            for (int i = 0; i < intensityMultipliers.Length; i++)
            {
                if (i != 1) //missing the peak at scan index 1
                    peakList5.Add(new IndexedMassSpectralPeak(intensity: 1e6 * intensityMultipliers[i], retentionTime: 1 + i / 10, zeroBasedScanIndex: i, mz: 503.0));
            }
            var xic5 = new ExtractedIonChromatogram(peakList5);
            corr = PrecursorFragmentsGroup.CalculateXicCorrelation(xic1, xic5);
            Assert.That(corr, Is.LessThan(1.0));
            Assert.That(corr, Is.GreaterThan(0.0)); //Correlation should still be positive but less than 1 due to missing point
        }

        [Test]
        public static void TestCorrelationCalculationWithSpline()
        {
            //When spline is available, the correlation is calculated with the spline data
            //xic1 and xic2 are perfectly aligned; xic3 does not have enough overlap points so corr will be NaN
            var peakList1 = new List<IIndexedPeak>();
            var peakList2 = new List<IIndexedPeak>();
            var peakList3 = new List<IIndexedPeak>();
            double[] intensityMultipliers = { 1, 2, 3, 2, 1 };
            for (int i = 0; i < intensityMultipliers.Length; i++)
            {
                peakList1.Add(new IndexedMassSpectralPeak(intensity: 1e5 * intensityMultipliers[i], retentionTime: 1 + i / 10, zeroBasedScanIndex: i, mz: 500.0));
                peakList2.Add(new IndexedMassSpectralPeak(intensity: 1e6 * intensityMultipliers[i], retentionTime: 1 + i / 10, zeroBasedScanIndex: i, mz: 501.0));
                peakList3.Add(new IndexedMassSpectralPeak(intensity: 1e6 * intensityMultipliers[i], retentionTime: 1.5 + i / 10, zeroBasedScanIndex: i + 5, mz: 501.0));
            }
            var xic1 = new ExtractedIonChromatogram(peakList1);
            var xic2 = new ExtractedIonChromatogram(peakList2);
            var xic3 = new ExtractedIonChromatogram(peakList3);

            //It should still return 1.0 for two perfectly aligned XICs (xic1 and xic2)
            var cubicSpline = new XicCubicSpline(scanIndexBased: true);
            cubicSpline.SetXicSplineXYData(xic1);
            cubicSpline.SetXicSplineXYData(xic2);
            var corr = PrecursorFragmentsGroup.CalculateXicCorrelation(xic1, xic2);
            Assert.That(corr, Is.EqualTo(1.0).Within(1e-6));
            var linearSpline = new XicLinearSpline(scanIndexBased: true);
            linearSpline.SetXicSplineXYData(xic1);
            linearSpline.SetXicSplineXYData(xic2);
            corr = PrecursorFragmentsGroup.CalculateXicCorrelation(xic1, xic2);
            Assert.That(corr, Is.EqualTo(1.0).Within(1e-6));

            //For the XICs with insufficient overlap points (xic1 and xic3), it should return NaN
            linearSpline.SetXicSplineXYData(xic3);
            corr = PrecursorFragmentsGroup.CalculateXicCorrelation(xic1, xic3);
            Assert.That(corr, Is.NaN);
            //Test on umpire correlation calculation, should also return NaN
            var bSpline = new Bspline(2, 2);
            bSpline.SetXicSplineXYData(xic1);
            bSpline.SetXicSplineXYData(xic3);
            corr = PrecursorFragmentsGroup.CalculateXicCorrXYData_Umpire(xic1, xic3, 2);
            Assert.That(corr, Is.NaN);
        }

        [Test]
        public static void TestCalculateOverlapRatio()
        {
            //perfectly aligned XICs, should return 1
            var peakList1 = new List<IIndexedPeak>();
            double[] intensityMultipliers = { 1, 2, 3, 2, 1 };
            for (int i = 0; i < intensityMultipliers.Length; i++)
            {
                peakList1.Add(new IndexedMassSpectralPeak(intensity: 1e5 * intensityMultipliers[i], retentionTime: 1 + i / 10, zeroBasedScanIndex: i, mz: 500.0));
            }
            var xic1 = new ExtractedIonChromatogram(peakList1);
            var peakList2 = new List<IIndexedPeak>();
            for (int i = 0; i < intensityMultipliers.Length; i++)
            {
                peakList2.Add(new IndexedMassSpectralPeak(intensity: 1e6 * intensityMultipliers[i], retentionTime: 1 + i / 10, zeroBasedScanIndex: i, mz: 501.0));
            }
            var xic2 = new ExtractedIonChromatogram(peakList2);
            double overlap = PrecursorFragmentsGroup.CalculateXicOverlapRatio(xic1, xic2);
            Assert.That(overlap, Is.EqualTo(1.0).Within(1e-6));

            //XICs with no overlap, should return 0
            var peakList3 = new List<IIndexedPeak>();
            for (int i = 0; i < intensityMultipliers.Length; i++)
            {
                peakList3.Add(new IndexedMassSpectralPeak(intensity: 1e6 * intensityMultipliers[i], retentionTime: 100 + i / 10, zeroBasedScanIndex: i + 100, mz: 501.0));
            }
            var xic3 = new ExtractedIonChromatogram(peakList3);
            overlap = PrecursorFragmentsGroup.CalculateXicOverlapRatio(xic1, xic3);
            Assert.That(overlap, Is.EqualTo(0));

            //XICs with partial overlap
            //xic1 covers scan 0 to 4, xic4 covers scan 1 to 5, so the overlap should be 0.6 (3 out of 5 scans)
            var peakList4 = new List<IIndexedPeak>();
            for (int i = 0; i < intensityMultipliers.Length; i++)
            {
                peakList4.Add(new IndexedMassSpectralPeak(intensity: 1e6 * intensityMultipliers[i], retentionTime: 1 + i / 10, zeroBasedScanIndex: i + 1, mz: 502.0));
            }
            var xic4 = new ExtractedIonChromatogram(peakList4);
            overlap = PrecursorFragmentsGroup.CalculateXicOverlapRatio(xic1, xic4);
            Assert.That(overlap, Is.EqualTo(0.6).Within(1e-6));
        }

        [Test]
        public static void TestFragmentAndPrecursorRank()
        {
            //This test creates a set of precursor and fragment XICs with random intensities and group them.
            //The resulting precursor-fragment pairs will have random correlations depending on their generated intensities. However, the rankings should be done as expected regardless of the correlation values.
            var random = new Random(42);

            // Create a set of precursor XICs with random intensities
            int numberOfPrecursors = 10;
            var precursorXics = new List<ExtractedIonChromatogram>();
            for (int i = 0; i < numberOfPrecursors; i++)
            {
                var peakList = new List<IIndexedPeak>();
                for (int j = 0; j < 10; j++)
                {
                    peakList.Add(new IndexedMassSpectralPeak(500.0 + i, random.NextDouble() * 1e5,  j, 1 + j / 10.0));
                }
                precursorXics.Add(new ExtractedIonChromatogram(peakList));
            }

            // Create a set of fragment XICs with random intensities
            int numberOfFragments = 20;
            var fragmentXics = new List<ExtractedIonChromatogram>();
            for (int i = 0; i < numberOfFragments; i++)
            {
                var peakList = new List<IIndexedPeak>();
                for (int j = 0; j < 10; j++)
                {
                    peakList.Add(new IndexedMassSpectralPeak(300.0 + i, random.NextDouble() * 1e5, j, 1.01 + j / 10.0));
                }
                fragmentXics.Add(new ExtractedIonChromatogram(peakList));
            }

            //Try to group fragment Xics with each precursor Xic. A set of reasonable threshold values are used for grouping.
            var pfGroups = new List<PrecursorFragmentsGroup>();
            foreach (var precursorXic in precursorXics)
            {
                var pfGroup = XicGroupingEngine.GroupFragmentsForOnePrecursor(precursorXic, fragmentXics, 0.3f, 0.3, 0, 1);
                if (pfGroup != null)
                {
                    pfGroups.Add(pfGroup);
                    //In this precursor-fragments grouping method, all pfPairs should have a valid Correlation value higher than the threshold
                    Assert.That(pfGroup.PFpairs.All(pf => pf.Correlation.HasValue));
                }
            }

            //Test fragment ranking. Fragment ranking is performed within each pfGroup.
            foreach (var pfGroup in pfGroups)
            {
                pfGroup.SetFragmentRankForPfPairs();
                //All PFpairs should be ordered with decreasing correlation values.
                for (int i = 0; i < pfGroup.PFpairs.Count - 1; i++)
                {
                    Assert.That(pfGroup.PFpairs[i].Correlation, Is.GreaterThanOrEqualTo(pfGroup.PFpairs[i + 1].Correlation));
                    Assert.That(pfGroup.PFpairs[i].FragmentRank, Is.EqualTo(i + 1));
                }
            }

            //Test precursor ranking. Precursor ranking needs to be performed across all groups, based on each fragment Xic. 
            PrecursorFragmentPair.SetPrecursorRankForPfPairs(pfGroups.SelectMany(g => g.PFpairs));
            //pfPairs are grouped by fragment Xics to rank precursor Xics, so precursor ranks are tested by each fragment Xic
            foreach(var fragment in fragmentXics)
            {
                //Get all pfPairs associated with this fragment Xic. Skip if no pfPairs are found.
                var fragmentPairs = pfGroups.SelectMany(g => g.PFpairs.Where(pf => pf.FragmentXic == fragment)).ToList();
                if (fragmentPairs.IsNullOrEmpty()) continue;
                //precursor ranking is performed across all pfGroups so any pfPair should have a valid PrecursorRank value
                Assert.That(fragmentPairs.All(pf => pf.PrecursorRank.HasValue));
                //The largest PrecursorRank of these pfPairs should equal to the total number of pfPairs associated with this fragment Xic
                Assert.That(fragmentPairs.Max(pf => pf.PrecursorRank), Is.EqualTo(fragmentPairs.Count()));
                //If we rank these pfPairs based on their PrecursorRank, they should be ordered in decreasing Correlation values.
                fragmentPairs.Sort((a, b) => a.PrecursorRank.Value.CompareTo(b.PrecursorRank.Value));
                for (int i = 0; i < fragmentPairs.Count - 1; i++)
                {
                    Assert.That(fragmentPairs[i].Correlation, Is.GreaterThanOrEqualTo(fragmentPairs[i + 1].Correlation));
                }
            }
        }

        [Test]
        public static void TestPfPairRankingExceptionHandling()
        {
            var peakList = new List<IIndexedPeak> { new IndexedMassSpectralPeak(500.0 , 1, 1, 1) };
            var emptyXic = new ExtractedIonChromatogram(peakList);
            var pfGroup = new PrecursorFragmentsGroup(emptyXic, null);
            var ex = Assert.Throws<MetaMorpheusException>(() => pfGroup.SetFragmentRankForPfPairs());
            Assert.That(ex.Message, Is.EqualTo("The PfGroup does not contain any PFpairs."));
        }
    }
}
