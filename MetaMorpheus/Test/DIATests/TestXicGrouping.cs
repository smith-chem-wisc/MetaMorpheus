using EngineLayer.DIA;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using EngineLayer;
using TaskLayer;
using Chemistry;
using FlashLFQ;
using Omics;
using Proteomics.ProteolyticDigestion;
using Omics.Modifications;
using Omics.Fragmentation;

namespace Test.DIATests
{
    public class TestXicGrouping
    { 
        [Test]
        public static void TestCorrelationCalculation()
        {
            //perfectly aligned XICs
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
            double corr = PrecursorFragmentsGroup.CalculateXicCorrelationXYData(xic1, xic2);
            Assert.That(corr, Is.EqualTo(1.0).Within(1e-6));

            //XICs with completely opposite trends
            var peakList3 = new List<IIndexedPeak>();
            double[] intensityMultipliers2 = { 3, 2, 1, 2, 3 };
            for (int i = 0; i < intensityMultipliers2.Length; i++)
            {
                peakList3.Add(new IndexedMassSpectralPeak(intensity: 1e6 * intensityMultipliers2[i], retentionTime: 1 + i / 10, zeroBasedScanIndex: i, mz: 501.0));
            }
            var xic3 = new ExtractedIonChromatogram(peakList3);
            corr = PrecursorFragmentsGroup.CalculateXicCorrelationXYData(xic1, xic3);
            Assert.That(corr, Is.EqualTo(-1.0).Within(1e-6));

            //XICs with insufficient overlap points
            var peakList4 = new List<IIndexedPeak>();
            for (int i = 0; i < intensityMultipliers.Length; i++)
            {
                peakList4.Add(new IndexedMassSpectralPeak(intensity: 1e6 * intensityMultipliers[i], retentionTime: 1 + i / 10, zeroBasedScanIndex: i + 4, mz: 502.0));
            }
            var xic4 = new ExtractedIonChromatogram(peakList4);
            corr = PrecursorFragmentsGroup.CalculateXicCorrelationXYData(xic1, xic4);
            Assert.That(corr, Is.EqualTo(double.NaN).Within(1e-6));

            //XICs with spline
            var cubicSpline = new XicCubicSpline();
            cubicSpline.SetXicSplineXYData(xic1, true);
            cubicSpline.SetXicSplineXYData(xic2, true);
            corr = PrecursorFragmentsGroup.CalculateXicCorrelationXYData(xic1, xic2);
            Assert.That(corr, Is.EqualTo(1.0).Within(1e-6));
            var linearSpline = new XicLinearSpline();
            linearSpline.SetXicSplineXYData(xic1, true);
            linearSpline.SetXicSplineXYData(xic2, true);
            corr = PrecursorFragmentsGroup.CalculateXicCorrelationXYData(xic1, xic2);
            Assert.That(corr, Is.EqualTo(1.0).Within(1e-6));
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
        public static void TestXicPfGrouping()
        {
            //Create two fake precursor peaks with no overlap in their retention time
            var ms1Peaks = new List<IIndexedPeak>();
            double[] intensityMultipliers = { 1, 2, 3, 2, 1 };
            for (int i = 0; i < intensityMultipliers.Length; i++)
            {
                ms1Peaks.Add(new IndexedMassSpectralPeak(intensity: 1e5 * intensityMultipliers[i], retentionTime: 1 + i / 10, zeroBasedScanIndex: i, mz: 500.0));
            }
            var ms1Xic = new ExtractedIonChromatogram(ms1Peaks);

            //Create a list of fragment XICs that align with the positive precursor XIC but have no overlap with the negative precursor XIC
            var matchedFragXics = new List<ExtractedIonChromatogram>();
            int numberOfmatchedFragmentXics = 5;
            for (int j = 0; j < numberOfmatchedFragmentXics; j++)
            {
                var fragPeaks = new List<IIndexedPeak>();
                for (int i = 0; i < intensityMultipliers.Length; i++)
                    fragPeaks.Add(new IndexedMassSpectralPeak(intensity: 1e5 * intensityMultipliers[i], retentionTime: 1.1 + i / 10, zeroBasedScanIndex: i, mz: 500.0));
                var fragXic = new ExtractedIonChromatogram(fragPeaks);
                matchedFragXics.Add(fragXic);
            }

            //Create a list of fragment Xics that fail to pass apexRtTolerance
            var apexXics = new List<ExtractedIonChromatogram>();
            for (int j = 0; j < 4; j++)
            {
                var fragPeaks = new List<IIndexedPeak>();
                for (int i = 0; i < intensityMultipliers.Length; i++)
                    fragPeaks.Add(new IndexedMassSpectralPeak(intensity: 1e5 * intensityMultipliers[i], retentionTime: 1.5 + i / 10, zeroBasedScanIndex: i + 50, mz: 500.0));
                var fragXic = new ExtractedIonChromatogram(fragPeaks);
                apexXics.Add(fragXic);
            }

            //Create a list of fragment Xics that passes apexRtTolerance and overlap threshold but not correlation threshold
            var corrXics = new List<ExtractedIonChromatogram>();
            double[] intensityMultipliers2 = { 1, 1, 1, 1.5, 1 };
            for (int j = 0; j < 3; j++)
            {
                var fragPeaks = new List<IIndexedPeak>();
                for (int i = 0; i < intensityMultipliers.Length; i++)
                    fragPeaks.Add(new IndexedMassSpectralPeak(intensity: 1e5 * intensityMultipliers2[i], retentionTime: 1.1 + i / 10, zeroBasedScanIndex: i, mz: 500.0));
                var fragXic = new ExtractedIonChromatogram(fragPeaks);
                corrXics.Add(fragXic);
            }

            //perform XicGrouping on the ms1 XIC and all fragment XICs with apexTolerance 0.2f, overlap threshold 0.5, correlation threshold 0.5, and minimum number of fragment XICs 2
            var allFragXics = matchedFragXics.Concat(apexXics).Concat(corrXics).ToList();
            var pfGroup = XicGroupingEngine.GroupFragmentsForOnePrecursor(ms1Xic, allFragXics, 0.2f, 0.5, 0.5, 2);
            //all and only the Xics in matchedFragXics should be grouped with the ms1 XIC, and there are 5 of them
            Assert.That(pfGroup, Is.Not.Null);
            Assert.That(pfGroup.PFpairs.Count, Is.EqualTo(5));
            //check that all fragment XICs grouped are from matchedFragXics, and they all should have perfect correlation and overlap with the precursor XIC
            foreach (var pfPair in pfGroup.PFpairs)
            {
                Assert.That(matchedFragXics.Contains(pfPair.FragmentXic));
                Assert.That(pfPair.Correlation, Is.EqualTo(1.0).Within(1e-6));
                Assert.That(pfPair.Overlap, Is.EqualTo(1.0).Within(1e-6));
            }

            //perform XicGrouping with apexTolerance 0.2f, overlap threshold 0.5, correlation threshold 0.5, and minimum number of fragment XICs 10; now there should be no groups returned because only 5 fragment XICs are in the group
            pfGroup = XicGroupingEngine.GroupFragmentsForOnePrecursor(ms1Xic, allFragXics, 0.2f, 0.5, 0.5, 10);
            Assert.That(pfGroup, Is.Null);

            //perform XicGrouping with apexTolerance 0.2f, overlap threshold 0.5, correlation threshold 0, and minimum number of fragment XICs 2; now the tolerance on correlation is looser and all fragment XICs in corrXics should be grouped
            pfGroup = XicGroupingEngine.GroupFragmentsForOnePrecursor(ms1Xic, allFragXics, 0.2f, 0.5, 0, 2);
            Assert.That(pfGroup, Is.Not.Null);
            Assert.That(pfGroup.PFpairs.Count, Is.EqualTo(5 + 3)); //5 from matchedFragXics and 3 from corrXics
            //check that all fragment XICs grouped are from matchedFragXics and corrXics
            foreach (var pfPair in pfGroup.PFpairs)
            {
                Assert.That(matchedFragXics.Contains(pfPair.FragmentXic) || corrXics.Contains(pfPair.FragmentXic));
            }
        }

        [Test]
        public static void TestXicGroupingWithRankFilter()
        {
            // Create two fake precursor Xics with slightly different elution profiles
            var intensityMultipliers1 = new[] { 1, 2, 10, 2, 1 };
            var intensityMultipliers2 = new[] { 1, 2, 5, 2, 1 };
            var intensityMultipliers3 = new[] { 1, 2, 2, 2, 1 }; 
            var preList1 = new List<IIndexedPeak>();
            var preList2 = new List<IIndexedPeak>();
            var preList3 = new List<IIndexedPeak>();
            for (int i = 0; i < intensityMultipliers1.Length; i++)
            {
                preList1.Add(new IndexedMassSpectralPeak(500.0 + i, intensityMultipliers1[i] * 1e6, i, 1 + i / 10.0));
                preList2.Add(new IndexedMassSpectralPeak(501.0 + i, intensityMultipliers2[i] * 1e6, i, 1 + i / 10.0));
                preList3.Add(new IndexedMassSpectralPeak(502.0 + i, intensityMultipliers3[i] * 1e6, i, 1 + i / 10.0));
            }
            var precursorXic1 = new ExtractedIonChromatogram(preList1);
            var precursorXic2 = new ExtractedIonChromatogram(preList2);
            var precursorXic3 = new ExtractedIonChromatogram(preList3);

            // Create two sets of fragment XICs that closely align with precursor 1 and 2 respectively
            //fragment XICs in set 1 correlate better with precursorXic1, and those in set 2 correlate better with precursorXic2
            int numberOfFragments = 5;
            var fragmentXicSet1 = new List<ExtractedIonChromatogram>();
            var fragmentXicSet2 = new List<ExtractedIonChromatogram>();
            for (int j = 0; j < numberOfFragments; j++)
            {
                var fragList1 = new List<IIndexedPeak>();
                var fragList2 = new List<IIndexedPeak>();
                for (int i = 0; i < intensityMultipliers1.Length; i++)
                {
                    //- 0.05 * j is just to introduce a little variation so the correlation values are not exactly the same but still close
                    fragList1.Add(new IndexedMassSpectralPeak(300.0 + j, (intensityMultipliers1[i] - 0.05 * j) * 1e5, i, 1 + i / 10.0));
                    fragList2.Add(new IndexedMassSpectralPeak(400.0 + j, (intensityMultipliers2[i] - 0.05 * j) * 1e5, i, 1 + i / 10.0));
                }
                var fragXic1 = new ExtractedIonChromatogram(fragList1);
                var fragXic2 = new ExtractedIonChromatogram(fragList2);
                fragmentXicSet1.Add(fragXic1);
                fragmentXicSet2.Add(fragXic2);
            }
            var allPrecursorXics = new List<ExtractedIonChromatogram> { precursorXic1, precursorXic2, precursorXic3};
            var allFragmentXics = fragmentXicSet1.Concat(fragmentXicSet2).ToList();

            //Create XicGroupingEngine and do precursor-fragment grouping for all fake Xics with a relatively loose threshold for correlation
            //When there is no rank filter, all fragment XICs should be grouped with each precursor XIC
            var xicGroupingEngine = new XicGroupingEngine(0.5f, 0.5, 0.25, 2, 1, null, null);
            var groups = xicGroupingEngine.PrecursorFragmentGrouping(allPrecursorXics, allFragmentXics);
            //So there will be 3 pfGroups, each with 10 grouped fragments
            Assert.That(groups.Count, Is.EqualTo(3));
            Assert.That(groups.All(g => g.PFpairs.Count == 10));

            //When we set a threshold of 4 for fragment rank but no precursor rank threshold
            xicGroupingEngine = new XicGroupingEngine(0.5f, 0.5, 0.25, 2, 1, null, 4);
            groups = xicGroupingEngine.PrecursorFragmentGrouping(allPrecursorXics, allFragmentXics);
            Assert.That(groups.Count, Is.EqualTo(3));
            //The 4 fragment Xics that correlate better with each precursor Xic should be retained
            Assert.That(groups.All(g => g.PFpairs.Count == 4));
            //all fragment Xics paired with precursor 1 should from fragment set 1, and all pfPairs paired with precursor 2 should from fragment set 2
            Assert.That(groups[0].PFpairs.All(pf => fragmentXicSet1.Contains(pf.FragmentXic)));
            Assert.That(groups[1].PFpairs.All(pf => fragmentXicSet2.Contains(pf.FragmentXic)));
            //Since fragment xics in set 2 are closer to precursor3, the pfPairs paired with precursor 3 should be from fragment set 2
            Assert.That(groups[2].PFpairs.All(pf => fragmentXicSet2.Contains(pf.FragmentXic)));

            //When we set a threshold of 1 for precursor rank but no fragment rank threshold
            xicGroupingEngine = new XicGroupingEngine(0.5f, 0.5, 0.25, 2, 1, 1, null);
            groups = xicGroupingEngine.PrecursorFragmentGrouping(allPrecursorXics, allFragmentXics);
            //Precursor 3 should not have any fragments grouped with it since each fragment Xic is limited to be paired with only one precursor
            Assert.That(groups.Count, Is.EqualTo(2));
            //The precursor Xic that correlate best with each fragment Xic should be retained
            //so the 5 fragment Xics in set 1 should go with precursor 1 and the 5 in set 2 should go with precursor 2
            Assert.That(groups.All(g => g.PFpairs.Count == 5));
            Assert.That(groups[0].PFpairs.All(pf => fragmentXicSet1.Contains(pf.FragmentXic)));
            Assert.That(groups[1].PFpairs.All(pf => fragmentXicSet2.Contains(pf.FragmentXic)));

            //When we set a threshold of 2 for precursor rank and 6 for fragment rank
            xicGroupingEngine = new XicGroupingEngine(0.5f, 0.5, 0.25, 2, 1, 2, 6);
            groups = xicGroupingEngine.PrecursorFragmentGrouping(allPrecursorXics, allFragmentXics);
            //Precursor 3 should still not have any fragments so there are two groups
            Assert.That(groups.Count, Is.EqualTo(2));
            //Each group should have 6 fragments
            Assert.That(groups.All(g => g.PFpairs.Count == 6));
            //precursor 1 should have all fragments from fragment set 1 plus one from fragment set 2
            Assert.That(groups[0].PFpairs.Count(pf => fragmentXicSet1.Contains(pf.FragmentXic)), Is.EqualTo(5));
        }
    }
}
