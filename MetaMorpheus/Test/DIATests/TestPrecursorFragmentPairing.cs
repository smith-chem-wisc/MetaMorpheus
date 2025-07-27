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

namespace Test
{
    public class TestPrecursorFragmentPairing
    {
        public static MsDataScan[] FakeMs1Scans { get; set; }
        public static IsotopicDistribution Dist { get; set; }
        public static double Intensity { get; set; } 

        [OneTimeSetUp]
        public void OneTimeSetup()
        {
            string peptide = "PEPTIDE";
            Intensity = 1e6;
            ChemicalFormula cf = new Proteomics.AminoAcidPolymer.Peptide(peptide).GetChemicalFormula();
            Dist = IsotopicDistribution.GetDistribution(cf, 0.125, 1e-8);

            //Create MS1 scans with two precursors that do not have any overlap in their retention time; one of the precursors has an apex RT at the third scan, the other at the eighth scan
            FakeMs1Scans = new MsDataScan[10];
            double[] intensityMultipliers = { 1, 2, 3, 2, 1, 1, 2, 3, 2, 1 };
            var massArray1 = Dist.Masses.SelectMany(v => new List<double> { v.ToMz(1) }).ToArray();
            var massArray2 = Dist.Masses.SelectMany(v => new List<double> { (v + 10).ToMz(1) }).ToArray();
            double[][] masses = Enumerable.Repeat(massArray1, 5).Concat(Enumerable.Repeat(massArray2, 5)).ToArray();

            for (int s = 0; s < FakeMs1Scans.Length; s++)
            {
                double[] intensities = Dist.Intensities.SelectMany(v => new List<double> { v * Intensity * intensityMultipliers[s] }).ToArray();
                FakeMs1Scans[s] = new MsDataScan(massSpectrum: new MzSpectrum(masses[s], intensities, false), oneBasedScanNumber: s + 1, msnOrder: 1, isCentroid: true,
                    polarity: Polarity.Positive, retentionTime: 1.0 + s / 10.0, scanWindowRange: new MzRange(400, 1600), scanFilter: "f",
                    mzAnalyzer: MZAnalyzerType.Orbitrap, totalIonCurrent: intensities.Sum(), injectionTime: 1.0, noiseData: null, nativeId: "scan=" + (s + 1));
            }
        }

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
        public static void PrecursorFragmentsGroupTest()
        {
            //Create MS2 scans with two sets of fragments each aligning perfectly with one of the precursors in FakeMs1Scans
            var fakeMs2Scans = new MsDataScan[10];
            double[] intensityMultipliers2 = { 0.1, 0.2, 0.3, 0.2, 0.1, 0.1, 0.2, 0.3, 0.2, 0.1 };
            var fragArray1 = Dist.Masses.SelectMany(v => new List<double> { v/2.0.ToMz(1) }).ToArray();
            var fragArray2 = Dist.Masses.SelectMany(v => new List<double> { (v + 10)/2.0.ToMz(1) }).ToArray();
            double[][] frags = Enumerable.Repeat(fragArray1, 5).Concat(Enumerable.Repeat(fragArray2, 5)).ToArray();
            for (int s = 0; s < fakeMs2Scans.Length; s++)
            {
                double[] intensities = Dist.Intensities.SelectMany(v => new List<double> { v * Intensity * intensityMultipliers2[s]}).ToArray();
                fakeMs2Scans[s] = new MsDataScan(massSpectrum: new MzSpectrum(frags[s], intensities, false), oneBasedScanNumber: s + 1, msnOrder: 2, isCentroid: true,
                    polarity: Polarity.Positive, retentionTime: 1.01 + s / 10.0, scanWindowRange: new MzRange(400, 1600), scanFilter: "f",
                    mzAnalyzer: MZAnalyzerType.Orbitrap, totalIonCurrent: intensities.Sum(), injectionTime: 1.0, noiseData: null, nativeId: "scan=" + (s + 1));
            }

            //Find XICs of the ms1 and ms2 scans
            var deconParameters = new ClassicDeconvolutionParameters(1, 20, 4, 3);
            var ms1XicConstructor = new NeutralMassXicConstructor(new PpmTolerance(20), 2, 1, 3, deconParameters);
            var ms2XicConstructor = new MzPeakXicConstructor(new PpmTolerance(10), 2, 1, 3);
            var ms1Xics = ms1XicConstructor.GetAllXics(FakeMs1Scans);
            var ms2Xics = ms2XicConstructor.GetAllXics(fakeMs2Scans);
            var xicGroupingEngine = new XicGrouping(0.1f, 0.5, 0.5);
            var allGroups = xicGroupingEngine.PrecursorFragmentGrouping(ms1Xics, ms2Xics);
            //There should be two pfGroups, one for each precursor in FakeMs1Scans
            Assert.That(allGroups.Count, Is.EqualTo(2));
            foreach (var pfGroup in allGroups)
            {
                //For each precursor XIC, half of the fragment XICs should be paired with it (which is 10 in this case)
                Assert.That(pfGroup.PFpairs.Count, Is.EqualTo(10));
                //All the fragment XICs should have the same apex scan index as the precursor XIC
                Assert.That(pfGroup.PFpairs.All(pf => pf.FragmentXic.ApexScanIndex == pf.PrecursorXic.ApexScanIndex));
            }

            //If we set the fragments paired with the second precursor to have an XIC with an opposite trend, they would not be paired with the second precursor and the pfGrouping would return null
            var fakeMs2Scans2 = new MsDataScan[10];
            double[] intensityMultipliers3 = { 0.1, 0.2, 0.3, 0.2, 0.1, 0.3, 0.2, 0.1, 0.2, 0.3 };
            for (int s = 0; s < fakeMs2Scans2.Length; s++)
            {
                double[] intensities = Dist.Intensities.SelectMany(v => new List<double> { v * Intensity * intensityMultipliers3[s]}).ToArray();
                fakeMs2Scans2[s] = new MsDataScan(massSpectrum: new MzSpectrum(frags[s], intensities, false), oneBasedScanNumber: s + 1, msnOrder: 2, isCentroid: true,
                    polarity: Polarity.Positive, retentionTime: 1.01 + s / 10.0, scanWindowRange: new MzRange(400, 1600), scanFilter: "f",
                    mzAnalyzer: MZAnalyzerType.Orbitrap, totalIonCurrent: intensities.Sum(), injectionTime: 1.0, noiseData: null, nativeId: "scan=" + (s + 1));
            }
            var ms2Xics2 = ms2XicConstructor.GetAllXics(fakeMs2Scans2);
            var allGroups2 = xicGroupingEngine.PrecursorFragmentGrouping(ms1Xics, ms2Xics2);
            //Only the first precursor should have a pfGroup returned
            Assert.That(allGroups2.Count, Is.EqualTo(1));
            Assert.That(allGroups2[0].PFpairs.Count, Is.EqualTo(10));
            Assert.That(allGroups2[0].PrecursorXic.ApexScanIndex, Is.EqualTo(2));
        }
    }
}
