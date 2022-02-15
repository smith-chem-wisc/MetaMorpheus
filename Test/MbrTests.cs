using EngineLayer;
using EngineLayer.FdrAnalysis;
using EngineLayer.HistogramAnalysis;
using EngineLayer.Localization;
using EngineLayer.ModificationAnalysis;
using FlashLFQ;
using MassSpectrometry;
using Chemistry;
using MzLibUtil;
using Proteomics;
using MathNet.Numerics.Distributions;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.IO.Compression;
using System.Linq;
using System.Text;
using UsefulProteomicsDatabases;
using NUnit.Framework;

namespace Test
{
    [TestFixture]
    public static class MbrTest
    {
        [Test]
        public static void TestFlashLfqMatchBetweenRuns()
        {
            List<string> filesToWrite = new List<string> { "mzml_1", "mzml_2" };
            List<string> pepSequences = new List<string> { "PEPTIDE", "PEPTIDEV", "PEPTIDEVV", "PEPTIDEVVV", "PEPTIDEVVVV" };
            double intensity = 1e6;

            double[] file1Rt = new double[] { 1.01, 1.02, 1.03, 1.04, 1.05 };
            double[] file2Rt = new double[] { 1.00, 1.025, 1.04, 1.055, 1.070 };

            Loaders.LoadElements();

            // generate mzml files (5 peptides each)
            for (int f = 0; f < filesToWrite.Count; f++)
            {
                // 1 MS1 scan per peptide
                MsDataScan[] scans = new MsDataScan[5];

                for (int p = 0; p < pepSequences.Count; p++)
                {
                    ChemicalFormula cf = new Proteomics.AminoAcidPolymer.Peptide(pepSequences[p]).GetChemicalFormula();
                    IsotopicDistribution dist = IsotopicDistribution.GetDistribution(cf, 0.125, 1e-8);
                    double[] mz = dist.Masses.Select(v => v.ToMz(1)).ToArray();
                    double[] intensities = dist.Intensities.Select(v => v * intensity).ToArray();
                    double rt;
                    if (f == 0)
                    {
                        rt = file1Rt[p];
                    }
                    else
                    {
                        rt = file2Rt[p];
                    }

                    // add the scan
                    scans[p] = new MsDataScan(massSpectrum: new MzSpectrum(mz, intensities, false), oneBasedScanNumber: p + 1, msnOrder: 1, isCentroid: true,
                        polarity: Polarity.Positive, retentionTime: rt, scanWindowRange: new MzRange(400, 1600), scanFilter: "f",
                        mzAnalyzer: MZAnalyzerType.Orbitrap, totalIonCurrent: intensities.Sum(), injectionTime: 1.0, noiseData: null, nativeId: "scan=" + (p + 1));
                }

                // write the .mzML
                IO.MzML.MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(new FakeMsDataFile(scans),
                    Path.Combine(TestContext.CurrentContext.TestDirectory, filesToWrite[f] + ".mzML"), false);
            }

            // set up spectra file info
            SpectraFileInfo file1 = new SpectraFileInfo(Path.Combine(TestContext.CurrentContext.TestDirectory, filesToWrite[0] + ".mzML"), "a", 0, 0, 0);
            SpectraFileInfo file2 = new SpectraFileInfo(Path.Combine(TestContext.CurrentContext.TestDirectory, filesToWrite[1] + ".mzML"), "a", 1, 0, 0);

            // create some PSMs
            var pg = new FlashLFQ.ProteinGroup("MyProtein", "gene", "org");
            Identification id1 = new Identification(file1, "PEPTIDE", "PEPTIDE",
                new Proteomics.AminoAcidPolymer.Peptide("PEPTIDE").MonoisotopicMass, file1Rt[0] + 0.001, 1, new List<FlashLFQ.ProteinGroup> { pg });
            Identification id2 = new Identification(file1, "PEPTIDEV", "PEPTIDEV",
                new Proteomics.AminoAcidPolymer.Peptide("PEPTIDEV").MonoisotopicMass, file1Rt[1] + 0.001, 1, new List<FlashLFQ.ProteinGroup> { pg });
            Identification id3 = new Identification(file1, "PEPTIDEVV", "PEPTIDEVV",
                new Proteomics.AminoAcidPolymer.Peptide("PEPTIDEVV").MonoisotopicMass, file1Rt[2] + 0.001, 1, new List<FlashLFQ.ProteinGroup> { pg });
            Identification id4 = new Identification(file1, "PEPTIDEVVV", "PEPTIDEVVV",
                new Proteomics.AminoAcidPolymer.Peptide("PEPTIDEVVV").MonoisotopicMass, file1Rt[3] + 0.001, 1, new List<FlashLFQ.ProteinGroup> { pg });
            Identification id5 = new Identification(file1, "PEPTIDEVVVV", "PEPTIDEVVVV",
                new Proteomics.AminoAcidPolymer.Peptide("PEPTIDEVVVV").MonoisotopicMass, file1Rt[4] + 0.001, 1, new List<FlashLFQ.ProteinGroup> { pg });

            Identification id6 = new Identification(file2, "PEPTIDE", "PEPTIDE",
                new Proteomics.AminoAcidPolymer.Peptide("PEPTIDE").MonoisotopicMass, file2Rt[0] + 0.001, 1, new List<FlashLFQ.ProteinGroup> { pg });
            Identification id7 = new Identification(file2, "PEPTIDEV", "PEPTIDEV",
                new Proteomics.AminoAcidPolymer.Peptide("PEPTIDEV").MonoisotopicMass, file2Rt[1] + 0.001, 1, new List<FlashLFQ.ProteinGroup> { pg });
            // missing ID 8 - MBR feature
            Identification id9 = new Identification(file2, "PEPTIDEVVV", "PEPTIDEVVV",
                new Proteomics.AminoAcidPolymer.Peptide("PEPTIDEVVV").MonoisotopicMass, file2Rt[3] + 0.001, 1, new List<FlashLFQ.ProteinGroup> { pg });
            Identification id10 = new Identification(file2, "PEPTIDEVVVV", "PEPTIDEVVVV",
                new Proteomics.AminoAcidPolymer.Peptide("PEPTIDEVVVV").MonoisotopicMass, file2Rt[4] + 0.001, 1, new List<FlashLFQ.ProteinGroup> { pg });

            // create the FlashLFQ engine
            FlashLfqEngine engine = new FlashLfqEngine(new List<Identification> { id1, id2, id3, id4, id5, id6, id7, id9, id10 }, matchBetweenRuns: true);

            // run the engine
            var results = engine.Run();

            // Accesses the zeroBasedMS1Index of the ApexScan

            // TODO: Figure out how to recover OneBasedIndex using the ZeroBasedMs1ScanIndex
            // Literally need to write a method that allows intercoversion OR
            // Modify IndexedMassSpectralPeak such that the onebased scan number is included
            int apexIndex = results.Peaks[file1][0].Apex.IndexedPeak.ZeroBasedMs1ScanIndex;

            Assert.That(results.Peaks[file2].Count == 5);
            Assert.That(results.Peaks[file2].Where(p => p.IsMbrPeak).Count() == 1);

            var peak = results.Peaks[file2].Where(p => p.IsMbrPeak).First();
            var otherFilePeak = results.Peaks[file1].Where(p => p.Identifications.First().BaseSequence ==
                peak.Identifications.First().BaseSequence).First();

            Assert.That(peak.Intensity > 0);
            Assert.That(peak.Intensity == otherFilePeak.Intensity);

            Assert.That(results.Peaks[file1].Count == 5);
            Assert.That(results.Peaks[file1].Where(p => p.IsMbrPeak).Count() == 0);

            Assert.That(results.ProteinGroups["MyProtein"].GetIntensity(file1) > 0);
            Assert.That(results.ProteinGroups["MyProtein"].GetIntensity(file2) > 0);
        }



        // TODO: Move these classes to their own files!! This is probably not great practice
        public class FakeMsDataFile : MsDataFile
        {
            public FakeMsDataFile(MsDataScan[] FakeScans) : base(FakeScans, new SourceFile(@"scan number only nativeID format", "mzML format", null, "SHA-1", @"C:\fake.mzML", null))
            {
                this.Scans = FakeScans;
            }

        }

    }



}
