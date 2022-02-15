using Chemistry;
using FlashLFQ;
using MassSpectrometry;
using MathNet.Numerics.Distributions;
using MathNet.Numerics.Statistics;
using MzLibUtil;
using NUnit.Framework;
using Proteomics.AminoAcidPolymer;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using UsefulProteomicsDatabases;
using ChromatographicPeak = FlashLFQ.ChromatographicPeak;
using Stopwatch = System.Diagnostics.Stopwatch;

namespace Test
{
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    internal class TestFlashLFQ
    {
        private static Stopwatch Stopwatch { get; set; }

        [SetUp]
        public static void Setuppp()
        {
            Stopwatch = new Stopwatch();
            Stopwatch.Start();
        }

        [TearDown]
        public static void TearDown()
        {
            Console.WriteLine($"Analysis time: {Stopwatch.Elapsed.Hours}h {Stopwatch.Elapsed.Minutes}m {Stopwatch.Elapsed.Seconds}s");
        }

        [Test]
        public static void TestFlashLfq()
        {
            // get the raw file paths
            SpectraFileInfo raw = new SpectraFileInfo(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", @"sliced-raw.raw"), "a", 0, 0, 0);
            SpectraFileInfo mzml = new SpectraFileInfo(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", @"sliced-mzml.mzml"), "a", 1, 0, 0);

            // create some PSMs
            var pg = new ProteinGroup("MyProtein", "gene", "org");
            Identification id1 = new Identification(raw, "EGFQVADGPLYR", "EGFQVADGPLYR", 1350.65681, 94.12193, 2, new List<ProteinGroup> { pg });
            Identification id2 = new Identification(raw, "EGFQVADGPLYR", "EGFQVADGPLYR", 1350.65681, 94.05811, 2, new List<ProteinGroup> { pg });
            Identification id3 = new Identification(mzml, "EGFQVADGPLYR", "EGFQVADGPLYR", 1350.65681, 94.12193, 2, new List<ProteinGroup> { pg });
            Identification id4 = new Identification(mzml, "EGFQVADGPLYR", "EGFQVADGPLYR", 1350.65681, 94.05811, 2, new List<ProteinGroup> { pg });

            // create the FlashLFQ engine
            FlashLfqEngine engine = new FlashLfqEngine(new List<Identification> { id1, id2, id3, id4 }, normalize: true, maxThreads: 1);

            // run the engine
            var results = engine.Run();

            // check raw results
            Assert.That(results.Peaks[raw].Count == 1);
            Assert.That(results.Peaks[raw].First().Intensity > 0);
            Assert.That(!results.Peaks[raw].First().IsMbrPeak);
            Assert.That(results.PeptideModifiedSequences["EGFQVADGPLYR"].GetIntensity(raw) > 0);
            Assert.That(results.ProteinGroups["MyProtein"].GetIntensity(raw) > 0);

            // check mzml results
            Assert.That(results.Peaks[mzml].Count == 1);
            Assert.That(results.Peaks[mzml].First().Intensity > 0);
            Assert.That(!results.Peaks[mzml].First().IsMbrPeak);
            Assert.That(results.PeptideModifiedSequences["EGFQVADGPLYR"].GetIntensity(mzml) > 0);
            Assert.That(results.ProteinGroups["MyProtein"].GetIntensity(mzml) > 0);

            // check that condition normalization worked
            int int1 = (int)System.Math.Round(results.Peaks[mzml].First().Intensity, 0);
            int int2 = (int)System.Math.Round(results.Peaks[raw].First().Intensity, 0);
            Assert.That(int1 == int2);

            // test peak output
            results.WriteResults(
                Path.Combine(TestContext.CurrentContext.TestDirectory, @"peaks.tsv"),
                Path.Combine(TestContext.CurrentContext.TestDirectory, @"modSeq.tsv"),
                Path.Combine(TestContext.CurrentContext.TestDirectory, @"protein.tsv"),
                null,
                true);
        }

        [Test]
        public static void TestEnvelopQuantification()
        {
            Loaders.LoadElements();

            double monoIsotopicMass = 1350.65681;
            double massOfAveragine = 111.1254;
            double numberOfAveragines = monoIsotopicMass / massOfAveragine;

            double averageC = 4.9384 * numberOfAveragines;
            double averageH = 7.7583 * numberOfAveragines;
            double averageO = 1.4773 * numberOfAveragines;
            double averageN = 1.3577 * numberOfAveragines;
            double averageS = 0.0417 * numberOfAveragines;

            ChemicalFormula myFormula = ChemicalFormula.ParseFormula(
                "C" + (int)Math.Round(averageC) +
                "H" + (int)Math.Round(averageH) +
                "O" + (int)Math.Round(averageO) +
                "N" + (int)Math.Round(averageN) +
                "S" + (int)Math.Round(averageS));


            // get the raw file paths
            SpectraFileInfo mzml = new SpectraFileInfo(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", @"sliced-mzml.mzml"), "a", 0, 0, 0);

            // create some PSMs
            var pg = new ProteinGroup("MyProtein", "gene", "org");

            Identification id3 = new Identification(mzml, "", "1", 1350.65681, 94.12193, 2, new List<ProteinGroup> { pg }, myFormula);
            Identification id4 = new Identification(mzml, "", "2", 1350.65681, 94.05811, 2, new List<ProteinGroup> { pg }, myFormula);

            // create the FlashLFQ engine
            FlashLfqEngine engine = new FlashLfqEngine(new List<Identification> { id3, id4 }, normalize: true);

            // run the engine
            var results = engine.Run();

            Assert.IsTrue(results.Peaks.First().Value.First().Intensity > 0);
        }

        [Test]
        public static void TestFlashLfqNormalization()
        {
            // ********************************* check biorep normalization *********************************
            // get the raw file paths
            SpectraFileInfo raw = new SpectraFileInfo(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", @"sliced-raw.raw"), "a", 0, 0, 0);
            SpectraFileInfo mzml = new SpectraFileInfo(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", @"sliced-mzml.mzml"), "a", 1, 0, 0);

            // create some PSMs
            var pg = new ProteinGroup("MyProtein", "gene", "org");
            Identification id1 = new Identification(raw, "EGFQVADGPLYR", "EGFQVADGPLYR", 1350.65681, 94.12193, 2, new List<ProteinGroup> { pg });
            Identification id2 = new Identification(mzml, "EGFQVADGPLYR", "EGFQVADGPLYR", 1350.65681, 94.12193, 2, new List<ProteinGroup> { pg });

            // create the FlashLFQ engine
            var results = new FlashLfqEngine(new List<Identification> { id1, id2 }, normalize: true).Run();

            // check that biorep normalization worked
            int int1 = (int)System.Math.Round(results.Peaks[mzml].First().Intensity, 0);
            int int2 = (int)System.Math.Round(results.Peaks[raw].First().Intensity, 0);
            Assert.That(int1 > 0);
            Assert.That(int1 == int2);

            // ********************************* check condition normalization *********************************
            raw = new SpectraFileInfo(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", @"sliced-raw.raw"), "a", 0, 0, 0);
            mzml = new SpectraFileInfo(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", @"sliced-mzml.mzml"), "b", 0, 0, 0);

            id1 = new Identification(raw, "EGFQVADGPLYR", "EGFQVADGPLYR", 1350.65681, 94.12193, 2, new List<ProteinGroup> { pg });
            id2 = new Identification(mzml, "EGFQVADGPLYR", "EGFQVADGPLYR", 1350.65681, 94.12193, 2, new List<ProteinGroup> { pg });

            results = new FlashLfqEngine(new List<Identification> { id1, id2 }, normalize: true).Run();

            int int3 = (int)System.Math.Round(results.Peaks[mzml].First().Intensity, 0);
            int int4 = (int)System.Math.Round(results.Peaks[raw].First().Intensity, 0);
            Assert.That(int3 > 0);
            Assert.That(int3 == int4);

            // ********************************* check techrep normalization *********************************
            raw = new SpectraFileInfo(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", @"sliced-raw.raw"), "a", 0, 0, 0);
            mzml = new SpectraFileInfo(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", @"sliced-mzml.mzml"), "a", 0, 1, 0);

            id1 = new Identification(raw, "EGFQVADGPLYR", "EGFQVADGPLYR", 1350.65681, 94.12193, 2, new List<ProteinGroup> { pg });
            id2 = new Identification(mzml, "EGFQVADGPLYR", "EGFQVADGPLYR", 1350.65681, 94.12193, 2, new List<ProteinGroup> { pg });

            results = new FlashLfqEngine(new List<Identification> { id1, id2 }, normalize: true).Run();

            int int5 = (int)System.Math.Round(results.Peaks[mzml].First().Intensity, 0);
            int int6 = (int)System.Math.Round(results.Peaks[raw].First().Intensity, 0);
            Assert.That(int5 > 0);
            Assert.That(int5 == int6);

            Assert.That(int1 == int3);
            Assert.That(int1 == int5);

            // ********************************* check fraction normalization *********************************
            raw = new SpectraFileInfo(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", @"sliced-raw.raw"), "a", 0, 0, 0);
            var raw2 = new SpectraFileInfo(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", @"sliced-raw.raw"), "a", 0, 0, 1);
            mzml = new SpectraFileInfo(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", @"sliced-mzml.mzml"), "a", 1, 0, 0);
            var mzml2 = new SpectraFileInfo(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", @"sliced-mzml.mzml"), "a", 1, 0, 1);

            id1 = new Identification(raw, "EGFQVADGPLYR", "EGFQVADGPLYR", 1350.65681, 94.12193, 2, new List<ProteinGroup> { pg });
            id2 = new Identification(raw2, "EGFQVADGPLYR", "EGFQVADGPLYR", 1350.65681, 94.12193, 2, new List<ProteinGroup> { pg });
            var id3 = new Identification(mzml, "EGFQVADGPLYR", "EGFQVADGPLYR", 1350.65681, 94.12193, 2, new List<ProteinGroup> { pg });
            var id4 = new Identification(mzml2, "EGFQVADGPLYR", "EGFQVADGPLYR", 1350.65681, 94.12193, 2, new List<ProteinGroup> { pg });

            results = new FlashLfqEngine(new List<Identification> { id1, id2, id3, id4 }, normalize: true).Run();

            int int7 = (int)System.Math.Round(results.PeptideModifiedSequences["EGFQVADGPLYR"].GetIntensity(raw) + results.PeptideModifiedSequences["EGFQVADGPLYR"].GetIntensity(raw2));
            int int8 = (int)System.Math.Round(results.PeptideModifiedSequences["EGFQVADGPLYR"].GetIntensity(mzml) + results.PeptideModifiedSequences["EGFQVADGPLYR"].GetIntensity(mzml2));
            Assert.That(int7 > 0);
            Assert.That(int7 == int8);
        }

        [Test]
        public static void TestFlashLfqMergeResults()
        {
            SpectraFileInfo rawA = new SpectraFileInfo(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", @"sliced-raw.raw"), "a", 0, 0, 0);
            SpectraFileInfo mzmlA = new SpectraFileInfo(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", @"sliced-mzml.mzml"), "a", 0, 1, 0);

            // create some PSMs
            var pgA = new ProteinGroup("MyProtein", "gene", "org");
            Identification id1A = new Identification(rawA, "EGFQVADGPLYR", "EGFQVADGPLYR", 1350.65681, 94.12193, 2, new List<ProteinGroup> { pgA });
            Identification id2A = new Identification(rawA, "EGFQVADGPLYR", "EGFQVADGPLYR", 1350.65681, 94.05811, 2, new List<ProteinGroup> { pgA });
            Identification id3A = new Identification(mzmlA, "EGFQVADGPLYR", "EGFQVADGPLYR", 1350.65681, 94.12193, 2, new List<ProteinGroup> { pgA });
            Identification id4A = new Identification(mzmlA, "EGFQVADGPLYR", "EGFQVADGPLYR", 1350.65681, 94.05811, 2, new List<ProteinGroup> { pgA });

            // create the FlashLFQ engine
            FlashLfqEngine engineA = new FlashLfqEngine(new List<Identification> { id1A, id2A, id3A, id4A });

            // run the engine
            var resultsA = engineA.Run();

            SpectraFileInfo rawB = new SpectraFileInfo(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", @"sliced-raw.raw"), "b", 0, 0, 0);
            SpectraFileInfo mzmlB = new SpectraFileInfo(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", @"sliced-mzml.mzml"), "b", 0, 1, 0);

            // create some PSMs
            var pgB = new ProteinGroup("MyProtein", "gene", "org");
            Identification id1 = new Identification(rawB, "EGFQVADGPLYR", "EGFQVADGPLYR", 1350.65681, 94.12193, 2, new List<ProteinGroup> { pgB });
            Identification id2 = new Identification(rawB, "EGFQVADGPLYR", "EGFQVADGPLYR", 1350.65681, 94.05811, 2, new List<ProteinGroup> { pgB });
            Identification id3 = new Identification(mzmlB, "EGFQVADGPLYR", "EGFQVADGPLYR", 1350.65681, 94.12193, 2, new List<ProteinGroup> { pgB });
            Identification id4 = new Identification(mzmlB, "EGFQVADGPLYR", "EGFQVADGPLYR", 1350.65681, 94.05811, 2, new List<ProteinGroup> { pgB });

            // create the FlashLFQ engine
            FlashLfqEngine engineB = new FlashLfqEngine(new List<Identification> { id1, id2, id3, id4 });

            // run the engine
            var resultsB = engineB.Run();

            resultsA.MergeResultsWith(resultsB);
            Assert.AreEqual(4, resultsA.Peaks.Count);
            Assert.AreEqual(1, resultsA.PeptideModifiedSequences.Count);
            Assert.AreEqual(1, resultsA.ProteinGroups.Count);
            Assert.AreEqual(4, resultsA.SpectraFiles.Count);
        }

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
            var pg = new ProteinGroup("MyProtein", "gene", "org");
            Identification id1 = new Identification(file1, "PEPTIDE", "PEPTIDE",
                new Proteomics.AminoAcidPolymer.Peptide("PEPTIDE").MonoisotopicMass, file1Rt[0] + 0.001, 1, new List<ProteinGroup> { pg });
            Identification id2 = new Identification(file1, "PEPTIDEV", "PEPTIDEV",
                new Proteomics.AminoAcidPolymer.Peptide("PEPTIDEV").MonoisotopicMass, file1Rt[1] + 0.001, 1, new List<ProteinGroup> { pg });
            Identification id3 = new Identification(file1, "PEPTIDEVV", "PEPTIDEVV",
                new Proteomics.AminoAcidPolymer.Peptide("PEPTIDEVV").MonoisotopicMass, file1Rt[2] + 0.001, 1, new List<ProteinGroup> { pg });
            Identification id4 = new Identification(file1, "PEPTIDEVVV", "PEPTIDEVVV",
                new Proteomics.AminoAcidPolymer.Peptide("PEPTIDEVVV").MonoisotopicMass, file1Rt[3] + 0.001, 1, new List<ProteinGroup> { pg });
            Identification id5 = new Identification(file1, "PEPTIDEVVVV", "PEPTIDEVVVV",
                new Proteomics.AminoAcidPolymer.Peptide("PEPTIDEVVVV").MonoisotopicMass, file1Rt[4] + 0.001, 1, new List<ProteinGroup> { pg });

            Identification id6 = new Identification(file2, "PEPTIDE", "PEPTIDE",
                new Proteomics.AminoAcidPolymer.Peptide("PEPTIDE").MonoisotopicMass, file2Rt[0] + 0.001, 1, new List<ProteinGroup> { pg });
            Identification id7 = new Identification(file2, "PEPTIDEV", "PEPTIDEV",
                new Proteomics.AminoAcidPolymer.Peptide("PEPTIDEV").MonoisotopicMass, file2Rt[1] + 0.001, 1, new List<ProteinGroup> { pg });
            // missing ID 8 - MBR feature
            Identification id9 = new Identification(file2, "PEPTIDEVVV", "PEPTIDEVVV",
                new Proteomics.AminoAcidPolymer.Peptide("PEPTIDEVVV").MonoisotopicMass, file2Rt[3] + 0.001, 1, new List<ProteinGroup> { pg });
            Identification id10 = new Identification(file2, "PEPTIDEVVVV", "PEPTIDEVVVV",
                new Proteomics.AminoAcidPolymer.Peptide("PEPTIDEVVVV").MonoisotopicMass, file2Rt[4] + 0.001, 1, new List<ProteinGroup> { pg });

            // create the FlashLFQ engine
            FlashLfqEngine engine = new FlashLfqEngine(new List<Identification> { id1, id2, id3, id4, id5, id6, id7, id9, id10 }, matchBetweenRuns: true);

            // run the engine
            var results = engine.Run();

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

        [Test]
        public static void TestPeakSplittingLeft()
        {
            string fileToWrite = "myMzml.mzML";
            string peptide = "PEPTIDE";
            double intensity = 1e6;

            Loaders.LoadElements();

            // generate mzml file

            // 1 MS1 scan per peptide
            MsDataScan[] scans = new MsDataScan[10];
            double[] intensityMultipliers = { 1, 3, 1, 1, 3, 5, 10, 5, 3, 1 };

            for (int s = 0; s < scans.Length; s++)
            {
                ChemicalFormula cf = new Proteomics.AminoAcidPolymer.Peptide(peptide).GetChemicalFormula();
                IsotopicDistribution dist = IsotopicDistribution.GetDistribution(cf, 0.125, 1e-8);
                double[] mz = dist.Masses.Select(v => v.ToMz(1)).ToArray();
                double[] intensities = dist.Intensities.Select(v => v * intensity * intensityMultipliers[s]).ToArray();

                // add the scan
                scans[s] = new MsDataScan(massSpectrum: new MzSpectrum(mz, intensities, false), oneBasedScanNumber: s + 1, msnOrder: 1, isCentroid: true,
                    polarity: Polarity.Positive, retentionTime: 1.0 + s / 10.0, scanWindowRange: new MzRange(400, 1600), scanFilter: "f",
                    mzAnalyzer: MZAnalyzerType.Orbitrap, totalIonCurrent: intensities.Sum(), injectionTime: 1.0, noiseData: null, nativeId: "scan=" + (s + 1));
            }

            // write the .mzML
            IO.MzML.MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(new FakeMsDataFile(scans),
                Path.Combine(TestContext.CurrentContext.TestDirectory, fileToWrite), false);

            // set up spectra file info
            SpectraFileInfo file1 = new SpectraFileInfo(Path.Combine(TestContext.CurrentContext.TestDirectory, fileToWrite), "", 0, 0, 0);

            // create some PSMs
            var pg = new ProteinGroup("MyProtein", "gene", "org");

            Identification id1 = new Identification(file1, peptide, peptide,
                new Proteomics.AminoAcidPolymer.Peptide(peptide).MonoisotopicMass, 1.7 + 0.001, 1, new List<ProteinGroup> { pg });

            // create the FlashLFQ engine
            FlashLfqEngine engine = new FlashLfqEngine(new List<Identification> { id1 });

            // run the engine
            var results = engine.Run();
            ChromatographicPeak peak = results.Peaks.First().Value.First();

            Assert.That(peak.Apex.IndexedPeak.RetentionTime == 1.6);
            Assert.That(peak.SplitRT == 1.3);
            Assert.That(!peak.IsotopicEnvelopes.Any(p => p.IndexedPeak.RetentionTime < 1.3));
            Assert.That(peak.IsotopicEnvelopes.Count == 6);
        }

        [Test]
        public static void TestPeakSplittingRight()
        {
            string fileToWrite = "myMzml.mzML";
            string peptide = "PEPTIDE";
            double intensity = 1e6;

            Loaders.LoadElements();

            // generate mzml file

            // 1 MS1 scan per peptide
            MsDataScan[] scans = new MsDataScan[10];
            double[] intensityMultipliers = { 1, 3, 5, 10, 5, 3, 1, 1, 3, 1 };

            for (int s = 0; s < scans.Length; s++)
            {
                ChemicalFormula cf = new Proteomics.AminoAcidPolymer.Peptide(peptide).GetChemicalFormula();
                IsotopicDistribution dist = IsotopicDistribution.GetDistribution(cf, 0.125, 1e-8);
                double[] mz = dist.Masses.Select(v => v.ToMz(1)).ToArray();
                double[] intensities = dist.Intensities.Select(v => v * intensity * intensityMultipliers[s]).ToArray();

                // add the scan
                scans[s] = new MsDataScan(massSpectrum: new MzSpectrum(mz, intensities, false), oneBasedScanNumber: s + 1, msnOrder: 1, isCentroid: true,
                    polarity: Polarity.Positive, retentionTime: 1.0 + s / 10.0, scanWindowRange: new MzRange(400, 1600), scanFilter: "f",
                    mzAnalyzer: MZAnalyzerType.Orbitrap, totalIonCurrent: intensities.Sum(), injectionTime: 1.0, noiseData: null, nativeId: "scan=" + (s + 1));
            }

            // write the .mzML
            IO.MzML.MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(new FakeMsDataFile(scans),
                Path.Combine(TestContext.CurrentContext.TestDirectory, fileToWrite), false);

            // set up spectra file info
            SpectraFileInfo file1 = new SpectraFileInfo(Path.Combine(TestContext.CurrentContext.TestDirectory, fileToWrite), "", 0, 0, 0);

            // create some PSMs
            var pg = new ProteinGroup("MyProtein", "gene", "org");

            Identification id1 = new Identification(file1, peptide, peptide,
                new Proteomics.AminoAcidPolymer.Peptide(peptide).MonoisotopicMass, 1.3 + 0.001, 1, new List<ProteinGroup> { pg });

            // create the FlashLFQ engine
            FlashLfqEngine engine = new FlashLfqEngine(new List<Identification> { id1 });

            // run the engine
            var results = engine.Run();
            ChromatographicPeak peak = results.Peaks.First().Value.First();

            Assert.That(peak.Apex.IndexedPeak.RetentionTime == 1.3);
            Assert.That(peak.SplitRT == 1.6);
            Assert.That(!peak.IsotopicEnvelopes.Any(p => p.IndexedPeak.RetentionTime > 1.6));
            Assert.That(peak.IsotopicEnvelopes.Count == 6);
        }

        [Test]
        public static void TestPeakSplittingRightWithEmptyScan()
        {
            string fileToWrite = "myMzml.mzML";
            string peptide = "PEPTIDE";
            double intensity = 1e6;

            Loaders.LoadElements();

            // generate mzml file

            // 1 MS1 scan per peptide
            MsDataScan[] scans = new MsDataScan[10];
            double[] intensityMultipliers = { 1, 3, 5, 10, 5, 3, 1, 1, 3, 1 };

            for (int s = 0; s < scans.Length; s++)
            {
                ChemicalFormula cf = new Proteomics.AminoAcidPolymer.Peptide(peptide).GetChemicalFormula();
                IsotopicDistribution dist = IsotopicDistribution.GetDistribution(cf, 0.125, 1e-8);
                double[] mz = dist.Masses.Select(v => v.ToMz(1)).ToArray();
                double[] intensities = dist.Intensities.Select(v => v * intensity * intensityMultipliers[s]).ToArray();

                if (s == 7)
                {
                    mz = new[] { 401.0 };
                    intensities = new[] { 1000.0 };
                }

                // add the scan
                scans[s] = new MsDataScan(massSpectrum: new MzSpectrum(mz, intensities, false), oneBasedScanNumber: s + 1, msnOrder: 1, isCentroid: true,
                    polarity: Polarity.Positive, retentionTime: 1.0 + s / 10.0, scanWindowRange: new MzRange(400, 1600), scanFilter: "f",
                    mzAnalyzer: MZAnalyzerType.Orbitrap, totalIonCurrent: intensities.Sum(), injectionTime: 1.0, noiseData: null, nativeId: "scan=" + (s + 1));
            }

            // write the .mzML
            IO.MzML.MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(new FakeMsDataFile(scans),
                Path.Combine(TestContext.CurrentContext.TestDirectory, fileToWrite), false);

            // set up spectra file info
            SpectraFileInfo file1 = new SpectraFileInfo(Path.Combine(TestContext.CurrentContext.TestDirectory, fileToWrite), "", 0, 0, 0);

            // create some PSMs
            var pg = new ProteinGroup("MyProtein", "gene", "org");

            Identification id1 = new Identification(file1, peptide, peptide,
                new Proteomics.AminoAcidPolymer.Peptide(peptide).MonoisotopicMass, 1.3 + 0.001, 1, new List<ProteinGroup> { pg });

            // create the FlashLFQ engine
            FlashLfqEngine engine = new FlashLfqEngine(new List<Identification> { id1 });

            // run the engine
            var results = engine.Run();
            ChromatographicPeak peak = results.Peaks.First().Value.First();

            Assert.That(peak.Apex.IndexedPeak.RetentionTime == 1.3);
            Assert.That(peak.SplitRT == 1.6);
            Assert.That(!peak.IsotopicEnvelopes.Any(p => p.IndexedPeak.RetentionTime > 1.6));
            Assert.That(peak.IsotopicEnvelopes.Count == 6);
        }

        [Test]
        public static void TestPeakSplittingLeftWithEmptyScan()
        {
            string fileToWrite = "myMzml.mzML";
            string peptide = "PEPTIDE";
            double intensity = 1e6;

            Loaders.LoadElements();

            // generate mzml file

            // 1 MS1 scan per peptide
            MsDataScan[] scans = new MsDataScan[10];
            double[] intensityMultipliers = { 1, 3, 1, 1, 3, 5, 10, 5, 3, 1 };

            for (int s = 0; s < scans.Length; s++)
            {
                ChemicalFormula cf = new Proteomics.AminoAcidPolymer.Peptide(peptide).GetChemicalFormula();
                IsotopicDistribution dist = IsotopicDistribution.GetDistribution(cf, 0.125, 1e-8);
                double[] mz = dist.Masses.Select(v => v.ToMz(1)).ToArray();
                double[] intensities = dist.Intensities.Select(v => v * intensity * intensityMultipliers[s]).ToArray();

                if (s == 2)
                {
                    mz = new[] { 401.0 };
                    intensities = new[] { 1000.0 };
                }

                // add the scan
                scans[s] = new MsDataScan(massSpectrum: new MzSpectrum(mz, intensities, false), oneBasedScanNumber: s + 1, msnOrder: 1, isCentroid: true,
                    polarity: Polarity.Positive, retentionTime: 1.0 + s / 10.0, scanWindowRange: new MzRange(400, 1600), scanFilter: "f",
                    mzAnalyzer: MZAnalyzerType.Orbitrap, totalIonCurrent: intensities.Sum(), injectionTime: 1.0, noiseData: null, nativeId: "scan=" + (s + 1));
            }

            // write the .mzML
            IO.MzML.MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(new FakeMsDataFile(scans),
                Path.Combine(TestContext.CurrentContext.TestDirectory, fileToWrite), false);

            // set up spectra file info
            SpectraFileInfo file1 = new SpectraFileInfo(Path.Combine(TestContext.CurrentContext.TestDirectory, fileToWrite), "", 0, 0, 0);

            // create some PSMs
            var pg = new ProteinGroup("MyProtein", "gene", "org");

            Identification id1 = new Identification(file1, peptide, peptide,
                new Proteomics.AminoAcidPolymer.Peptide(peptide).MonoisotopicMass, 1.3 + 0.001, 1, new List<ProteinGroup> { pg });

            // create the FlashLFQ engine
            FlashLfqEngine engine = new FlashLfqEngine(new List<Identification> { id1 });

            // run the engine
            var results = engine.Run();
            ChromatographicPeak peak = results.Peaks.First().Value.First();

            Assert.That(peak.Apex.IndexedPeak.RetentionTime == 1.6);
            Assert.That(peak.SplitRT == 1.3);
            Assert.That(!peak.IsotopicEnvelopes.Any(p => p.IndexedPeak.RetentionTime < 1.3));
            Assert.That(peak.IsotopicEnvelopes.Count == 6);
        }

        [Test]
        public static void TestToString()
        {
            // many of these are just to check that the ToString methods don't cause crashes
            var indexedPeak = new IndexedMassSpectralPeak(1.0, 2.0, 4, 5.0);
            Assert.That(indexedPeak.ToString().Equals("1.000; 4"));

            var spectraFile = new SpectraFileInfo("myFullPath", "", 0, 0, 0);
            string spectraString = spectraFile.ToString();

            var proteinGroup = new ProteinGroup("Accession", "Gene", "Organism");
            string pgString = proteinGroup.ToString(new List<SpectraFileInfo> { spectraFile });

            var identification = new Identification(
                spectraFile, "PEPTIDE", "PEPTIDE", 1.0, 2.0, 3,
                new List<ProteinGroup> { proteinGroup });
            string idString = identification.ToString();

            var chromPeak = new ChromatographicPeak(identification, false, spectraFile);
            string chromPeakString = chromPeak.ToString();
            chromPeak.CalculateIntensityForThisFeature(true);
            string peakAfterCalculatingIntensity = chromPeak.ToString();

            var peptide = new FlashLFQ.Peptide("PEPTIDE", "PEPTIDE", true, new HashSet<ProteinGroup>());
            var peptideString = peptide.ToString(new List<SpectraFileInfo> { spectraFile });
            Assert.That(peptideString == "PEPTIDE\tPEPTIDE\t\t\t\t0\tNotDetected\t");

            peptide = new FlashLFQ.Peptide("PEPTIDE", "PEPTIDE", true, new HashSet<ProteinGroup> { proteinGroup });
            peptideString = peptide.ToString(new List<SpectraFileInfo> { spectraFile });
            Assert.That(peptideString == "PEPTIDE\tPEPTIDE\tAccession\tGene\tOrganism\t0\tNotDetected\t");

            peptide = new FlashLFQ.Peptide("PEPTIDE", "PEPTIDE", true, new HashSet<ProteinGroup> { proteinGroup, new ProteinGroup("Accession2", "Gene2", "Organism2") });
            peptideString = peptide.ToString(new List<SpectraFileInfo> { spectraFile });
            Assert.That(peptideString == "PEPTIDE\tPEPTIDE\tAccession;Accession2\tGene;Gene2\tOrganism;Organism2\t0\tNotDetected\t");
        }

        [Test]
        public static void TestNotFound()
        {
            FlashLFQ.Peptide p = new FlashLFQ.Peptide("Seq", "SEQ", true, new HashSet<ProteinGroup>());
            var notFound = p.GetDetectionType(new SpectraFileInfo("", "", 0, 0, 0));
            Assert.That(notFound == DetectionType.NotDetected);
        }

        [Test]
        public static void TestAmbiguous()
        {
            // get the raw file paths
            SpectraFileInfo mzml = new SpectraFileInfo(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", @"sliced-mzml.mzml"), "a", 0, 1, 0);

            // create some PSMs
            var pg = new ProteinGroup("MyProtein", "gene", "org");
            Identification id3 = new Identification(mzml, "EGFQVADGPLRY", "EGFQVADGPLRY", 1350.65681, 94.12193, 2, new List<ProteinGroup> { pg });
            Identification id4 = new Identification(mzml, "EGFQVADGPLYR", "EGFQVADGPLYR", 1350.65681, 94.05811, 2, new List<ProteinGroup> { pg });

            // create the FlashLFQ engine
            FlashLfqEngine engine = new FlashLfqEngine(new List<Identification> { id3, id4 });

            // run the engine
            var results = engine.Run();

            Assert.That(results.Peaks[mzml].Count == 1);
            Assert.That(results.Peaks[mzml].First().Intensity > 0);
            Assert.That(!results.Peaks[mzml].First().IsMbrPeak);
            Assert.That(results.Peaks[mzml].First().NumIdentificationsByFullSeq == 2);
            Assert.That(results.PeptideModifiedSequences["EGFQVADGPLYR"].GetIntensity(mzml) == 0);
            Assert.That(results.PeptideModifiedSequences["EGFQVADGPLRY"].GetIntensity(mzml) == 0);
            Assert.That(results.ProteinGroups["MyProtein"].GetIntensity(mzml) == 0);

            // test peak output
            results.WriteResults(
                Path.Combine(TestContext.CurrentContext.TestDirectory, @"peaks.tsv"),
                Path.Combine(TestContext.CurrentContext.TestDirectory, @"modSeq.tsv"),
                Path.Combine(TestContext.CurrentContext.TestDirectory, @"protein.tsv"),
                null, true);
        }

        [Test]
        public static void TestMatchBetweenRunsWithNoIdsInCommon()
        {
            List<string> filesToWrite = new List<string> { "mzml_1", "mzml_2" };
            List<string> pepSequences = new List<string> { "PEPTIDE", "PEPTIDEV", "PEPTIDEVV", "PEPTIDEVVV", "PEPTIDEVVVV" };
            double intensity = 1e6;

            double[] file1Rt = new double[] { 1.01, 1.02, 1.03, 1.04, 1.05 };
            double[] file2Rt = new double[] { 1.015, 1.030, 1.036, 1.050, 1.065 };

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
            var pg = new ProteinGroup("MyProtein", "gene", "org");
            var myMbrProteinGroup = new ProteinGroup("MyMbrProtein", "MbrGene", "org");

            Identification id1 = new Identification(file1, "PEPTIDE", "PEPTIDE",
                new Proteomics.AminoAcidPolymer.Peptide("PEPTIDE").MonoisotopicMass, file1Rt[0] + 0.001, 1, new List<ProteinGroup> { pg });
            Identification id2 = new Identification(file1, "PEPTIDEV", "PEPTIDEV",
                new Proteomics.AminoAcidPolymer.Peptide("PEPTIDEV").MonoisotopicMass, file1Rt[1] + 0.001, 1, new List<ProteinGroup> { pg });
            Identification id3 = new Identification(file1, "PEPTIDEVV", "PEPTIDEVV",
                new Proteomics.AminoAcidPolymer.Peptide("PEPTIDEVV").MonoisotopicMass, file1Rt[2] + 0.001, 1, new List<ProteinGroup> { myMbrProteinGroup });
            Identification id4 = new Identification(file1, "PEPTIDEVVV", "PEPTIDEVVV",
                new Proteomics.AminoAcidPolymer.Peptide("PEPTIDEVVV").MonoisotopicMass, file1Rt[3] + 0.001, 1, new List<ProteinGroup> { pg });
            Identification id5 = new Identification(file1, "PEPTIDEVVVV", "PEPTIDEVVVV",
                new Proteomics.AminoAcidPolymer.Peptide("PEPTIDEVVVV").MonoisotopicMass, file1Rt[4] + 0.001, 1, new List<ProteinGroup> { pg });

            Identification id6 = new Identification(file2, "PEPTIED", "PEPTIED",
                new Proteomics.AminoAcidPolymer.Peptide("PEPTIED").MonoisotopicMass, file2Rt[0] + 0.001, 1, new List<ProteinGroup> { pg });
            Identification id7 = new Identification(file2, "PEPTIEDV", "PEPTIEDV",
                new Proteomics.AminoAcidPolymer.Peptide("PEPTIEDV").MonoisotopicMass, file2Rt[1] + 0.001, 1, new List<ProteinGroup> { pg });
            // missing ID 8 - MBR feature
            Identification id9 = new Identification(file2, "PEPTIEDVVV", "PEPTIEDVVV",
                new Proteomics.AminoAcidPolymer.Peptide("PEPTIEDVVV").MonoisotopicMass, file2Rt[3] + 0.001, 1, new List<ProteinGroup> { pg });
            Identification id10 = new Identification(file2, "PEPTIEDVVVV", "PEPTIEDVVVV",
                new Proteomics.AminoAcidPolymer.Peptide("PEPTIEDVVVV").MonoisotopicMass, file2Rt[4] + 0.001, 1, new List<ProteinGroup> { pg });

            FlashLfqEngine engine = new FlashLfqEngine(new List<Identification> { id1, id2, id3, id4, id5, id6, id7, id9, id10 }, matchBetweenRuns: true);
            var results = engine.Run();

            // no assertions - just don't crash
        }

        [Test]
        public static void TestFlashLfqDoesNotRemovePeptides()
        {
            Loaders.LoadElements();

            Residue x = new Residue("a", 'a', "a", Chemistry.ChemicalFormula.ParseFormula("C{13}6H12N{15}2O"), ModificationSites.All); //+8 lysine
            Residue lightLysine = Residue.GetResidue('K');

            Residue.AddNewResiduesToDictionary(new List<Residue> { new Residue("heavyLysine", 'a', "a", x.ThisChemicalFormula, ModificationSites.All) });

            SpectraFileInfo fileInfo = new SpectraFileInfo(Path.Combine(TestContext.CurrentContext.TestDirectory, @"SilacTest.mzML"), "", 0, 0, 0);
            FlashLfqEngine engine = new FlashLfqEngine(
                new List<Identification>
                {
                    new Identification(fileInfo,"RDILSSNNQHGILPLSWNIPELVNMGQWK","RDILSSNNQHGILPLSWNIPELVNM[Common Variable:Oxidation on M]GQWK",3374.7193792,98.814005,3,new List<FlashLFQ.ProteinGroup>{new FlashLFQ.ProteinGroup("P01027","C3","Mus") },null, true),
                    new Identification(fileInfo,"RDILSSNNQHGILPLSWNIPELVNMGQWa","RDILSSNNQHGILPLSWNIPELVNM[Common Variable:Oxidation on M]GQWa",3382.733578,98.814005,3,new List<FlashLFQ.ProteinGroup>{new FlashLFQ.ProteinGroup("P01027+8.014","C3","Mus") },null, true),
                    new Identification(fileInfo,"RDILSSNNQHGILPLSWNIPELVNMGQWK","RDILSSNNQHGILPLSWNIPELVNM[Common Variable:Oxidation on M]GQWK",3374.7193792,98.7193782,4,new List<FlashLFQ.ProteinGroup>{new FlashLFQ.ProteinGroup("P01027","C3","Mus") },null, true),
                    new Identification(fileInfo,"RDILSSNNQHGILPLSWNIPELVNMGQWa","RDILSSNNQHGILPLSWNIPELVNM[Common Variable:Oxidation on M]GQWa",3382.733578,98.7193782,4,new List<FlashLFQ.ProteinGroup>{new FlashLFQ.ProteinGroup("P01027+8.014","C3","Mus") },null, true),
                },
                ppmTolerance: 5,
                silent: true,
                maxThreads: 7
                );
            var results = engine.Run();
            Assert.IsTrue(results.PeptideModifiedSequences.Count == 2);
        }

        [Test]
        public static void TestBayesianProteinQuantification()
        {
            // this mostly just tests that the Bayesian quant algorithm produces "reasonable"
            // estimates of the posterior error probability (PEP) and mean differences.

            // the idea here is to create a protein, assign some peptides with intensities to it,
            // and compute the probability that the protein is changing between the conditions "a" and "b".
            // the intensities in condition "b" are about double that of condition "a", and the results
            // of the Bayesian estimation should reflect that.
            ProteinGroup pg = new ProteinGroup("Accession", "Gene", "Organism");
            var p = new FlashLFQ.Peptide("PEPTIDE", "PEPTIDE", true, new HashSet<ProteinGroup> { pg });

            var files = new List<SpectraFileInfo>
            {
                new SpectraFileInfo("a1", "a", 0, 0, 0),
                new SpectraFileInfo("a2", "a", 1, 0, 0),
                new SpectraFileInfo("a3", "a", 2, 0, 0),
                new SpectraFileInfo("b1", "b", 0, 0, 0),
                new SpectraFileInfo("b2", "b", 1, 0, 0),
                new SpectraFileInfo("b3", "b", 2, 0, 0)
            };

            var res = new FlashLfqResults(files, new List<Identification>
            {
                new Identification(null, "SEQUENCE", "SEQUENCE", 0, 0, 0, new List<ProteinGroup>{ new ProteinGroup("Accession", "Gene", "Organism") })
            });

            FlashLFQ.Peptide peptide = res.PeptideModifiedSequences.First().Value;
            ProteinGroup proteinGroup = res.ProteinGroups.First().Value;

            peptide.SetIntensity(files[0], 900);
            peptide.SetIntensity(files[1], 1000);
            peptide.SetIntensity(files[2], 1100);

            peptide.SetIntensity(files[3], 1950);
            peptide.SetIntensity(files[4], 2000);
            peptide.SetIntensity(files[5], 2050);

            var engine = new ProteinQuantificationEngine(res, maxThreads: 1, controlCondition: "a", randomSeed: 0, foldChangeCutoff: 0.1);
            engine.Run();

            var quantResult = (UnpairedProteinQuantResult)proteinGroup.ConditionToQuantificationResults["b"];

            Assert.That(Math.Round(quantResult.NullHypothesisInterval.Value, 3) == 0.100);
            Assert.That(Math.Round(quantResult.PosteriorErrorProbability, 3) == 0.000);
            Assert.That(Math.Round(quantResult.FoldChangePointEstimate, 3) == 1.000);
            Assert.That(quantResult.ConditionsWithPeptideSampleQuantities["a"].Count == 3);
            Assert.That(quantResult.ConditionsWithPeptideSampleQuantities["b"].Count == 3);

            string filepath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"bayesianProteinQuant.tsv");
            res.WriteResults(null, null, null, filepath, true);

            var textResults = File.ReadAllLines(filepath);
            Assert.That(textResults.Length == 2);
            var line = textResults[1].Split(new char[] { '\t' });
            Assert.That(Math.Round(double.Parse(line[17]), 3) == 0.000);
            File.Delete(filepath);

            // try with some missing values
            peptide.SetIntensity(files[1], 0);
            peptide.SetIntensity(files[5], 0);

            proteinGroup.ConditionToQuantificationResults.Clear();
            engine = new ProteinQuantificationEngine(res, maxThreads: 1, controlCondition: "a", randomSeed: 2, foldChangeCutoff: 0.5);
            engine.Run();

            quantResult = (UnpairedProteinQuantResult)proteinGroup.ConditionToQuantificationResults["b"];

            Assert.That(Math.Round(quantResult.NullHypothesisInterval.Value, 3) == 0.500);
            Assert.That(Math.Round(quantResult.PosteriorErrorProbability, 3) == 1.000);
            Assert.That(Math.Round(quantResult.FoldChangePointEstimate, 3) == 0.922);
            Assert.That(quantResult.ConditionsWithPeptideSampleQuantities["a"].Count == 2);
            Assert.That(quantResult.ConditionsWithPeptideSampleQuantities["b"].Count == 2);

            //// try with paired samples
            //peptide.SetIntensity(files[0], 100);
            //peptide.SetIntensity(files[1], 1000);
            //peptide.SetIntensity(files[2], 10000);

            //peptide.SetIntensity(files[3], 210);
            //peptide.SetIntensity(files[4], 2200);
            //peptide.SetIntensity(files[5], 21500);

            //proteinGroup.ConditionToQuantificationResults.Clear();
            //engine = new ProteinQuantificationEngine(res, maxThreads: 1, controlCondition: "a", randomSeed: 3, pairedSamples: true);
            //engine.Run();

            //quantResult = proteinGroup.ConditionToQuantificationResults["b"];

            //Assert.That(Math.Round(quantResult.PosteriorErrorProbability, 3) == 0.098);
            //Assert.That(Math.Round(quantResult.FoldChangePointEstimate, 3) == 1.103);
            //Assert.That(quantResult.PeptideFoldChangeMeasurements.Count == 1);
            //Assert.That(quantResult.PeptideFoldChangeMeasurements.SelectMany(v => v.foldChanges).Count() == 3);
        }

        [Test]
        public static void RealDataMbrTest()
        {
            string psmFile = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", @"PSMsForMbrTest.psmtsv");

            SpectraFileInfo f1r1 = new SpectraFileInfo(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", @"f1r1_sliced_mbr.raw"), "a", 0, 0, 0);
            SpectraFileInfo f1r2 = new SpectraFileInfo(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", @"f1r2_sliced_mbr.raw"), "a", 1, 0, 0);

            List<Identification> ids = new List<Identification>();
            Dictionary<string, ProteinGroup> allProteinGroups = new Dictionary<string, ProteinGroup>();
            foreach (string line in File.ReadAllLines(psmFile))
            {
                var split = line.Split(new char[] { '\t' });

                if (split.Contains("File Name") || string.IsNullOrWhiteSpace(line))
                {
                    continue;
                }

                SpectraFileInfo file = null;

                if (split[0].Contains("f1r1"))
                {
                    file = f1r1;
                }
                else if (split[0].Contains("f1r2"))
                {
                    file = f1r2;
                }

                string baseSequence = split[12];
                string fullSequence = split[13];
                double monoMass = double.Parse(split[21]);
                double rt = double.Parse(split[2]);
                int z = (int)double.Parse(split[6]);
                var proteins = split[24].Split(new char[] { '|' });
                List<ProteinGroup> proteinGroups = new List<ProteinGroup>();
                foreach (var protein in proteins)
                {
                    if (allProteinGroups.TryGetValue(protein, out var proteinGroup))
                    {
                        proteinGroups.Add(proteinGroup);
                    }
                    else
                    {
                        allProteinGroups.Add(protein, new ProteinGroup(protein, "", ""));
                        proteinGroups.Add(allProteinGroups[protein]);
                    }
                }

                Identification id = new Identification(file, baseSequence, fullSequence, monoMass, rt, z, proteinGroups);
                ids.Add(id);
            }

            var engine = new FlashLfqEngine(ids, matchBetweenRuns: true, requireMsmsIdInCondition: false, maxThreads: 1);
            var results = engine.Run();

            var f1r1MbrResults = results
                .PeptideModifiedSequences
                .Where(p => p.Value.GetDetectionType(f1r1) == DetectionType.MBR && p.Value.GetDetectionType(f1r2) == DetectionType.MSMS).ToList();

            Assert.That(f1r1MbrResults.Count >= 132);

            var f1r2MbrResults = results.PeptideModifiedSequences
                .Where(p => p.Value.GetDetectionType(f1r1) == DetectionType.MSMS && p.Value.GetDetectionType(f1r2) == DetectionType.MBR).ToList();

            Assert.That(f1r2MbrResults.Count >= 77);

            List<(double, double)> peptideIntensities = new List<(double, double)>();

            foreach (var peptide in f1r1MbrResults)
            {
                double mbrIntensity = Math.Log(peptide.Value.GetIntensity(f1r1));
                double msmsIntensity = Math.Log(peptide.Value.GetIntensity(f1r2));
                peptideIntensities.Add((mbrIntensity, msmsIntensity));
            }

            double corr = Correlation.Pearson(peptideIntensities.Select(p => p.Item1), peptideIntensities.Select(p => p.Item2));
            Assert.That(corr > 0.8);

            peptideIntensities.Clear();
            foreach (var peptide in f1r2MbrResults)
            {
                double mbrIntensity = Math.Log(peptide.Value.GetIntensity(f1r2));
                double msmsIntensity = Math.Log(peptide.Value.GetIntensity(f1r1));
                peptideIntensities.Add((mbrIntensity, msmsIntensity));
            }

            corr = Correlation.Pearson(peptideIntensities.Select(p => p.Item1), peptideIntensities.Select(p => p.Item2));

            Assert.That(corr > 0.7);

            // the "requireMsmsIdInCondition" field requires that at least one MS/MS identification from a protein
            // has to be observed in a condition for match-between-runs
            f1r1.Condition = "b";
            engine = new FlashLfqEngine(ids, matchBetweenRuns: true, requireMsmsIdInCondition: true, maxThreads: 1);
            results = engine.Run();
            var proteinsObservedInF1 = ids.Where(p => p.FileInfo == f1r1).SelectMany(p => p.ProteinGroups).Distinct().ToList();
            var proteinsObservedInF2 = ids.Where(p => p.FileInfo == f1r2).SelectMany(p => p.ProteinGroups).Distinct().ToList();
            var proteinsObservedInF1ButNotF2 = proteinsObservedInF1.Except(proteinsObservedInF2).ToList();
            foreach (ProteinGroup protein in proteinsObservedInF1ButNotF2)
            {
                Assert.That(results.ProteinGroups[protein.ProteinGroupName].GetIntensity(f1r2) == 0);
            }
        }

        [Test]
        public static void ProteoformPeakfindingTest()
        {
            Loaders.LoadElements();
            string sequence =
                "PEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDE" +
                "PEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDE" +
                "PEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDE" +
                "PEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDE" +
                "PEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDE" +
                "PEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDEPEPTIDE";

            int charge = 7;

            MsDataScan[] scans = new MsDataScan[10];

            ChemicalFormula cf = new Proteomics.AminoAcidPolymer.Peptide(sequence).GetChemicalFormula();

            IsotopicDistribution dist = IsotopicDistribution.GetDistribution(cf, 0.125, 1e-8);
            double[] mz = dist.Masses.Select(v => v.ToMz(charge)).ToArray();
            double[] intensities = dist.Intensities.Select(v => v * 1e7).ToArray();

            List<double> filteredMzs = new List<double>();
            List<double> filteredIntensities = new List<double>();

            for (int i = 0; i < mz.Length; i++)
            {
                if (intensities[i] < 1000)
                {
                    continue;
                }

                filteredMzs.Add(mz[i]);
                filteredIntensities.Add(intensities[i]);
            }

            for (int s = 0; s < scans.Length; s++)
            {
                double rt = s + 0.1;

                // add the scan
                scans[s] = new MsDataScan(massSpectrum: new MzSpectrum(filteredMzs.ToArray(), filteredIntensities.ToArray(), false),
                    oneBasedScanNumber: s + 1, msnOrder: 1, isCentroid: true,
                    polarity: Polarity.Positive, retentionTime: rt, scanWindowRange: new MzRange(400, 1600), scanFilter: "f",
                    mzAnalyzer: MZAnalyzerType.Orbitrap, totalIonCurrent: intensities.Sum(), injectionTime: 1.0, noiseData: null, nativeId: "scan=" + (s + 1));
            }

            // write the .mzML
            string path = Path.Combine(TestContext.CurrentContext.TestDirectory, "testFile.mzML");
            IO.MzML.MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(new FakeMsDataFile(scans), path, false);

            // set up spectra file info
            SpectraFileInfo file1 = new SpectraFileInfo(path, "a", 0, 0, 0);

            // create the PSM
            var pg = new ProteinGroup("MyProtein", "gene", "org");
            Identification id1 = new Identification(file1, sequence, sequence,
                cf.MonoisotopicMass, scans[0].RetentionTime + 0.001, charge, new List<ProteinGroup> { pg });

            // create the FlashLFQ engine
            FlashLfqEngine engine = new FlashLfqEngine(new List<Identification> { id1 });

            // run the engine
            var results = engine.Run();

            Assert.That((int)results.PeptideModifiedSequences[sequence].GetIntensity(file1) == 1386491);
            ChromatographicPeak peak = results.Peaks[file1].First(p => p.Identifications.First().ModifiedSequence == sequence);
            Assert.That(Math.Round(peak.MassError, 3) == 0);
            Assert.That(peak.IsotopicEnvelopes.Count == 10);
        }

        [Test]
        public static void TestUseSharedPeptidesForQuant()
        {
            ProteinGroup pg1 = new ProteinGroup("Accession1", "Gene1", "Organism");
            ProteinGroup pg2 = new ProteinGroup("Accession2", "Gene2", "Organism");

            var files = new List<SpectraFileInfo>
            {
                new SpectraFileInfo("a1", "a", 0, 0, 0)
            };

            // this peptide is shared between protein groups
            var id = new Identification(files[0], "PEPTIDE", "PEPTIDE", 0, 0, 1, new List<ProteinGroup> { pg1, pg2 });

            // this is unique to a single protein group
            var id2 = new Identification(files[0], "PEPTIDEA", "PEPTIDEA", 0, 0, 1, new List<ProteinGroup> { pg1 });

            var res = new FlashLfqResults(files, new List<Identification> { id, id2 });

            res.PeptideModifiedSequences["PEPTIDE"].SetIntensity(files[0], 1000); // shared peptide
            res.PeptideModifiedSequences["PEPTIDEA"].SetIntensity(files[0], 2000); // unique peptide
            res.PeptideModifiedSequences["PEPTIDE"].SetDetectionType(files[0], DetectionType.MSMS);
            res.PeptideModifiedSequences["PEPTIDEA"].SetDetectionType(files[0], DetectionType.MSMS);

            res.CalculateProteinResultsTop3(useSharedPeptides: true);
            Assert.That(res.ProteinGroups["Accession1"].GetIntensity(files[0]) == 3000); // protein intensity should be the sum of shared+unique peptides

            res.CalculateProteinResultsTop3(useSharedPeptides: false);
            Assert.That(res.ProteinGroups["Accession1"].GetIntensity(files[0]) == 2000); // protein intensity should be from the unique peptide only
        }

        [Test]
        public static void TestIntensityDependentProteinQuant()
        {
            List<double> diffAbundantFractions = new List<double> { 0.1 };
            int peptidesPerProtein = 4;

            foreach (var differentiallyAbundantFraction in diffAbundantFractions)
            {
                // create the files, peptides, and proteins
                int numTotalProteins = 300;
                int numDifferentiallyAbundant = (int)(numTotalProteins * differentiallyAbundantFraction);

                Random randomSource = new Random(0);
                List<ProteinGroup> pgs = new List<ProteinGroup>();
                List<Identification> ids = new List<Identification>();

                var files = new List<SpectraFileInfo>
            {
                new SpectraFileInfo("a1", "a", 0, 0, 0),
                new SpectraFileInfo("a2", "a", 1, 0, 0),
                new SpectraFileInfo("a3", "a", 2, 0, 0),
                new SpectraFileInfo("b1", "b", 0, 0, 0),
                new SpectraFileInfo("b2", "b", 1, 0, 0),
                new SpectraFileInfo("b3", "b", 2, 0, 0)
            };

                for (int i = 0; i < numTotalProteins; i++)
                {
                    string organism = i < numDifferentiallyAbundant ? "changing" : "not_changing";

                    var pg = new ProteinGroup("protein_" + i, "", organism);

                    for (int j = 0; j < peptidesPerProtein; j++)
                    {
                        string peptideName = "peptide_" + i + "_" + j;

                        foreach (var file in files)
                        {
                            var id = new Identification(file, peptideName, peptideName, 0, 0, 0, new List<ProteinGroup> { pg });
                            ids.Add(id);
                        }
                    }

                    pgs.Add(pg);
                }

                var res = new FlashLfqResults(files, ids);

                // set intensity values
                var differentiallyAbundantProteins = new HashSet<ProteinGroup>(pgs.Take(numDifferentiallyAbundant));
                Normal n = new Normal(20, 5, randomSource);

                foreach (var peptide in res.PeptideModifiedSequences)
                {
                    var protein = peptide.Value.ProteinGroups.First();

                    bool differentiallyAbundant = differentiallyAbundantProteins.Contains(protein);

                    double groupAIntensity = 0;
                    double groupBIntensity = 0;
                    while (groupAIntensity <= 0)
                    {
                        groupAIntensity = n.Sample();
                    }

                    if (differentiallyAbundant)
                    {
                        groupBIntensity = groupAIntensity + 1;
                    }
                    else
                    {
                        groupBIntensity = groupAIntensity;
                    }

                    foreach (var file in files)
                    {
                        if (file.Condition == "a")
                        {
                            double noise = Normal.Sample(randomSource, 0, 1 / (groupAIntensity * 0.1));
                            double fileIntensity = Math.Pow(2, groupAIntensity + noise);
                            peptide.Value.SetIntensity(file, fileIntensity);
                        }
                        else
                        {
                            double noise = Normal.Sample(randomSource, 0, 1 / (groupBIntensity * 0.1));
                            double fileIntensity = Math.Pow(2, groupBIntensity + noise);
                            peptide.Value.SetIntensity(file, fileIntensity);
                        }
                    }
                }

                // run the protein quant engine
                ProteinQuantificationEngine engine = new ProteinQuantificationEngine(res, -1, "a", false, 0.1, 0, 1000, 1000, false);
                engine.Run();

                // test the quant engine results
                var proteinQuantResults = res.ProteinGroups.Values.Select(p => (UnpairedProteinQuantResult)p.ConditionToQuantificationResults["b"]).ToList();

                var proteinsBelow5percentFdr = proteinQuantResults.Where(p => p.FalseDiscoveryRate < 0.05).ToList();
                var fdp = proteinsBelow5percentFdr.Count(p => p.Protein.Organism == "not_changing") / (double)proteinsBelow5percentFdr.Count;

                // require FDR control
                Assert.That(fdp < 0.05);

                // require 70% sensitivity
                // in this unit test there are 4 peptides per protein
                // if the number of peptides is increased, the classifier will be more sensitive
                Assert.That(proteinsBelow5percentFdr.Count >= proteinQuantResults.Count * differentiallyAbundantFraction * 0.7);
            }
        }

        [Test]
        public static void TestAmbiguousFraction()
        {
            SpectraFileInfo fraction1 = new SpectraFileInfo("", "", 0, 0, fraction: 0);
            SpectraFileInfo fraction2 = new SpectraFileInfo("", "", 0, 0, fraction: 1);
            Identification id1 = new Identification(fraction1, "peptide1", "peptide1", 0, 0, 0, new List<ProteinGroup>());

            Identification id2 = new Identification(fraction2, "peptide1", "peptide1", 0, 0, 0, new List<ProteinGroup>());
            Identification id3 = new Identification(fraction2, "peptide2", "peptide2", 0, 0, 0, new List<ProteinGroup>());

            ChromatographicPeak peak1 = new ChromatographicPeak(id1, false, fraction1);
            ChromatographicPeak peak2 = new ChromatographicPeak(id2, false, fraction1);
            peak2.Identifications.Add(id3);

            peak1.ResolveIdentifications();
            peak2.ResolveIdentifications();

            peak1.IsotopicEnvelopes.Add(new FlashLFQ.IsotopicEnvelope(new IndexedMassSpectralPeak(0, 0, 0, 0), 1, 1000));
            peak2.IsotopicEnvelopes.Add(new FlashLFQ.IsotopicEnvelope(new IndexedMassSpectralPeak(0, 0, 0, 0), 1, 10000));

            peak1.CalculateIntensityForThisFeature(false);
            peak2.CalculateIntensityForThisFeature(false);

            FlashLfqResults res = new FlashLfqResults(new List<SpectraFileInfo> { fraction1, fraction2 }, new List<Identification> { id1, id2, id3 });
            res.Peaks[fraction1].Add(peak1);
            res.Peaks[fraction2].Add(peak2);
            res.CalculatePeptideResults();

            var peptides = res.PeptideModifiedSequences;
            Assert.That(peptides["peptide1"].GetIntensity(fraction1) == 0);
            Assert.That(peptides["peptide1"].GetIntensity(fraction2) == 0);
            Assert.That(peptides["peptide2"].GetIntensity(fraction1) == 0);
            Assert.That(peptides["peptide2"].GetIntensity(fraction2) == 0);

            Assert.That(peptides["peptide1"].GetDetectionType(fraction1) == DetectionType.MSMS);
            Assert.That(peptides["peptide1"].GetDetectionType(fraction2) == DetectionType.MSMSAmbiguousPeakfinding);
            Assert.That(peptides["peptide2"].GetDetectionType(fraction1) == DetectionType.NotDetected);
            Assert.That(peptides["peptide2"].GetDetectionType(fraction2) == DetectionType.MSMSAmbiguousPeakfinding);
        }

        [Test]
        public static void TestMedianPolishProteinQuant()
        {
            // 2 groups, 2 samples each, 3 fractions per sample
            List<SpectraFileInfo> spectraFileInfos = new List<SpectraFileInfo>
            {
                new SpectraFileInfo("", "group1", 0, 0, 0),
                new SpectraFileInfo("", "group1", 0, 0, 1),
                new SpectraFileInfo("", "group1", 0, 0, 2),

                new SpectraFileInfo("", "group1", 1, 0, 0),
                new SpectraFileInfo("", "group1", 1, 0, 1),
                new SpectraFileInfo("", "group1", 1, 0, 2),

                new SpectraFileInfo("", "group2", 0, 0, 0),
                new SpectraFileInfo("", "group2", 0, 0, 1),
                new SpectraFileInfo("", "group2", 0, 0, 2),

                new SpectraFileInfo("", "group2", 1, 0, 0),
                new SpectraFileInfo("", "group2", 1, 0, 1),
                new SpectraFileInfo("", "group2", 1, 0, 2),
            };

            // 2 proteins
            ProteinGroup pg1 = new ProteinGroup("accession1", "gene1", "organism1");
            ProteinGroup pg2 = new ProteinGroup("accession2", "gene2", "organism1");

            // 3 peptides, 1 peptide is shared b/w protein1 and protein2
            FlashLFQ.Peptide pep1 = new FlashLFQ.Peptide("PEPTIDE", "PEPTIDE1", true, new HashSet<ProteinGroup> { pg1 });
            FlashLFQ.Peptide pep2 = new FlashLFQ.Peptide("PEPTIDEE", "PEPTIDE2", true, new HashSet<ProteinGroup> { pg1 });
            FlashLFQ.Peptide pep3 = new FlashLFQ.Peptide("PEPTIDEEE", "PEPTIDE3", true, new HashSet<ProteinGroup> { pg1, pg2 }); // shared peptide

            FlashLfqResults res = new FlashLfqResults(spectraFileInfos, new List<Identification>());
            res.PeptideModifiedSequences.Add(pep1.Sequence, pep1);
            res.PeptideModifiedSequences.Add(pep2.Sequence, pep2);
            res.PeptideModifiedSequences.Add(pep3.Sequence, pep3);

            res.ProteinGroups.Add(pg1.ProteinGroupName, pg1);
            res.ProteinGroups.Add(pg2.ProteinGroupName, pg2);

            Random r = new Random(1);
            List<FlashLFQ.Peptide> peptides = new List<FlashLFQ.Peptide> { pep1, pep2, pep3 };
            Normal randomIonizationEfficiencyGenerator = new Normal(20, 2, r);

            // create peptide quantities
            foreach (var peptide in peptides)
            {
                double logIonizationEfficiency = randomIonizationEfficiencyGenerator.Sample();

                foreach (var group in spectraFileInfos.GroupBy(p => p.Condition))
                {
                    foreach (var sample in group.GroupBy(p => p.BiologicalReplicate))
                    {
                        foreach (var fraction in sample.OrderBy(p => p.Fraction))
                        {
                            // peptide 1 will elute in fraction 1, peptide 2 elutes in fraction 2...
                            int peptideNumber = peptides.IndexOf(peptide);

                            if (peptideNumber == fraction.Fraction)
                            {
                                if (group.Key == "group2" && sample.Key == 1 && peptideNumber == 1)
                                {
                                    // create a missing peptide value
                                    // the protein should still get quantified because its other peptide had a valid value
                                }
                                else
                                {
                                    peptide.SetIntensity(fraction, Math.Pow(2, logIonizationEfficiency));
                                    peptide.SetDetectionType(fraction, DetectionType.MSMS);
                                }
                            }
                        }
                    }
                }
            }

            // do the protein quant, skipping shared peptides
            res.CalculateProteinResultsMedianPolish(useSharedPeptides: false);

            // write/read the protein quantification output
            string filepath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"proteinQuant.tsv");
            res.WriteResults(null, null, filepath, null, true);

            var textResults = File.ReadAllLines(filepath);
            Assert.That(textResults.Length == 3);

            // the header should show the names of the samples, not the fractionated file names
            var header = textResults[0].Split(new char[] { '\t' });
            Assert.That(header[3] == "Intensity_group1_1");
            Assert.That(header[4] == "Intensity_group1_2");
            Assert.That(header[5] == "Intensity_group2_1");
            Assert.That(header[6] == "Intensity_group2_2");

            // the quantities reported for protein1 should have no missing values and should be identical
            var protein1Results = textResults[1].Split(new char[] { '\t' });
            Assert.That((int)double.Parse(protein1Results[3]) == 1501270);
            Assert.That((int)double.Parse(protein1Results[4]) == 1501270);
            Assert.That((int)double.Parse(protein1Results[5]) == 1501270);
            Assert.That((int)double.Parse(protein1Results[6]) == 1501270);

            // protein2 doesn't get quantified because it only has 1 peptide and it's shared,
            // and we said to not quantified shared peptides
            var protein2Results = textResults[2].Split(new char[] { '\t' });
            Assert.That(double.Parse(protein2Results[3]) == 0);
            Assert.That(double.Parse(protein2Results[4]) == 0);
            Assert.That(double.Parse(protein2Results[5]) == 0);
            Assert.That(double.Parse(protein2Results[6]) == 0);

            File.Delete(filepath);
        }
    }
}