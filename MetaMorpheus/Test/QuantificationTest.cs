using Chemistry;
using EngineLayer;
using FlashLFQ;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework; using Assert = NUnit.Framework.Legacy.ClassicAssert;
using Proteomics;
using Omics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using Readers;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using EngineLayer.DatabaseLoading;
using Omics.Modifications;
using TaskLayer;

namespace Test
{
    internal static class QuantificationTest
    {
        [Test]
        public static void WriteExperimentalDesignTest()
        {
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"ExperimentalDesignTest");
            Directory.CreateDirectory(outputFolder);

            List<SpectraFileInfo> spectraFiles = new List<SpectraFileInfo>();

            spectraFiles.Add(new SpectraFileInfo(Path.Combine(outputFolder, @"myFile.1.raw"), "condition1", 0, 0, 0));
            spectraFiles.Add(new SpectraFileInfo(Path.Combine(outputFolder, @"myFile.2.raw"), "condition1", 1, 0, 0));
            spectraFiles.Add(new SpectraFileInfo(Path.Combine(outputFolder, @"myFile.3.raw"), "condition1", 2, 0, 0));

            ExperimentalDesign.WriteExperimentalDesignToFile(spectraFiles);

            var readIn = ExperimentalDesign.ReadExperimentalDesign(
                Path.Combine(outputFolder, @"ExperimentalDesign.tsv"),
                spectraFiles.Select(p => p.FullFilePathWithExtension).ToList(),
                out var errors);

            Assert.That(!errors.Any());
            Assert.That(readIn.Count == 3);

            Directory.Delete(outputFolder, true);
        }

        [Test]
        public static void TestExperimentalDesignErrors()
        {
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"ExperimentalDesignTest");
            Directory.CreateDirectory(outputFolder);

            // test non-consecutive bioreps (a warning, not an error: the numbers are kept as given)
            List<SpectraFileInfo> spectraFiles = new List<SpectraFileInfo>();

            spectraFiles.Add(new SpectraFileInfo(Path.Combine(outputFolder, @"myFile.1.raw"), "condition1", 0, 0, 0));
            spectraFiles.Add(new SpectraFileInfo(Path.Combine(outputFolder, @"myFile.2.raw"), "condition1", 1, 0, 0));
            spectraFiles.Add(new SpectraFileInfo(Path.Combine(outputFolder, @"myFile.3.raw"), "condition1", 3, 0, 0));

            ExperimentalDesign.WriteExperimentalDesignToFile(spectraFiles);

            var readIn = ExperimentalDesign.ReadExperimentalDesign(
                Path.Combine(outputFolder, @"ExperimentalDesign.tsv"),
                spectraFiles.Select(p => p.FullFilePathWithExtension).ToList(),
                out var errors);

            Assert.That(errors, Is.Empty);
            Assert.That(readIn.Select(p => p.BiologicalReplicate), Is.EquivalentTo(new[] { 0, 1, 3 }));
            Assert.That(ExperimentalDesign.GetWarningsInExperimentalDesign(readIn), Has.Count.EqualTo(1));

            // test non-consecutive fractions, here 2, 3, 4 (a warning, not an error)
            spectraFiles.Clear();

            spectraFiles.Add(new SpectraFileInfo(Path.Combine(outputFolder, @"myFile.1.raw"), "condition1", 0, 0, 1));
            spectraFiles.Add(new SpectraFileInfo(Path.Combine(outputFolder, @"myFile.2.raw"), "condition1", 0, 0, 2));
            spectraFiles.Add(new SpectraFileInfo(Path.Combine(outputFolder, @"myFile.3.raw"), "condition1", 0, 0, 3));

            ExperimentalDesign.WriteExperimentalDesignToFile(spectraFiles);

            readIn = ExperimentalDesign.ReadExperimentalDesign(
                Path.Combine(outputFolder, @"ExperimentalDesign.tsv"),
                spectraFiles.Select(p => p.FullFilePathWithExtension).ToList(),
                out errors);

            Assert.That(errors, Is.Empty);
            Assert.That(ExperimentalDesign.GetWarningsInExperimentalDesign(readIn), Has.Count.EqualTo(1));

            // test non-consecutive techreps, here 1, 3, 4 (a warning, not an error)
            spectraFiles.Clear();

            spectraFiles.Add(new SpectraFileInfo(Path.Combine(outputFolder, @"myFile.1.raw"), "condition1", 0, 0, 0));
            spectraFiles.Add(new SpectraFileInfo(Path.Combine(outputFolder, @"myFile.2.raw"), "condition1", 0, 2, 0));
            spectraFiles.Add(new SpectraFileInfo(Path.Combine(outputFolder, @"myFile.3.raw"), "condition1", 0, 3, 0));

            ExperimentalDesign.WriteExperimentalDesignToFile(spectraFiles);

            readIn = ExperimentalDesign.ReadExperimentalDesign(
                Path.Combine(outputFolder, @"ExperimentalDesign.tsv"),
                spectraFiles.Select(p => p.FullFilePathWithExtension).ToList(),
                out errors);

            Assert.That(errors, Is.Empty);
            Assert.That(ExperimentalDesign.GetWarningsInExperimentalDesign(readIn), Has.Count.EqualTo(1));

            // test duplicates (should produce an error)
            spectraFiles.Clear();

            spectraFiles.Add(new SpectraFileInfo(Path.Combine(outputFolder, @"myFile.1.raw"), "condition1", 0, 0, 0));
            spectraFiles.Add(new SpectraFileInfo(Path.Combine(outputFolder, @"myFile.2.raw"), "condition1", 0, 0, 0));

            ExperimentalDesign.WriteExperimentalDesignToFile(spectraFiles);

            readIn = ExperimentalDesign.ReadExperimentalDesign(
                Path.Combine(outputFolder, @"ExperimentalDesign.tsv"),
                spectraFiles.Select(p => p.FullFilePathWithExtension).ToList(),
                out errors);

            Assert.That(errors.Any());

            // test situation where experimental design does not contain a file (should produce an error)
            spectraFiles.Clear();
            spectraFiles.Add(new SpectraFileInfo(Path.Combine(outputFolder, @"myFile.1.raw"), "condition1", 0, 0, 0));
            ExperimentalDesign.WriteExperimentalDesignToFile(spectraFiles);
            spectraFiles.Add(new SpectraFileInfo(Path.Combine(outputFolder, @"myFile.2.raw"), "condition1", 1, 0, 0));

            readIn = ExperimentalDesign.ReadExperimentalDesign(Path.Combine(outputFolder, @"ExperimentalDesign.tsv"),
                spectraFiles.Select(p => p.FullFilePathWithExtension).ToList(), out errors);

            Assert.That(errors.Any());

            // tests for error messages if experimental design file has non-integer bioreps/techreps/fractions
            List<string> output = new List<string> { "FileName\tCondition\tBiorep\tFraction\tTechrep" };

            output.Add("myFile.1.raw\tcondition1\ta\t0\t0");
            File.WriteAllLines(Path.Combine(outputFolder, @"ExperimentalDesign.tsv"), output);
            readIn = ExperimentalDesign.ReadExperimentalDesign(Path.Combine(outputFolder, @"ExperimentalDesign.tsv"),
                spectraFiles.Select(p => p.FullFilePathWithExtension).ToList(), out errors);
            Assert.That(errors.Any());

            output[1] = "myFile.1.raw\tcondition1\t0\ta\t0";
            File.WriteAllLines(Path.Combine(outputFolder, @"ExperimentalDesign.tsv"), output);
            readIn = ExperimentalDesign.ReadExperimentalDesign(Path.Combine(outputFolder, @"ExperimentalDesign.tsv"),
                spectraFiles.Select(p => p.FullFilePathWithExtension).ToList(), out errors);
            Assert.That(errors.Any());

            output[1] = "myFile.1.raw\tcondition1\t0\t0\ta";
            File.WriteAllLines(Path.Combine(outputFolder, @"ExperimentalDesign.tsv"), output);
            readIn = ExperimentalDesign.ReadExperimentalDesign(Path.Combine(outputFolder, @"ExperimentalDesign.tsv"),
                spectraFiles.Select(p => p.FullFilePathWithExtension).ToList(), out errors);
            Assert.That(errors.Any());

            // tests complete. delete output directory
            Directory.Delete(outputFolder, true);
        }

        /// <summary>
        /// A biorep number can name a subject across conditions (study-wide numbering, as an SDRF often
        /// gives it) or mark a lost sample, so a gap is kept and warned about rather than refused.
        /// Fractions and techreps of the bioreps that are present are still checked.
        /// </summary>
        [Test]
        public static void TestBiorepGapIsAWarningNotAnError()
        {
            static SpectraFileInfo Info(string name, string condition, int oneBasedBiorep, int oneBasedFraction = 1) =>
                new(name + ".raw", condition, oneBasedBiorep - 1, 0, oneBasedFraction - 1);

            // 1, 2, 4: biorep 3 is named as missing
            var gap = new List<SpectraFileInfo> { Info("a1", "A", 1), Info("a2", "A", 2), Info("a4", "A", 4) };
            Assert.That(ExperimentalDesign.GetErrorsInExperimentalDesign(gap), Is.Null);
            Assert.That(ExperimentalDesign.GetWarningsInExperimentalDesign(gap), Is.EqualTo(new[]
            {
                "Condition \"A\": biological replicates 1, 2, 4, quantified as numbered; biological replicate 3 is not in the design. " +
                "A missing number may be a sample or file that was lost or not searched."
            }));

            // study-wide numbering: each condition is judged on its own numbers
            var studyWide = new List<SpectraFileInfo> { Info("c1", "control", 1), Info("c2", "control", 2), Info("t3", "treated", 3), Info("t4", "treated", 4) };
            Assert.That(ExperimentalDesign.GetErrorsInExperimentalDesign(studyWide), Is.Null);
            var warnings = ExperimentalDesign.GetWarningsInExperimentalDesign(studyWide);
            Assert.That(warnings, Has.Count.EqualTo(1));
            Assert.That(warnings[0], Does.StartWith("Condition \"treated\": biological replicates 3, 4, quantified as numbered; biological replicates 1, 2 are not in the design."));

            // 1..N says nothing
            var complete = new List<SpectraFileInfo> { Info("a1", "A", 1), Info("a2", "A", 2) };
            Assert.That(ExperimentalDesign.GetWarningsInExperimentalDesign(complete), Is.Empty);

            // fractions and techreps: a gap is a warning too, named down to the biorep and fraction
            var fractionGap = new List<SpectraFileInfo> { Info("a1f1", "A", 1), Info("a1f3", "A", 1, oneBasedFraction: 3) };
            Assert.That(ExperimentalDesign.GetErrorsInExperimentalDesign(fractionGap), Is.Null);
            Assert.That(ExperimentalDesign.GetWarningsInExperimentalDesign(fractionGap), Is.EqualTo(new[]
            {
                "Condition \"A\" biorep 1: fractions 1, 3, quantified as numbered; fraction 2 is not in the design. " +
                "A missing number may be a sample or file that was lost or not searched."
            }));
            var techrepGap = new List<SpectraFileInfo> { new("t2.raw", "A", 0, 1, 0) };
            Assert.That(ExperimentalDesign.GetErrorsInExperimentalDesign(techrepGap), Is.Null);
            Assert.That(ExperimentalDesign.GetWarningsInExperimentalDesign(techrepGap),
                Has.One.StartWith("Condition \"A\" biorep 1 fraction 1: technical replicates 2, quantified as numbered; technical replicate 1 is not in the design."));

            // a duplicate is still refused
            var duplicate = new List<SpectraFileInfo> { Info("x", "A", 4), Info("y", "A", 4) };
            Assert.That(ExperimentalDesign.GetErrorsInExperimentalDesign(duplicate),
                Is.EqualTo("Duplicates are not allowed:\nCondition \"A\" biorep 4 fraction 1 techrep 1"));
        }

        /// <summary>
        /// End to end, with normalization on: a design numbered study-wide (A: 2, 3; B: 4, 5) used to make
        /// the search skip quantification with one warning. It is now quantified under the numbers given,
        /// and the gaps are warned about, in the log and in results.txt. The first condition has no
        /// biorep 1, which is the normalization reference before mzLib #1422; without #1422 biorep
        /// normalization silently does nothing here.
        /// </summary>
        [Test]
        public static void TestSearchQuantifiesADesignWithABiorepGap()
        {
            var run = RunOnePeptideSearch("TestSearchQuantifiesADesignWithABiorepGap", new[]
            {
                ("A2", "A", 2, 1, 1, 2e6), ("A3", "A", 3, 1, 1, 3e6), ("B4", "B", 4, 1, 1, 4e6), ("B5", "B", 5, 1, 1, 5e6)
            });

            Assert.That(run.Warnings, Has.None.Contain("Skipping quantification"));
            foreach (string expected in new[]
            {
                "Condition \"A\": biological replicates 2, 3, quantified as numbered; biological replicate 1 is not in the design.",
                "Condition \"B\": biological replicates 4, 5, quantified as numbered; biological replicates 1, 2, 3 are not in the design."
            })
            {
                Assert.That(run.Warnings, Has.Some.StartWith(expected));
                Assert.That(run.ResultsTxt, Does.Contain(expected), "the warning is recorded in results.txt too");
            }

            // the protein table's sample columns carry the numbers the design gives
            Assert.That(run.ProteinHeader.Where(h => h.StartsWith("Intensity_")),
                Is.EquivalentTo(new[] { "Intensity_A_2", "Intensity_A_3", "Intensity_B_4", "Intensity_B_5" }));

            // and every file's peptide was quantified and normalized: the files were written at 2, 3, 4
            // and 5 times one intensity, and normalization brings this one peptide back to one value
            var intensities = run.PeptideIntensityByFile.Values.ToList();
            Assert.That(intensities, Has.Count.EqualTo(4));
            Assert.That(intensities, Has.All.GreaterThan(0));
            Assert.That(intensities.Max() / intensities.Min(), Is.LessThan(1.01));
        }

        /// <summary>
        /// End to end, with normalization on: the same five files numbered with a gap in fractions
        /// (1, 3) and in techreps (1, 3), and numbered without one. Both are quantified, the gapped one
        /// with warnings, and every file's normalized peptide intensity is the same in both: a gap
        /// changes nothing but the numbers.
        /// </summary>
        [Test]
        public static void TestSearchQuantifiesFractionAndTechrepGapsAsNumbered()
        {
            var gapped = RunOnePeptideSearch("TestFractionTechrepGapsGapped", new[]
            {
                ("Af1t1", "A", 1, 1, 1, 1e6), ("Af1tB", "A", 1, 1, 3, 2e6), ("AfBt1", "A", 1, 3, 1, 3e6),
                ("Bf1t1", "B", 1, 1, 1, 4e6), ("BfBt1", "B", 1, 3, 1, 5e6)
            });
            var contiguous = RunOnePeptideSearch("TestFractionTechrepGapsContiguous", new[]
            {
                ("Af1t1", "A", 1, 1, 1, 1e6), ("Af1tB", "A", 1, 1, 2, 2e6), ("AfBt1", "A", 1, 2, 1, 3e6),
                ("Bf1t1", "B", 1, 1, 1, 4e6), ("BfBt1", "B", 1, 2, 1, 5e6)
            });

            Assert.That(gapped.Warnings, Has.None.Contain("Skipping quantification"));
            foreach (string expected in new[]
            {
                "Condition \"A\" biorep 1: fractions 1, 3, quantified as numbered; fraction 2 is not in the design.",
                "Condition \"A\" biorep 1 fraction 1: technical replicates 1, 3, quantified as numbered; technical replicate 2 is not in the design.",
                "Condition \"B\" biorep 1: fractions 1, 3, quantified as numbered; fraction 2 is not in the design."
            })
            {
                Assert.That(gapped.Warnings, Has.Some.StartWith(expected));
                Assert.That(gapped.ResultsTxt, Does.Contain(expected));
            }
            Assert.That(contiguous.Warnings, Has.None.Contain("not in the design"));

            Assert.That(gapped.PeptideIntensityByFile.Keys, Is.EquivalentTo(contiguous.PeptideIntensityByFile.Keys));
            Assert.That(gapped.PeptideIntensityByFile.Values, Has.All.GreaterThan(0));
            foreach (var (file, intensity) in contiguous.PeptideIntensityByFile)
            {
                Assert.That(gapped.PeptideIntensityByFile[file], Is.EqualTo(intensity).Within(1e-9).Percent, file);
            }
        }

        /// <summary>
        /// Runs a search with normalization on over one-peptide files written at the given intensities,
        /// under a design with the given one-based numbers, and returns what it warned, results.txt, the
        /// protein table's header and the peptide table's intensity per file.
        /// </summary>
        private static (List<string> Warnings, string ResultsTxt, string[] ProteinHeader, Dictionary<string, double> PeptideIntensityByFile)
            RunOnePeptideSearch(string folderName, (string Name, string Condition, int Biorep, int Fraction, int Techrep, double Intensity)[] files)
        {
            string unitTestFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, folderName);
            if (Directory.Exists(unitTestFolder))
            {
                Directory.Delete(unitTestFolder, true);
            }
            _ = Directory.CreateDirectory(unitTestFolder);

            string peptide = "PEPTIDE";
            string dbName = Path.Combine(unitTestFolder, "testDB.fasta");
            UsefulProteomicsDatabases.ProteinDbWriter.WriteFastaDatabase(new List<Protein> { new(peptide, @"test") }, dbName, ">");

            var fileInfos = new List<SpectraFileInfo>();
            foreach (var f in files)
            {
                string fullPath = Path.Combine(unitTestFolder, f.Name + ".mzML");
                WriteOnePeptideMzml(peptide, f.Intensity, fullPath);
                fileInfos.Add(new SpectraFileInfo(fullPath, f.Condition, f.Biorep - 1, f.Techrep - 1, f.Fraction - 1));
            }
            _ = ExperimentalDesign.WriteExperimentalDesignToFile(fileInfos);

            var warnings = new List<string>();
            EventHandler<StringEventArgs> handler = (o, e) => warnings.Add(e.S);
            MetaMorpheusTask.WarnHandler += handler;
            try
            {
                SearchTask task = new SearchTask();
                task.SearchParameters.Normalize = true;
                task.RunTask(unitTestFolder, new List<DbForTask> { new DbForTask(dbName, false) }, fileInfos.Select(p => p.FullFilePathWithExtension).ToList(), "");
            }
            finally
            {
                MetaMorpheusTask.WarnHandler -= handler;
            }

            string resultsTxt = File.ReadAllText(Path.Combine(unitTestFolder, "results.txt"));
            string[] proteinHeader = File.ReadAllLines(Path.Combine(unitTestFolder, "AllQuantifiedProteinGroups.tsv"))[0].Split('\t');

            string[] peptideLines = File.ReadAllLines(Path.Combine(unitTestFolder, "AllQuantifiedPeptides.tsv"));
            var peptideHeader = peptideLines[0].Split('\t').ToList();
            var peptideIntensityByFile = files.ToDictionary(f => f.Name,
                f => double.Parse(peptideLines[1].Split('\t')[peptideHeader.IndexOf("Intensity_" + f.Name)]));

            Directory.Delete(unitTestFolder, true);
            return (warnings, resultsTxt, proteinHeader, peptideIntensityByFile);
        }

        /// <summary>
        /// An mzML holding one MS1 scan of the peptide's isotopic envelope at the given intensity, and one
        /// MS2 scan of its b and y ions. The same scans TestProteinQuantFileHeaders writes.
        /// </summary>
        private static void WriteOnePeptideMzml(string peptide, double ionIntensity, string fullPath)
        {
            MsDataScan[] scans = new MsDataScan[2];

            ChemicalFormula cf = new Proteomics.AminoAcidPolymer.Peptide(peptide).GetChemicalFormula();
            IsotopicDistribution dist = IsotopicDistribution.GetDistribution(cf, 0.125, 1e-8);
            double[] mz = dist.Masses.Select(v => v.ToMz(1)).ToArray();
            double[] intensities = dist.Intensities.Select(v => v * ionIntensity).ToArray();

            scans[0] = new MsDataScan(massSpectrum: new MzSpectrum(mz, intensities, false), oneBasedScanNumber: 1, msnOrder: 1, isCentroid: true,
                polarity: Polarity.Positive, retentionTime: 1.0, scanWindowRange: new MzRange(400, 1600), scanFilter: "f",
                mzAnalyzer: MZAnalyzerType.Orbitrap, totalIonCurrent: intensities.Sum(), injectionTime: 1.0, noiseData: null, nativeId: "scan=1");

            var pep = new PeptideWithSetModifications(peptide, new Dictionary<string, Modification>());
            List<Product> frags = new List<Product>();
            pep.Fragment(DissociationType.HCD, FragmentationTerminus.Both, frags);
            double[] mz2 = frags.Select(v => v.NeutralMass.ToMz(1)).ToArray();
            double[] intensities2 = frags.Select(v => 1e6).ToArray();

            scans[1] = new MsDataScan(massSpectrum: new MzSpectrum(mz2, intensities2, false), oneBasedScanNumber: 2, msnOrder: 2, isCentroid: true,
                polarity: Polarity.Positive, retentionTime: 1.01, scanWindowRange: new MzRange(100, 1600), scanFilter: "f",
                mzAnalyzer: MZAnalyzerType.Orbitrap, totalIonCurrent: intensities.Sum(), injectionTime: 1.0, noiseData: null, nativeId: "scan=2", selectedIonMz: pep.MonoisotopicMass.ToMz(1),
                selectedIonChargeStateGuess: 1, selectedIonIntensity: 1e6, isolationMZ: pep.MonoisotopicMass.ToMz(1), isolationWidth: 1.5, dissociationType: DissociationType.HCD,
                oneBasedPrecursorScanNumber: 1, selectedIonMonoisotopicGuessMz: pep.MonoisotopicMass.ToMz(1), hcdEnergy: "35");

            Readers.MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(
                new GenericMsDataFile(scans, new SourceFile(@"scan number only nativeID format", "mzML format", null, "SHA-1", @"C:\fake.mzML", null)),
                fullPath, false);
        }

        [Test]
        [TestCase(false, 2, 1, 1)]
        [TestCase(true, 2, 3, 1)]
        [TestCase(true, 2, 3, 2)]
        public static void TestProteinQuantFileHeaders(bool hasDefinedExperimentalDesign, int bioreps, int fractions, int techreps)
        {
            // create the unit test directory
            string unitTestFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestProteinQuantFileHeaders");
            _ = Directory.CreateDirectory(unitTestFolder);

            List<SpectraFileInfo> fileInfos = new();
            string peptide = "PEPTIDE";
            double ionIntensity = 1e6;
            string condition = hasDefinedExperimentalDesign ? "TestCondition" : "";

            // create the protein database
            Protein prot = new(peptide, @"test"); // necessary to pass name to protein. otherwise dbloader will do crazy things
            string dbName = Path.Combine(unitTestFolder, "testDB.fasta");
            UsefulProteomicsDatabases.ProteinDbWriter.WriteFastaDatabase(new List<Protein> { prot }, dbName, ">");

            // create the .mzML files to search/quantify
            for (int b = 0; b < bioreps; b++)
            {
                for (int f = 0; f < fractions; f++)
                {
                    for (int r = 0; r < techreps; r++)
                    {
                        string fileToWrite = "file_" + "b" + b + "f" + f + "r" + r + ".mzML";

                        // generate mzml file
                        MsDataScan[] scans = new MsDataScan[2];

                        // create the MS1 scan
                        ChemicalFormula cf = new Proteomics.AminoAcidPolymer.Peptide(peptide).GetChemicalFormula();
                        IsotopicDistribution dist = IsotopicDistribution.GetDistribution(cf, 0.125, 1e-8);
                        double[] mz = dist.Masses.Select(v => v.ToMz(1)).ToArray();
                        double[] intensities = dist.Intensities.Select(v => v * ionIntensity * (b + 1)).ToArray();

                        scans[0] = new MsDataScan(massSpectrum: new MzSpectrum(mz, intensities, false), oneBasedScanNumber: 1, msnOrder: 1, isCentroid: true,
                            polarity: Polarity.Positive, retentionTime: 1.0, scanWindowRange: new MzRange(400, 1600), scanFilter: "f",
                            mzAnalyzer: MZAnalyzerType.Orbitrap, totalIonCurrent: intensities.Sum(), injectionTime: 1.0, noiseData: null, nativeId: "scan=1");

                        // create the MS2 scan
                        var pep = new PeptideWithSetModifications(peptide, new Dictionary<string, Modification>());
                        List<Product> frags = new List<Product>();
                        pep.Fragment(DissociationType.HCD, FragmentationTerminus.Both, frags);
                        double[] mz2 = frags.Select(v => v.NeutralMass.ToMz(1)).ToArray();
                        double[] intensities2 = frags.Select(v => 1e6).ToArray();

                        scans[1] = new MsDataScan(massSpectrum: new MzSpectrum(mz2, intensities2, false), oneBasedScanNumber: 2, msnOrder: 2, isCentroid: true,
                            polarity: Polarity.Positive, retentionTime: 1.01, scanWindowRange: new MzRange(100, 1600), scanFilter: "f",
                            mzAnalyzer: MZAnalyzerType.Orbitrap, totalIonCurrent: intensities.Sum(), injectionTime: 1.0, noiseData: null, nativeId: "scan=2", selectedIonMz: pep.MonoisotopicMass.ToMz(1),
                            selectedIonChargeStateGuess: 1, selectedIonIntensity: 1e6, isolationMZ: pep.MonoisotopicMass.ToMz(1), isolationWidth: 1.5, dissociationType: DissociationType.HCD,
                            oneBasedPrecursorScanNumber: 1, selectedIonMonoisotopicGuessMz: pep.MonoisotopicMass.ToMz(1), hcdEnergy: "35");

                        // write the .mzML
                        string fullPath = Path.Combine(unitTestFolder, fileToWrite);
                        Readers.MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(
                            new GenericMsDataFile(scans, new SourceFile(@"scan number only nativeID format", "mzML format", null, "SHA-1", @"C:\fake.mzML", null)),
                            fullPath, false);

                        SpectraFileInfo spectraFileInfo = new(fullPath, condition, b, r, f);
                        fileInfos.Add(spectraFileInfo);
                    }
                }
            }

            // write the experimental design for this quantification test
            if (hasDefinedExperimentalDesign)
            {
                _ = ExperimentalDesign.WriteExperimentalDesignToFile(fileInfos);
            }

            // run the search/quantification
            SearchTask task = new SearchTask();
            task.RunTask(unitTestFolder, new List<DbForTask> { new DbForTask(dbName, false) }, fileInfos.Select(p => p.FullFilePathWithExtension).ToList(), "");

            // read in the protein quant results
            Assert.That(File.Exists(Path.Combine(unitTestFolder, "AllQuantifiedProteinGroups.tsv")));
            string[] lines = File.ReadAllLines(Path.Combine(unitTestFolder, "AllQuantifiedProteinGroups.tsv"));

            // check the intensity column headers
            List<string> splitHeader = lines[0].Split(new char[] { '\t' }).ToList();
            List<string> intensityColumnHeaders = splitHeader.Where(p => p.Contains("Intensity_", StringComparison.OrdinalIgnoreCase)).ToList();

            Assert.That(intensityColumnHeaders.Count == 2);

            if (!hasDefinedExperimentalDesign)
            {
                Assert.That(intensityColumnHeaders[0] == "Intensity_file_b0f0r0");
                Assert.That(intensityColumnHeaders[1] == "Intensity_file_b1f0r0");
            }
            else
            {
                Assert.That(intensityColumnHeaders[0] == "Intensity_TestCondition_1");
                Assert.That(intensityColumnHeaders[1] == "Intensity_TestCondition_2");
            }

            // check the protein intensity values
            int ind1 = splitHeader.IndexOf(intensityColumnHeaders[0]);
            int ind2 = splitHeader.IndexOf(intensityColumnHeaders[1]);
            double intensity1 = double.Parse(lines[1].Split(new char[] { '\t' })[ind1]);
            double intensity2 = double.Parse(lines[1].Split(new char[] { '\t' })[ind2]);

            Assert.That(intensity1 > 0);
            Assert.That(intensity2 > 0);
            Assert.That(intensity1 < intensity2);

            Directory.Delete(unitTestFolder, true);
        }
    }
}
