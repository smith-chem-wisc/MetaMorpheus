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

            // test non-consecutive fractions (should produce an error)
            spectraFiles.Clear();

            spectraFiles.Add(new SpectraFileInfo(Path.Combine(outputFolder, @"myFile.1.raw"), "condition1", 0, 0, 1));
            spectraFiles.Add(new SpectraFileInfo(Path.Combine(outputFolder, @"myFile.2.raw"), "condition1", 0, 0, 2));
            spectraFiles.Add(new SpectraFileInfo(Path.Combine(outputFolder, @"myFile.3.raw"), "condition1", 0, 0, 3));

            ExperimentalDesign.WriteExperimentalDesignToFile(spectraFiles);

            readIn = ExperimentalDesign.ReadExperimentalDesign(
                Path.Combine(outputFolder, @"ExperimentalDesign.tsv"),
                spectraFiles.Select(p => p.FullFilePathWithExtension).ToList(),
                out errors);

            Assert.That(errors.Any());

            // test non-consecutive techreps (should produce an error)
            spectraFiles.Clear();

            spectraFiles.Add(new SpectraFileInfo(Path.Combine(outputFolder, @"myFile.1.raw"), "condition1", 0, 0, 0));
            spectraFiles.Add(new SpectraFileInfo(Path.Combine(outputFolder, @"myFile.2.raw"), "condition1", 0, 2, 0));
            spectraFiles.Add(new SpectraFileInfo(Path.Combine(outputFolder, @"myFile.3.raw"), "condition1", 0, 3, 0));

            ExperimentalDesign.WriteExperimentalDesignToFile(spectraFiles);

            readIn = ExperimentalDesign.ReadExperimentalDesign(
                Path.Combine(outputFolder, @"ExperimentalDesign.tsv"),
                spectraFiles.Select(p => p.FullFilePathWithExtension).ToList(),
                out errors);

            Assert.That(errors.Any());

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
                "Condition \"A\" has biorep 1, 2, 4 but not biorep 3. The bioreps are quantified as numbered; " +
                "a missing number may be a sample that was lost or not searched."
            }));

            // study-wide numbering: each condition is judged on its own numbers
            var studyWide = new List<SpectraFileInfo> { Info("c1", "control", 1), Info("c2", "control", 2), Info("t3", "treated", 3), Info("t4", "treated", 4) };
            Assert.That(ExperimentalDesign.GetErrorsInExperimentalDesign(studyWide), Is.Null);
            var warnings = ExperimentalDesign.GetWarningsInExperimentalDesign(studyWide);
            Assert.That(warnings, Has.Count.EqualTo(1));
            Assert.That(warnings[0], Does.StartWith("Condition \"treated\" has biorep 3, 4 but not biorep 1, 2."));

            // 1..N says nothing
            var complete = new List<SpectraFileInfo> { Info("a1", "A", 1), Info("a2", "A", 2) };
            Assert.That(ExperimentalDesign.GetWarningsInExperimentalDesign(complete), Is.Empty);

            // a biorep after a gap still has its fractions checked
            var fractionGap = new List<SpectraFileInfo> { Info("a1", "A", 1), Info("a4f2", "A", 4, oneBasedFraction: 2) };
            Assert.That(ExperimentalDesign.GetErrorsInExperimentalDesign(fractionGap), Is.EqualTo("Condition \"A\" biorep 4 fraction 1 is missing!"));
        }

        /// <summary>
        /// End to end, with normalization on: a design numbered study-wide (A: 2, 3; B: 4, 5) used to make
        /// the search skip quantification with one warning. It is now quantified under the numbers given,
        /// and the gaps are warned about. The first condition has no biorep 1, which is the normalization
        /// reference before mzLib #1422; without #1422 biorep normalization silently does nothing here.
        /// </summary>
        [Test]
        public static void TestSearchQuantifiesADesignWithABiorepGap()
        {
            string unitTestFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestSearchQuantifiesADesignWithABiorepGap");
            if (Directory.Exists(unitTestFolder))
            {
                Directory.Delete(unitTestFolder, true);
            }
            _ = Directory.CreateDirectory(unitTestFolder);

            string peptide = "PEPTIDE";
            Protein prot = new(peptide, @"test");
            string dbName = Path.Combine(unitTestFolder, "testDB.fasta");
            UsefulProteomicsDatabases.ProteinDbWriter.WriteFastaDatabase(new List<Protein> { prot }, dbName, ">");

            var fileInfos = new List<SpectraFileInfo>();
            foreach (var (condition, oneBasedBiorep) in new[] { ("A", 2), ("A", 3), ("B", 4), ("B", 5) })
            {
                string fullPath = Path.Combine(unitTestFolder, condition + oneBasedBiorep + ".mzML");
                WriteOnePeptideMzml(peptide, 1e6 * oneBasedBiorep, fullPath);
                fileInfos.Add(new SpectraFileInfo(fullPath, condition, oneBasedBiorep - 1, 0, 0));
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

            Assert.That(warnings, Has.None.Contain("Skipping quantification"));
            Assert.That(warnings, Has.Some.StartWith("Condition \"A\" has biorep 2, 3 but not biorep 1."));
            Assert.That(warnings, Has.Some.StartWith("Condition \"B\" has biorep 4, 5 but not biorep 1, 2, 3."));

            // the protein table's sample columns carry the numbers the design gives
            string proteinTable = Path.Combine(unitTestFolder, "AllQuantifiedProteinGroups.tsv");
            Assert.That(File.Exists(proteinTable));
            var proteinHeader = File.ReadAllLines(proteinTable)[0].Split('\t');
            Assert.That(proteinHeader.Where(h => h.StartsWith("Intensity_")),
                Is.EquivalentTo(new[] { "Intensity_A_2", "Intensity_A_3", "Intensity_B_4", "Intensity_B_5" }));

            // and every file's peptide was quantified and normalized: the files were written at 2, 3, 4
            // and 5 times one intensity, and normalization brings this one peptide back to one value
            string[] lines = File.ReadAllLines(Path.Combine(unitTestFolder, "AllQuantifiedPeptides.tsv"));
            var header = lines[0].Split('\t').ToList();
            var intensities = new[] { "Intensity_A2", "Intensity_A3", "Intensity_B4", "Intensity_B5" }
                .Select(column => double.Parse(lines[1].Split('\t')[header.IndexOf(column)]))
                .ToList();
            Assert.That(intensities, Has.All.GreaterThan(0));
            Assert.That(intensities.Max() / intensities.Min(), Is.LessThan(1.01));

            Directory.Delete(unitTestFolder, true);
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
