using Chemistry;
using EngineLayer;
using FlashLFQ;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Proteomics;
using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
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

            List<FlashLFQ.SpectraFileInfo> spectraFiles = new List<FlashLFQ.SpectraFileInfo>();

            spectraFiles.Add(new FlashLFQ.SpectraFileInfo(Path.Combine(outputFolder, @"myFile.1.raw"), "condition1", 0, 0, 0));
            spectraFiles.Add(new FlashLFQ.SpectraFileInfo(Path.Combine(outputFolder, @"myFile.2.raw"), "condition1", 1, 0, 0));
            spectraFiles.Add(new FlashLFQ.SpectraFileInfo(Path.Combine(outputFolder, @"myFile.3.raw"), "condition1", 2, 0, 0));

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

            // test non-consecutive bioreps (should produce an error)
            List<FlashLFQ.SpectraFileInfo> spectraFiles = new List<FlashLFQ.SpectraFileInfo>();

            spectraFiles.Add(new FlashLFQ.SpectraFileInfo(Path.Combine(outputFolder, @"myFile.1.raw"), "condition1", 0, 0, 0));
            spectraFiles.Add(new FlashLFQ.SpectraFileInfo(Path.Combine(outputFolder, @"myFile.2.raw"), "condition1", 1, 0, 0));
            spectraFiles.Add(new FlashLFQ.SpectraFileInfo(Path.Combine(outputFolder, @"myFile.3.raw"), "condition1", 3, 0, 0));

            ExperimentalDesign.WriteExperimentalDesignToFile(spectraFiles);

            var readIn = ExperimentalDesign.ReadExperimentalDesign(
                Path.Combine(outputFolder, @"ExperimentalDesign.tsv"),
                spectraFiles.Select(p => p.FullFilePathWithExtension).ToList(),
                out var errors);

            Assert.That(errors.Any());

            // test non-consecutive fractions (should produce an error)
            spectraFiles.Clear();

            spectraFiles.Add(new FlashLFQ.SpectraFileInfo(Path.Combine(outputFolder, @"myFile.1.raw"), "condition1", 0, 0, 1));
            spectraFiles.Add(new FlashLFQ.SpectraFileInfo(Path.Combine(outputFolder, @"myFile.2.raw"), "condition1", 0, 0, 2));
            spectraFiles.Add(new FlashLFQ.SpectraFileInfo(Path.Combine(outputFolder, @"myFile.3.raw"), "condition1", 0, 0, 3));

            ExperimentalDesign.WriteExperimentalDesignToFile(spectraFiles);

            readIn = ExperimentalDesign.ReadExperimentalDesign(
                Path.Combine(outputFolder, @"ExperimentalDesign.tsv"),
                spectraFiles.Select(p => p.FullFilePathWithExtension).ToList(),
                out errors);

            Assert.That(errors.Any());

            // test non-consecutive techreps (should produce an error)
            spectraFiles.Clear();

            spectraFiles.Add(new FlashLFQ.SpectraFileInfo(Path.Combine(outputFolder, @"myFile.1.raw"), "condition1", 0, 0, 0));
            spectraFiles.Add(new FlashLFQ.SpectraFileInfo(Path.Combine(outputFolder, @"myFile.2.raw"), "condition1", 0, 2, 0));
            spectraFiles.Add(new FlashLFQ.SpectraFileInfo(Path.Combine(outputFolder, @"myFile.3.raw"), "condition1", 0, 3, 0));

            ExperimentalDesign.WriteExperimentalDesignToFile(spectraFiles);

            readIn = ExperimentalDesign.ReadExperimentalDesign(
                Path.Combine(outputFolder, @"ExperimentalDesign.tsv"),
                spectraFiles.Select(p => p.FullFilePathWithExtension).ToList(),
                out errors);

            Assert.That(errors.Any());

            // test duplicates (should produce an error)
            spectraFiles.Clear();

            spectraFiles.Add(new FlashLFQ.SpectraFileInfo(Path.Combine(outputFolder, @"myFile.1.raw"), "condition1", 0, 0, 0));
            spectraFiles.Add(new FlashLFQ.SpectraFileInfo(Path.Combine(outputFolder, @"myFile.2.raw"), "condition1", 0, 0, 0));

            ExperimentalDesign.WriteExperimentalDesignToFile(spectraFiles);

            readIn = ExperimentalDesign.ReadExperimentalDesign(
                Path.Combine(outputFolder, @"ExperimentalDesign.tsv"),
                spectraFiles.Select(p => p.FullFilePathWithExtension).ToList(),
                out errors);

            Assert.That(errors.Any());

            // test situation where experimental design does not contain a file (should produce an error)
            spectraFiles.Clear();
            spectraFiles.Add(new FlashLFQ.SpectraFileInfo(Path.Combine(outputFolder, @"myFile.1.raw"), "condition1", 0, 0, 0));
            ExperimentalDesign.WriteExperimentalDesignToFile(spectraFiles);
            spectraFiles.Add(new FlashLFQ.SpectraFileInfo(Path.Combine(outputFolder, @"myFile.2.raw"), "condition1", 1, 0, 0));

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
            Protein prot = new(peptide, @"");
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
                        PeptideWithSetModifications pep = new(peptide, new Dictionary<string, Proteomics.Modification>());
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
                        IO.MzML.MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(
                            new MsDataFile(scans, new SourceFile(@"scan number only nativeID format", "mzML format", null, "SHA-1", @"C:\fake.mzML", null)),
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
            _ = task.RunTask(unitTestFolder, new List<DbForTask> { new DbForTask(dbName, false) }, fileInfos.Select(p => p.FullFilePathWithExtension).ToList(), "");

            // read in the protein quant results
            Assert.That(File.Exists(Path.Combine(unitTestFolder, "AllQuantifiedProteinGroups.tsv")));
            string[] lines = File.ReadAllLines(Path.Combine(unitTestFolder, "AllQuantifiedProteinGroups.tsv"));

            // check the intensity column headers
            List<string> splitHeader = lines[0].Split(new char[] { '\t' }).ToList();
            List<string> intensityColumnHeaders = splitHeader.Where(p => p.Contains("Intensity", StringComparison.OrdinalIgnoreCase)).ToList();

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
