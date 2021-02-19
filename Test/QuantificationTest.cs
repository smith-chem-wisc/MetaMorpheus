using EngineLayer;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;

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
    }
}
