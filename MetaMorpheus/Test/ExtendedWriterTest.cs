using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using FlashLFQ;
using Nett;
using NUnit.Framework;
using TaskLayer;

namespace Test
{
    [TestFixture]
    public class ExtendedWriterTest
    {
        [Test]
        public static void TestExtendedWriter()
        {
            SpectraFileInfo fraction1 = new SpectraFileInfo("", "", 0, 0, fraction: 0);
            SpectraFileInfo fraction2 = new SpectraFileInfo("", "", 0, 0, fraction: 1);

            Identification id1 = new Identification(fraction1, "peptide1", "peptide1", 0, 0, 0, new List<ProteinGroup>());
            Identification id2 = new Identification(fraction2, "peptide1", "peptide1", 0, 0, 0, new List<ProteinGroup>());

            ChromatographicPeak peak1 = new ChromatographicPeak(id1, false, fraction1);
            ChromatographicPeak peak2 = new ChromatographicPeak(id2, false, fraction1);
            List<ChromatographicPeak> peakList = new List<ChromatographicPeak> {peak1, peak2};

            List<(string, int)> newField = new List<(string, int)> { ("Spectral Contrast Angle", 16) };

            ExtendedWriter writer = new ExtendedWriter(peakList, ChromatographicPeak.TabSeparatedHeader, newField);
            string extendedHeader = writer.ExtendedHeader;

            string[] peakHeaderSplit = ChromatographicPeak.TabSeparatedHeader.Split('\t');
            StringBuilder sb = new();
            sb.Append(string.Join('\t', peakHeaderSplit[0..16]));
            sb.Append('\t');
            sb.Append("Spectral Contrast Angle");
            sb.Append('\t');
            sb.Append(string.Join('\t', peakHeaderSplit[16..]));
            string expectedHeader = sb.ToString().Trim();

            Assert.AreEqual(extendedHeader, expectedHeader);

            writer.AddInfo(peak1, "Spectral Contrast Angle", "1");
            writer.AddInfo(peak2, "Spectral Contrast Angle", "0");

            Assert.AreEqual(writer.WriteExtendedString(peak1).Split('\t')[16], "1");
            Assert.AreEqual(writer.WriteExtendedString(peak2).Split('\t')[16], "0");
        }

        [Test]
        public static void TestExtendedWriterAlternateConstructor()
        {
            SpectraFileInfo fraction1 = new SpectraFileInfo("", "", 0, 0, fraction: 0);
            SpectraFileInfo fraction2 = new SpectraFileInfo("", "", 0, 0, fraction: 1);

            Identification id1 = new Identification(fraction1, "peptide1", "peptide1", 0, 0, 0, new List<ProteinGroup>());
            Identification id2 = new Identification(fraction2, "peptide1", "peptide1", 0, 0, 0, new List<ProteinGroup>());

            ChromatographicPeak peak1 = new ChromatographicPeak(id1, false, fraction1);
            ChromatographicPeak peak2 = new ChromatographicPeak(id2, false, fraction1);
            List<ChromatographicPeak> peakList = new List<ChromatographicPeak> { peak1, peak2 };

            List<string> newField = new List<string> { "Spectral Contrast Angle" };

            ExtendedWriter writer = new ExtendedWriter(peakList, ChromatographicPeak.TabSeparatedHeader, newField, insertionIndex: 16);
            string extendedHeader = writer.ExtendedHeader;

            string[] peakHeaderSplit = ChromatographicPeak.TabSeparatedHeader.Split('\t');
            StringBuilder sb = new();
            sb.Append(string.Join('\t', peakHeaderSplit[0..16]));
            sb.Append('\t');
            sb.Append("Spectral Contrast Angle");
            sb.Append('\t');
            sb.Append(string.Join('\t', peakHeaderSplit[16..]));
            string expectedHeader = sb.ToString().Trim();

            Assert.AreEqual(extendedHeader, expectedHeader);

            writer.AddInfo(peak1, "Spectral Contrast Angle", "1");
            writer.AddInfo(peak2, "Spectral Contrast Angle", "0");

            Assert.AreEqual(writer.WriteExtendedString(peak1).Split('\t')[16], "1");
            Assert.AreEqual(writer.WriteExtendedString(peak2).Split('\t')[16], "0");
        }

        [Test]
        public static void TestWriteExtendedPeakQuant()
        {
            var task = Toml.ReadFile<SearchTask>(Path.Combine(TestContext.CurrentContext.TestDirectory, @"SlicedSearchTaskConfig.toml"),
                MetaMorpheusTask.tomlConfig);

            DbForTask db = new DbForTask(Path.Combine(TestContext.CurrentContext.TestDirectory, @"sliced-db.fasta"), false);
            string raw = Path.Combine(TestContext.CurrentContext.TestDirectory, @"sliced-raw.mzML");
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestSlicedTest1");
            EverythingRunnerEngine a = new EverythingRunnerEngine(
                new List<(string, MetaMorpheusTask)> { ("Task", task) },
                new List<string> { raw },
                new List<DbForTask> { db }, outputFolder);

            a.Run();
        }
    }
}
