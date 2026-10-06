using EngineLayer;
using EngineLayer.DatabaseLoading;
using MassSpectrometry;
using NUnit.Framework;
using Readers;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using TaskLayer;

namespace Test
{
    /// <summary>
    /// A site a search saw modified in one file is reported as 0/N in a file that covered it without ever seeing it
    /// modified, rather than left out (the user's ruling on gap AB1; mzLib #1411). The second file is the first with
    /// its one oxidized spectrum removed, so it covers the oxidation site of P10591 only unmodified.
    /// </summary>
    [TestFixture]
    public static class OccupancyZeroOfNTests
    {
        [Test]
        public static void AFileThatCoversASiteModifiedOnlyInAnotherFileReportsZeroOfN()
        {
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, "OccupancyZeroOfN");
            if (Directory.Exists(outputFolder)) Directory.Delete(outputFolder, true);
            string inputFolder = Path.Combine(outputFolder, "inputs");
            Directory.CreateDirectory(inputFolder);
            try
            {
                string fasta = Path.Combine(inputFolder, "DbForPrunedDb.fasta");
                File.Copy(Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\DbForPrunedDb.fasta"), fasta, true);
                string modifiedFile = Path.Combine(inputFolder, "withMod.mzml");
                File.Copy(Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestData\PrunedDbSpectra.mzml"), modifiedFile, true);

                // find the oxidized spectrum with one search, then write the same file without it
                var firstPass = RunGptmdThenSearch(Path.Combine(outputFolder, "firstPass"), fasta, new List<string> { modifiedFile });
                var oxidizedScans = File.ReadAllLines(Path.Combine(firstPass, "AllPSMs.psmtsv")).Skip(1)
                    .Select(line => line.Split('\t'))
                    .Where(fields => fields[FullSequenceIndex(firstPass)].Contains("Oxidation on S"))
                    .Select(fields => int.Parse(fields[ScanNumberIndex(firstPass)]))
                    .ToHashSet();
                Assert.That(oxidizedScans, Is.Not.Empty, "premise: the first file has an oxidized PSM");

                // the oxidized scans keep their place (mzML precursor references are by position) but lose their
                // fragments, so nothing can be identified in them
                string unmodifiedFile = Path.Combine(inputFolder, "withoutMod.mzml");
                var source = MsDataFileReader.GetDataFile(modifiedFile).LoadAllStaticData();
                var scans = source.GetAllScansList().Select(scan => !oxidizedScans.Contains(scan.OneBasedScanNumber) ? scan
                    : new MsDataScan(new MzSpectrum(new[] { 100.0 }, new[] { 1.0 }, false), scan.OneBasedScanNumber, scan.MsnOrder,
                        scan.IsCentroid, scan.Polarity, scan.RetentionTime, scan.ScanWindowRange, scan.ScanFilter, scan.MzAnalyzer,
                        1.0, scan.InjectionTime, scan.NoiseData, scan.NativeId, scan.SelectedIonMZ, scan.SelectedIonChargeStateGuess,
                        scan.SelectedIonIntensity, scan.IsolationMz, scan.IsolationWidth, scan.DissociationType,
                        scan.OneBasedPrecursorScanNumber, scan.SelectedIonMonoisotopicGuessMz, scan.HcdEnergy))
                    .ToArray();
                MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(new GenericMsDataFile(scans, source.SourceFile), unmodifiedFile, false);

                File.WriteAllLines(Path.Combine(inputFolder, GlobalVariables.ExperimentalDesignFileName), new[]
                {
                    "FileName\tCondition\tBiorep\tFraction\tTechrep",
                    "withMod.mzml\tA\t1\t1\t1",
                    "withoutMod.mzml\tA\t2\t1\t1",
                });

                string search = RunGptmdThenSearch(Path.Combine(outputFolder, "both"), fasta, new List<string> { modifiedFile, unmodifiedFile });

                string withMod = CountOccupancy(Path.Combine(search, "AllQuantifiedProteinGroups.tsv"), "CountOccupancy_A_1");
                string withoutMod = CountOccupancy(Path.Combine(search, "AllQuantifiedProteinGroups.tsv"), "CountOccupancy_A_2");
                Assert.That(withMod, Does.Contain("Oxidation on S").And.Not.Contain("(0/"), "premise: the first file sees the site modified");
                Assert.That(withoutMod, Does.Match(@"pos\d+\[Oxidation on S,info:fraction=0\.00\(0/\d+\)\]"),
                    "the second file covers the site unmodified, so it reports 0/N instead of nothing");

                // the per-file table of the second file says the same: its group inherits the search-wide sites
                string individual = Path.Combine(search, "Individual File Results", "withoutMod_" + GlobalVariables.AnalyteType.GetBioPolymerLabel() + "Groups.tsv");
                Assert.That(File.Exists(individual), "premise: the per-file protein table is written");
                string perFile = CountOccupancy(individual, File.ReadLines(individual).First().Split('\t').First(h => h.StartsWith("CountOccupancy_")));
                Assert.That(perFile, Does.Match(@"pos\d+\[Oxidation on S,info:fraction=0\.00\(0/\d+\)\]"));
            }
            finally
            {
                if (Directory.Exists(outputFolder)) Directory.Delete(outputFolder, true);
            }
        }

        private static string RunGptmdThenSearch(string outputFolder, string fasta, List<string> spectra)
        {
            var gptmd = new GptmdTask
            {
                CommonParameters = new CommonParameters(),
                GptmdParameters = new GptmdParameters
                {
                    ListOfModsGptmd = GlobalVariables.AllModsKnown.Where(b =>
                        b.ModificationType.Equals("Common Artifact") || b.ModificationType.Equals("Common Biological")
                        || b.ModificationType.Equals("Metal") || b.ModificationType.Equals("Less Common"))
                        .Select(b => (b.ModificationType, b.IdWithMotif)).ToList()
                }
            };
            var search = new SearchTask
            {
                CommonParameters = new CommonParameters(),
                SearchParameters = new SearchParameters { DoParsimony = true, SearchTarget = true, SearchType = SearchType.Classic }
            };
            new EverythingRunnerEngine(new List<(string, MetaMorpheusTask)> { ("gptmd", gptmd), ("search", search) },
                spectra, new List<DbForTask> { new DbForTask(fasta, false) }, outputFolder).Run();
            return Path.Combine(outputFolder, "search");
        }

        private static int FullSequenceIndex(string folder) => Header(Path.Combine(folder, "AllPSMs.psmtsv")).IndexOf("Full Sequence");
        private static int ScanNumberIndex(string folder) => Header(Path.Combine(folder, "AllPSMs.psmtsv")).IndexOf("Scan Number");
        private static List<string> Header(string path) => File.ReadLines(path).First().Split('\t').ToList();

        /// <summary>The named occupancy column of P10591's row.</summary>
        private static string CountOccupancy(string proteinTable, string column)
        {
            int index = Header(proteinTable).IndexOf(column);
            Assert.That(index, Is.GreaterThanOrEqualTo(0), "no column " + column + " in " + proteinTable);
            return File.ReadLines(proteinTable).Skip(1).First(line => line.StartsWith("P10591")).Split('\t')[index];
        }
    }
}
