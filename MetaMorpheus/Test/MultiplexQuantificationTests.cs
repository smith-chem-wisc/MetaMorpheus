using EngineLayer;
using EngineLayer.DatabaseLoading;
using MassSpectrometry;
using Nett;
using NUnit.Framework;
using Quantification;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using TaskLayer;

namespace Test
{
    /// <summary>
    /// Covers the handoff from MetaMorpheus to mzLib's <see cref="QuantificationEngine"/> for isobaric
    /// (TMT/iTRAQ) searches: a TMT search that finds a TmtDesign.txt beside its spectra files must roll
    /// reporter ion intensities up to peptides and protein groups, not merely report them per PSM.
    /// </summary>
    [TestFixture]
    public static class MultiplexQuantificationTests
    {
        private static readonly string[] Tmt11Channels =
            { "126", "127N", "127C", "128N", "128C", "129N", "129C", "130N", "130C", "131N", "131C" };

        /// <summary>
        /// Writes a TmtDesign.txt beside <paramref name="stagedMzmlPath"/> annotating every channel of
        /// the TMT11 plex as a study sample, alternating between two conditions.
        /// </summary>
        private static void WriteTmtDesign(string dataFolder, string stagedMzmlPath)
        {
            var rows = Tmt11Channels.Select((tag, i) =>
                $"{stagedMzmlPath}\tPlex1\tSample{i + 1}\t{tag}\tCond{(i % 2 == 0 ? "A" : "B")}\t{i / 2 + 1}\t1\t1\tstudy sample");

            File.WriteAllLines(
                Path.Combine(dataFolder, GlobalVariables.TmtExperimentalDesignFileName),
                new[] { TmtExperimentalDesign.Header }.Concat(rows));
        }

        private static void RunTmtSearch(string stagedMzmlPath, string outputFolder, bool doParsimony = true)
        {
            var searchTask = Toml.ReadFile<SearchTask>(
                Path.Combine(TestContext.CurrentContext.TestDirectory, @"TMT_test\TMT-Task1-SearchTaskconfig.toml"),
                MetaMorpheusTask.tomlConfig);

            // The stock TMT toml leaves parsimony off, and multiplex quantification rolls up to protein
            // groups, so there would be nothing to quantify onto. Turning it on is what a user doing
            // protein-level TMT would do.
            searchTask.SearchParameters.DoParsimony = doParsimony;

            var taskList = new List<(string, MetaMorpheusTask)> { ("search", searchTask) };
            string fasta = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TMT_test\mouseTmt.fasta");

            new EverythingRunnerEngine(
                taskList,
                new List<string> { stagedMzmlPath },
                new List<DbForTask> { new DbForTask(fasta, false) },
                outputFolder).Run();
        }

        /// <summary>
        /// Stages the TMT mzML in a private folder so a design file can be written beside it without
        /// touching the shared TMT_test data directory.
        /// </summary>
        private static string StageMzml(string dataFolder)
        {
            Directory.CreateDirectory(dataFolder);
            string staged = Path.Combine(dataFolder, "VA084TQ_6.mzML");
            File.Copy(Path.Combine(TestContext.CurrentContext.TestDirectory, @"TMT_test\VA084TQ_6.mzML"), staged);
            return staged;
        }

        /// <summary>
        /// The step this whole pipeline existed to reach: a TMT search now writes per-channel peptide and
        /// protein tables, instead of stopping at reporter ion columns in the .psmtsv.
        /// </summary>
        [Test]
        public static void TmtSearchWithDesignFile_WritesPerChannelPeptideAndProteinTables()
        {
            string root = Path.Combine(TestContext.CurrentContext.TestDirectory, "TmtMultiplexQuant");
            string dataFolder = Path.Combine(root, "data");
            string outputFolder = Path.Combine(root, "out");
            if (Directory.Exists(root)) Directory.Delete(root, true);

            try
            {
                string stagedMzml = StageMzml(dataFolder);
                WriteTmtDesign(dataFolder, stagedMzml);

                RunTmtSearch(stagedMzml, outputFolder);

                string searchOut = Path.Combine(outputFolder, "search");

                Assert.That(File.Exists(Path.Combine(searchOut, QuantificationWriter.RawFileName)), Is.True,
                    "the reloadable PSM-level snapshot must be written");
                Assert.That(File.Exists(Path.Combine(searchOut, QuantificationWriter.PeptideFileName)), Is.True);

                string proteinPath = Path.Combine(searchOut, QuantificationWriter.ProteinGroupFileName);
                Assert.That(File.Exists(proteinPath), Is.True, "the per-channel protein table must be written");

                var proteinLines = File.ReadAllLines(proteinPath);
                Assert.That(proteinLines, Has.Length.GreaterThan(1), "at least one protein group must be quantified");

                var header = proteinLines[0].Split('\t');
                Assert.That(header[0], Is.EqualTo("Protein Group"));
                Assert.That(header.Length, Is.EqualTo(1 + Tmt11Channels.Length),
                    "one column per TMT11 channel, even for channels the design left unannotated");
                Assert.That(header.Skip(1).Distinct().Count(), Is.EqualTo(Tmt11Channels.Length),
                    "sample columns must be distinguishable");

                var firstRow = proteinLines[1].Split('\t').Skip(1).Select(double.Parse).ToList();
                Assert.That(firstRow, Has.Count.EqualTo(Tmt11Channels.Length));
                Assert.That(firstRow.Any(v => v > 0), Is.True, "a quantified protein group must carry a value");

                // The design the results were produced under travels with them.
                Assert.That(File.Exists(Path.Combine(searchOut, GlobalVariables.TmtExperimentalDesignFileName)), Is.True);

                AssertPeptideValuesAreTheSumOfTheirRawValues(searchOut);
            }
            finally
            {
                if (Directory.Exists(root)) Directory.Delete(root, true);
            }
        }

        /// <summary>
        /// Checks the roll-up arithmetic AND the channel alignment that the whole handoff rests on: each
        /// peptide's value in a given column must be the sum of that peptide's PSM values in the matching
        /// reporter column of the raw snapshot.
        /// </summary>
        /// <remarks>
        /// A positional comparison is valid here because the run has one spectra file, so the peptide
        /// matrix's columns are that file's channels in ascending reporter m/z -- the same order
        /// <c>Reporter_1..Reporter_n</c> are written in. If MetaMorpheus's intensity array and the
        /// design's ISampleInfo array ever stopped agreeing on that order, the totals would still look
        /// plausible but would land in the wrong columns, and this is what would catch it.
        /// </remarks>
        private static void AssertPeptideValuesAreTheSumOfTheirRawValues(string searchOut)
        {
            var rawLines = File.ReadAllLines(Path.Combine(searchOut, QuantificationWriter.RawFileName));
            var rawHeader = rawLines[0].Split('\t');
            int fullSequenceColumn = System.Array.IndexOf(rawHeader, "FullSequence");
            int firstReporterColumn = System.Array.IndexOf(rawHeader, "Reporter_1");
            Assert.That(firstReporterColumn, Is.GreaterThan(0), "the raw file must carry reporter columns");

            var expectedByPeptide = new Dictionary<string, double[]>();
            foreach (var line in rawLines.Skip(1))
            {
                var cells = line.Split('\t');
                string sequence = cells[fullSequenceColumn];
                if (!expectedByPeptide.TryGetValue(sequence, out var totals))
                {
                    totals = new double[Tmt11Channels.Length];
                    expectedByPeptide[sequence] = totals;
                }
                for (int channel = 0; channel < Tmt11Channels.Length; channel++)
                {
                    totals[channel] += double.Parse(cells[firstReporterColumn + channel]);
                }
            }

            var peptideLines = File.ReadAllLines(Path.Combine(searchOut, QuantificationWriter.PeptideFileName));
            Assert.That(peptideLines, Has.Length.GreaterThan(1), "at least one peptide must be quantified");
            Assert.That(peptideLines[0].Split('\t'), Has.Length.EqualTo(1 + Tmt11Channels.Length));

            int compared = 0;
            foreach (var line in peptideLines.Skip(1))
            {
                var cells = line.Split('\t');
                if (!expectedByPeptide.TryGetValue(cells[0], out var expected)) continue;

                for (int channel = 0; channel < Tmt11Channels.Length; channel++)
                {
                    Assert.That(double.Parse(cells[1 + channel]), Is.EqualTo(expected[channel]).Within(1e-6),
                        $"peptide {cells[0]}, channel {Tmt11Channels[channel]}");
                }
                compared++;
            }

            Assert.That(compared, Is.GreaterThan(0), "no peptide in the peptide table was found in the raw table");
        }

        /// <summary>
        /// Quantification is additive, not a replacement: the reporter ion columns that were the only
        /// TMT output before must still be in the .psmtsv, and the protein table's channel columns must
        /// agree with them in number.
        /// </summary>
        [Test]
        public static void TmtSearchWithDesignFile_LeavesThePsmtsvReporterColumnsIntact()
        {
            string root = Path.Combine(TestContext.CurrentContext.TestDirectory, "TmtMultiplexQuantPsmtsv");
            string dataFolder = Path.Combine(root, "data");
            string outputFolder = Path.Combine(root, "out");
            if (Directory.Exists(root)) Directory.Delete(root, true);

            try
            {
                string stagedMzml = StageMzml(dataFolder);
                WriteTmtDesign(dataFolder, stagedMzml);

                RunTmtSearch(stagedMzml, outputFolder);

                string searchOut = Path.Combine(outputFolder, "search");
                var psmtsvHeader = File.ReadLines(Path.Combine(searchOut, "AllPeptides.psmtsv")).First().Trim().Split('\t');

                Assert.That(psmtsvHeader[^Tmt11Channels.Length..], Is.EquivalentTo(Tmt11Channels),
                    "the per-PSM reporter ion columns must survive");
            }
            finally
            {
                if (Directory.Exists(root)) Directory.Delete(root, true);
            }
        }

        /// <summary>
        /// A TMT search with no design file is still a legitimate search — there is simply no
        /// channel-to-sample mapping to quantify against. It must complete and still write the
        /// reporter ion columns, without any quantification files appearing.
        /// </summary>
        [Test]
        public static void TmtSearchWithoutDesignFile_StillSucceedsAndWritesNoQuantTables()
        {
            string root = Path.Combine(TestContext.CurrentContext.TestDirectory, "TmtMultiplexNoDesign");
            string dataFolder = Path.Combine(root, "data");
            string outputFolder = Path.Combine(root, "out");
            if (Directory.Exists(root)) Directory.Delete(root, true);

            try
            {
                string stagedMzml = StageMzml(dataFolder);   // deliberately no TmtDesign.txt beside it

                RunTmtSearch(stagedMzml, outputFolder);

                string searchOut = Path.Combine(outputFolder, "search");
                Assert.That(File.Exists(Path.Combine(searchOut, "AllPeptides.psmtsv")), Is.True,
                    "the search itself must still complete");
                Assert.That(File.Exists(Path.Combine(searchOut, QuantificationWriter.ProteinGroupFileName)), Is.False);
                Assert.That(File.Exists(Path.Combine(searchOut, QuantificationWriter.RawFileName)), Is.False);
            }
            finally
            {
                if (Directory.Exists(root)) Directory.Delete(root, true);
            }
        }
    }
}
