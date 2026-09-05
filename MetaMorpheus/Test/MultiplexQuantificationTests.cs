using EngineLayer;
using EngineLayer.DatabaseLoading;
using MassSpectrometry;
using Nett;
using NUnit.Framework;
using Omics;
using Omics.Fragmentation;
using Omics.Modifications;
using Omics.SpectralMatch;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using Quantification;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Reflection;
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

        /// <summary>The staged spectra file's name without extension, which prefixes every sample column.</summary>
        private const string StagedMzmlName = "VA084TQ_6";

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

                // By identity, not by count and distinctness. Counting and de-duplicating passes just as
                // happily on a design rotated by one channel, which is the desync worth catching: the
                // values would still look plausible and would be under the wrong sample.
                for (int channel = 0; channel < Tmt11Channels.Length; channel++)
                {
                    Assert.That(header[1 + channel], Is.EqualTo($"{StagedMzmlName}_{Tmt11Channels[channel]}"),
                        "column order must follow the plex's channel order, not merely be distinct");
                }

                var firstRow = proteinLines[1].Split('\t').Skip(1).Select(double.Parse).ToList();
                Assert.That(firstRow, Has.Count.EqualTo(Tmt11Channels.Length));
                Assert.That(firstRow.Any(v => v > 0), Is.True, "a quantified protein group must carry a value");

                // The design the results were produced under travels with them.
                Assert.That(File.Exists(Path.Combine(searchOut, GlobalVariables.TmtExperimentalDesignFileName)), Is.True);

                AssertPeptideValuesAreTheSumOfTheirRawValues(searchOut);
                AssertProteinTableKeepsItsColumnsAndGainsChannelIntensities(searchOut);
                AssertQuantifiedGroupsAreTheGroupsTheProteinTableShows(searchOut, proteinLines);
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
        /// <c>Reporter_1..Reporter_n</c> are written in.
        ///
        /// It does NOT catch a desync between MetaMorpheus's intensity array and the design's
        /// ISampleInfo array: this compares index j on one side to index j on the other, so a design
        /// rotated by a channel agrees with itself and the totals still balance. The header assertion
        /// in the caller, which names the channel each column must belong to, is what catches that.
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
        /// AllProteinGroups.tsv must come out of a quantified run with MORE columns than an unquantified
        /// one, not fewer.
        /// </summary>
        /// <remarks>
        /// The engine delivers its results by assigning IntensitiesBySample and SamplesForQuantification
        /// on each group, and both setters clear the cached SampleGroupResults. MetaMorpheus's copy of the
        /// protein header reads that cache without rebuilding it, unlike mzLib's base header, so before
        /// the analysis repopulated it the file silently lost its SpectralCount_ and CountOccupancy_
        /// columns and never gained the Intensity_ ones -- with every remaining column unchanged, which is
        /// why nothing else here noticed.
        /// </remarks>
        private static void AssertProteinTableKeepsItsColumnsAndGainsChannelIntensities(string searchOut)
        {
            var lines = File.ReadAllLines(Path.Combine(searchOut, "AllProteinGroups.tsv"));
            Assert.That(lines, Has.Length.GreaterThan(1), "at least one protein group must be written");

            var header = lines[0].Split('\t');

            foreach (string channel in Tmt11Channels)
            {
                string label = $"{StagedMzmlName}_{channel}";
                Assert.That(header, Contains.Item($"SpectralCount_{label}"),
                    "the count columns must survive quantification");
                Assert.That(header, Contains.Item($"CountOccupancy_{label}"),
                    "the occupancy columns must survive quantification");
                Assert.That(header, Contains.Item($"Intensity_{label}"),
                    "the per-channel intensity columns the parameters comment promises must appear");
            }

            foreach (var line in lines.Skip(1))
            {
                Assert.That(line.Split('\t'), Has.Length.EqualTo(header.Length),
                    "every row must carry the columns the header advertises");
            }
        }

        /// <summary>
        /// The quantification tables and the protein table must agree about which groups exist. The PSMs
        /// handed to the engine are filtered to targets below threshold, so handing it every group put
        /// rows in ProteinGroupQuantification.tsv -- decoys, groups above the q-value threshold -- that
        /// AllProteinGroups.tsv suppresses.
        /// </summary>
        private static void AssertQuantifiedGroupsAreTheGroupsTheProteinTableShows(string searchOut, string[] proteinQuantLines)
        {
            var written = File.ReadAllLines(Path.Combine(searchOut, "AllProteinGroups.tsv"))
                .Skip(1)
                .Select(line => line.Split('\t')[0])
                .ToList();

            var quantified = proteinQuantLines.Skip(1).Select(line => line.Split('\t')[0]).ToList();

            Assert.That(quantified, Is.EquivalentTo(written),
                "the quantified groups must be exactly the groups AllProteinGroups.tsv shows");
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

        #region Declining to quantify

        // The seven guards in MultiplexQuantificationAnalysis all end the same way -- a warning and an
        // untouched search -- and all but one are reached before the engine is constructed. Driving the
        // method directly, the way PostSearchAnalysisTaskTests drives the other private members of this
        // task, is what makes them testable at all: reaching them through EverythingRunnerEngine would
        // cost a full search each to observe an early return, and three of them (an unreadable design,
        // an unknown multiplex label, a design naming no file in the run) cannot be provoked through
        // the GUI or the toml at all, because both front ends validate the design before the search.

        /// <summary>
        /// Invokes <c>MultiplexQuantificationAnalysis</c> against a hand-built parameter set and returns
        /// the warnings it raised together with the parameters it wrote back into.
        /// </summary>
        private static (List<string> Warnings, PostSearchAnalysisParameters Parameters) RunMultiplexAnalysis(
            List<string> currentRawFileList,
            string outputFolder,
            List<SpectralMatch> allSpectralMatches = null,
            List<ProteinGroup> proteinGroups = null,
            string multiplexModId = "TMT11")
        {
            var task = new PostSearchAnalysisTask { CommonParameters = new CommonParameters() };
            var parameters = new PostSearchAnalysisParameters
            {
                SearchParameters = new SearchParameters { MultiplexModId = multiplexModId },
                CurrentRawFileList = currentRawFileList,
                AllSpectralMatches = allSpectralMatches ?? new List<SpectralMatch>(),
                OutputFolder = outputFolder,
                SearchTaskId = "MultiplexGuardTest"
            };

            task.GetType().GetProperty("Parameters").SetValue(task, parameters);
            task.GetType().GetProperty("ProteinGroups", BindingFlags.NonPublic | BindingFlags.Instance)
                .SetValue(task, proteinGroups);

            var warnings = new List<string>();
            EventHandler<StringEventArgs> handler = (o, e) => warnings.Add(e.S);
            MetaMorpheusTask.WarnHandler += handler;
            try
            {
                task.GetType()
                    .GetMethod("MultiplexQuantificationAnalysis", BindingFlags.NonPublic | BindingFlags.Instance)
                    .Invoke(task, null);
            }
            finally
            {
                MetaMorpheusTask.WarnHandler -= handler;
            }

            return (warnings, parameters);
        }

        /// <summary>Creates an empty scratch folder, replacing any left by an earlier run.</summary>
        private static string StageFolder(string name)
        {
            string folder = Path.Combine(TestContext.CurrentContext.TestDirectory, name);
            if (Directory.Exists(folder)) Directory.Delete(folder, true);
            Directory.CreateDirectory(folder);
            return folder;
        }

        /// <summary>
        /// The task output folder these guards write into. It must not be the folder the design lives
        /// in: the analysis archives the design beside the results, and a copy onto itself throws.
        /// </summary>
        private static string StageOutput(string folder)
        {
            string output = Path.Combine(folder, "out");
            Directory.CreateDirectory(output);
            return output;
        }

        /// <summary>
        /// The spectra file path these guards run against. The mzML is never staged, because every guard
        /// below fires before anything opens it -- the design is read from beside the path, not from it.
        /// </summary>
        private static string RawPathIn(string folder) => Path.Combine(folder, "VA084TQ_6.mzML");

        private static void WriteDesign(string folder, IEnumerable<string> rows) =>
            File.WriteAllLines(
                Path.Combine(folder, GlobalVariables.TmtExperimentalDesignFileName),
                new[] { TmtExperimentalDesign.Header }.Concat(rows));

        /// <summary>Rows annotating every TMT11 channel of <paramref name="rawPath"/> as a study sample.</summary>
        private static IEnumerable<string> ValidDesignRows(string rawPath) =>
            Tmt11Channels.Select((tag, i) =>
                $"{rawPath}\tPlex1\tSample{i + 1}\t{tag}\tCond{(i % 2 == 0 ? "A" : "B")}\t{i / 2 + 1}\t1\t1\tstudy sample");

        /// <summary>
        /// With no spectra files there is nothing to look beside for a design and nothing to quantify.
        /// That is not a misconfiguration, so it returns without warning anyone about anything.
        /// </summary>
        [Test]
        public static void NoSpectraFiles_DeclinesWithoutWarning()
        {
            var (warnings, parameters) = RunMultiplexAnalysis(
                new List<string>(), TestContext.CurrentContext.TestDirectory);

            Assert.That(warnings, Is.Empty);
            Assert.That(parameters.MultiplexQuantificationResults, Is.Null);
        }

        /// <summary>
        /// The warning for a missing design has to name the file the user is expected to create. The
        /// end-to-end test covers that no tables appear; this covers that the user is told why.
        /// </summary>
        [Test]
        public static void NoDesignFile_WarnsNamingTheFileItLookedFor()
        {
            string folder = StageFolder("TmtGuardNoDesign");
            try
            {
                var (warnings, parameters) = RunMultiplexAnalysis(
                    new List<string> { RawPathIn(folder) }, folder);

                Assert.That(warnings, Has.Exactly(1).Contains(GlobalVariables.TmtExperimentalDesignFileName));
                Assert.That(parameters.MultiplexQuantificationResults, Is.Null);
            }
            finally
            {
                Directory.Delete(folder, true);
            }
        }

        /// <summary>
        /// The path this guard exists for. Every row parses and every value is legal, but they all name
        /// some other file, so the design describes nothing about this run. Quantifying it would produce
        /// a table of nothing rather than a table of wrong numbers -- which is why it has to be reported
        /// as an error rather than accepted as an empty design.
        /// </summary>
        [Test]
        public static void DesignNamingNoFileInThisRun_WarnsAndSkips()
        {
            string folder = StageFolder("TmtGuardNoMatch");
            try
            {
                WriteDesign(folder, ValidDesignRows(Path.Combine(folder, "SomeOtherRun.mzML")));

                var (warnings, parameters) = RunMultiplexAnalysis(
                    new List<string> { RawPathIn(folder) }, StageOutput(folder));

                Assert.That(warnings, Has.Count.EqualTo(1),
                    "the design itself archived cleanly, so declining is the only thing to report");
                Assert.That(warnings, Has.Exactly(1).Contains("Error reading TMT design file"));
                Assert.That(warnings, Has.Exactly(1).Contains("name a file in this run"),
                    "the user has to be told the design is pointed at the wrong place");
                Assert.That(parameters.MultiplexQuantificationResults, Is.Null);
            }
            finally
            {
                Directory.Delete(folder, true);
            }
        }

        /// <summary>
        /// A design that does not parse is reported rather than partially applied. A Biological Replicate
        /// of zero is the example because it is a plausible typo rather than a corrupted file.
        /// </summary>
        [Test]
        public static void UnreadableDesign_WarnsAndSkips()
        {
            string folder = StageFolder("TmtGuardBadDesign");
            try
            {
                WriteDesign(folder, new[]
                {
                    $"{RawPathIn(folder)}\tPlex1\tSample1\t126\tCondA\t0\t1\t1\tstudy sample"
                });

                var (warnings, parameters) = RunMultiplexAnalysis(
                    new List<string> { RawPathIn(folder) }, StageOutput(folder));

                Assert.That(warnings, Has.Count.EqualTo(1));
                Assert.That(warnings, Has.Exactly(1).Contains("Error reading TMT design file"));
                Assert.That(warnings, Has.Exactly(1).Contains("Biological Replicate"));
                Assert.That(parameters.MultiplexQuantificationResults, Is.Null);
            }
            finally
            {
                Directory.Delete(folder, true);
            }
        }

        /// <summary>
        /// The design can be perfect and still not be projectable: without a recognised multiplex label
        /// there are no reporter ion m/z values to order the channels by, and that order is the whole
        /// contract between this design and <c>ISpectralMatch.Intensities</c>. Guessing an order here
        /// would mislabel every channel while leaving the totals plausible, so it declines instead.
        /// </summary>
        [Test]
        public static void UnknownMultiplexLabel_WarnsAndSkips()
        {
            string folder = StageFolder("TmtGuardNoTag");
            try
            {
                WriteDesign(folder, ValidDesignRows(RawPathIn(folder)));

                var (warnings, parameters) = RunMultiplexAnalysis(
                    new List<string> { RawPathIn(folder) }, StageOutput(folder),
                    multiplexModId: "NotATagWeKnow");

                Assert.That(warnings, Has.Count.EqualTo(1));
                Assert.That(warnings, Has.Exactly(1).Contains("Could not build a quantification design"));
                Assert.That(warnings, Has.Exactly(1).Contains("multiplex label"));
                Assert.That(parameters.MultiplexQuantificationResults, Is.Null);
            }
            finally
            {
                Directory.Delete(folder, true);
            }
        }

        /// <summary>
        /// A readable, projectable design over a search whose PSMs carry no reporter ions. This is what a
        /// design file left behind from an earlier TMT experiment looks like to a search that found none,
        /// and it must decline rather than write a table of zeroes.
        /// </summary>
        [Test]
        public static void NoPsmCarriesReporterIons_WarnsAndSkips()
        {
            string folder = StageFolder("TmtGuardNoReporterIons");
            try
            {
                WriteDesign(folder, ValidDesignRows(RawPathIn(folder)));

                var (warnings, parameters) = RunMultiplexAnalysis(
                    new List<string> { RawPathIn(folder) }, StageOutput(folder),
                    allSpectralMatches: new List<SpectralMatch>());

                Assert.That(warnings, Has.Count.EqualTo(1));
                Assert.That(warnings, Has.Exactly(1).Contains("No spectral matches carried reporter ion intensities"));
                Assert.That(File.Exists(Path.Combine(StageOutput(folder), GlobalVariables.TmtExperimentalDesignFileName)),
                    Is.True, "a design it could read is archived beside the results even when it declines to quantify");
                Assert.That(parameters.MultiplexQuantificationResults, Is.Null);
            }
            finally
            {
                Directory.Delete(folder, true);
            }
        }

        /// <summary>
        /// Failing to archive the design beside the results is an inconvenience, not a reason to abandon
        /// quantification, so the copy is guarded on its own and the method carries on. Asserting the
        /// second warning is the point of the test: it is what proves the copy did not return early.
        /// </summary>
        [Test]
        public static void DesignCopyFailure_WarnsAndKeepsGoing()
        {
            string folder = StageFolder("TmtGuardCopyFails");
            try
            {
                WriteDesign(folder, ValidDesignRows(RawPathIn(folder)));

                // An output folder that was never created: the copy throws, nothing else does.
                string missingOutput = Path.Combine(folder, "no", "such", "folder");

                var (warnings, parameters) = RunMultiplexAnalysis(
                    new List<string> { RawPathIn(folder) }, missingOutput,
                    allSpectralMatches: new List<SpectralMatch>());

                Assert.That(warnings, Has.Exactly(1).Contains("Could not copy the TMT design file"));
                Assert.That(warnings, Has.Exactly(1).Contains("No spectral matches carried reporter ion intensities"),
                    "a failed copy must not stop the analysis from reaching the next guard");
                Assert.That(parameters.MultiplexQuantificationResults, Is.Null);
            }
            finally
            {
                Directory.Delete(folder, true);
            }
        }

        /// <summary>
        /// The one guard that cannot be reached without a real search, and the one a user is most likely
        /// to hit: multiplex quantification rolls up to protein groups, and the stock
        /// TMT-Task1-SearchTaskconfig.toml leaves parsimony off. A design file and reporter ions are
        /// both present here -- there is simply nothing to roll up onto -- so the search must finish
        /// normally and write no quantification tables.
        /// </summary>
        [Test]
        public static void TmtSearchWithoutParsimony_WarnsAndWritesNoQuantTables()
        {
            string root = Path.Combine(TestContext.CurrentContext.TestDirectory, "TmtMultiplexNoParsimony");
            string dataFolder = Path.Combine(root, "data");
            string outputFolder = Path.Combine(root, "out");
            if (Directory.Exists(root)) Directory.Delete(root, true);

            var warnings = new List<string>();
            EventHandler<StringEventArgs> handler = (o, e) => warnings.Add(e.S);
            MetaMorpheusTask.WarnHandler += handler;
            try
            {
                string stagedMzml = StageMzml(dataFolder);
                WriteTmtDesign(dataFolder, stagedMzml);

                RunTmtSearch(stagedMzml, outputFolder, doParsimony: false);

                string searchOut = Path.Combine(outputFolder, "search");
                Assert.That(File.Exists(Path.Combine(searchOut, "AllPeptides.psmtsv")), Is.True,
                    "the search itself must still complete");
                Assert.That(warnings, Has.Some.Contains("Multiplex quantification needs protein groups"));
                Assert.That(File.Exists(Path.Combine(searchOut, QuantificationWriter.ProteinGroupFileName)), Is.False);
                Assert.That(File.Exists(Path.Combine(searchOut, QuantificationWriter.PeptideFileName)), Is.False);
                Assert.That(File.Exists(Path.Combine(searchOut, QuantificationWriter.RawFileName)), Is.False);
            }
            finally
            {
                MetaMorpheusTask.WarnHandler -= handler;
                if (Directory.Exists(root)) Directory.Delete(root, true);
            }
        }

        /// <summary>
        /// The case neither TMT fixture can reach, and the reason the excluded-match warning exists.
        /// mzLib quantifies only a match that identifies exactly one biopolymer, and one sequence found in
        /// two proteins is two unequal objects, so a peptide shared between two protein groups is dropped.
        /// Every "Shared Peptides" cell in both fixture databases is empty, so nothing else here exercises
        /// it -- and the failure is quiet: the raw table stays intact while the roll-ups total less, or,
        /// when every peptide is shared, nothing at all.
        /// </summary>
        /// <remarks>
        /// The exclusion itself is mzLib's to fix. What belongs here is saying that it happened: the engine
        /// already counts the dropped matches and reports the count on a SUCCESSFUL run, which this
        /// analysis previously read only in the failure branch -- so a total loss produced a clean search
        /// and no warning at all.
        /// </remarks>
        [Test]
        public static void PeptideSharedBetweenTwoProteins_WarnsThatMatchesWereExcluded()
        {
            string folder = StageFolder("TmtGuardSharedPeptide");
            try
            {
                string rawPath = RawPathIn(folder);
                WriteDesign(folder, ValidDesignRows(rawPath));

                var proteinA = new Protein("PEPTIDEKPEPTIDER", "PROTEINA", "ORGANISM");
                var proteinB = new Protein("PEPTIDEKPEPTIDER", "PROTEINB", "ORGANISM");
                var digestionParams = new DigestionParams();
                var fromA = proteinA.Digest(digestionParams, new List<Modification>(), new List<Modification>()).First();
                var fromB = proteinB.Digest(digestionParams, new List<Modification>(), new List<Modification>()).First();

                var psm = SharedPeptidePsm(rawPath, fromA, fromB);

                // Unique peptides is deliberately empty: the two proteins are indistinguishable, so the
                // sequence belongs to the group without being unique to either protein in it.
                var group = new ProteinGroup(
                    new HashSet<IBioPolymer> { proteinA, proteinB },
                    new HashSet<IBioPolymerWithSetMods> { fromA, fromB },
                    new HashSet<IBioPolymerWithSetMods>())
                {
                    AllPsmsBelowOnePercentFDR = new HashSet<ISpectralMatch> { psm }
                };

                var (warnings, parameters) = RunMultiplexAnalysis(
                    new List<string> { rawPath }, StageOutput(folder),
                    allSpectralMatches: new List<SpectralMatch> { psm },
                    proteinGroups: new List<ProteinGroup> { group });

                Assert.That(parameters.MultiplexQuantificationResults, Is.Not.Null,
                    "the run succeeds, which is exactly why the loss has to be said out loud");
                Assert.That(parameters.MultiplexQuantificationResults.AmbiguousSpectralMatchesExcluded, Is.EqualTo(1),
                    "the one match, shared between two proteins, is the one the engine drops");
                Assert.That(warnings, Has.Exactly(1).Contains("did not identify exactly one biopolymer"));
                Assert.That(warnings, Has.Exactly(1).Contains("shared between two protein groups"),
                    "the warning has to name the cause, or a user cannot act on it");
            }
            finally
            {
                Directory.Delete(folder, true);
            }
        }

        /// <summary>
        /// One PSM carrying TMT11 reporter ions whose sequence was found in two proteins. Both candidates
        /// are kept -- equal scores, reportAllAmbiguity -- which is what makes the match ambiguous to
        /// mzLib's filter without making its BaseSequence ambiguous to MetaMorpheus's, so it passes
        /// includeAmbiguous: false and is dropped later.
        /// </summary>
        private static SpectralMatch SharedPeptidePsm(string rawPath, IBioPolymerWithSetMods fromA, IBioPolymerWithSetMods fromB)
        {
            var tag = IsobaricMassTag.GetIsobaricMassTag("TMT11");
            double[] reporterMzs = tag.ReporterIonMzs.ToArray();
            double[] reporterIntensities = Enumerable.Range(1, reporterMzs.Length).Select(i => 1000.0 * i).ToArray();

            var dataScan = new MsDataScan(
                new MzSpectrum(reporterMzs, reporterIntensities, false),
                oneBasedScanNumber: 1, msnOrder: 2, isCentroid: true,
                polarity: Polarity.Positive, retentionTime: 10.0,
                scanWindowRange: null, scanFilter: "f",
                mzAnalyzer: MZAnalyzerType.Orbitrap, totalIonCurrent: reporterIntensities.Sum(),
                injectionTime: 1.0, noiseData: null, nativeId: "scan=1",
                selectedIonMz: 500.0, selectedIonChargeStateGuess: 2,
                selectedIonIntensity: 1, isolationMZ: 500.0, isolationWidth: 2,
                dissociationType: DissociationType.HCD, oneBasedPrecursorScanNumber: null,
                selectedIonMonoisotopicGuessMz: 500.0);

            var scan = new Ms2ScanWithSpecificMass(dataScan, 500.0, 2, rawPath, new CommonParameters());

            // Before the PSM is built: the spectral match copies the array off the scan in its constructor.
            scan.SetIsobaricMassTagReporterIonIntensities(tag);

            SpectralMatch psm = new PeptideSpectralMatch(fromA, 0, 10, 0, scan, new CommonParameters(), new List<MatchedFragmentIon>());
            psm.AddOrReplace(fromB, 10, 0, true, new List<MatchedFragmentIon>());
            psm.SetFdrValues(1, 0, 0.0, 1, 0, 0.0, 0, 0.0);
            psm.ResolveAllAmbiguities();
            return psm;
        }

        #endregion

    }
}
