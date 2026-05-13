using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Threading;
using GuiFunctions;
using MassSpectrometry;
using NUnit.Framework;
using Omics.Fragmentation;
using Omics.SpectrumMatch;
using OxyPlot;
using OxyPlot.Wpf;
using Readers;

namespace Test.MetaDraw
{
    [TestFixture, Apartment(ApartmentState.STA)]
    public class RnaMirrorPlotGeneration
    {
        private const string OsmTsvPath = @"D:\Projects\SingleOligoPaperGit\results\burkemers\Search for M-Ions and Adducted Oligos\Task1-Burkemers\AllOSMs.osmtsv";
        private const string MzmlPathA = @"B:\Users\Burke\RNA\25-06-12_standard_runs_1-5_mix\MS_files\renamed_MS_files_02-26-26\22mer#1-Am.raw";
        private const string MzmlPathB = @"B:\Users\Burke\RNA\25-06-12_standard_runs_1-5_mix\MS_files\renamed_MS_files_02-26-26\22mer#1-m6A.raw";
        private const string OutputDir = @"D:\Projects\SingleOligoPaperGit\results\burkemers\Search for M-Ions and Adducted Oligos\MirrorPlots";

        private const string TargetM6A = "UUCAAGUAA[Biological:N6-methyladenosine on A]UCCAGGAUAGGCU";
        private const string Target2OMe = "UUCAAGUAA[Biological:2'-O-Methyladenosine on A]UCCAGGAUAGGCU";

        [SetUp]
        public void SetUpMetaDrawSettings()
        {
            MetaDrawSettings.ExportType = "png";
            MetaDrawSettings.CanvasPdfExportDpi = 300;
            MetaDrawSettings.AnnotationBold = true;
            MetaDrawSettings.StrokeThicknessAnnotated = 2;
            MetaDrawSettings.UseShortIonAnnotationsWhenPossible = true;

            MetaDrawSettingsViewModel.Instance.DeselectAllSpectrumDescriptorsCommand.Execute(null);
            MetaDrawSettings.SpectrumDescription["Protein: "] = true;
            MetaDrawSettings.SpectrumDescription["Score: "] = true;
            MetaDrawSettings.SpectrumDescription["Retention Time: "] = true;
        }

        [Test]
        public void GenerateMirrorPlots()
        {
            // Step 1: Load .osmtsv
            var parsingParams = new SpectrumMatchParsingParameters
            {
                FragmentIonsHavePlaceholderForEnvelope = true,
                ParseMatchedFragmentIons = true
            };
            var allPsms = SpectrumMatchTsvReader.ReadOsmTsv(OsmTsvPath, out var warnings, parsingParams);
            if (warnings.Count > 0)
                TestContext.WriteLine($"Warnings loading {OsmTsvPath}: {string.Join("; ", warnings)}");
            TestContext.WriteLine($"Loaded {allPsms.Count} PSMs from {OsmTsvPath}");

            // Step 2: Filter by target sequences, sort by Score descending
            var m6aPsms = allPsms.Where(p => p.FullSequence == TargetM6A)
                .OrderByDescending(p => p.Score).ToList();
            var twoOMePsms = allPsms.Where(p => p.FullSequence == Target2OMe)
                .OrderByDescending(p => p.Score).ToList();

            TestContext.WriteLine($"Found {m6aPsms.Count} PSMs for m6A sequence");
            TestContext.WriteLine($"Found {twoOMePsms.Count} PSMs for 2'-OMe sequence");

            Assert.That(m6aPsms, Is.Not.Empty, $"No PSMs found for sequence: {TargetM6A}");
            Assert.That(twoOMePsms, Is.Not.Empty, $"No PSMs found for sequence: {Target2OMe}");

            // Step 3: Load spectra (one file per sequence)
            var msDataFileA = MsDataFileReader.GetDataFile(MzmlPathA).LoadAllStaticData();
            var msDataFileB = MsDataFileReader.GetDataFile(MzmlPathB).LoadAllStaticData();
            TestContext.WriteLine($"Loaded spectra from A: {MzmlPathA}");
            TestContext.WriteLine($"Loaded spectra from B: {MzmlPathB}");

            // Step 4: Pair by index (sorted by score descending)
            int pairCount = Math.Min(m6aPsms.Count, twoOMePsms.Count);
            if (m6aPsms.Count != twoOMePsms.Count)
            {
                TestContext.WriteLine($"Warning: PSM count mismatch — m6A has {m6aPsms.Count}, " +
                    $"2'-OMe has {twoOMePsms.Count}. Generating {pairCount} pairs " +
                    $"(unmatched: {Math.Abs(m6aPsms.Count - twoOMePsms.Count)}).");
            }

            // Step 5: Generate mirror plots
            Directory.CreateDirectory(OutputDir);
            int successfulExports = 0;

            for (int i = 0; i < pairCount; i++)
            {
                var psmM6A = m6aPsms[i];
                var psm2OMe = twoOMePsms[i];

                if (psmM6A.MatchedIons == null || psmM6A.MatchedIons.Count == 0)
                {
                    TestContext.WriteLine($"Warning: No matched ions for m6A PSM #{i + 1} (Scan {psmM6A.Ms2ScanNumber}), skipping");
                    continue;
                }
                if (psm2OMe.MatchedIons == null || psm2OMe.MatchedIons.Count == 0)
                {
                    TestContext.WriteLine($"Warning: No matched ions for 2'-OMe PSM #{i + 1} (Scan {psm2OMe.Ms2ScanNumber}), skipping");
                    continue;
                }

                var scanA = msDataFileA.GetOneBasedScan(psmM6A.Ms2ScanNumber);
                var scanB = msDataFileB.GetOneBasedScan(psm2OMe.Ms2ScanNumber);

                var plotView = new PlotView();
                var mirrorPlot = new MirrorSpectrumMatchPlot(
                    plotView,
                    psmM6A, scanA, psmM6A.MatchedIons,
                    psm2OMe, scanB, psm2OMe.MatchedIons);

                var outputPath = Path.Combine(OutputDir,
                    $"Mirror_Pair{i + 1:D3}_Scan{psmM6A.Ms2ScanNumber}_Score{psmM6A.Score:F1}" +
                    $"_vs_Scan{psm2OMe.Ms2ScanNumber}_Score{psm2OMe.Score:F1}.png");

                mirrorPlot.ExportToPng(outputPath, 700, 370);

                TestContext.WriteLine($"  [{i + 1}/{pairCount}] Exported: {Path.GetFileName(outputPath)}");
                successfulExports++;
            }

            TestContext.WriteLine($"Done. Successfully exported {successfulExports} / {pairCount} mirror plots to: {OutputDir}");
        }
    }
}
