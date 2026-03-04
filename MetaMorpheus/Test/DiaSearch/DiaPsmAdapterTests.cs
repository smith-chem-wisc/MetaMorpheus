// Copyright 2026 MetaMorpheus Contributors
// Licensed under the MIT License

using EngineLayer.DiaSearch;
using MassSpectrometry.Dia;
using NUnit.Framework;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace Test.DiaSearch
{
    /// <summary>
    /// Unit tests for:
    ///   - <see cref="DiaPsmAdapter"/>     (field mapping, fallback logic)
    ///   - <see cref="DiaSequenceParser"/> (base sequence extraction, mod summary)
    ///   - <see cref="DiaPsmTsvWriter"/>   (schema, filtering, sorting, aggregate)
    ///
    /// All tests are self-contained — no raw files, no spectral libraries.
    /// Temp files are cleaned up in [TearDown].
    /// </summary>
    [TestFixture]
    public class DiaPsmAdapterTests
    {
        private List<string> _tempFiles;

        [SetUp]
        public void SetUp() => _tempFiles = new List<string>();

        [TearDown]
        public void TearDown()
        {
            foreach (var f in _tempFiles)
                if (File.Exists(f)) File.Delete(f);
        }

        // =====================================================================
        // DiaPsmAdapter — field mapping (8 tests)
        // =====================================================================

        [Test]
        public void Adapter_FileName_IsPassedThrough()
        {
            var adapter = MakeAdapter("HeLa_01");
            Assert.That(adapter.FileName, Is.EqualTo("HeLa_01"));
        }

        [Test]
        public void Adapter_ScanNumber_IsWindowId()
        {
            var result = MakeResult(windowId: 42);
            var adapter = new DiaPsmAdapter(result, "file");
            Assert.That(adapter.ScanNumber, Is.EqualTo(42));
        }

        [Test]
        public void Adapter_ScanRetentionTime_UsesObservedApexRt_WhenValid()
        {
            // rtWindowStart=10, rtWindowEnd=14 → midpoint=12; apex=12.5 overrides
            var result = MakeResult(rtWindowStart: 10f, rtWindowEnd: 14f);
            result.ObservedApexRt = 12.5f;

            var adapter = new DiaPsmAdapter(result, "f");
            Assert.That(adapter.ScanRetentionTime, Is.EqualTo(12.5).Within(1e-6));
        }

        [Test]
        public void Adapter_ScanRetentionTime_FallsBackToWindowMidpoint_WhenApexNaN()
        {
            var result = MakeResult(rtWindowStart: 10f, rtWindowEnd: 14f);
            result.ObservedApexRt = float.NaN;

            var adapter = new DiaPsmAdapter(result, "f");
            Assert.That(adapter.ScanRetentionTime, Is.EqualTo(12.0).Within(1e-6));
        }

        [Test]
        public void Adapter_Score_UsesClassifierScore_WhenAvailable()
        {
            var result = MakeResult();
            result.ClassifierScore = 0.95f;
            result.SpectralAngleScore = 0.70f;

            var adapter = new DiaPsmAdapter(result, "f");
            Assert.That(adapter.Score, Is.EqualTo(0.95).Within(1e-6));
        }

        [Test]
        public void Adapter_Score_FallsBackToSpectralAngle_WhenClassifierScoreNaN()
        {
            var result = MakeResult();
            result.ClassifierScore = float.NaN;
            result.SpectralAngleScore = 0.82f;

            var adapter = new DiaPsmAdapter(result, "f");
            Assert.That(adapter.Score, Is.EqualTo(0.82).Within(1e-6));
        }

        [Test]
        public void Adapter_QValue_ReturnsOne_WhenFdrInfoNull()
        {
            var result = MakeResult();
            result.FdrInfo = null;

            var adapter = new DiaPsmAdapter(result, "f");
            Assert.That(adapter.QValue, Is.EqualTo(1.0).Within(1e-6));
        }

        [Test]
        public void Adapter_TargetDecoyLabel_IsT_ForTarget_AndD_ForDecoy()
        {
            var target = new DiaPsmAdapter(MakeResult(isDecoy: false), "f");
            var decoy = new DiaPsmAdapter(MakeResult(isDecoy: true), "f");

            Assert.That(target.TargetDecoyLabel, Is.EqualTo("T"));
            Assert.That(decoy.TargetDecoyLabel, Is.EqualTo("D"));
        }

        // =====================================================================
        // DiaSequenceParser (7 tests)
        // =====================================================================

        [Test]
        public void SequenceParser_GetBaseSequence_ReturnsUnmodified_WhenNoMods()
        {
            Assert.That(DiaSequenceParser.GetBaseSequence("PEPTIDEK"), Is.EqualTo("PEPTIDEK"));
        }

        [Test]
        public void SequenceParser_GetBaseSequence_StripsSingleBracketMod()
        {
            Assert.That(DiaSequenceParser.GetBaseSequence("PEPTM[Oxidation]IDEK"),
                Is.EqualTo("PEPTMIDEK"));
        }

        [Test]
        public void SequenceParser_GetBaseSequence_StripsMultipleMods()
        {
            Assert.That(DiaSequenceParser.GetBaseSequence("S[Phospho]PEPTM[Oxidation]IDE"),
                Is.EqualTo("SPEPTMIDE"));
        }

        [Test]
        public void SequenceParser_Parse_ExtractsModNameAndPosition()
        {
            DiaSequenceParser.Parse("PEPTM[Oxidation]IDEK",
                out string baseSeq, out var mods);

            Assert.That(baseSeq, Is.EqualTo("PEPTMIDEK"));
            Assert.That(mods.Count, Is.EqualTo(1));
            var (position, name) = mods[0];
            Assert.That(position, Is.EqualTo(5)); // 1-based: M is 5th residue
            Assert.That(name, Is.EqualTo("Oxidation"));
        }

        [Test]
        public void SequenceParser_Parse_HandlesNoMods()
        {
            DiaSequenceParser.Parse("PEPTIDEK", out string baseSeq, out var mods);
            Assert.That(baseSeq, Is.EqualTo("PEPTIDEK"));
            Assert.That(mods.Count, Is.EqualTo(0));
        }

        [Test]
        public void SequenceParser_FormatModificationSummary_EmptyString_WhenNoMods()
        {
            Assert.That(DiaSequenceParser.FormatModificationSummary("PEPTIDEK"),
                Is.EqualTo(string.Empty));
        }

        [Test]
        public void SequenceParser_FormatModificationSummary_ContainsModNameAndPosition()
        {
            string summary = DiaSequenceParser.FormatModificationSummary("PEPTM[Oxidation]IDEK");
            // Expected format: "Oxidation@5"
            Assert.That(summary, Does.Contain("Oxidation"));
            Assert.That(summary, Does.Contain("5"));
        }

        // =====================================================================
        // DiaPsmTsvWriter (8 tests)
        // =====================================================================

        [Test]
        public void TsvWriter_ColumnHeaders_Has27Columns()
        {
            Assert.That(DiaPsmTsvWriter.ColumnHeaders.Count, Is.EqualTo(27));
        }

        [Test]
        public void TsvWriter_WriteToFile_ProducesHeaderAndOneDataRow()
        {
            var adapters = new List<DiaPsmAdapter> { MakeAdapter("file1") };
            string path = TempPath("single.psmtsv");

            DiaPsmTsvWriter.WriteToFile(adapters, path,
                qValueThreshold: 1.0, includeDecoys: false);

            var lines = NonEmptyLines(path);
            Assert.That(lines.Length, Is.EqualTo(2));                 // header + 1 row
            Assert.That(lines[0].Split('\t').Length, Is.EqualTo(27)); // correct column count
        }

        [Test]
        public void TsvWriter_WriteToFile_FiltersDecoys_WhenIncludeDecoysFalse()
        {
            var adapters = new List<DiaPsmAdapter>
            {
                new DiaPsmAdapter(MakeResult(isDecoy: false), "f"),
                new DiaPsmAdapter(MakeResult(isDecoy: true),  "f"),
            };
            string path = TempPath("nodecoy.psmtsv");

            DiaPsmTsvWriter.WriteToFile(adapters, path,
                qValueThreshold: 1.0, includeDecoys: false);

            Assert.That(NonEmptyLines(path).Length, Is.EqualTo(2)); // header + 1 target only
        }

        [Test]
        public void TsvWriter_WriteToFile_IncludesDecoys_WhenIncludeDecoysTrue()
        {
            var adapters = new List<DiaPsmAdapter>
            {
                new DiaPsmAdapter(MakeResult(isDecoy: false), "f"),
                new DiaPsmAdapter(MakeResult(isDecoy: true),  "f"),
            };
            string path = TempPath("withdecoy.psmtsv");

            DiaPsmTsvWriter.WriteToFile(adapters, path,
                qValueThreshold: 1.0, includeDecoys: true);

            Assert.That(NonEmptyLines(path).Length, Is.EqualTo(3)); // header + 2 rows
        }

        [Test]
        public void TsvWriter_WriteToFile_FiltersHighQValueResults()
        {
            var r1 = MakeResult(); r1.FdrInfo = new DiaFdrInfo { QValue = 0.005 };
            var r2 = MakeResult(); r2.FdrInfo = new DiaFdrInfo { QValue = 0.050 };
            var adapters = new List<DiaPsmAdapter>
            {
                new DiaPsmAdapter(r1, "f"),
                new DiaPsmAdapter(r2, "f"),
            };
            string path = TempPath("qvalfilter.psmtsv");

            DiaPsmTsvWriter.WriteToFile(adapters, path,
                qValueThreshold: 0.01, includeDecoys: false);

            Assert.That(NonEmptyLines(path).Length, Is.EqualTo(2)); // header + r1 only
        }

        [Test]
        public void TsvWriter_BuildAdapters_SortsByScoreDescending()
        {
            var r1 = MakeResult(); r1.ClassifierScore = 0.9f;
            var r2 = MakeResult(); r2.ClassifierScore = 0.6f;
            var r3 = MakeResult(); r3.ClassifierScore = 0.8f;

            var adapters = DiaPsmTsvWriter.BuildAdapters(new[] { r1, r2, r3 }, "f");

            Assert.That(adapters[0].Score, Is.EqualTo(0.9).Within(0.001));
            Assert.That(adapters[1].Score, Is.EqualTo(0.8).Within(0.001));
            Assert.That(adapters[2].Score, Is.EqualTo(0.6).Within(0.001));
        }

        [Test]
        public void TsvWriter_WriteAggregate_MergesMultipleFiles()
        {
            var file1 = new List<DiaPsmAdapter> { MakeAdapter("file1") };
            var file2 = new List<DiaPsmAdapter> { MakeAdapter("file2") };
            string path = TempPath("agg.psmtsv");

            DiaPsmTsvWriter.WriteAggregate(
                new[] { file1, file2 }, path,
                qValueThreshold: 1.0, includeDecoys: false);

            Assert.That(NonEmptyLines(path).Length, Is.EqualTo(3)); // header + 2 rows
        }

        [Test]
        public void TsvWriter_WriteToFile_HandlesEmptyList_WithoutThrowing()
        {
            string path = TempPath("empty.psmtsv");

            Assert.DoesNotThrow(() =>
                DiaPsmTsvWriter.WriteToFile(
                    new List<DiaPsmAdapter>(), path,
                    qValueThreshold: 0.01, includeDecoys: false));

            Assert.That(NonEmptyLines(path).Length, Is.EqualTo(1)); // header only
        }

        // =====================================================================
        // Helpers
        // =====================================================================

        private string TempPath(string name)
        {
            string path = Path.Combine(
                Path.GetTempPath(), $"DiaPsmAdapterTests_{name}");
            _tempFiles.Add(path);
            return path;
        }

        private static string[] NonEmptyLines(string path) =>
            File.ReadAllLines(path).Where(l => !string.IsNullOrWhiteSpace(l)).ToArray();

        /// <summary>
        /// Creates a minimal valid DiaSearchResult with FdrInfo populated.
        /// Constructor args match the DiaSearchResult(string, int, double, int, bool,
        /// int, double?, float, float) signature exactly.
        /// </summary>
        private static DiaSearchResult MakeResult(
            int windowId = 1,
            bool isDecoy = false,
            float rtWindowStart = 21.0f,
            float rtWindowEnd = 24.0f)
        {
            var r = new DiaSearchResult(
                sequence: "PEPTIDEK",
                chargeState: 2,
                precursorMz: 453.2567,
                windowId: windowId,
                isDecoy: isDecoy,
                fragmentsQueried: 6,
                libraryRetentionTime: 22.5,
                rtWindowStart: rtWindowStart,
                rtWindowEnd: rtWindowEnd);

            r.SpectralAngleScore = 0.75f;
            r.ClassifierScore = float.NaN;
            r.ObservedApexRt = 22.3f;
            r.FdrInfo = new DiaFdrInfo
            {
                QValue = 0.005,
                PEP = 0.002,
                PEP_QValue = 0.003,
                PeptideQValue = 0.004,
            };
            return r;
        }

        private static DiaPsmAdapter MakeAdapter(string fileName = "testFile") =>
            new DiaPsmAdapter(MakeResult(), fileName);
    }
}