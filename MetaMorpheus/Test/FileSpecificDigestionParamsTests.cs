using EngineLayer;
using GuiFunctions.Util;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Omics.Digestion;
using Omics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System;
using System.IO;
using TaskLayer;

namespace Test
{
    /// <summary>
    /// Locks down that file-specific settings change only the digestion settings they name, and carry every other digestion
    /// setting of the task through unchanged.
    /// </summary>
    /// <remarks>
    /// <para><see cref="MetaMorpheusTask.SetAllFileSpecificCommonParams"/> rebuilds the task's <see cref="DigestionParams"/>
    /// from an explicit list of fields whenever a spectra file has its own settings. That list left out three fields, which
    /// then silently fell back to their constructor defaults for every file with file-specific settings, even a file whose
    /// settings only change a mass tolerance:</para>
    /// <list type="bullet">
    /// <item><description><see cref="DigestionParams.KeepNGlycopeptide"/> and <see cref="DigestionParams.KeepOGlycopeptide"/>
    /// (default false): a glycopeptide digest restricted to candidates carrying a glycosylation motif became an unrestricted
    /// digest for those files.</description></item>
    /// <item><description><see cref="DigestionParams.GeneratehUnlabeledProteinsForSilac"/> (default true): a SILAC search set
    /// not to quantify unlabeled peptides quantified them anyway for those files (PostSearchAnalysisTask reads this flag).</description></item>
    /// </list>
    /// <para>The second test compares the whole object with <see cref="DigestionParams.Equals(DigestionParams)"/>, which checks
    /// every field, so a digestion setting added in the future that is forgotten in the rebuild fails here too --
    /// <b>but only if this fixture sets that setting away from its constructor default</b>. The rebuild's fallback IS the
    /// constructor default, so a field left at its default agrees by accident and the guard says nothing.
    /// <see cref="DigestionParams.RespectCleavageBlockingModifications"/> was added to mzLib defaulting to false and was
    /// forgotten in the rebuild, and this test passed anyway until the line below set it true. Set every new digestion
    /// setting here, not just the one being added.</para>
    /// </remarks>
    [TestFixture]
    [System.Diagnostics.CodeAnalysis.ExcludeFromCodeCoverage]
    public static class FileSpecificDigestionParamsTests
    {
        /// <summary>
        /// With a file-specific override of one digestion setting (missed cleavages), the glycopeptide and SILAC flags keep
        /// the task's values in every combination, including the non-default ones.
        /// </summary>
        [Test]
        [TestCase(false, false, true)]
        [TestCase(true, false, true)]
        [TestCase(false, true, true)]
        [TestCase(true, true, false)]
        [TestCase(false, false, false)]
        public static void FileSpecificRebuild_KeepsGlycopeptideAndSilacFlags(bool keepNGlycopeptide, bool keepOGlycopeptide, bool generateUnlabeledProteinsForSilac)
        {
            var task = new CommonParameters(digestionParams: new DigestionParams("trypsin",
                generateUnlabeledProteinsForSilac: generateUnlabeledProteinsForSilac,
                keepNGlycopeptide: keepNGlycopeptide, keepOGlycopeptide: keepOGlycopeptide));

            var combined = MetaMorpheusTask.SetAllFileSpecificCommonParams(task, new FileSpecificParameters { MaxMissedCleavages = 5 });

            var digestionParams = (DigestionParams)combined.DigestionParams;
            Assert.That(digestionParams.MaxMissedCleavages, Is.EqualTo(5), "the file-specific override still applies");
            Assert.That(digestionParams.KeepNGlycopeptide, Is.EqualTo(keepNGlycopeptide), "KeepNGlycopeptide was reset by the file-specific rebuild");
            Assert.That(digestionParams.KeepOGlycopeptide, Is.EqualTo(keepOGlycopeptide), "KeepOGlycopeptide was reset by the file-specific rebuild");
            Assert.That(digestionParams.GeneratehUnlabeledProteinsForSilac, Is.EqualTo(generateUnlabeledProteinsForSilac), "GeneratehUnlabeledProteinsForSilac was reset by the file-specific rebuild");
        }

        /// <summary>
        /// File-specific settings that change nothing about digestion (here a precursor tolerance) must leave the digestion
        /// parameters equal to the task's, field for field. Every field is set away from its default so a dropped one shows.
        /// </summary>
        [Test]
        public static void FileSpecificRebuild_WithoutDigestionOverrides_KeepsEveryDigestionSetting()
        {
            var original = new DigestionParams("trypsin", maxMissedCleavages: 4, minPeptideLength: 6, maxPeptideLength: 40,
                maxModificationIsoforms: 256, initiatorMethionineBehavior: InitiatorMethionineBehavior.Cleave, maxModsForPeptides: 3,
                searchModeType: CleavageSpecificity.Semi, fragmentationTerminus: FragmentationTerminus.C,
                generateUnlabeledProteinsForSilac: false, keepNGlycopeptide: true, keepOGlycopeptide: true,
                respectCleavageBlockingModifications: true);
            var task = new CommonParameters(digestionParams: original);

            var combined = MetaMorpheusTask.SetAllFileSpecificCommonParams(task, new FileSpecificParameters { PrecursorMassTolerance = new PpmTolerance(7) });

            Assert.That(combined.PrecursorMassTolerance.Value, Is.EqualTo(7), "the file-specific override still applies");
            Assert.That(combined.DigestionParams, Is.EqualTo(original),
                $"file-specific rebuild changed digestion settings it was not asked to change: {original} became {combined.DigestionParams}");
        }

        /// <summary>
        /// The cleavage-blocking flag keeps the task's value through a file-specific rebuild, in both states.
        /// </summary>
        /// <remarks>
        /// Named separately from the whole-object test because the flag decides whether digestion is modification-aware at
        /// all: silently falling back to false gives a file the historical modification-blind digest while the task reports
        /// the corrected one, and the flag also feeds the index cache fingerprint, so the two runs are not distinguishable
        /// afterwards. The true case is the one that failed before the rebuild carried the field.
        /// </remarks>
        [Test]
        [TestCase(true)]
        [TestCase(false)]
        public static void FileSpecificRebuild_KeepsRespectCleavageBlockingModifications(bool respectCleavageBlockingModifications)
        {
            var task = new CommonParameters(digestionParams: new DigestionParams("trypsin",
                respectCleavageBlockingModifications: respectCleavageBlockingModifications));

            var combined = MetaMorpheusTask.SetAllFileSpecificCommonParams(task, new FileSpecificParameters { MaxMissedCleavages = 5 });

            var digestionParams = (DigestionParams)combined.DigestionParams;
            Assert.That(digestionParams.MaxMissedCleavages, Is.EqualTo(5), "the file-specific override still applies");
            Assert.That(digestionParams.RespectCleavageBlockingModifications, Is.EqualTo(respectCleavageBlockingModifications),
                "RespectCleavageBlockingModifications was reset by the file-specific rebuild");
        }

        [Test]
        public static void FileSpecificRebuild_PreservesRetentionTimeRange()
        {
            var task = new CommonParameters(
                retentionTimeRange: new DoubleRange(12.5, 48.75));

            var combined = MetaMorpheusTask.SetAllFileSpecificCommonParams(task,
                new FileSpecificParameters { PrecursorMassTolerance = new PpmTolerance(7) });

            Assert.That(combined.RetentionTimeRange.Minimum, Is.EqualTo(12.5));
            Assert.That(combined.RetentionTimeRange.Maximum, Is.EqualTo(48.75));
        }

        [Test]
        [TestCase(-1, 10)]
        [TestCase(double.NaN, 10)]
        public static void CommonParameters_RejectsInvalidRetentionTimeMinimum(double minimum, double maximum)
        {
            var exception = Assert.Throws<ArgumentOutOfRangeException>(() => new CommonParameters(
                retentionTimeRange: new DoubleRange(minimum, maximum)));

            Assert.That(exception.ParamName, Is.EqualTo("retentionTimeRange"));
        }

        [Test]
        [TestCase(0, double.NaN)]
        [TestCase(10, 9)]
        public static void CommonParameters_RejectsInvalidRetentionTimeMaximum(double minimum, double maximum)
        {
            var exception = Assert.Throws<ArgumentOutOfRangeException>(() => new CommonParameters(
                retentionTimeRange: new DoubleRange(minimum, maximum)));

            Assert.That(exception.ParamName, Is.EqualTo("retentionTimeRange"));
        }

        [Test]
        public static void RetentionTimeRange_RoundTripsThroughToml()
        {
            string path = Path.Combine(TestContext.CurrentContext.TestDirectory, "RetentionTimeRangeRoundTrip.toml");
            var task = new SearchTask
            {
                CommonParameters = new CommonParameters(
                    retentionTimeRange: new DoubleRange(12.5, 48.75))
            };

            try
            {
                Nett.Toml.WriteFile(task, path, MetaMorpheusTask.tomlConfig);
                SearchTask loaded = Nett.Toml.ReadFile<SearchTask>(path, MetaMorpheusTask.tomlConfig);

                Assert.That(loaded.CommonParameters.RetentionTimeRange.Minimum, Is.EqualTo(12.5));
                Assert.That(loaded.CommonParameters.RetentionTimeRange.Maximum, Is.EqualTo(48.75));
            }
            finally
            {
                if (File.Exists(path))
                    File.Delete(path);
            }
        }

        [Test]
        public static void RetentionTimeRangeParser_UsesBlankDefaults()
        {
            bool parsed = RetentionTimeRangeParser.TryParse(string.Empty, string.Empty,
                out double minimum, out double maximum, out string errorMessage);

            Assert.That(parsed, Is.True);
            Assert.That(minimum, Is.EqualTo(0));
            Assert.That(maximum, Is.EqualTo(double.MaxValue));
            Assert.That(errorMessage, Is.Null);
        }

        [TestCase("NaN", "10")]
        [TestCase("0", "NaN")]
        [TestCase("Infinity", "10")]
        [TestCase("10", "9")]
        public static void RetentionTimeRangeParser_RejectsInvalidBounds(string minimumText, string maximumText)
        {
            bool parsed = RetentionTimeRangeParser.TryParse(minimumText, maximumText,
                out _, out _, out string errorMessage);

            Assert.That(parsed, Is.False);
            Assert.That(errorMessage, Is.Not.Null.And.Not.Empty);
        }

        [TestCase("12.5;9")]
        [TestCase("NaN;10")]
        [TestCase("12.5;not-a-number")]
        public static void RetentionTimeRange_TomlRejectsInvalidValues(string invalidRange)
        {
            var task = new SearchTask
            {
                CommonParameters = new CommonParameters(retentionTimeRange: new DoubleRange(12.5, 48.75))
            };
            string toml = Nett.Toml.WriteString(task, MetaMorpheusTask.tomlConfig);
            Assert.That(toml, Does.Contain("12.5;48.75"));
            toml = toml.Replace("12.5;48.75", invalidRange, StringComparison.Ordinal);

            var exception = Assert.Throws<InvalidOperationException>(() => Nett.Toml.ReadString<SearchTask>(toml, MetaMorpheusTask.tomlConfig));
            Assert.That(exception.InnerException?.InnerException, Is.TypeOf<MetaMorpheusException>());
        }
    }
}
