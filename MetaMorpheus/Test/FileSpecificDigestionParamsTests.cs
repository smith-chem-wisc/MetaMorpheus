using EngineLayer;
using MzLibUtil;
using NUnit.Framework;
using Omics.Digestion;
using Omics.Fragmentation;
using Proteomics.ProteolyticDigestion;
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
    /// every field, so a digestion setting added in the future that is forgotten in the rebuild fails here too.</para>
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
                generateUnlabeledProteinsForSilac: false, keepNGlycopeptide: true, keepOGlycopeptide: true);
            var task = new CommonParameters(digestionParams: original);

            var combined = MetaMorpheusTask.SetAllFileSpecificCommonParams(task, new FileSpecificParameters { PrecursorMassTolerance = new PpmTolerance(7) });

            Assert.That(combined.PrecursorMassTolerance.Value, Is.EqualTo(7), "the file-specific override still applies");
            Assert.That(combined.DigestionParams, Is.EqualTo(original),
                $"file-specific rebuild changed digestion settings it was not asked to change: {original} became {combined.DigestionParams}");
        }

        [Test]
        public static void FileSpecificRebuild_PreservesRetentionTimeRange()
        {
            var task = new CommonParameters(
                minRetentionTimeToSearch: 12.5,
                maxRetentionTimeToSearch: 48.75);

            var combined = MetaMorpheusTask.SetAllFileSpecificCommonParams(task,
                new FileSpecificParameters { PrecursorMassTolerance = new PpmTolerance(7) });

            Assert.That(combined.MinRetentionTimeToSearch, Is.EqualTo(12.5));
            Assert.That(combined.MaxRetentionTimeToSearch, Is.EqualTo(48.75));
        }
    }
}
