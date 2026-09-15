using EngineLayer;
using Nett;
using NUnit.Framework;
using Omics.Digestion;
using Omics.Fragmentation;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using TaskLayer;

namespace Test
{
    /// <summary>
    /// Task settings files that name the protease "semi-trypsin" must keep loading, and must keep meaning a semi-specific
    /// tryptic digest.
    /// </summary>
    /// <remarks>
    /// <para><b>Background.</b> mzLib shipped a protease called "semi-trypsin" (cleave after K or R, no proline rule,
    /// specificity semi) until mzLib #1005 removed it in February 2026. The replacement it intended is trypsin with
    /// <c>SearchModeType = Semi</c>. Settings files written before that name it, including the O-Pair Search paper's glyco
    /// settings (StcE + trypsin mucin digests). Since #1005, loading one throws KeyNotFoundException: the command line
    /// crashes and the GUI refuses the file.</para>
    ///
    /// <para><b>What loading now does.</b> The digestion parameters are read as the shipped fully specific protease with
    /// exactly semi-trypsin's cleavage motifs (found by motif, not by name, because the name of the no-proline-rule trypsin
    /// is changing in mzLib #1186) and <c>SearchModeType = Semi</c>. Every other digestion setting is kept, and the user is
    /// warned through <see cref="MetaMorpheusTask.WarnHandler"/>.</para>
    ///
    /// <para><b>Why SearchModeType Semi with FragmentationTerminus Both is the right translation.</b> With terminus Both,
    /// mzLib returns the real semi-specific peptides (since mzLib #1303), which is what a Semi protease returned. With N or C
    /// it would return seeds that only the non-specific search engine can use; see
    /// <c>DigestionParams.SearchModeType</c> in mzLib. So a file that said Full becomes Semi with its terminus unchanged
    /// (Both in every glyco file), and a non-specific search file that already said Semi or None keeps that.</para>
    ///
    /// <para>These tests need an mzLib that contains #1303; before it, Semi + Both returned C-terminal seeds.</para>
    /// </remarks>
    [TestFixture]
    [NonParallelizable] // subscribes to the static MetaMorpheusTask.WarnHandler
    public static class LegacySemiTrypsinSettingsTests
    {
        private static string GlycoTestData(string fileName) => Path.Combine(TestContext.CurrentContext.TestDirectory, "GlycoTestData", fileName);

        /// <summary>The motifs "semi-trypsin" had in mzLib's proteases.tsv before #1005.</summary>
        private const string SemiTrypsinMotifs = "K|,R|";

        private static string MotifSignature(IEnumerable<DigestionMotif> motifs) => string.Join(";", motifs
            .Select(m => $"{m.InducingCleavage}|{m.PreventingCleavage}|{m.CutIndex}|{m.ExcludeFromWildcard}")
            .OrderBy(s => s, StringComparer.Ordinal));

        /// <summary>Loads a task file the way the GUI or command line does, collecting the warnings raised while loading.</summary>
        private static (GlycoSearchTask Task, List<string> Warnings) LoadGlycoTask(string path, bool throughLowResFallbackReader)
        {
            var warnings = new List<string>();
            EventHandler<StringEventArgs> collect = (_, e) => warnings.Add(e.S);
            MetaMorpheusTask.WarnHandler += collect;
            try
            {
                GlycoSearchTask task = throughLowResFallbackReader
                    ? MetaMorpheusTask.ReadTaskTomlWithLowResFallback<GlycoSearchTask>(path)
                    : Toml.ReadFile<GlycoSearchTask>(path, MetaMorpheusTask.tomlConfig);
                return (task, warnings);
            }
            finally
            {
                MetaMorpheusTask.WarnHandler -= collect;
            }
        }

        /// <summary>Writes a copy of a settings file with some lines replaced, for the cases no shipped file covers.</summary>
        private static string WriteModifiedCopy(string sourcePath, string testName, params (string OldLine, string NewLine)[] replacements)
        {
            string text = File.ReadAllText(sourcePath);
            foreach (var (oldLine, newLine) in replacements)
            {
                Assert.That(text, Does.Contain(oldLine), $"premise: {Path.GetFileName(sourcePath)} contains the line being replaced");
                text = text.Replace(oldLine, newLine);
            }
            string path = Path.Combine(TestContext.CurrentContext.TestDirectory, $"{nameof(LegacySemiTrypsinSettingsTests)}_{testName}.toml");
            File.WriteAllText(path, text);
            return path;
        }

        /// <summary>
        /// The two shipped glyco settings files that name "semi-trypsin" (EThcD mucin searches) load, through both readers
        /// the GUI and command line use, as the no-proline-rule fully specific trypsin plus SearchModeType Semi, with every
        /// other digestion setting kept and a warning that names the old protease.
        /// </summary>
        [Test]
        [TestCase("GlycoSearchTaskconfig_ETD.toml", true, true)]
        [TestCase("GlycoSearchTaskconfig_ETD.toml", false, true)]
        [TestCase("GlycoSearchTaskconfig_ETD_Run3.toml", true, false)]
        [TestCase("GlycoSearchTaskconfig_ETD_Run3.toml", false, false)]
        public static void LegacyGlycoSettingsNamingSemiTrypsin_LoadAsTrypsinWithSemiSearchMode(string fileName, bool throughLowResFallbackReader, bool expectedUnlabeledSilac)
        {
            string path = GlycoTestData(fileName);
            Assert.That(File.ReadAllText(path), Does.Contain("Protease = \"semi-trypsin\""), "premise: the file names semi-trypsin");

            var (task, warnings) = LoadGlycoTask(path, throughLowResFallbackReader);
            var digestionParams = (DigestionParams)task.CommonParameters.DigestionParams;

            Assert.That(MotifSignature(digestionParams.Protease.DigestionMotifs), Is.EqualTo(MotifSignature(DigestionMotif.ParseDigestionMotifsFromString(SemiTrypsinMotifs))),
                "the protease must cut exactly where semi-trypsin cut: after every K and R, with no proline rule");
            Assert.That(digestionParams.Protease.CleavageSpecificity, Is.EqualTo(CleavageSpecificity.Full), "the semi-specificity now comes from the search mode, not the protease");
            Assert.That(digestionParams.SpecificProtease, Is.SameAs(digestionParams.Protease));
            Assert.That(digestionParams.SearchModeType, Is.EqualTo(CleavageSpecificity.Semi));
            Assert.That(digestionParams.FragmentationTerminus, Is.EqualTo(FragmentationTerminus.Both), "Both asks for peptides; N or C would ask for seeds a glyco search cannot use");

            // every other setting in the file survives the translation
            Assert.That(digestionParams.MaxMissedCleavages, Is.EqualTo(12));
            Assert.That(digestionParams.MinPeptideLength, Is.EqualTo(5));
            Assert.That(digestionParams.MaxPeptideLength, Is.EqualTo(60));
            Assert.That(digestionParams.MaxModificationIsoforms, Is.EqualTo(1024));
            Assert.That(digestionParams.MaxModsForPeptide, Is.EqualTo(2));
            Assert.That(digestionParams.InitiatorMethionineBehavior, Is.EqualTo(InitiatorMethionineBehavior.Variable));
            Assert.That(digestionParams.GeneratehUnlabeledProteinsForSilac, Is.EqualTo(expectedUnlabeledSilac));
            Assert.That(digestionParams.KeepNGlycopeptide, Is.False);
            Assert.That(digestionParams.KeepOGlycopeptide, Is.False);

            Assert.That(warnings.Count(w => w.Contains("semi-trypsin")), Is.EqualTo(1), "the user must be told, once, that the file was translated");
            Assert.That(task.CommonParameters.DissociationType, Is.EqualTo(MassSpectrometry.DissociationType.Custom), "premise: the rest of the task still loaded");
        }

        /// <summary>
        /// The translation changes the settings but not the search: the loaded parameters digest a protein into exactly
        /// the peptides the old Semi protease gave. Uses the mucin leukosialin (CD43, P16150) that these EThcD settings
        /// were written for, with the old protease re-created under a name of its own.
        /// </summary>
        [Test]
        public static void LegacyGlycoSettingsNamingSemiTrypsin_DigestTheSamePeptidesTheOldSemiProteaseDid()
        {
            const string oldSemiTrypsinName = "LegacySemiTrypsinSettingsTests-semi-trypsin-as-it-was";
            ProteaseDictionary.Dictionary[oldSemiTrypsinName] = new Protease(oldSemiTrypsinName, CleavageSpecificity.Semi, null, null,
                DigestionMotif.ParseDigestionMotifsFromString(SemiTrypsinMotifs));
            try
            {
                var migrated = (DigestionParams)LoadGlycoTask(GlycoTestData("GlycoSearchTaskconfig_ETD.toml"), true).Task.CommonParameters.DigestionParams;
                var asItWas = new DigestionParams(oldSemiTrypsinName, migrated.MaxMissedCleavages, migrated.MinPeptideLength, migrated.MaxPeptideLength,
                    migrated.MaxModificationIsoforms, migrated.InitiatorMethionineBehavior, migrated.MaxModsForPeptide);

                Protein leukosialin = UsefulProteomicsDatabases.ProteinDbLoader.LoadProteinFasta(GlycoTestData("P16150.fasta"), true, UsefulProteomicsDatabases.DecoyType.None, false, out _).Single();
                List<(int, int)> Peptides(DigestionParams dp) => leukosialin.Digest(dp, new List<Omics.Modifications.Modification>(), new List<Omics.Modifications.Modification>())
                    .Select(p => (p.OneBasedStartResidueInProtein, p.OneBasedEndResidueInProtein)).OrderBy(p => p).ToList();

                List<(int, int)> expected = Peptides(asItWas);
                Assert.That(expected, Has.Count.GreaterThan(1000), "premise: a 12-missed-cleavage semi-tryptic digest of leukosialin is large");
                Assert.That(Peptides(migrated), Is.EqualTo(expected), "trypsin + SearchModeType Semi must give exactly semi-trypsin's peptides (needs mzLib #1303)");
            }
            finally
            {
                ProteaseDictionary.Dictionary.Remove(oldSemiTrypsinName);
            }
        }

        /// <summary>
        /// A non-specific search file already asks for seeds through its search mode, so only the protease is translated:
        /// Semi stays Semi and None stays None, each with the terminus the file gave.
        /// </summary>
        [Test]
        [TestCase("Semi", "N", "semi-trypsin")]
        [TestCase("Semi", "C", "semi-trypsin")]
        [TestCase("None", "N", "singleN")]
        [TestCase("None", "C", "singleC")]
        public static void SettingsNamingSemiTrypsinWithASeedSearchMode_KeepTheirSearchModeAndTerminus(string searchModeType, string terminus, string savedProtease)
        {
            string path = WriteModifiedCopy(GlycoTestData("GlycoSearchTaskconfig_ETD.toml"), $"{searchModeType}{terminus}",
                // anchored at the line start so SpecificProtease (which keeps naming semi-trypsin) is not also replaced
                ("\nProtease = \"semi-trypsin\"", $"\nProtease = \"{savedProtease}\""),
                ("SearchModeType = \"Full\"", $"SearchModeType = \"{searchModeType}\""),
                ("FragmentationTerminus = \"Both\"", $"FragmentationTerminus = \"{terminus}\""));
            try
            {
                var (task, warnings) = LoadGlycoTask(path, true);
                var digestionParams = (DigestionParams)task.CommonParameters.DigestionParams;

                Assert.That(digestionParams.SearchModeType.ToString(), Is.EqualTo(searchModeType));
                Assert.That(digestionParams.FragmentationTerminus.ToString(), Is.EqualTo(terminus));
                Assert.That(MotifSignature(digestionParams.SpecificProtease.DigestionMotifs), Is.EqualTo(MotifSignature(DigestionMotif.ParseDigestionMotifsFromString(SemiTrypsinMotifs))));
                if (searchModeType == "None")
                    Assert.That(digestionParams.Protease.Name, Is.EqualTo(savedProtease), "a non-specific search still digests with singleN or singleC");
                Assert.That(warnings.Count(w => w.Contains("semi-trypsin")), Is.EqualTo(1));
            }
            finally
            {
                File.Delete(path);
            }
        }

        /// <summary>
        /// A user who defined their own protease called "semi-trypsin" (in a custom proteases file) gets it exactly as
        /// defined: nothing is translated and nobody is warned.
        /// </summary>
        [Test]
        public static void SettingsNamingSemiTrypsin_WhenTheUserDefinedAProteaseWithThatName_UseIt()
        {
            var userDefined = new Protease("semi-trypsin", CleavageSpecificity.Semi, null, null, DigestionMotif.ParseDigestionMotifsFromString("K|"));
            ProteaseDictionary.Dictionary["semi-trypsin"] = userDefined;
            try
            {
                var (task, warnings) = LoadGlycoTask(GlycoTestData("GlycoSearchTaskconfig_ETD.toml"), true);
                var digestionParams = (DigestionParams)task.CommonParameters.DigestionParams;

                Assert.That(digestionParams.Protease, Is.SameAs(userDefined));
                Assert.That(digestionParams.SearchModeType, Is.EqualTo(CleavageSpecificity.Full));
                Assert.That(warnings.Where(w => w.Contains("semi-trypsin")), Is.Empty);
            }
            finally
            {
                ProteaseDictionary.Dictionary.Remove("semi-trypsin");
            }
        }

        /// <summary>
        /// Only names mzLib actually removed are translated. Any other unknown protease still fails to load, as it always
        /// has, rather than silently becoming something else.
        /// </summary>
        [Test]
        public static void SettingsNamingAnUnknownProtease_StillFailToLoad()
        {
            string path = WriteModifiedCopy(GlycoTestData("GlycoSearchTaskconfig_ETD.toml"), "UnknownProtease",
                ("\"semi-trypsin\"", "\"no-such-protease\"")); // both Protease and SpecificProtease
            try
            {
                Assert.That(() => MetaMorpheusTask.ReadTaskTomlWithLowResFallback<GlycoSearchTask>(path), Throws.Exception);
            }
            finally
            {
                File.Delete(path);
            }
        }
    }
}
