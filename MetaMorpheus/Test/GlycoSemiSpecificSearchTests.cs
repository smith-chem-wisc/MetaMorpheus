using EngineLayer;
using EngineLayer.DatabaseLoading;
using Nett;
using NUnit.Framework;
using Omics.Digestion;
using Proteomics.ProteolyticDigestion;
using Readers;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using TaskLayer;

namespace Test
{
    /// <summary>
    /// A glyco search with trypsin and <c>SearchModeType = Semi</c> must identify glycopeptides whose N-terminus is not
    /// tryptic.
    /// </summary>
    /// <remarks>
    /// <para><b>The bug this guards.</b> Until mzLib #1303, SearchModeType Semi with FragmentationTerminus Both (what a glyco
    /// task always has) made mzLib return C-terminal "seed" peptides meant for the non-specific search engine, which trims
    /// them after the search. The glyco search does no trimming: it takes precursor mass minus peptide mass as the glycan
    /// mass, which is meaningless for a seed. So a semi-specific glyco search silently lost most of its identifications,
    /// with no error. This is the configuration that replaces the removed "semi-trypsin" protease (see
    /// <see cref="LegacySemiTrypsinSettingsTests"/>), so it has to work end to end.</para>
    /// <para><b>The data.</b> One EThcD spectrum of leukosialin (CD43, P16150) from a StcE + trypsin mucin digest,
    /// identified by <c>TestNOGlyco.NandO_GlycoSearch_onlyO</c> as TTGSLEPSSGASGPQVSSVK with three O-glycans. In
    /// leukosialin that peptide follows a Met (...GPPVTM|TTGS...), so its N-terminus was made by StcE, not trypsin: it is
    /// semi-tryptic. A fully tryptic search cannot find it; a semi-tryptic one must.</para>
    /// </remarks>
    [TestFixture]
    public static class GlycoSemiSpecificSearchTests
    {
        private const string SemiTrypticGlycopeptide = "TTGSLEPSSGASGPQVSSVK";

        private static string GlycoTestData(string fileName) => Path.Combine(TestContext.CurrentContext.TestDirectory, "GlycoTestData", fileName);

        /// <summary>
        /// Runs the O-glyco search from GlycoSearchTaskconfigOGlycoTest_Run.toml on the leukosialin spectrum, with trypsin
        /// and the given search mode in place of StcE-trypsin, and returns the base sequences of the PSMs it reports.
        /// </summary>
        private static List<string> RunOGlycoSearchWithTrypsin(CleavageSpecificity searchModeType, string testName)
        {
            string settingsText = File.ReadAllText(GlycoTestData("GlycoSearchTaskconfigOGlycoTest_Run.toml"));
            Assert.That(settingsText, Does.Contain("Protease = \"StcE-trypsin\"").And.Contain("SearchModeType = \"Full\"").And.Contain("FragmentationTerminus = \"Both\""),
                "premise: the settings file this test edits");
            settingsText = settingsText
                .Replace("Protease = \"StcE-trypsin\"", "Protease = \"trypsin\"") // also SpecificProtease
                .Replace("SearchModeType = \"Full\"", $"SearchModeType = \"{searchModeType}\"");

            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, $"{nameof(GlycoSemiSpecificSearchTests)}_{testName}");
            string settingsPath = outputFolder + ".toml";
            try
            {
                if (Directory.Exists(outputFolder))
                    Directory.Delete(outputFolder, true);
                Directory.CreateDirectory(outputFolder);
                File.WriteAllText(settingsPath, settingsText);
                var task = MetaMorpheusTask.ReadTaskTomlWithLowResFallback<GlycoSearchTask>(settingsPath);

                var digestionParams = (DigestionParams)task.CommonParameters.DigestionParams;
                Assert.That(digestionParams.Protease.Name, Is.EqualTo("trypsin"), "premise: trypsin");
                Assert.That(digestionParams.SearchModeType, Is.EqualTo(searchModeType), "premise: the search mode under test");

                new EverythingRunnerEngine(new List<(string, MetaMorpheusTask)> { ("Task", task) },
                    new List<string> { GlycoTestData("2019_09_16_StcEmix_35trig_EThcD25_rep1_9906.mgf") },
                    new List<DbForTask> { new DbForTask(GlycoTestData("P16150.fasta"), false) },
                    outputFolder).Run();

                string resultsPath = Path.Combine(outputFolder, "Task", "results.txt");
                Assert.That(File.Exists(resultsPath), Is.True, "the task must run");
                Assert.That(File.ReadAllText(resultsPath), Does.Not.Contain("Exception"), "the task must finish without an error");

                // WriteGlycoFile.WritePsmGlycoToTsv writes no file at all when nothing was identified
                string allPsmsPath = Path.Combine(outputFolder, "Task", "AllPSMs.psmtsv");
                if (!File.Exists(allPsmsPath))
                    return new List<string>();
                List<PsmFromTsv> psms = SpectrumMatchTsvReader.ReadPsmTsv(allPsmsPath, out List<string> errors);
                Assert.That(errors, Is.Empty);
                return psms.Select(p => p.BaseSeq).ToList();
            }
            finally
            {
                if (Directory.Exists(outputFolder))
                    Directory.Delete(outputFolder, true);
                if (File.Exists(settingsPath))
                    File.Delete(settingsPath);
            }
        }

        /// <summary>
        /// Trypsin with SearchModeType Semi identifies the semi-tryptic leukosialin glycopeptide. Before mzLib #1303 the
        /// peptide was never a candidate, because digestion returned seeds.
        /// </summary>
        [Test]
        public static void OGlycoSearch_TrypsinWithSemiSearchMode_IdentifiesTheSemiTrypticGlycopeptide()
        {
            List<string> identified = RunOGlycoSearchWithTrypsin(CleavageSpecificity.Semi, "Semi");
            Assert.That(identified, Does.Contain(SemiTrypticGlycopeptide),
                "a semi-specific glyco search must consider peptides with one non-tryptic terminus");
        }

        /// <summary>
        /// The control: the same search fully tryptic does not identify it, so the test above really depends on
        /// semi-specific digestion and not on something else in the setup.
        /// </summary>
        [Test]
        public static void OGlycoSearch_TrypsinFullySpecific_DoesNotIdentifyTheSemiTrypticGlycopeptide()
        {
            List<string> identified = RunOGlycoSearchWithTrypsin(CleavageSpecificity.Full, "Full");
            Assert.That(identified, Does.Not.Contain(SemiTrypticGlycopeptide),
                "premise: the peptide's N-terminus follows a Met, so fully tryptic digestion cannot produce it");
        }
    }
}
