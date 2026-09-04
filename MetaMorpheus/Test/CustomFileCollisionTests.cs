using EngineLayer;
using NUnit.Framework;
using Proteomics.ProteolyticDigestion;
using System;
using System.IO;
using System.Linq;
using System.Reflection;

namespace Test
{
    /// <summary>
    /// mzLib refuses to let a custom digestion agent shadow one of its own and reports it through
    /// CustomDigestionAgentLoadResult.Skipped rather than throwing, so that the caller can tell the user.
    /// Nothing consumed that before, so a user who named a custom protease after a built-in got silence
    /// and a protease that was not theirs. These cover the reporting, and the embedded-resource half of
    /// the template builder.
    /// </summary>
    [TestFixture]
    public static class CustomFileCollisionTests
    {
        private static string _backup;
        private static bool _hadFile;

        [SetUp]
        public static void SetUp()
        {
            _hadFile = File.Exists(GlobalVariables.CustomProteasePath);
            if (_hadFile)
            {
                _backup = GlobalVariables.CustomProteasePath + ".bak";
                File.Copy(GlobalVariables.CustomProteasePath, _backup, true);
            }
        }

        [TearDown]
        public static void TearDown()
        {
            if (_hadFile && _backup != null && File.Exists(_backup))
            {
                File.Copy(_backup, GlobalVariables.CustomProteasePath, true);
                File.Delete(_backup);
            }
            else if (!_hadFile && File.Exists(GlobalVariables.CustomProteasePath))
            {
                File.Delete(GlobalVariables.CustomProteasePath);
            }

            _backup = null;
            GlobalVariables.SetUpGlobalVariables();
        }

        /// <summary>
        /// A custom protease that reuses a built-in's name is refused by mzLib -- silently, until now. The
        /// user has to be told, because the protease they get is not the one they wrote.
        /// </summary>
        [Test]
        public static void ACustomProteaseThatShadowsABuiltInIsReportedToTheUser()
        {
            string shadowed = ProteaseDictionary.Dictionary.Keys.First();

            File.WriteAllLines(GlobalVariables.CustomProteasePath, new[]
            {
                "Name\tMotif\tSpecificity\tPSI-MS Accession\tPSI-MS Name\tCleavage Modification",
                $"{shadowed}\tX|\tfull\t\t\t",
            });

            GlobalVariables.SetUpGlobalVariables();

            Assert.Multiple(() =>
            {
                Assert.That(GlobalVariables.StartupWarnings.Any(w => w.Contains(shadowed)),
                    Is.True, "the shadowed name has to be named: " + string.Join(" | ", GlobalVariables.StartupWarnings));
                Assert.That(GlobalVariables.StartupWarnings.Any(w => w.Contains("proteases_custom.tsv")),
                    Is.True, "and the file to fix it in");
                // the built-in definition is the one that survives, which is the reason to warn at all
                Assert.That(ProteaseDictionary.Dictionary[shadowed].DigestionMotifs.Any(m => m.InducingCleavage == "X"),
                    Is.False, "the custom entry must not have replaced the built-in");
            });
        }

        /// <summary>
        /// The counterpart: a custom protease with a name of its own is not reported. Without this, the
        /// assertion above would be satisfied by warning unconditionally, which is worse than silence.
        /// </summary>
        [Test]
        public static void ACustomProteaseWithItsOwnNameIsNotReported()
        {
            File.WriteAllLines(GlobalVariables.CustomProteasePath, new[]
            {
                "Name\tMotif\tSpecificity\tPSI-MS Accession\tPSI-MS Name\tCleavage Modification",
                "NotAShadowingName_ForTest\tX|\tfull\t\t\t",
            });

            GlobalVariables.SetUpGlobalVariables();

            Assert.That(GlobalVariables.StartupWarnings.Any(w => w.Contains("could not be") || w.Contains("ignored")),
                Is.False, "nothing to report: " + string.Join(" | ", GlobalVariables.StartupWarnings));
        }

        /// <summary>
        /// The seeded custom-protease template is derived from mzLib's embedded proteases.tsv, which is
        /// where the "banner plus header, no data rows" shape this whole policy copies comes from.
        /// </summary>
        [Test]
        public static void TheProteaseTemplateComesFromMzLibsEmbeddedFile_AndCarriesNoProteases()
        {
            string template = CustomDataFile.BannerAndHeaderFrom(typeof(ProteaseDictionary).Assembly,
                "Proteomics.ProteolyticDigestion.proteases.tsv", "Name\t");

            var lines = template.Split(new[] { "\r\n", "\n" }, StringSplitOptions.RemoveEmptyEntries);
            Assert.Multiple(() =>
            {
                Assert.That(lines.Last(), Does.StartWith("Name\t"), "the header is the last line kept");
                Assert.That(lines.Count(l => l.StartsWith("#")), Is.GreaterThan(0),
                    "the banner is the part that makes the template worth seeding");
                Assert.That(lines.Any(l => l.StartsWith("Arg-C\t")), Is.False,
                    "a shipped protease in the custom file would look like the user's own");
            });
        }

        [Test]
        public static void AMissingEmbeddedTemplateIsReportedByName()
        {
            var ex = Assert.Throws<MetaMorpheusException>(() =>
                CustomDataFile.BannerAndHeaderFrom(typeof(ProteaseDictionary).Assembly,
                    "Proteomics.ProteolyticDigestion.this.resource.does.not.exist.tsv", "Name\t"));

            Assert.That(ex!.Message, Does.Contain("this.resource.does.not.exist.tsv"));
        }
    }
}
