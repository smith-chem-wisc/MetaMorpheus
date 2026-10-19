using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text.RegularExpressions;
using EngineLayer;
using EngineLayer.DatabaseLoading;
using Nett;
using NUnit.Framework;
using Omics.Digestion;
using Proteomics.ProteolyticDigestion;
using Readers;
using TaskLayer;


namespace Test
{
    [TestFixture]
    public static class ProteaseTests
    {
        [Test]
        public static void TestChymotrypsinCleaveAfterLeucine()
        {
            var myTomlPath = Path.Combine(TestContext.CurrentContext.TestDirectory, @"DatabaseTests\Task1-SearchTaskconfig.toml");
            var searchTaskLoaded = Toml.ReadFile<SearchTask>(myTomlPath, MetaMorpheusTask.tomlConfig);
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"DatabaseTests\TestProtease");
            Directory.CreateDirectory(outputFolder);
            string myFile = Path.Combine(TestContext.CurrentContext.TestDirectory, @"DatabaseTests\Q9UHB6_Chym_snip.mzML");
            string myDatabase = Path.Combine(TestContext.CurrentContext.TestDirectory, @"DatabaseTests\Q9UHB6.fasta");

            var engineToml = new EverythingRunnerEngine(new List<(string, MetaMorpheusTask)> { ("SearchTOML", searchTaskLoaded) }, new List<string> { myFile }, new List<DbForTask> { new DbForTask(myDatabase, false) }, outputFolder);
            engineToml.Run();

            string psmFile = Path.Combine(outputFolder, @"SearchTOML\AllPSMs.psmtsv");

            List<PsmFromTsv> parsedPsms = SpectrumMatchTsvReader.ReadPsmTsv(psmFile, out var warnings);
            PsmFromTsv psm = parsedPsms.First();
            Assert.That(psm.BaseSeq, Is.EqualTo("TTQNQKSQDVELWEGEVVKEL")); //base sequence ends in leucine as expected
            Directory.Delete(outputFolder,true);
        }

        private static string MotifSignature(IEnumerable<DigestionMotif> motifs) => string.Join(";", motifs
            .Select(m => $"{m.InducingCleavage}|{m.PreventingCleavage}|{m.CutIndex}|{m.ExcludeFromWildcard}")
            .OrderBy(s => s, StringComparer.Ordinal));

        /// <summary>
        /// A task file that spells its protease the way an older release wrote it has to load as the enzyme it meant.
        /// The converter used the dictionary's raw indexer, which knows only today's names, so these files threw
        /// KeyNotFoundException. The expected enzyme is pinned by motif rather than by name, because the name of the
        /// proline-restricted variant itself changes in mzLib #1186.
        /// </summary>
        [Test]
        [TestCase("chymotrypsin (don't cleave before proline)", "F[P]|,W[P]|,Y[P]|,L[P]|")]
        [TestCase("trypsin (don't cleave before proline)", "K[P]|,R[P]|")]
        public static void ATaskFileWithAHistoricalProteaseNameLoads(string historicalName, string motifs)
        {
            string path = Path.Combine(TestContext.CurrentContext.WorkDirectory, nameof(ATaskFileWithAHistoricalProteaseNameLoads) + ".toml");
            Toml.WriteFile(new SearchTask(), path, MetaMorpheusTask.tomlConfig);
            string text = File.ReadAllText(path);
            Assert.That(Regex.Matches(text, "^(Specific)?Protease = \"[^\"]*\"", RegexOptions.Multiline).Count, Is.EqualTo(2),
                "premise: the saved task names its protease on exactly two lines");
            File.WriteAllText(path, Regex.Replace(text, "^((?:Specific)?Protease) = \"[^\"]*\"", $"$1 = \"{historicalName}\"", RegexOptions.Multiline));

            var digestionParams = (DigestionParams)Toml.ReadFile<SearchTask>(path, MetaMorpheusTask.tomlConfig).CommonParameters.DigestionParams;
            File.Delete(path);

            Assert.That(MotifSignature(digestionParams.Protease.DigestionMotifs),
                Is.EqualTo(MotifSignature(DigestionMotif.ParseDigestionMotifsFromString(motifs))));
            Assert.That(digestionParams.SpecificProtease, Is.SameAs(digestionParams.Protease));
        }
    }
}
