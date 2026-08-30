using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using EngineLayer;
using EngineLayer.DatabaseLoading;
using Nett;
using Proteomics;
using NUnit.Framework;
using TaskLayer;
using UsefulProteomicsDatabases;

namespace Test
{
    [TestFixture]
    public static class FastaHeaderParsingTest
    {
        // A legacy NCBI defline: the pipes are field separators, and the accession is the fourth field.
        private const string LegacyNcbiHeader =
            ">gi|16128008|ref|NP_414555.1| thr operon leader peptide [Escherichia coli str. K-12 substr. MG1655]";

        // Taken from mzLib's Test/Transcriptomics/TestData/ModomicsUnmodifiedTrimmed.fasta.
        private const string ModomicsHeader =
            ">id:1|Name:tdbR00000010|SOterm:SO:0000254|Type:tRNA|Subtype:Ala|Feature:VGC|Cellular_Localization:prokaryotic cytosol|Species:Escherichia coli";

        private const string UniProtHeader =
            ">sp|P12345|AATM_RABIT Aspartate aminotransferase, mitochondrial OS=Oryctolagus cuniculus OX=9986 GN=GOT2 PE=1 SV=2";

        private static string WriteFasta(string fileName, string header, string sequence)
        {
            string path = Path.Combine(TestContext.CurrentContext.TestDirectory, fileName);
            Directory.CreateDirectory(Path.GetDirectoryName(path));
            File.WriteAllLines(path, new[] { header, sequence });
            return path;
        }

        private static List<Protein> Load(string path, FastaHeaderParsingParameters headerParsing)
        {
            var commonParameters = new CommonParameters(fastaHeaderParsing: headerParsing);
            return DatabaseLoadingEngine.LoadProteinDb(path, true, DecoyType.None, new List<string>(), false,
                out _, out _, commonParameters).ToList();
        }

        [Test]
        public static void DefaultIsUniProtAndKeepsTheRegexesMetaMorpheusHasAlwaysPassed()
        {
            var regexes = new FastaHeaderParsingParameters().GetFieldRegexes();

            Assert.That(new FastaHeaderParsingParameters().HeaderFormat, Is.EqualTo(FastaHeaderFormat.UniProt));
            Assert.That(regexes.Accession, Is.SameAs(ProteinDbLoader.UniprotAccessionRegex));
            Assert.That(regexes.FullName, Is.SameAs(ProteinDbLoader.UniprotFullNameRegex));
            // Name has always been parsed with the full-name regex, and organism id has never been
            // parsed at all. Both are preserved so the default produces identical output.
            Assert.That(regexes.Name, Is.SameAs(ProteinDbLoader.UniprotFullNameRegex));
            Assert.That(regexes.GeneName, Is.SameAs(ProteinDbLoader.UniprotGeneNameRegex));
            Assert.That(regexes.Organism, Is.SameAs(ProteinDbLoader.UniprotOrganismRegex));
            Assert.That(regexes.OrganismId, Is.Null);
        }

        [Test]
        public static void UniProtDefaultStillParsesAUniProtHeader()
        {
            string path = WriteFasta(Path.Combine("FastaHeaderParsing", "uniprot.fasta"), UniProtHeader, "PEPTIDEK");
            var protein = Load(path, new FastaHeaderParsingParameters()).Single();

            Assert.That(protein.Accession, Is.EqualTo("P12345"));
            Assert.That(protein.Organism, Is.EqualTo("Oryctolagus cuniculus"));
            Assert.That(protein.GeneNames.Single().Item2, Is.EqualTo("GOT2"));
        }

        [Test]
        public static void UniProtDefaultGluesLegacyNcbiFieldsTogetherAndACustomRegexDoesNot()
        {
            string path = WriteFasta(Path.Combine("FastaHeaderParsing", "legacyNcbi.fasta"), LegacyNcbiHeader, "PEPTIDEK");

            var withDefault = Load(path, new FastaHeaderParsingParameters()).Single();
            Assert.That(withDefault.Accession, Is.EqualTo("16128008|ref|NP_414555.1"));
            Assert.That(withDefault.Accession, Does.Contain("|"),
                "the accession contains the character MetaMorpheus uses to separate multiple values");

            var withCustom = Load(path, new FastaHeaderParsingParameters(FastaHeaderFormat.Custom,
                accessionRegex: @">.*\|(.*)\|")).Single();
            Assert.That(withCustom.Accession, Is.EqualTo("NP_414555.1"));
            Assert.That(withCustom.Organism, Is.Null.Or.Empty, "no organism regex was given, so the field is skipped");
        }

        [Test]
        public static void NcbiPresetParsesALegacyDefline()
        {
            string path = WriteFasta(Path.Combine("FastaHeaderParsing", "ncbiPreset.fasta"), LegacyNcbiHeader, "PEPTIDEK");
            var protein = Load(path, new FastaHeaderParsingParameters(FastaHeaderFormat.Ncbi)).Single();

            Assert.That(protein.Accession, Is.EqualTo("NP_414555.1"));
            Assert.That(protein.FullName, Is.EqualTo("thr operon leader peptide"));
            Assert.That(protein.Organism, Is.EqualTo("Escherichia coli str. K-12 substr. MG1655"));
        }

        [Test]
        public static void ModomicsPresetParsesAModomicsDefline()
        {
            string path = WriteFasta(Path.Combine("FastaHeaderParsing", "modomicsPreset.fasta"), ModomicsHeader, "PEPTIDEK");
            var protein = Load(path, new FastaHeaderParsingParameters(FastaHeaderFormat.Modomics)).Single();

            // Name, not SOterm: SOterm repeats across entries, so it collides as an accession.
            Assert.That(protein.Accession, Is.EqualTo("tdbR00000010"));
            Assert.That(protein.Name, Is.EqualTo("SO:0000254"));
            Assert.That(protein.FullName, Is.EqualTo("tRNA"));
            Assert.That(protein.Organism, Is.EqualTo("Escherichia coli"));
        }

        [Test]
        public static void ModomicsHeaderIsUnusableUnderTheUniProtDefault()
        {
            string path = WriteFasta(Path.Combine("FastaHeaderParsing", "modomicsDefault.fasta"), ModomicsHeader, "PEPTIDEK");
            var protein = Load(path, new FastaHeaderParsingParameters()).Single();

            Assert.That(protein.Accession, Does.StartWith("Name:tdbR00000010|"));
            Assert.That(protein.Accession.Count(c => c == '|'), Is.GreaterThan(1));
        }

        [Test]
        public static void ATaskTomlRoundTripsTheCustomHeaderRegexes()
        {
            var task = new SearchTask
            {
                CommonParameters = new CommonParameters(fastaHeaderParsing:
                    new FastaHeaderParsingParameters(FastaHeaderFormat.Custom,
                        accessionRegex: @">.*\|(.*)\|",
                        fullNameRegex: @">[^ ]*\s(.*)$",
                        organismRegex: @"\[([^\]]+)\]\s*$"))
            };

            string tomlPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "FastaHeaderParsing", "roundTrip.toml");
            Directory.CreateDirectory(Path.GetDirectoryName(tomlPath));
            Toml.WriteFile(task, tomlPath, MetaMorpheusTask.tomlConfig);

            string written = File.ReadAllText(tomlPath);
            Assert.That(written, Does.Contain("HeaderFormat"), "the setting must be visible in the toml");

            var readBack = Toml.ReadFile<SearchTask>(tomlPath, MetaMorpheusTask.tomlConfig);
            var headerParsing = readBack.CommonParameters.FastaHeaderParsing;

            Assert.That(headerParsing.HeaderFormat, Is.EqualTo(FastaHeaderFormat.Custom));
            Assert.That(headerParsing.CustomAccessionRegex, Is.EqualTo(@">.*\|(.*)\|"));
            Assert.That(headerParsing.CustomFullNameRegex, Is.EqualTo(@">[^ ]*\s(.*)$"));
            Assert.That(headerParsing.CustomOrganismRegex, Is.EqualTo(@"\[([^\]]+)\]\s*$"));
            Assert.That(headerParsing.CustomNameRegex, Is.Empty);
            Assert.That(headerParsing.CustomGeneNameRegex, Is.Empty);
            Assert.That(headerParsing.CustomOrganismIdRegex, Is.Empty);
        }

        [Test]
        public static void ADefaultTaskTomlRoundTripsToTheUniProtDefault()
        {
            string tomlPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "FastaHeaderParsing", "defaultRoundTrip.toml");
            Directory.CreateDirectory(Path.GetDirectoryName(tomlPath));
            Toml.WriteFile(new SearchTask(), tomlPath, MetaMorpheusTask.tomlConfig);

            var readBack = Toml.ReadFile<SearchTask>(tomlPath, MetaMorpheusTask.tomlConfig);
            Assert.That(readBack.CommonParameters.FastaHeaderParsing.HeaderFormat, Is.EqualTo(FastaHeaderFormat.UniProt));
            Assert.That(readBack.CommonParameters.FastaHeaderParsing.CustomAccessionRegex, Is.Empty);
        }

        [Test]
        public static void ATomlWrittenBeforeThisSettingExistedStillReadsAsUniProt()
        {
            string tomlPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "FastaHeaderParsing", "legacy.toml");
            Directory.CreateDirectory(Path.GetDirectoryName(tomlPath));
            Toml.WriteFile(new SearchTask(), tomlPath, MetaMorpheusTask.tomlConfig);

            // Strip the new table to imitate a toml saved by an earlier version.
            var kept = new List<string>();
            bool inSection = false;
            foreach (string line in File.ReadAllLines(tomlPath))
            {
                if (line.StartsWith("["))
                    inSection = line.Trim() == "[CommonParameters.FastaHeaderParsing]";
                if (!inSection)
                    kept.Add(line);
            }
            File.WriteAllLines(tomlPath, kept);
            Assert.That(File.ReadAllText(tomlPath), Does.Not.Contain("FastaHeaderParsing"));

            var readBack = Toml.ReadFile<SearchTask>(tomlPath, MetaMorpheusTask.tomlConfig);
            Assert.That(readBack.CommonParameters.FastaHeaderParsing, Is.Not.Null);
            Assert.That(readBack.CommonParameters.FastaHeaderParsing.HeaderFormat, Is.EqualTo(FastaHeaderFormat.UniProt));
        }

        [Test]
        public static void PresetsNeedNoValidationAndAMalformedCustomRegexFailsIt()
        {
            foreach (FastaHeaderFormat format in Enum.GetValues<FastaHeaderFormat>().Where(f => f != FastaHeaderFormat.Custom))
                Assert.That(new FastaHeaderParsingParameters(format).Validate(out _), Is.True, format.ToString());

            var unbalanced = new FastaHeaderParsingParameters(FastaHeaderFormat.Custom, accessionRegex: @">(.*\|");
            Assert.That(unbalanced.Validate(out var errors), Is.False);
            Assert.That(errors.Single(), Does.Contain("accession").And.Contain("not usable"));
        }

        [Test]
        public static void ACustomFormatWithoutAnAccessionRegexIsRefused()
        {
            var noAccession = new FastaHeaderParsingParameters(FastaHeaderFormat.Custom, organismRegex: @"OS=(\S+)");
            Assert.That(noAccession.Validate(out var errors), Is.False);
            Assert.That(errors.Single(), Does.Contain("requires an accession regex"));
        }

        [Test]
        public static void ARegexWithNoCaptureGroupIsRefused()
        {
            var noGroup = new FastaHeaderParsingParameters(FastaHeaderFormat.Custom, accessionRegex: @">\S+");
            Assert.That(noGroup.Validate(out var errors), Is.False);
            Assert.That(errors.Single(), Does.Contain("no capture group"));
        }

        [Test]
        public static void ACatastrophicallyBacktrackingRegexIsRefused()
        {
            var runaway = new FastaHeaderParsingParameters(FastaHeaderFormat.Custom,
                accessionRegex: @"^(\s*\w+|\W)*(z9q)$");
            Assert.That(runaway.Validate(out var errors), Is.False);
            Assert.That(errors.Single(), Does.Contain("longer than"));
        }

        [Test]
        public static void AMalformedRegexStopsTheRunWithAWarningAndStartsNoEngine()
        {
            var startedEngines = new List<string>();
            void RecordEngine(object sender, SingleEngineEventArgs e) => startedEngines.Add(e.MyEngine.GetType().Name);

            string fastaPath = WriteFasta(Path.Combine("FastaHeaderParsing", "gateDb.fasta"), UniProtHeader, "PEPTIDEK");
            string mzmlPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "PrunedDbSpectra.mzml");
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestFastaHeaderRegexGate");

            var task = new SearchTask
            {
                CommonParameters = new CommonParameters(fastaHeaderParsing:
                    new FastaHeaderParsingParameters(FastaHeaderFormat.Custom, accessionRegex: @">(.*\|"))
            };

            MetaMorpheusEngine.StartingSingleEngineHander += RecordEngine;
            EverythingRunnerEngine engine;
            try
            {
                engine = new EverythingRunnerEngine(
                    new List<(string, MetaMorpheusTask)> { ("Task1-SearchTask", task) },
                    new List<string> { mzmlPath },
                    new List<DbForTask> { new DbForTask(fastaPath, false) },
                    outputFolder);
                engine.Run();
            }
            finally
            {
                MetaMorpheusEngine.StartingSingleEngineHander -= RecordEngine;
            }

            Assert.That(startedEngines, Is.Empty, "no engine should start: " + string.Join(", ", startedEngines));
            Assert.That(engine.Warnings, Is.Not.Empty);
            Assert.That(engine.Warnings.Single(), Does.StartWith("Cannot proceed. Task1-SearchTask:")
                .And.Contain("not usable"));
        }
    }
}
