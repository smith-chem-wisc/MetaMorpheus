using System;
using System.Collections.Generic;
using System.IO;
using System.IO.Compression;
using System.Linq;
using System.Text.RegularExpressions;
using EngineLayer;
using EngineLayer.DatabaseLoading;
using MassSpectrometry;
using MzLibUtil;
using Nett;
using Proteomics;
using NUnit.Framework;
using TaskLayer;
using Transcriptomics.Digestion;
using UsefulProteomicsDatabases;

namespace Test
{
    [TestFixture]
    public static class FastaHeaderParsingTest
    {
        // Taken from mzLib's Test/Transcriptomics/TestData/ModomicsUnmodifiedTrimmed.fasta.
        private const string ModomicsHeader =
            ">id:1|Name:tdbR00000010|SOterm:SO:0000254|Type:tRNA|Subtype:Ala|Feature:VGC|Cellular_Localization:prokaryotic cytosol|Species:Escherichia coli";

        // Gencode: exactly seven pipes, which is what DetectFastaHeaderFormat keys on.
        private const string GencodeHeader =
            ">ENSMUSP00000034487.2|ENSMUST00000034487.3|ENSMUSG00000018620.3|OTTMUSG00000062304.1|OTTMUST00000151979.1|Mmp20-201|Mmp20|482";

        private const string EnsemblHeader =
            ">ENSP00000493376.2 pep chromosome:GRCh38:1:65565:71585:1 gene:ENSG00000186092.7 transcript:ENST00000641515.2 gene_biotype:protein_coding transcript_biotype:protein_coding gene_symbol:OR4F5 description:olfactory receptor family 4 subfamily F member 5 [Source:HGNC Symbol;Acc:HGNC:14825]";

        private const string UniProtHeader =
            ">sp|P12345|AATM_RABIT Aspartate aminotransferase, mitochondrial OS=Oryctolagus cuniculus OX=9986 GN=GOT2 PE=1 SV=2";

        // NCBI has emitted pipe-less deflines since GI numbers were retired.
        private const string ModernNcbiHeader = ">NP_414555.1 thr operon leader peptide [Escherichia coli str. K-12 substr. MG1655]";

        [TearDown]
        public static void ResetAnalyteType()
        {
            // An RNA run below sets this process-wide, and LoadBioPolymers branches on it.
            GlobalVariables.AnalyteType = AnalyteType.Peptide;
        }

        private static string TestPath(params string[] parts) =>
            Path.Combine(new[] { TestContext.CurrentContext.TestDirectory, "FastaHeaderParsing" }.Concat(parts).ToArray());

        private static string WriteFasta(string fileName, params string[] lines)
        {
            string path = TestPath(fileName);
            Directory.CreateDirectory(Path.GetDirectoryName(path));
            File.WriteAllLines(path, lines);
            return path;
        }

        private static List<Protein> Load(string path, FastaHeaderParsingParameters headerParsing)
        {
            var commonParameters = new CommonParameters(fastaHeaderParsing: headerParsing);
            return DatabaseLoadingEngine.LoadProteinDb(path, true, DecoyType.None, new List<string>(), false,
                out _, out _, commonParameters).ToList();
        }

        private static EverythingRunnerEngine RunAndRecordEngines(List<(string, MetaMorpheusTask)> tasks, string spectraPath,
            List<DbForTask> databases, string outputFolderName, List<string> startedEngines)
        {
            void RecordEngine(object sender, SingleEngineEventArgs e) => startedEngines.Add(e.MyEngine.GetType().Name);
            MetaMorpheusEngine.StartingSingleEngineHander += RecordEngine;
            try
            {
                var engine = new EverythingRunnerEngine(tasks, new List<string> { spectraPath }, databases,
                    Path.Combine(TestContext.CurrentContext.TestDirectory, outputFolderName));
                engine.Run();
                return engine;
            }
            finally
            {
                MetaMorpheusEngine.StartingSingleEngineHander -= RecordEngine;
            }
        }

        private static string ProteinSpectra => Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "PrunedDbSpectra.mzml");

        #region Presets

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
            string path = WriteFasta("uniprot.fasta", UniProtHeader, "PEPTIDEK");
            var protein = Load(path, new FastaHeaderParsingParameters()).Single();

            Assert.That(protein.Accession, Is.EqualTo("P12345"));
            Assert.That(protein.Organism, Is.EqualTo("Oryctolagus cuniculus"));
            Assert.That(protein.GeneNames.Single().Item2, Is.EqualTo("GOT2"));
        }

        [TestCase(FastaHeaderFormat.Ensembl, EnsemblHeader, "ENSP00000493376.2", "ENSG00000186092.7")]
        [TestCase(FastaHeaderFormat.Gencode, GencodeHeader, "ENSMUSP00000034487.2", "Mmp20")]
        public static void EachPresetParsesItsOwnHeader(FastaHeaderFormat format, string header, string accession, string gene)
        {
            string path = WriteFasta($"preset{format}.fasta", header, "PEPTIDEK");
            var protein = Load(path, new FastaHeaderParsingParameters(format)).Single();

            Assert.That(protein.Accession, Is.EqualTo(accession));
            Assert.That(protein.GeneNames.Single().Item2, Is.EqualTo(gene));
        }

        [Test]
        public static void ACustomAccessionRegexReachesTheLoader()
        {
            string path = WriteFasta("customNcbi.fasta", ModernNcbiHeader, "PEPTIDEK");

            // The default has no pipes to match here, so mzLib falls back to the whole defline.
            var withDefault = Load(path, new FastaHeaderParsingParameters()).Single();
            Assert.That(withDefault.Accession, Is.Not.EqualTo("NP_414555.1"));

            var withCustom = Load(path, new FastaHeaderParsingParameters(FastaHeaderFormat.Custom,
                accessionRegex: @"^>(\S+)", organismRegex: @"\[([^\[\]]+)\]\s*$")).Single();
            Assert.That(withCustom.Accession, Is.EqualTo("NP_414555.1"));
            Assert.That(withCustom.Organism, Is.EqualTo("Escherichia coli str. K-12 substr. MG1655"));
            Assert.That(withCustom.GeneNames, Is.Empty, "no gene regex was given, so the field is skipped");
        }

        [Test]
        public static void ModomicsHeaderIsUnusableUnderTheUniProtDefault()
        {
            // An invariant rather than a value: mzLib's UniProt accession regex has changed under this
            // test once already, and neither version yields the Modomics id.
            string path = WriteFasta("modomicsDefault.fasta", ModomicsHeader, "PEPTIDEK");
            var protein = Load(path, new FastaHeaderParsingParameters()).Single();

            Assert.That(protein.Accession, Is.Not.EqualTo("tdbR00000010"));
        }

        [Test]
        public static void AutoDetectsGencodeWhereTheUniProtPresetManglesIt()
        {
            string path = WriteFasta("gencodeAuto.fasta", GencodeHeader, "PEPTIDEK");

            var mangled = Load(path, new FastaHeaderParsingParameters(FastaHeaderFormat.UniProt)).Single();
            Assert.That(mangled.Accession, Is.Not.EqualTo("ENSMUSP00000034487.2"));

            var detected = Load(path, new FastaHeaderParsingParameters(FastaHeaderFormat.Auto)).Single();
            Assert.That(detected.Accession, Is.EqualTo("ENSMUSP00000034487.2"));
            Assert.That(detected.GeneNames.Single().Item2, Is.EqualTo("Mmp20"));
        }

        [Test]
        public static void AutoChangesNameAndTaxonomyIdOnAUniProtHeader()
        {
            string path = WriteFasta("uniprotAuto.fasta", UniProtHeader, "PEPTIDEK");

            var preset = Load(path, new FastaHeaderParsingParameters(FastaHeaderFormat.UniProt)).Single();
            var auto = Load(path, new FastaHeaderParsingParameters(FastaHeaderFormat.Auto)).Single();

            // Same accession, but Auto is not a synonym for the UniProt preset: mzLib's detection
            // uses UniprotNameRegex for Name and defaults an organism id regex, neither of which
            // MetaMorpheus has ever passed. Pinned so the difference is documented, not discovered.
            Assert.That(auto.Accession, Is.EqualTo(preset.Accession));
            Assert.That(preset.Name, Is.EqualTo("Aspartate aminotransferase, mitochondrial"));
            Assert.That(auto.Name, Is.EqualTo("AATM_RABIT"));
            Assert.That(preset.NcbiTaxonomyId, Is.Null.Or.Empty);
            Assert.That(auto.NcbiTaxonomyId, Is.EqualTo("9986"));
        }

        [Test]
        public static void EveryCopyOfCommonParametersKeepsTheHeaderSetting()
        {
            // Each of these rebuilds CommonParameters from an argument list kept by hand; one that drops
            // the setting silently reverts a partitioned, semi-specific or per-file search to UniProt.
            var headerParsing = new FastaHeaderParsingParameters(FastaHeaderFormat.Gencode);
            var original = new CommonParameters(fastaHeaderParsing: headerParsing);

            Assert.That(original.CloneWithNewTotalPartitions(3).FastaHeaderParsing, Is.SameAs(headerParsing));
            Assert.That(original.CloneWithNewMaxThreadsToUsePerFile(1).FastaHeaderParsing, Is.SameAs(headerParsing));
            Assert.That(original.CloneWithNewTerminus().FastaHeaderParsing, Is.SameAs(headerParsing));
            Assert.That(MetaMorpheusTask.SetAllFileSpecificCommonParams(original, new FileSpecificParameters { MinPeptideLength = 6 }, "a.mzML")
                .FastaHeaderParsing, Is.SameAs(headerParsing));
        }

        #endregion

        #region TOML

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

            string tomlPath = TestPath("roundTrip.toml");
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
            string tomlPath = TestPath("defaultRoundTrip.toml");
            Directory.CreateDirectory(Path.GetDirectoryName(tomlPath));
            Toml.WriteFile(new SearchTask(), tomlPath, MetaMorpheusTask.tomlConfig);

            var readBack = Toml.ReadFile<SearchTask>(tomlPath, MetaMorpheusTask.tomlConfig);
            Assert.That(readBack.CommonParameters.FastaHeaderParsing.HeaderFormat, Is.EqualTo(FastaHeaderFormat.UniProt));
            Assert.That(readBack.CommonParameters.FastaHeaderParsing.CustomAccessionRegex, Is.Empty);
        }

        [Test]
        public static void ATomlWrittenBeforeThisSettingExistedStillReadsAsUniProt()
        {
            string tomlPath = TestPath("legacy.toml");
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

        // "7" because Enum.TryParse accepts any number; "Ncbi" because an earlier revision of this branch had it.
        [TestCase("Uniprot2")]
        [TestCase("7")]
        [TestCase("Ncbi")]
        public static void AnUnrecognizedHeaderFormatNamesTheValidValues(string badValue)
        {
            string tomlPath = TestPath($"badFormatName{badValue}.toml");
            Directory.CreateDirectory(Path.GetDirectoryName(tomlPath));
            Toml.WriteFile(new SearchTask(), tomlPath, MetaMorpheusTask.tomlConfig);

            string original = File.ReadAllText(tomlPath);
            string corrupted = Regex.Replace(original, @"HeaderFormat\s*=\s*""[^""]*""",
                $"HeaderFormat = \"{badValue}\"");
            Assert.That(corrupted, Is.Not.EqualTo(original), "the written toml must carry a HeaderFormat to corrupt");
            File.WriteAllText(tomlPath, corrupted);

            // Catch, not Throws: Throws<T> is exact-type and Nett may wrap ours.
            var thrown = Assert.Catch<Exception>(
                () => Toml.ReadFile<SearchTask>(tomlPath, MetaMorpheusTask.tomlConfig),
                "an unknown format must not read as a silent default")!;

            // ToString() so the assertion holds whether or not Nett wraps ours.
            Assert.That(thrown.ToString(), Does.Contain($"'{badValue}'").And.Contain("Valid values are"));
        }

        #endregion

        #region Configuration checks, before any task runs

        [TestCase(FastaHeaderFormat.UniProt)]
        [TestCase(FastaHeaderFormat.Ensembl)]
        [TestCase(FastaHeaderFormat.Gencode)]
        public static void APresetNeedsNoConfiguration(FastaHeaderFormat preset)
        {
            Assert.That(new FastaHeaderParsingParameters(preset).ValidateConfiguration(out var errors), Is.True);
            Assert.That(errors, Is.Empty);
        }

        [Test]
        public static void AMalformedCustomRegexFailsConfiguration()
        {
            var unbalanced = new FastaHeaderParsingParameters(FastaHeaderFormat.Custom, accessionRegex: @">(.*\|");
            Assert.That(unbalanced.ValidateConfiguration(out var errors), Is.False);
            Assert.That(errors.Single(), Does.Contain("accession").And.Contain("not usable"));
        }

        [Test]
        public static void ACustomFormatWithoutAnAccessionRegexIsRefused()
        {
            // With no accession regex mzLib auto-detects and overwrites the other five regexes with the
            // detected preset, so this organism pattern would be silently discarded, not used.
            var noAccession = new FastaHeaderParsingParameters(FastaHeaderFormat.Custom, organismRegex: @"OS=(\S+)");
            Assert.That(noAccession.ValidateConfiguration(out var errors), Is.False);
            Assert.That(errors.Single(), Does.Contain("requires an accession regex"));
        }

        [Test]
        public static void ARegexWithNoCaptureGroupIsRefused()
        {
            var noGroup = new FastaHeaderParsingParameters(FastaHeaderFormat.Custom, accessionRegex: @">\S+");
            Assert.That(noGroup.ValidateConfiguration(out var errors), Is.False);
            Assert.That(errors.Single(), Does.Contain("no capture group"));
        }

        [Test]
        public static void ACatastrophicallyBacktrackingRegexIsRefused()
        {
            var runaway = new FastaHeaderParsingParameters(FastaHeaderFormat.Custom,
                accessionRegex: @"^(\s*\w+|\W)*(z9q)$");
            Assert.That(runaway.ValidateConfiguration(out var errors), Is.False);
            Assert.That(errors.Single(), Does.Contain("longer than"));
        }

        [Test]
        public static void AMalformedRegexInALaterTaskStopsTheRunBeforeTheFirstTaskStarts()
        {
            string fastaPath = WriteFasta("gateDb.fasta", UniProtHeader, "PEPTIDEK");
            var valid = new SearchTask();
            var malformed = new SearchTask
            {
                CommonParameters = new CommonParameters(fastaHeaderParsing:
                    new FastaHeaderParsingParameters(FastaHeaderFormat.Custom, accessionRegex: @">(.*\|"))
            };

            var startedEngines = new List<string>();
            var engine = RunAndRecordEngines(
                new List<(string, MetaMorpheusTask)> { ("Task1-SearchTask", valid), ("Task2-SearchTask", malformed) },
                ProteinSpectra, new List<DbForTask> { new DbForTask(fastaPath, false) }, "TestFastaHeaderRegexGate", startedEngines);

            Assert.That(startedEngines, Is.Empty, "Task1 must not run: " + string.Join(", ", startedEngines));
            Assert.That(engine.Warnings.Single(), Does.StartWith("Cannot proceed. Task2-SearchTask:").And.Contain("not usable"));
        }

        #endregion

        #region Header checks, per database

        [Test]
        public static void TheDefaultIsNeverRefusedOnItsHeaders()
        {
            // Plain deflines have always loaded under the default, as whole-defline accessions.
            var ok = new FastaHeaderParsingParameters().ValidateAgainstHeaders("db.fasta", new[] { ">rat_reductase" },
                out var errors, out var warnings);

            Assert.That(ok, Is.True);
            Assert.That(errors, Is.Empty);
            Assert.That(warnings, Is.Empty);
        }

        [Test]
        public static void APresetThatMatchesNoSampledHeaderIsRefused()
        {
            var gencode = new FastaHeaderParsingParameters(FastaHeaderFormat.Gencode);

            Assert.That(gencode.ValidateAgainstHeaders("ncbi.fasta", new[] { ModernNcbiHeader }, out var errors, out _), Is.False);
            Assert.That(errors.Single(), Does.Contain("Gencode accession regex matches none of the first 1 headers of ncbi.fasta"));

            Assert.That(gencode.ValidateAgainstHeaders("gencode.fasta", new[] { GencodeHeader }, out errors, out _), Is.True);
        }

        [Test]
        public static void ACustomRegexThatMatchesNoSampledHeaderIsRefused()
        {
            var custom = new FastaHeaderParsingParameters(FastaHeaderFormat.Custom, accessionRegex: @"^>gi\|(\d+)");

            Assert.That(custom.ValidateConfiguration(out _), Is.True, "nothing is wrong with the pattern itself");
            Assert.That(custom.ValidateAgainstHeaders("ncbi.fasta", new[] { ModernNcbiHeader }, out var errors, out _), Is.False);
            Assert.That(errors.Single(), Does.Contain("custom accession regex matches none"));
        }

        [Test]
        public static void AnEmptyCaptureCountsAsNoMatch()
        {
            // mzLib falls back to the whole defline on an empty capture exactly as on a miss.
            var custom = new FastaHeaderParsingParameters(FastaHeaderFormat.Custom, accessionRegex: @"^>(\d*)");

            Assert.That(custom.ValidateAgainstHeaders("ncbi.fasta", new[] { ModernNcbiHeader }, out var errors, out _), Is.False);
            Assert.That(errors.Single(), Does.Contain("matches none"));
        }

        [Test]
        public static void ACustomRegexIsTimedOnTheDatabaseHeadersNotOnlyTheCannedOnes()
        {
            // Fails at the second character of every canned defline, so it is fast on those, and
            // backtracks catastrophically on the header supplied below.
            var parameters = new FastaHeaderParsingParameters(FastaHeaderFormat.Custom, accessionRegex: @"^>(a+)+X$");
            Assert.That(parameters.ValidateConfiguration(out _), Is.True, "the canned headers alone cannot show this pattern is slow");

            string ownHeader = ">" + new string('a', 40);
            Assert.That(parameters.ValidateAgainstHeaders("db.fasta", new[] { ownHeader }, out var errors, out _), Is.False);
            Assert.That(errors.Single(), Does.Contain("longer than").And.Contain("db.fasta"));
        }

        [Test]
        public static void AutoPassesARecognizableHeaderAndRefusesAnUnrecognizableOne()
        {
            var auto = new FastaHeaderParsingParameters(FastaHeaderFormat.Auto);

            Assert.That(auto.ValidateAgainstHeaders("gencode.fasta", new[] { GencodeHeader }, out var errors, out _), Is.True);
            Assert.That(errors, Is.Empty);

            Assert.That(auto.ValidateAgainstHeaders("plain.fasta", new[] { ">rat_reductase" }, out errors, out _), Is.False);
            Assert.That(errors.Single(), Does.Contain("could not be detected automatically"));
        }

        [Test]
        public static void AutoRefusesAHeaderItDetectsButWhoseAccessionItCannotParse()
        {
            // Detected as UniProt from ">sp|", but the UniProt accession needs a second pipe.
            var auto = new FastaHeaderParsingParameters(FastaHeaderFormat.Auto);

            Assert.That(auto.ValidateAgainstHeaders("onePipe.fasta", new[] { ">sp|P12345 Aspartate aminotransferase OS=Oryctolagus cuniculus" },
                out var errors, out _), Is.False);
            Assert.That(errors.Single(), Does.Contain("Auto accession regex matches none"));
        }

        [TestCase(FastaHeaderFormat.Auto, false)]
        [TestCase(FastaHeaderFormat.Custom, false)]
        [TestCase(FastaHeaderFormat.Gencode, true)]
        public static void ADatabaseWithNoReadableHeaderIsRefusedOnlyWhereTheFormatDependsOnIt(FastaHeaderFormat format, bool passes)
        {
            var parameters = new FastaHeaderParsingParameters(format, accessionRegex: @"^>(\S+)");

            Assert.That(parameters.ValidateAgainstHeaders("empty.fasta", null, out var errors, out _), Is.EqualTo(passes));
            Assert.That(parameters.ValidateAgainstHeaders("empty.fasta", new List<string>(), out _, out _), Is.EqualTo(passes));
            if (!passes)
                Assert.That(errors.Single(), Does.Contain("No FASTA header could be read from empty.fasta"));
        }

        [Test]
        public static void APresetOnADatabaseOfAnotherFormatWarnsButProceeds()
        {
            // The contaminant case: one task-wide setting also applies to a UniProt-style cRAP FASTA.
            var gencode = new FastaHeaderParsingParameters(FastaHeaderFormat.Gencode);

            Assert.That(gencode.ValidateAgainstHeaders("crap.fasta", new[] { ">sp|P00761|TRYP_PIG Trypsin OS=Sus scrofa OX=9823 PE=1 SV=1" },
                out var errors, out var warnings), Is.True);
            Assert.That(errors, Is.Empty);
            Assert.That(warnings.Single(), Does.Contain("crap.fasta looks like a UniProt database").And.Contain("Gencode"));

            gencode.ValidateAgainstHeaders("gencode.fasta", new[] { GencodeHeader }, out _, out warnings);
            Assert.That(warnings, Is.Empty);
        }

        [Test]
        public static void TheSamplerReadsTheFirstHeadersOfAGzippedFasta()
        {
            string path = TestPath("sampled.fasta.gz");
            Directory.CreateDirectory(Path.GetDirectoryName(path));
            using (var gz = new GZipStream(File.Create(path), CompressionLevel.Fastest))
            using (var writer = new StreamWriter(gz))
            {
                for (int i = 0; i < 7; i++)
                    writer.Write($">h{i}\nPEPTIDEK\n");
            }

            Assert.That(FastaHeaderSampler.IsFasta(path), Is.True);
            Assert.That(FastaHeaderSampler.ReadFirstHeaders(path), Is.EqualTo(new[] { ">h0", ">h1", ">h2", ">h3", ">h4" }));
            Assert.That(FastaHeaderSampler.ReadFirstHeaders(path, 2), Is.EqualTo(new[] { ">h0", ">h1" }));
        }

        [Test]
        public static void TheSamplerReturnsNullForAFileItCannotRead()
        {
            Assert.That(FastaHeaderSampler.ReadFirstHeaders(TestPath("doesNotExist.fasta")), Is.Null);
        }

        [TestCase("db.fasta", true)]
        [TestCase("db.fa", true)]
        [TestCase("db.FASTA.gz", true)]
        [TestCase("db.xml", false)]
        [TestCase("db.xml.gz", false)]
        public static void TheSamplerUsesTheLoadersExtensionRule(string fileName, bool isFasta)
        {
            Assert.That(FastaHeaderSampler.IsFasta(fileName), Is.EqualTo(isFasta));
        }

        [Test]
        public static void AutoRefusesAnUndetectableGzippedDatabaseWithAWarningInsteadOfThrowing()
        {
            // Not UniProt (no two-letter code and pipe), no Ensembl field markers, not seven pipes.
            // mzLib would throw MzLibException on this, which reaches a GUI user as a crash report.
            string fastaPath = TestPath("undetectable.fasta.gz");
            Directory.CreateDirectory(Path.GetDirectoryName(fastaPath));
            using (var gz = new GZipStream(File.Create(fastaPath), CompressionLevel.Fastest))
            using (var writer = new StreamWriter(gz))
                writer.Write(">rat_reductase\nPEPTIDEK\n");

            var task = new SearchTask
            {
                CommonParameters = new CommonParameters(fastaHeaderParsing: new FastaHeaderParsingParameters(FastaHeaderFormat.Auto))
            };

            var startedEngines = new List<string>();
            var engine = RunAndRecordEngines(new List<(string, MetaMorpheusTask)> { ("Task1-SearchTask", task) },
                ProteinSpectra, new List<DbForTask> { new DbForTask(fastaPath, false) }, "TestFastaHeaderAutoGate", startedEngines);

            Assert.That(startedEngines, Is.Empty, "no engine should start: " + string.Join(", ", startedEngines));
            Assert.That(engine.Warnings.Single(), Does.StartWith("Cannot proceed. Task1-SearchTask:")
                .And.Contain("undetectable.fasta.gz").And.Contain("could not be detected automatically"));
        }

        #endregion

        #region Nucleic acid tasks

        [Test]
        public static void ANonDefaultFormatOnAnRnaTaskWarnsOnceAndIsNotCheckedAgainstTheRnaHeaders()
        {
            // RnaDbLoader takes no header regexes, so the setting is inert on an RNA task. Checked against
            // this database, the custom regex would be refused, blocking a search it cannot affect.
            string db = Path.Combine(TestContext.CurrentContext.TestDirectory, "Transcriptomics", "TestData", "6mer.fasta");
            string spectra = Path.Combine(TestContext.CurrentContext.TestDirectory, "Transcriptomics", "TestData", "GUACUG_NegativeMode_Sliced.mzML");
            var matchesNothing = new FastaHeaderParsingParameters(FastaHeaderFormat.Custom, accessionRegex: @"^>gi\|(\d+)");
            Assert.That(matchesNothing.ValidateAgainstHeaders("6mer.fasta", FastaHeaderSampler.ReadFirstHeaders(db), out _, out _), Is.False);

            CommonParameters Rna(FastaHeaderParsingParameters headerParsing) => new(
                dissociationType: DissociationType.CID, deconvolutionMaxAssumedChargeState: -20, deconvolutionIntensityRatio: 3,
                deconvolutionMassTolerance: new PpmTolerance(20), precursorMassTolerance: new PpmTolerance(10),
                productMassTolerance: new PpmTolerance(5), scoreCutoff: 5, totalPartitions: 1, maxThreadsToUsePerFile: 1,
                digestionParams: new RnaDigestionParams(), fastaHeaderParsing: headerParsing);
            SearchTask RnaSearch(FastaHeaderParsingParameters headerParsing) => new()
            {
                CommonParameters = Rna(headerParsing),
                SearchParameters = new RnaSearchParameters { DecoyType = DecoyType.Reverse, MassDiffAcceptorType = MassDiffAcceptorType.Custom, CustomMdac = "Custom interval [-5,5]" }
            };

            var startedEngines = new List<string>();
            var engine = RunAndRecordEngines(new List<(string, MetaMorpheusTask)>
                {
                    ("Task1-Default", RnaSearch(new FastaHeaderParsingParameters())),
                    ("Task2-Custom", RnaSearch(matchesNothing)),
                },
                spectra, new List<DbForTask> { new DbForTask(db, false) }, "TestFastaHeaderRnaGate", startedEngines);

            Assert.That(engine.Warnings.Single(), Does.StartWith("Task2-Custom:").And.Contain("no effect on this nucleic acid search"));
            Assert.That(startedEngines, Does.Contain(nameof(DatabaseLoadingEngine)), "the header setting must not stop an RNA search");
        }

        [Test]
        public static void ANonDefaultFormatOnAProteinTaskDoesNotGetTheRnaWarning()
        {
            var task = new SearchTask
            {
                CommonParameters = new CommonParameters(fastaHeaderParsing: new FastaHeaderParsingParameters(FastaHeaderFormat.Custom))
            };

            var engine = RunAndRecordEngines(new List<(string, MetaMorpheusTask)> { ("Task1-SearchTask", task) },
                ProteinSpectra, new List<DbForTask> { new DbForTask(WriteFasta("proteinOnly.fasta", UniProtHeader, "PEPTIDEK"), false) },
                "TestFastaHeaderProteinGate", new List<string>());

            // Custom with no accession regex is refused before the run; that is the only warning.
            Assert.That(engine.Warnings.Single(), Does.Contain("requires an accession regex"));
        }

        #endregion

        #region Load errors and the index cache

        [Test]
        public static void TheLoadersOwnErrorsAreReportedRatherThanDiscarded()
        {
            string path = WriteFasta("emptyEntry.fasta", ">sp|P00001|EMPTY_HUMAN Empty OS=Homo sapiens", UniProtHeader, "PEPTIDEK");

            var proteins = DatabaseLoadingEngine.LoadBioPolymers("task", new List<DbForTask> { new DbForTask(path, false) }, true,
                DecoyType.None, new List<string>(), new CommonParameters(), out var errors);

            Assert.That(proteins.Single().Accession, Is.EqualTo("P12345"));
            Assert.That(errors.Single(), Does.Contain("Protein Length of 0").And.Contain("emptyEntry.fasta"));
        }

        [Test]
        public static void AnIndexCachedUnderOneFormatIsNotReusedUnderAnother()
        {
            // Under one key the second run reads the cached index, cannot find the old accessions, warns
            // and rebuilds. With the format in the key it builds its own without trying.
            string root = TestPath("indexCache");
            if (Directory.Exists(root))
                Directory.Delete(root, true);
            Directory.CreateDirectory(Path.Combine(root, "out", "Task Settings"));
            string db = Path.Combine(root, "db.fasta");
            File.Copy(Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", "DbForPrunedDb.fasta"), db);

            void RunModern(string name, FastaHeaderParsingParameters headerParsing)
            {
                string output = Path.Combine(root, "out", name);
                Directory.CreateDirectory(output);
                new SearchTask
                {
                    SearchParameters = new SearchParameters { SearchType = SearchType.Modern, MaxFragmentSize = 3000 },
                    CommonParameters = new CommonParameters(fastaHeaderParsing: headerParsing)
                }.RunTask(output, new List<DbForTask> { new DbForTask(db, false) }, new List<string> { ProteinSpectra }, name);
            }

            var taskWarnings = new List<string>();
            void RecordWarning(object sender, StringEventArgs e) => taskWarnings.Add(e.S);
            MetaMorpheusTask.WarnHandler += RecordWarning;
            try
            {
                RunModern("uniProt", new FastaHeaderParsingParameters());
                // The UniProt entry name as accession: same protein count, different accessions.
                RunModern("entryName", new FastaHeaderParsingParameters(FastaHeaderFormat.Custom, accessionRegex: @"^>\w+\|[^|]*\|(\S+)"));
            }
            finally
            {
                MetaMorpheusTask.WarnHandler -= RecordWarning;
            }

            Assert.That(taskWarnings, Has.None.Contains("Could not read the existing index"));

            var keys = Directory.GetDirectories(Path.Combine(root, MetaMorpheusTask.IndexFolderName))
                .Select(d => File.ReadAllText(Path.Combine(d, MetaMorpheusTask.IndexEngineParamsFileName))).ToList();
            Assert.That(keys, Has.Count.EqualTo(2));
            Assert.That(keys, Has.Exactly(1).Contains("FASTA header format: Custom; accession=^>\\w+\\|[^|]*\\|(\\S+);"));
            // The default adds nothing, so indexes cached by earlier versions still match.
            Assert.That(keys, Has.Exactly(1).Not.Contains("FASTA header format"));
        }

        #endregion
    }
}
