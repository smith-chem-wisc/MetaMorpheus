#nullable enable
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Net.Http;
using System.Text.RegularExpressions;
using System.Threading;
using System.Threading.Tasks;
using EngineLayer;
using EngineLayer.DatabaseLoading;
using MassSpectrometry;
using Newtonsoft.Json;
using Newtonsoft.Json.Linq;
using NUnit.Framework;
using Omics.Digestion;
using Omics.SpectralMatch.MslSpectralLibrary;
using PredictionClients.Koina.Client;
using PredictionClients.Koina.SupportedModels.FragmentIntensityModels;
using PredictionClients.Koina.SupportedModels.RetentionTimeModels;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using Readers;
using Readers.SpectralLibrary;
using TaskLayer;
using Transcriptomics.Digestion;

namespace Test.DiaLibrarySearch;

/// <summary>
/// A DIA library search given proteins and no library predicts one (Prosit 2025 through Koina, with decoys), writes it as
/// DiaPredictedLibrary.msl and searches it. Koina is answered from memory here; the live call is in
/// <see cref="DiaPredictedLibraryLiveTests"/>.
/// </summary>
[TestFixture]
[NonParallelizable] // GlobalVariables.AnalyteType is static
[ExcludeFromCodeCoverage]
public class DiaPredictedLibraryTests
{
    private string _directory = "";

    [OneTimeSetUp]
    public void SetUp()
    {
        _directory = Path.Combine(TestContext.CurrentContext.TestDirectory, "DiaPredictedLibraryTests");
        Directory.CreateDirectory(_directory);
    }

    [OneTimeTearDown]
    public void TearDown() => Directory.Delete(_directory, true);

    [TearDown]
    public void ResetAnalyte() => GlobalVariables.AnalyteType = AnalyteType.Peptide;

    /// <summary>
    /// Koina answered from memory, the same in every build: each peptide's fragments and iRT come from a seed its sequence
    /// gives, so a library built here and one the task builds are the same library.
    /// </summary>
    internal sealed class CannedProsit
    {
        public bool Unreachable { get; set; }
        public int Calls;
        public ConcurrentBag<string> Requested { get; } = [];
        public FakeIntensity Intensity { get; }
        public FakeIrt Irt { get; }

        public CannedProsit()
        {
            Intensity = new FakeIntensity(this);
            Irt = new FakeIrt(this);
        }

        /// <summary>The residues of a peptide as Koina is sent it, modifications and terminal marks left out.</summary>
        public static string Residues(string koinaSequence) => Regex.Replace(koinaSequence, @"\[[^\]]*\]|-", "");

        private static Random SeededBy(string text) => new(text.Aggregate(17, (h, c) => unchecked(h * 31 + c)));

        private static string Numbers(IEnumerable<double> values) => string.Join(",", values.Select(v => v.ToString("F4", CultureInfo.InvariantCulture)));

        private Task<string> Answer(Dictionary<string, object> request, Func<string[], string> outputs)
        {
            Interlocked.Increment(ref Calls);
            if (Unreachable)
                throw new HttpRequestException("Koina is unreachable");
            var peptides = JObject.Parse(JsonConvert.SerializeObject(request))["inputs"]![0]!["data"]!.Select(t => t.ToString()).ToArray();
            foreach (string peptide in peptides)
                Requested.Add(peptide);
            return Task.FromResult("{\"outputs\":[" + outputs(peptides) + "]}");
        }

        public sealed class FakeIntensity(CannedProsit owner) : Prosit2025Intensity40PTM
        {
            private static readonly string[] Ions = ["y3+1", "y4+1", "y5+1", "y6+1", "y7+1", "y8+1", "b3+1", "b4+1"];

            protected override Task<string> SendInferenceRequestAsync(string modelName, Dictionary<string, object> request, CancellationToken cancellationToken) =>
                owner.Answer(request, peptides =>
                {
                    int n = peptides.Length, k = Ions.Length;
                    string annotations = string.Join(",", Enumerable.Repeat(string.Join(",", Ions.Select(i => $"\"{i}\"")), n));
                    string mz = Numbers(peptides.SelectMany(p => { var r = SeededBy(p); return Ions.Select(_ => 200 + 1200 * r.NextDouble()); }));
                    string intensities = Numbers(peptides.SelectMany(p => { var r = SeededBy(p + "/intensity"); return Ions.Select(_ => 0.1 + 0.9 * r.NextDouble()); }));
                    return $"{{\"name\":\"annotation\",\"datatype\":\"BYTES\",\"shape\":[{n},{k}],\"data\":[{annotations}]}}," +
                        $"{{\"name\":\"mz\",\"datatype\":\"FP32\",\"shape\":[{n},{k}],\"data\":[{mz}]}}," +
                        $"{{\"name\":\"intensities\",\"datatype\":\"FP32\",\"shape\":[{n},{k}],\"data\":[{intensities}]}}";
                });
        }

        public sealed class FakeIrt(CannedProsit owner) : Prosit2025iRT40PTM
        {
            protected override Task<string> SendInferenceRequestAsync(string modelName, Dictionary<string, object> request, CancellationToken cancellationToken) =>
                owner.Answer(request, peptides =>
                    $"{{\"name\":\"irt\",\"datatype\":\"FP32\",\"shape\":[{peptides.Length}],\"data\":[{Numbers(peptides.Select(p => -15 + 130 * SeededBy(p + "/irt").NextDouble()))}]}}");
        }
    }

    /// <summary>
    /// Proteins of two tryptic pieces each, P00000, P00001, ...: random residues with no K, R, C or M inside (so the task's
    /// default carbamidomethyl and oxidation never apply), each piece ending in K.
    /// </summary>
    private static List<Protein> Proteins(int count, int seed = 7, int pieceLength = 12)
    {
        var random = new Random(seed);
        const string residues = "ADEFGHILNPQSTVWY";
        string Piece(int length) => new string(Enumerable.Range(0, length - 1).Select(_ => residues[random.Next(residues.Length)]).ToArray()) + "K";
        return Enumerable.Range(0, count).Select(i => new Protein(Piece(pieceLength) + Piece(pieceLength), $"P{i:D5}")).ToList();
    }

    private static string WriteFasta(string folder, IEnumerable<Protein> proteins)
    {
        string path = Path.Combine(folder, "proteins.fasta");
        File.WriteAllLines(path, proteins.SelectMany(p => new[] { $">sp|{p.Accession}|{p.Accession}_HUMAN Protein {p.Accession} OS=Homo sapiens OX=9606 GN=G{p.Accession} PE=1 SV=1", p.BaseSequence }));
        return path;
    }

    private static string WriteRun(string folder, SyntheticDiaRun run)
    {
        string mzml = Path.Combine(folder, "run.mzML");
        MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(new GenericMsDataFile(DiaLibrarySearchTaskTests.AsWritten(run.Scans),
            new SourceFile("no nativeID format", "mzML format", null, null, null)), mzml, false);
        return mzml;
    }

    private static List<MslLibraryEntry> Build(IReadOnlyList<Protein> proteins, CannedProsit models, out PredictionClients.SpectralLibraryGeneration.MslLibraryBuildResult result,
        CommonParameters? commonParameters = null, DiaLibrarySearchTaskParameters? parameters = null) =>
        DiaPredictedLibrary.Build(proteins, commonParameters ?? new CommonParameters(), parameters ?? new DiaLibrarySearchTaskParameters(), [], [],
            models.Intensity, models.Irt, out result);

    /// <summary>
    /// End to end: given a FASTA and no library, the task predicts one from the proteins (targets and permuted decoys),
    /// writes it as DiaPredictedLibrary.msl beside its tables, searches the run against it, infers protein groups from the
    /// decoys' accessions, and says in results.txt what it built and with which models.
    /// </summary>
    [Test]
    public void WithProteinsAndNoLibraryTheTaskPredictsOneAndSearchesIt()
    {
        string folder = Path.Combine(_directory, "endtoend");
        Directory.CreateDirectory(folder);
        var proteins = Proteins(150);
        string fasta = WriteFasta(folder, proteins);
        // The library the task will build, built here too, so the run can be planted from it
        var expected = Build(proteins, new CannedProsit(), out var expectedResult);
        var run = SyntheticDiaRun.Plant(expected, e => !e.IsDecoy && SyntheticDiaRun.Bucket(e, 4) != 0, withMs1: true);
        string mzml = WriteRun(folder, run);
        string output = Path.Combine(folder, "output", "Task1");
        Directory.CreateDirectory(output);
        var models = new CannedProsit();

        new DiaLibrarySearchTask { PredictionModels = () => (models.Intensity, models.Irt) }.RunTask(output, [new DbForTask(fasta, false)], [mzml], "Task1");

        string libraryPath = Path.Combine(output, DiaPredictedLibrary.FileName);
        Assert.That(File.Exists(libraryPath));
        using var library = MslLibrary.Load(libraryPath);
        var entries = Enumerable.Range(0, library.PrecursorCount).Select(i => library.GetEntry(i, includeFragments: false)!).ToList();
        Assert.That(entries.Select(e => (e.FullSequence, e.ChargeState, e.IsDecoy)), Is.EquivalentTo(expected.Select(e => (e.FullSequence, e.ChargeState, e.IsDecoy))),
            "the task builds the library built here");
        var accessions = proteins.Select(p => p.Accession).ToHashSet();
        Assert.That(entries.Where(e => !e.IsDecoy).SelectMany(e => e.ProteinAccession.Split('|')).Distinct(), Is.SubsetOf(accessions));
        var decoys = entries.Where(e => e.IsDecoy).ToList();
        Assert.That(decoys, Is.Not.Empty);
        Assert.That(decoys, Has.None.Matches<MslLibraryEntry>(e => e.IsEntrapment), "a permuted decoy is a decoy, not entrapment");
        Assert.That(decoys.SelectMany(e => e.ProteinAccession.Split('|')).Select(a => Regex.Match(a, @"^DECOY_(P\d{5})_f1$")),
            Has.All.Matches<Match>(m => m.Success && accessions.Contains(m.Groups[1].Value)), "each decoy names its owner's permuted decoy protein");
        Assert.That(models.Calls, Is.GreaterThan(0));

        string results = File.ReadAllText(Path.Combine(output, "results.txt"));
        Assert.That(results, Does.Contain($"Predicted library {DiaPredictedLibrary.FileName}: {expectedResult.TargetPrecursors} target and {expectedResult.DecoyPrecursors} decoy precursors from 150 proteins, permuted decoys"));
        Assert.That(results, Does.Contain(models.Intensity.ModelName).And.Contain(models.Irt.ModelName).And.Contain(HTTP.ModelsURL));
        int found = int.Parse(Regex.Match(results, @"All target precursors with q-value <= 0\.01: (\d+)").Groups[1].Value);
        int planted = run.Library.Count(e => !e.IsDecoy && run.PlantedSequences.Contains(e.FullSequence)
            && e.PrecursorMz >= SyntheticDiaRun.FirstWindowLowMz && e.PrecursorMz < SyntheticDiaRun.FirstWindowLowMz + SyntheticDiaRun.WindowCount * SyntheticDiaRun.WindowWidth);
        Assert.That(found, Is.InRange(planted * 7 / 10, planted * 11 / 10), "most planted targets are found, few others");
        Assert.That(File.Exists(Path.Combine(output, "AllDiaProteinGroups.tsv")), "the decoys name decoy proteins, so protein groups are inferred");
    }

    /// <summary>
    /// Prosit takes peptides of at most 30 residues, and MetaMorpheus's default maximum length is unlimited: the library's
    /// digestion stops at 30, so nothing longer is asked of Koina or counted as not predicted.
    /// </summary>
    [Test]
    public void NoPeptideLongerThanPrositTakesIsAskedFor()
    {
        var models = new CannedProsit();
        var proteins = Proteins(5, pieceLength: 20);

        var entries = Build(proteins, models, out var result);

        Assert.That(models.Requested.Select(CannedProsit.Residues).Max(p => p.Length), Is.EqualTo(20), "the 40-residue joins of two pieces are not digested");
        Assert.That(result.NotPredicted, Is.Zero);
        Assert.That(entries.Where(e => !e.IsDecoy).Select(e => e.BaseSequence.Length).Distinct(), Is.EqualTo(new[] { 20 }));
    }

    /// <summary>Reversed decoys, for a digestion the permuter refuses: mzLib's DECOY_ and the target's accession.</summary>
    [Test]
    public void ReversedDecoysCanBeAskedFor()
    {
        var entries = Build(Proteins(5), new CannedProsit(), out var result, parameters: new DiaLibrarySearchTaskParameters { LibraryDecoys = DiaLibraryDecoys.Reversed });

        var decoys = entries.Where(e => e.IsDecoy).ToList();
        Assert.That(decoys, Is.Not.Empty);
        Assert.That(decoys.Select(e => e.ProteinAccession), Has.All.Matches<string>(a => Regex.IsMatch(a, @"^DECOY_P\d{5}$")));
        Assert.That(result.DecoyPrecursors, Is.EqualTo(decoys.Count));
    }

    /// <summary>
    /// A contaminant protein is a target of the library too, and gets a permuted decoy like any other (the permuter itself
    /// refuses a contaminant, so the task hands it a plain copy).
    /// </summary>
    [Test]
    public void AContaminantGetsAPermutedDecoyToo()
    {
        var proteins = Proteins(3);
        var contaminant = Proteins(1, seed: 11)[0];
        proteins.Add(new Protein(contaminant.BaseSequence, "Q99999", isContaminant: true));

        var entries = Build(proteins, new CannedProsit(), out _);

        Assert.That(entries.Where(e => !e.IsDecoy).Select(e => e.ProteinAccession), Does.Contain("Q99999"));
        Assert.That(entries.Where(e => e.IsDecoy).Select(e => e.ProteinAccession), Does.Contain("DECOY_Q99999_f1"));
    }

    /// <summary>A digestion the permuter cannot pair (here a semi-specific one) is refused before Koina is asked, saying what to set instead.</summary>
    [Test]
    public void ADigestionThePermuterRefusesSaysToUseReversedDecoys()
    {
        var models = new CannedProsit();
        var semi = new CommonParameters(digestionParams: new DigestionParams(searchModeType: CleavageSpecificity.Semi));

        var e = Assert.Throws<MetaMorpheusException>(() => Build(Proteins(3), models, out _, commonParameters: semi))!;

        Assert.That(e.Message, Does.Contain("LibraryDecoys = \"Reversed\""));
        Assert.That(models.Calls, Is.Zero);
    }

    /// <summary>
    /// Koina unreachable: nothing is left to search, so the task fails, saying the library needs network access to Koina and
    /// that an .msl library can be given instead, and writes no library.
    /// </summary>
    [Test]
    public void UnreachableKoinaFailsTheTaskSayingWhatItNeeds()
    {
        string folder = Path.Combine(_directory, "unreachable");
        Directory.CreateDirectory(folder);
        string fasta = WriteFasta(folder, Proteins(5));
        string mzml = WriteRun(folder, SyntheticDiaRun.Build(20, _ => true));
        string output = Path.Combine(folder, "output", "Task1");
        Directory.CreateDirectory(output);
        var models = new CannedProsit { Unreachable = true };

        var e = Assert.Throws<MetaMorpheusException>(() =>
            new DiaLibrarySearchTask { PredictionModels = () => (models.Intensity, models.Irt) }.RunTask(output, [new DbForTask(fasta, false)], [mzml], "Task1"))!;

        Assert.That(e.Message, Does.Contain(HTTP.ModelsURL).And.Contain("network").And.Contain(".msl"));
        Assert.That(e.InnerException, Is.InstanceOf<HttpRequestException>());
        Assert.That(File.Exists(Path.Combine(output, DiaPredictedLibrary.FileName)), Is.False);
    }

    /// <summary>
    /// What Prosit 2025 cannot predict for is refused before Koina is asked: an analyte other than peptides, a digestion
    /// that is not a protease's, a fragmentation other than HCD or CID, charges outside 1-6, and an energy that is not one.
    /// </summary>
    [TestCase("ETD", "ETD")]
    [TestCase("Autodetect", "Autodetect")]
    [TestCase("rna", "protease")]
    [TestCase("oligo", "peptides")]
    [TestCase("nocharges", "charge")]
    [TestCase("charge7", "charge")]
    [TestCase("nce0", "NCE")]
    public void WhatPrositCannotPredictForIsRefused(string setting, string said)
    {
        var commonParameters = setting switch
        {
            "ETD" => new CommonParameters(dissociationType: DissociationType.ETD),
            "Autodetect" => new CommonParameters(dissociationType: DissociationType.Autodetect),
            "rna" => new CommonParameters(digestionParams: new RnaDigestionParams()),
            _ => new CommonParameters(),
        };
        var parameters = setting switch
        {
            "nocharges" => new DiaLibrarySearchTaskParameters { LibraryPrecursorCharges = [] },
            "charge7" => new DiaLibrarySearchTaskParameters { LibraryPrecursorCharges = [2, 7] },
            "nce0" => new DiaLibrarySearchTaskParameters { LibraryNce = 0 },
            _ => new DiaLibrarySearchTaskParameters(),
        };
        if (setting == "oligo")
            GlobalVariables.AnalyteType = AnalyteType.Oligo;

        string? reason = DiaPredictedLibrary.Unsupported(commonParameters, [], parameters);
        GlobalVariables.AnalyteType = AnalyteType.Peptide;

        Assert.That(reason, Does.Contain(said));
        Assert.That(DiaPredictedLibrary.Unsupported(new CommonParameters(), [], new DiaLibrarySearchTaskParameters()), Is.Null, "the defaults can be predicted");
        Assert.That(DiaPredictedLibrary.Unsupported(new CommonParameters(dissociationType: DissociationType.CID), [], new DiaLibrarySearchTaskParameters()), Is.Null);
        Assert.That(DiaPredictedLibrary.Unsupported(new CommonParameters(), [DissociationType.ETD], new DiaLibrarySearchTaskParameters()), Does.Contain("ETD"),
            "a file's own fragmentation counts too");
    }

    /// <summary>The task refuses an unpredictable search before asking Koina, and says to give it an .msl library instead.</summary>
    [Test]
    public void TheTaskRefusesAnUnpredictableSearchBeforeAskingKoina()
    {
        string folder = Path.Combine(_directory, "etd");
        Directory.CreateDirectory(folder);
        string fasta = WriteFasta(folder, Proteins(5));
        string mzml = WriteRun(folder, SyntheticDiaRun.Build(20, _ => true));
        string output = Path.Combine(folder, "output", "Task1");
        Directory.CreateDirectory(output);
        var models = new CannedProsit();
        var task = new DiaLibrarySearchTask { CommonParameters = new CommonParameters(dissociationType: DissociationType.ETD), PredictionModels = () => (models.Intensity, models.Irt) };

        var e = Assert.Throws<MetaMorpheusException>(() => task.RunTask(output, [new DbForTask(fasta, false)], [mzml], "Task1"))!;

        Assert.That(e.Message, Does.Contain("ETD").And.Contain(".msl"));
        Assert.That(models.Calls, Is.Zero);
    }

    /// <summary>With neither a library nor proteins to predict one from, the task says it needs one or the other.</summary>
    [Test]
    public void WithNeitherALibraryNorProteinsTheTaskSaysWhatItNeeds()
    {
        string folder = Path.Combine(_directory, "nothing");
        Directory.CreateDirectory(folder);
        string mzml = WriteRun(folder, SyntheticDiaRun.Build(20, _ => true));
        string output = Path.Combine(folder, "output", "Task1");
        Directory.CreateDirectory(output);

        var e = Assert.Throws<MetaMorpheusException>(() => new DiaLibrarySearchTask().RunTask(output, [], [mzml], "Task1"))!;

        Assert.That(e.Message, Does.Contain(".msl").And.Contain("protein database"));
    }

    /// <summary>Given a library, the task predicts nothing, whatever protein database comes with it.</summary>
    [Test]
    public void WithALibraryNothingIsPredicted()
    {
        string folder = Path.Combine(_directory, "given");
        Directory.CreateDirectory(folder);
        var run = SyntheticDiaRun.Build(100, e => !e.IsDecoy, withMs1: true);
        string mzml = WriteRun(folder, run);
        string library = run.WriteLibrary(folder);
        string fasta = WriteFasta(folder, Proteins(5));
        string output = Path.Combine(folder, "output", "Task1");
        Directory.CreateDirectory(output);
        var models = new CannedProsit { Unreachable = true };

        new DiaLibrarySearchTask { PredictionModels = () => (models.Intensity, models.Irt) }
            .RunTask(output, [new DbForTask(library, false), new DbForTask(fasta, false)], [mzml], "Task1");

        Assert.That(models.Calls, Is.Zero);
        Assert.That(File.Exists(Path.Combine(output, DiaPredictedLibrary.FileName)), Is.False);
        Assert.That(File.ReadAllText(Path.Combine(output, "results.txt")), Does.Not.Contain("Predicted library"));
    }
}

/// <summary>Live: a small library through the real Prosit 2025 models on Koina. Skipped when Koina cannot be reached.</summary>
[TestFixture]
[Category("ExternalService")]
[ExcludeFromCodeCoverage]
public class DiaPredictedLibraryLiveTests
{
    [Test]
    public void ASmallLibraryIsPredictedThroughKoina()
    {
        var (intensity, irt) = DiaPredictedLibrary.KoinaModels();
        List<MslLibraryEntry> entries;
        PredictionClients.SpectralLibraryGeneration.MslLibraryBuildResult result;
        try
        {
            entries = DiaPredictedLibrary.Build([new Protein("MSLPDKAGTCVLVEGDRPYVLNSR", "P1")], new CommonParameters(), new DiaLibrarySearchTaskParameters(), [], [],
                intensity, irt, out result);
        }
        catch (MetaMorpheusException e) when (e.InnerException is HttpRequestException or OperationCanceledException)
        {
            Assert.Ignore($"Koina unavailable: {e.Message}");
            return;
        }

        Assert.That(result.TargetPrecursors, Is.GreaterThan(0));
        Assert.That(entries.Where(e => !e.IsDecoy), Has.All.Matches<MslLibraryEntry>(e => e.MatchedFragmentIons.Count > 0 && double.IsFinite(e.RetentionTime)));
        Assert.That(entries.Where(e => e.IsDecoy).Select(e => e.ProteinAccession), Has.All.StartWith("DECOY_P1_f1"));
    }
}
