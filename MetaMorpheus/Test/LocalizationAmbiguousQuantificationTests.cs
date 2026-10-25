using Chemistry;
using EngineLayer.DatabaseLoading;
using EngineLayer;
using EngineLayer.FdrAnalysis;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Omics.Digestion;
using Omics.Fragmentation;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using Readers;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.IO;
using System.Linq;
using System.Reflection;
using TaskLayer;

namespace Test
{
    /// <summary>
    /// Label-free quantification of a match that is ambiguous only in where its modifications sit: two or more full
    /// sequences of one base sequence, carrying the same modifications at different sites, match the MS2 scan equally
    /// well. Such a match is unique to its protein group when all its candidates belong to that one group, and it is
    /// quantified for it, under a name made of its candidates (no site is claimed). A match whose candidates span
    /// protein groups, are different peptides, or carry different modifications stays out of quantification.
    /// </summary>
    [TestFixture]
    [ExcludeFromCodeCoverage]
    public static class LocalizationAmbiguousQuantificationTests
    {
        private static readonly DigestionParams Digestion = new(protease: "trypsin", minPeptideLength: 1);

        private static Modification Mod(string id, string residue, string formula)
        {
            ModificationMotif.TryGetMotif(residue, out var motif);
            return new Modification(_originalId: id, _modificationType: "Test", _target: motif,
                _locationRestriction: "Anywhere.", _chemicalFormula: ChemicalFormula.ParseFormula(formula));
        }

        private static readonly Modification OxidationOnM = Mod("Oxidation", "M", "O");
        private static readonly Modification HydroxylationOnD = Mod("Hydroxylation", "D", "O");
        private static readonly Modification DioxidationOnM = Mod("Dioxidation", "M", "O2");
        private static readonly Modification PhosphorylationOnS = Mod("Phosphorylation", "S", "HO3P");
        private static readonly Modification PhosphorylationOnT = Mod("Phosphorylation", "T", "HO3P");

        private static Protein Protein(string sequence, string accession) =>
            new(sequence, accession, null, new List<System.Tuple<string, string>>(), new Dictionary<int, List<Modification>>());

        /// <summary>A peptide of <paramref name="protein"/>, with <paramref name="mods"/> keyed by 1-based residue in the peptide.</summary>
        private static PeptideWithSetModifications Peptide(Protein protein, int start, int end, params (int Residue, Modification Mod)[] mods) =>
            new(protein: protein, digestionParams: Digestion, oneBasedStartResidueInProtein: start,
                oneBasedEndResidueInProtein: end, cleavageSpecificity: CleavageSpecificity.Full, peptideDescription: "",
                missedCleavages: 0, allModsOneIsNterminus: mods.ToDictionary(m => m.Residue + 1, m => m.Mod), numFixedMods: 0);

        private static Ms2ScanWithSpecificMass Scan(string filePath, int scanNumber)
        {
            MsDataScan scan = new(
                new MzSpectrum(new[] { 100.0, 200.0 }, new[] { 1.0, 1.0 }, false),
                oneBasedScanNumber: scanNumber, msnOrder: 2, isCentroid: true,
                polarity: Polarity.Positive, retentionTime: 10.0,
                scanWindowRange: new MzRange(0, 2000), scanFilter: "f",
                mzAnalyzer: MZAnalyzerType.Orbitrap, totalIonCurrent: 2,
                injectionTime: 1.0, noiseData: null, nativeId: $"scan={scanNumber}",
                selectedIonMz: 500.0, selectedIonChargeStateGuess: 2,
                selectedIonIntensity: 1, isolationMZ: 500.0, isolationWidth: 2,
                dissociationType: DissociationType.HCD, oneBasedPrecursorScanNumber: null,
                selectedIonMonoisotopicGuessMz: 500.0);
            return new Ms2ScanWithSpecificMass(scan, 500.0, 2, filePath, new CommonParameters());
        }

        /// <summary>One confident match; several peptides make it ambiguous, each matching the scan with the same score.</summary>
        private static SpectralMatch Psm(string filePath, int scanNumber, params PeptideWithSetModifications[] peptides) =>
            ScoredPsm(10, filePath, scanNumber, peptides);

        /// <summary>As <see cref="Psm"/>, with the score given: peptide-level q-values are computed in score order.</summary>
        private static SpectralMatch ScoredPsm(double score, string filePath, int scanNumber, params PeptideWithSetModifications[] peptides)
        {
            var common = new CommonParameters(digestionParams: Digestion);
            SpectralMatch psm = new PeptideSpectralMatch(peptides[0], 0, score, scanNumber, Scan(filePath, scanNumber), common, new List<MatchedFragmentIon>());
            foreach (var peptide in peptides.Skip(1))
            {
                psm.AddOrReplace(peptide, score, 0, true, new List<MatchedFragmentIon>());
            }
            psm.ResolveAllAmbiguities();
            psm.SetFdrValues(1, 0, 0, 1, 0, 0, 0, 0);
            psm.PeptideFdrInfo = new FdrInfo { QValue = 0, QValueNotch = 0, PEP = 0, PEP_QValue = 0 };
            return psm;
        }

        /// <summary>
        /// Three MS1 scans around 10 min holding the charge-2 isotope envelope of each peptide, so FlashLFQ finds a peak
        /// for it. Without the signal every intensity is zero, and a test of "was it quantified" would pass vacuously.
        /// </summary>
        private static void WriteMs1Fixture(string mzmlPath, params PeptideWithSetModifications[] peptides)
        {
            var envelopes = peptides.Select(p => IsotopicDistribution.GetDistribution(p.FullChemicalFormula, 0.125, 1e-8)).ToList();
            var scans = new[] { 9.8, 10.0, 10.2 }
                .Select((rt, i) =>
                {
                    double scale = i == 1 ? 1e6 : 5e5;
                    var peaks = envelopes
                        .SelectMany(e => e.Masses.Select(m => m.ToMz(2)).Zip(e.Intensities, (mz, intensity) => (mz, intensity: intensity * scale)))
                        .OrderBy(peak => peak.mz)
                        .ToList();
                    double[] mzs = peaks.Select(peak => peak.mz).ToArray();
                    double[] intensities = peaks.Select(peak => peak.intensity).ToArray();
                    return new MsDataScan(new MzSpectrum(mzs, intensities, false), i + 1, 1, true, Polarity.Positive, rt,
                        new MzRange(0, 2000), "f", MZAnalyzerType.Orbitrap, intensities.Sum(), 1.0, null, "scan=" + (i + 1));
                })
                .ToArray();
            Readers.MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(
                new GenericMsDataFile(scans, new SourceFile(null, null, null, null, null)), mzmlPath, false);
        }

        /// <summary>Parsimony, protein scoring, then label-free quantification, as a search runs them.</summary>
        private static PostSearchAnalysisParameters Quantify(string outputFolder, string mzmlPath, List<SpectralMatch> psms,
            System.Action<SearchParameters> configure = null, bool scrubContaminants = false) =>
            Quantify(outputFolder, new List<string> { mzmlPath }, psms, configure, scrubContaminants);

        private static PostSearchAnalysisParameters Quantify(string outputFolder, List<string> mzmlPaths, List<SpectralMatch> psms,
            System.Action<SearchParameters> configure = null, bool scrubContaminants = false)
        {
            var common = new CommonParameters(digestionParams: Digestion);
            PostSearchAnalysisParameters parameters = new()
            {
                SearchParameters = new SearchParameters
                {
                    DoParsimony = true,
                    DoLabelFreeQuantification = true,
                    DoMultiplexQuantification = false,
                    MatchBetweenRuns = false,
                    Normalize = false,
                },
                OutputFolder = outputFolder,
                IndividualResultsOutputFolder = outputFolder,
                SearchTaskId = "TestTask",
                AllSpectralMatches = psms,
                CurrentRawFileList = mzmlPaths,
                MyFileManager = new MyFileManager(true),
                FixedModifications = new List<Modification>(),
                VariableModifications = new List<Modification>(),
                ListOfDigestionParams = new HashSet<IDigestionParams> { Digestion },
                DatabaseFilenameList = new List<DbForTask>(),
                FileSettingsList = mzmlPaths.Select(_ => (FileSpecificParameters)null).ToArray(),
            };
            PostSearchAnalysisTask task = new()
            {
                Parameters = parameters,
                CommonParameters = common,
                FileSpecificParameters = mzmlPaths.Select(path => (path, common)).ToList(),
            };
            configure?.Invoke(parameters.SearchParameters);
            var steps = new List<string> { "ProteinAnalysis", "QuantificationAnalysis" };
            if (scrubContaminants)
            {
                steps.Add("RemoveContaminantsFromQuantificationTables");
            }
            foreach (string step in steps)
            {
                typeof(PostSearchAnalysisTask).GetMethod(step, BindingFlags.NonPublic | BindingFlags.Instance)!.Invoke(task, null);
            }
            Assert.That(File.Exists(Path.Combine(outputFolder, "Quantification_crash.txt")), Is.False);
            return parameters;
        }

        private static string NewOutputFolder(string name)
        {
            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, name);
            if (Directory.Exists(outputFolder))
            {
                Directory.Delete(outputFolder, true);
            }
            Directory.CreateDirectory(outputFolder);
            return outputFolder;
        }

        /// <summary>
        /// The user's case: one protein, one base sequence, oxidation on M3 or on M5 explains the spectrum equally well.
        /// The match is unique to the protein's group (#2899 made it unique for parsimony), so it is quantified for that
        /// group, as a peptide named by both candidates. Before, every match without one full sequence was left out of
        /// quantification, and the group got nothing from it.
        /// </summary>
        [Test]
        public static void AMatchAmbiguousOnlyInModificationSiteIsQuantifiedForItsProteinGroup()
        {
            string outputFolder = NewOutputFolder("LocalizationAmbiguousQuant_OneGroup");
            string mzmlPath = Path.Combine(outputFolder, "file.mzML");
            var protein = Protein("ACMDMEKGGGR", "P1");
            var oxidizedAt3 = Peptide(protein, 1, 7, (3, OxidationOnM));
            var oxidizedAt5 = Peptide(protein, 1, 7, (5, OxidationOnM));
            WriteMs1Fixture(mzmlPath, oxidizedAt3);

            try
            {
                var ambiguous = Psm(mzmlPath, 1, oxidizedAt3, oxidizedAt5);
                Assume.That(ambiguous.FullSequence, Is.Null);
                Assume.That(ambiguous.BaseSequence, Is.EqualTo("ACMDMEK"));

                var results = Quantify(outputFolder, mzmlPath, new List<SpectralMatch> { ambiguous }).FlashLfqResults;

                string name = string.Join("|", new[] { oxidizedAt3.FullSequence, oxidizedAt5.FullSequence }.OrderBy(s => s, System.StringComparer.Ordinal));
                Assert.That(results, Is.Not.Null, "the match reaches FlashLFQ");
                Assert.That(results.PeptideModifiedSequences.Keys, Is.EquivalentTo(new[] { name }),
                    "quantified as one peptide named by both candidates, not as either site");
                var file = results.SpectraFiles.Single();
                Assert.That(results.PeptideModifiedSequences[name].GetIntensity(file), Is.GreaterThan(0));
                Assert.That(results.PeptideModifiedSequences[name].ProteinGroups.Select(g => g.ProteinGroupName), Is.EqualTo(new[] { "P1" }));
                Assert.That(results.ProteinGroups["P1"].GetIntensity(file), Is.GreaterThan(0),
                    "unique to P1's group, so it counts toward P1 with shared peptides off");
            }
            finally
            {
                Directory.Delete(outputFolder, true);
            }
        }

        /// <summary>
        /// Several proteins, one group: P1 and P2 carry the same sequence and merge into the group P1|P2. The match is
        /// ambiguous in site and in protein, but every candidate belongs to that one group, so it counts for it.
        /// </summary>
        [Test]
        public static void AMatchAmbiguousInSiteAcrossProteinsOfOneGroupIsQuantifiedForThatGroup()
        {
            string outputFolder = NewOutputFolder("LocalizationAmbiguousQuant_OneGroupTwoProteins");
            string mzmlPath = Path.Combine(outputFolder, "file.mzML");
            var p1 = Protein("ACMDMEKGGGR", "P1");
            var p2 = Protein("ACMDMEKGGGR", "P2");
            var candidates = new[]
            {
                Peptide(p1, 1, 7, (3, OxidationOnM)), Peptide(p1, 1, 7, (5, OxidationOnM)),
                Peptide(p2, 1, 7, (3, OxidationOnM)), Peptide(p2, 1, 7, (5, OxidationOnM)),
            };
            WriteMs1Fixture(mzmlPath, candidates[0]);

            try
            {
                var results = Quantify(outputFolder, mzmlPath, new List<SpectralMatch> { Psm(mzmlPath, 1, candidates) }).FlashLfqResults;

                Assert.That(results, Is.Not.Null);
                var peptide = results.PeptideModifiedSequences.Values.Single();
                Assert.That(peptide.Sequence.Split('|'), Has.Length.EqualTo(2), "two sites; the two proteins do not add names");
                Assert.That(peptide.ProteinGroups.Select(g => g.ProteinGroupName), Is.EqualTo(new[] { "P1|P2" }));
                Assert.That(results.ProteinGroups["P1|P2"].GetIntensity(results.SpectraFiles.Single()), Is.GreaterThan(0));
            }
            finally
            {
                Directory.Delete(outputFolder, true);
            }
        }

        /// <summary>
        /// Control: ACMDMEK is in P1 and in P2, and each has a peptide of its own (GGGR, WWWR), so they are two groups.
        /// A match ambiguous in site could have come from either group; crediting it would be a guess, so it stays out.
        /// </summary>
        [Test]
        public static void AMatchAmbiguousInSiteAcrossTwoProteinGroupsIsLeftOut()
        {
            string outputFolder = NewOutputFolder("LocalizationAmbiguousQuant_TwoGroups");
            string mzmlPath = Path.Combine(outputFolder, "file.mzML");
            var p1 = Protein("ACMDMEKGGGR", "P1");
            var p2 = Protein("ACMDMEKWWWR", "P2");
            var ambiguous = Psm(mzmlPath, 1,
                Peptide(p1, 1, 7, (3, OxidationOnM)), Peptide(p1, 1, 7, (5, OxidationOnM)),
                Peptide(p2, 1, 7, (3, OxidationOnM)), Peptide(p2, 1, 7, (5, OxidationOnM)));
            var gggr = Peptide(p1, 8, 11);
            var wwwr = Peptide(p2, 8, 11);
            WriteMs1Fixture(mzmlPath, Peptide(p1, 1, 7, (3, OxidationOnM)), gggr, wwwr);

            try
            {
                var parameters = Quantify(outputFolder, mzmlPath,
                    new List<SpectralMatch> { ambiguous, Psm(mzmlPath, 2, gggr), Psm(mzmlPath, 3, wwwr) });

                Assert.That(parameters.FlashLfqResults.PeptideModifiedSequences.Keys, Is.EquivalentTo(new[] { "GGGR", "WWWR" }));
            }
            finally
            {
                Directory.Delete(outputFolder, true);
            }
        }

        /// <summary>
        /// Oxidation on M3 is identified with its site, and in the same file another scan is ambiguous between M3 and
        /// M5. The localized match measures the feature; the ambiguous one is left out there. Passing both put two
        /// peptides on one peak, and FlashLFQ zeroed both: on the islet set that cost 153 localized peptides 206 values.
        /// </summary>
        [Test]
        public static void AnAmbiguousMatchBesideALocalizedIsomerInItsFileIsLeftOutAndTheIsomerKeepsItsValue()
        {
            string outputFolder = NewOutputFolder("LocalizationAmbiguousQuant_BesideIsomer");
            string mzmlPath = Path.Combine(outputFolder, "file.mzML");
            var protein = Protein("ACMDMEKGGGR", "P1");
            var oxidizedAt3 = Peptide(protein, 1, 7, (3, OxidationOnM));
            var oxidizedAt5 = Peptide(protein, 1, 7, (5, OxidationOnM));
            WriteMs1Fixture(mzmlPath, oxidizedAt3);

            try
            {
                var results = Quantify(outputFolder, mzmlPath, new List<SpectralMatch>
                {
                    Psm(mzmlPath, 1, oxidizedAt3),
                    Psm(mzmlPath, 2, oxidizedAt3, oxidizedAt5),
                }).FlashLfqResults;

                var file = results.SpectraFiles.Single();
                Assert.That(results.PeptideModifiedSequences.Keys, Is.EquivalentTo(new[] { oxidizedAt3.FullSequence }));
                Assert.That(results.PeptideModifiedSequences[oxidizedAt3.FullSequence].GetDetectionType(file), Is.EqualTo(FlashLFQ.DetectionType.MSMS));
                Assert.That(results.PeptideModifiedSequences[oxidizedAt3.FullSequence].GetIntensity(file), Is.GreaterThan(0));
            }
            finally
            {
                Directory.Delete(outputFolder, true);
            }
        }

        /// <summary>
        /// The rule is per file: M3 is localized in the first file only, and the second file has only the ambiguous
        /// match, so the second file's feature is quantified, under the ambiguous name.
        /// </summary>
        [Test]
        public static void AnAmbiguousMatchIsQuantifiedInAFileWhereNoIsomerWasLocalized()
        {
            string outputFolder = NewOutputFolder("LocalizationAmbiguousQuant_OtherFile");
            string first = Path.Combine(outputFolder, "first.mzML");
            string second = Path.Combine(outputFolder, "second.mzML");
            var protein = Protein("ACMDMEKGGGR", "P1");
            var oxidizedAt3 = Peptide(protein, 1, 7, (3, OxidationOnM));
            var oxidizedAt5 = Peptide(protein, 1, 7, (5, OxidationOnM));
            WriteMs1Fixture(first, oxidizedAt3);
            WriteMs1Fixture(second, oxidizedAt3);

            try
            {
                var results = Quantify(outputFolder, new List<string> { first, second }, new List<SpectralMatch>
                {
                    Psm(first, 1, oxidizedAt3),
                    Psm(second, 1, oxidizedAt3, oxidizedAt5),
                }).FlashLfqResults;

                string name = string.Join("|", new[] { oxidizedAt3.FullSequence, oxidizedAt5.FullSequence }.OrderBy(s => s, System.StringComparer.Ordinal));
                var files = results.SpectraFiles.ToDictionary(f => Path.GetFileName(f.FullFilePathWithExtension));
                Assert.That(results.PeptideModifiedSequences.Keys, Is.EquivalentTo(new[] { oxidizedAt3.FullSequence, name }));
                Assert.That(results.PeptideModifiedSequences[oxidizedAt3.FullSequence].GetIntensity(files["first.mzML"]), Is.GreaterThan(0));
                Assert.That(results.PeptideModifiedSequences[name].GetIntensity(files["second.mzML"]), Is.GreaterThan(0));
                Assert.That(results.PeptideModifiedSequences[name].ProteinGroups.Select(g => g.ProteinGroupName), Is.EqualTo(new[] { "P1" }),
                    "unique to P1's group, so it counts toward P1's value in the second file");
            }
            finally
            {
                Directory.Delete(outputFolder, true);
            }
        }

        /// <summary>
        /// One scan is ambiguous between M1 and M3, another between M3 and M5. Both are one base sequence with one
        /// oxidation, so they are one peptide, named by all three sites: two names would put two peptides on one peak.
        /// </summary>
        [Test]
        public static void AmbiguousMatchesOfOneIsomerAreOnePeptideNamedByEverySiteTheyLeaveOpen()
        {
            string outputFolder = NewOutputFolder("LocalizationAmbiguousQuant_OnePeptide");
            string mzmlPath = Path.Combine(outputFolder, "file.mzML");
            var protein = Protein("MDMEMKGGGR", "P1");
            var at1 = Peptide(protein, 1, 6, (1, OxidationOnM));
            var at3 = Peptide(protein, 1, 6, (3, OxidationOnM));
            var at5 = Peptide(protein, 1, 6, (5, OxidationOnM));
            WriteMs1Fixture(mzmlPath, at1);

            try
            {
                var results = Quantify(outputFolder, mzmlPath, new List<SpectralMatch>
                {
                    Psm(mzmlPath, 1, at1, at3),
                    Psm(mzmlPath, 2, at3, at5),
                }).FlashLfqResults;

                string name = string.Join("|", new[] { at1.FullSequence, at3.FullSequence, at5.FullSequence }.OrderBy(s => s, System.StringComparer.Ordinal));
                Assert.That(results.PeptideModifiedSequences.Keys, Is.EquivalentTo(new[] { name }));
                Assert.That(results.PeptideModifiedSequences[name].GetIntensity(results.SpectraFiles.Single()), Is.GreaterThan(0));
            }
            finally
            {
                Directory.Delete(outputFolder, true);
            }
        }

        /// <summary>
        /// A modification with no chemical formula (a mass-only one, as some database annotations are) is compared by its
        /// mass, so the same mass-only modification at two sites is still one peptide whose site is open.
        /// </summary>
        [Test]
        public static void AMassOnlyModificationAtTwoSitesIsLocalizationAmbiguous()
        {
            ModificationMotif.TryGetMotif("M", out var motif);
            var massOnly = new Modification(_originalId: "MassShift", _modificationType: "Test", _target: motif,
                _locationRestriction: "Anywhere.", _monoisotopicMass: 21.9819);
            Assume.That(massOnly.ChemicalFormula, Is.Null);
            var protein = Protein("ACMDMEKGGGR", "P1");

            var psm = Psm("f", 1, Peptide(protein, 1, 7, (3, massOnly)), Peptide(protein, 1, 7, (5, massOnly)));

            Assert.That(psm.IsLocalizationAmbiguous, Is.True);
            Assert.That(psm.GetSiteFreeIdentity(), Does.EndWith("21.98190"));
        }

        /// <summary>
        /// With contaminants not written, a peptide only contaminants identify is withheld from the quantification
        /// tables, and that holds for a localization-ambiguous one under its sites' name. A match ambiguous between two
        /// peptides has no such name and is passed over.
        /// </summary>
        [Test]
        public static void AContaminantOnlyAmbiguousRowIsWithheldWhenContaminantsAreNotWritten()
        {
            string outputFolder = NewOutputFolder("LocalizationAmbiguousQuant_Contaminant");
            string mzmlPath = Path.Combine(outputFolder, "file.mzML");
            var contaminant = new Protein("ACMDMEKGGGR", "CONTAMINANT_P1", null, new List<System.Tuple<string, string>>(),
                new Dictionary<int, List<Modification>>(), isContaminant: true);
            var target = Protein("WWWWWKAAAAR", "P2");
            var at3 = Peptide(contaminant, 1, 7, (3, OxidationOnM));
            var at5 = Peptide(contaminant, 1, 7, (5, OxidationOnM));
            var targetPeptide = Peptide(target, 1, 6);
            WriteMs1Fixture(mzmlPath, at3, targetPeptide);

            try
            {
                var results = Quantify(outputFolder, mzmlPath, new List<SpectralMatch>
                    {
                        Psm(mzmlPath, 1, at3, at5),
                        Psm(mzmlPath, 2, targetPeptide),
                        Psm(mzmlPath, 3, Peptide(target, 1, 6), Peptide(target, 7, 11)),
                    },
                    configure: search => search.WriteContaminants = false, scrubContaminants: true).FlashLfqResults;

                Assert.That(results.PeptideModifiedSequences.Keys, Is.EquivalentTo(new[] { targetPeptide.FullSequence }),
                    "the contaminant's ambiguous row is withheld; the target peptide is kept");
            }
            finally
            {
                Directory.Delete(outputFolder, true);
            }
        }

        /// <summary>
        /// A match whose protein forms no group (here the two-peptide rule removes P1's group, which has one peptide) has
        /// no group to count for, so a localization-ambiguous one stays out rather than being quantified under UNDEFINED.
        /// </summary>
        [Test]
        public static void AnAmbiguousMatchWhoseProteinFormsNoGroupIsLeftOut()
        {
            string outputFolder = NewOutputFolder("LocalizationAmbiguousQuant_NoGroup");
            string mzmlPath = Path.Combine(outputFolder, "file.mzML");
            var p1 = Protein("ACMDMEKGGGR", "P1");
            var p2 = Protein("WWWWWKAAAAR", "P2");
            var first = Peptide(p2, 1, 6);
            var second = Peptide(p2, 7, 11);
            WriteMs1Fixture(mzmlPath, Peptide(p1, 1, 7, (3, OxidationOnM)), first, second);

            try
            {
                var results = Quantify(outputFolder, mzmlPath, new List<SpectralMatch>
                    {
                        Psm(mzmlPath, 1, Peptide(p1, 1, 7, (3, OxidationOnM)), Peptide(p1, 1, 7, (5, OxidationOnM))),
                        Psm(mzmlPath, 2, first),
                        Psm(mzmlPath, 3, second),
                    },
                    configure: search => search.NoOneHitWonders = true).FlashLfqResults;

                Assert.That(results.PeptideModifiedSequences.Keys, Is.EquivalentTo(new[] { first.FullSequence, second.FullSequence }));
            }
            finally
            {
                Directory.Delete(outputFolder, true);
            }
        }

        /// <summary>The sites are the candidates' full sequences, each once, in ordinal order, whatever order they were found in.</summary>
        [Test]
        public static void TheSitesDoNotDependOnTheOrderTheCandidatesWereFound()
        {
            var p1 = Protein("ACMDMEKGGGR", "P1");
            var p2 = Protein("ACMDMEKGGGR", "P2");
            var at3 = Peptide(p1, 1, 7, (3, OxidationOnM));
            var at5 = Peptide(p1, 1, 7, (5, OxidationOnM));
            var at3InP2 = Peptide(p2, 1, 7, (3, OxidationOnM));

            var forward = Psm("f", 1, at3, at5, at3InP2).GetLocalizationCandidates();
            var backward = Psm("f", 2, at5, at3InP2, at3).GetLocalizationCandidates();

            Assert.That(forward, Is.EqualTo(new[] { at3.FullSequence, at5.FullSequence }.OrderBy(s => s, System.StringComparer.Ordinal)));
            Assert.That(backward, Is.EqualTo(forward));
        }

        /// <summary>
        /// Which matches are ambiguous only in where their modifications sit. A modification is judged by what it is
        /// (type, name without its residue, formula), not by its formula alone: Phosphorylation on S or on T qualifies (one
        /// modification, two residues); Oxidation on M against Hydroxylation on D does not (one oxygen, but two
        /// modifications, so the match is ambiguous in which one it carries); nor do two oxidations against one
        /// dioxidation (one total formula, different modifications).
        /// </summary>
        [Test]
        public static void OnlyAMatchWhoseCandidatesCarryTheSameModificationsIsLocalizationAmbiguous()
        {
            var protein = Protein("ACMDMEKGGGR", "P1");
            var phosphoProtein = Protein("ASDTEKGGGR", "P2");
            var unambiguous = Psm("f", 1, Peptide(protein, 1, 7, (3, OxidationOnM)));
            var oxidationOrHydroxylation = Psm("f", 2, Peptide(protein, 1, 7, (3, OxidationOnM)), Peptide(protein, 1, 7, (4, HydroxylationOnD)));
            var twoOxidationsOrOneDioxidation = Psm("f", 3,
                Peptide(protein, 1, 7, (3, OxidationOnM), (5, OxidationOnM)), Peptide(protein, 1, 7, (3, DioxidationOnM)));
            var twoPeptides = Psm("f", 4, Peptide(protein, 1, 7), Peptide(protein, 8, 11));
            var phosphoOnSOrT = Psm("f", 5, Peptide(phosphoProtein, 1, 6, (2, PhosphorylationOnS)), Peptide(phosphoProtein, 1, 6, (4, PhosphorylationOnT)));
            var localizedPhospho = Psm("f", 6, Peptide(phosphoProtein, 1, 6, (4, PhosphorylationOnT)));

            Assume.That(oxidationOrHydroxylation.FullSequence, Is.Null);
            Assume.That(oxidationOrHydroxylation.BioPolymerWithSetModsMonoisotopicMass, Is.Not.Null, "one mass: both are one oxygen");
            Assume.That(twoOxidationsOrOneDioxidation.BioPolymerWithSetModsMonoisotopicMass, Is.Not.Null, "one mass, one total formula");
            Assert.That(unambiguous.IsLocalizationAmbiguous, Is.False, "one full sequence: quantified under it, as before");
            Assert.That(unambiguous.GetLocalizationCandidates(), Is.Empty);
            Assert.That(oxidationOrHydroxylation.IsLocalizationAmbiguous, Is.False, "ambiguous in which modification, not only where");
            Assert.That(oxidationOrHydroxylation.GetSiteFreeIdentity(), Is.Null);
            Assert.That(twoOxidationsOrOneDioxidation.IsLocalizationAmbiguous, Is.False);
            Assert.That(twoPeptides.IsLocalizationAmbiguous, Is.False, "two base sequences");
            Assert.That(twoPeptides.GetSiteFreeIdentity(), Is.Null);
            Assert.That(phosphoOnSOrT.IsLocalizationAmbiguous, Is.True, "one modification on S or on T");
            Assert.That(phosphoOnSOrT.GetSiteFreeIdentity(), Is.EqualTo(localizedPhospho.GetSiteFreeIdentity()), "an isomer of the localized match");
        }

        /// <summary>
        /// The filter's includeLocalizationAmbiguousMods decides for matches ambiguous only in their sites, and
        /// includeAmbiguousMods for the other unresolved ones. Left null, it follows includeAmbiguousMods, so callers that
        /// do not pass it filter as before.
        /// </summary>
        [Test]
        public static void TheFilterKeepsLocalizationAmbiguousMatchesOnlyWhenAsked()
        {
            var protein = Protein("ACMDMEKGGGR", "P1");
            var resolved = Psm("f", 1, Peptide(protein, 1, 7, (3, OxidationOnM)));
            var siteOnly = Psm("f", 2, Peptide(protein, 1, 7, (3, OxidationOnM)), Peptide(protein, 1, 7, (5, OxidationOnM)));
            var whichModification = Psm("f", 3, Peptide(protein, 1, 7, (3, OxidationOnM)), Peptide(protein, 1, 7, (4, HydroxylationOnD)));
            var psms = new List<SpectralMatch> { resolved, siteOnly, whichModification };
            var common = new CommonParameters(digestionParams: Digestion);

            List<SpectralMatch> Kept(bool includeAmbiguousMods, bool? includeLocalizationAmbiguousMods) =>
                EngineLayer.SpectrumMatch.FilteredPsms.Filter(psms, common, includeAmbiguousMods: includeAmbiguousMods,
                    includeLocalizationAmbiguousMods: includeLocalizationAmbiguousMods).ToList();

            Assert.That(Kept(true, null), Is.EquivalentTo(psms), "null follows includeAmbiguousMods: true");
            Assert.That(Kept(false, null), Is.EquivalentTo(new[] { resolved }), "null follows includeAmbiguousMods: false");
            Assert.That(Kept(false, true), Is.EquivalentTo(new[] { resolved, siteOnly }));
            Assert.That(Kept(true, false), Is.EquivalentTo(new[] { resolved, whichModification }));
        }

        /// <summary>
        /// End to end: a match that cannot tell Oxidation on M from Hydroxylation on D (one oxygen either way) is ambiguous
        /// in which modification it carries, so it stays out of quantification, as every ambiguous match did before.
        /// </summary>
        [Test]
        public static void AMatchAmbiguousInWhichModificationItCarriesIsLeftOut()
        {
            string outputFolder = NewOutputFolder("LocalizationAmbiguousQuant_TwoModifications");
            string mzmlPath = Path.Combine(outputFolder, "file.mzML");
            var protein = Protein("ACMDMEKGGGR", "P1");
            var oxidized = Peptide(protein, 1, 7, (3, OxidationOnM));
            var hydroxylated = Peptide(protein, 1, 7, (4, HydroxylationOnD));
            var gggr = Peptide(protein, 8, 11);
            WriteMs1Fixture(mzmlPath, oxidized, gggr);

            try
            {
                var results = Quantify(outputFolder, mzmlPath, new List<SpectralMatch>
                {
                    Psm(mzmlPath, 1, oxidized, hydroxylated),
                    Psm(mzmlPath, 2, gggr),
                }).FlashLfqResults;

                Assert.That(results.PeptideModifiedSequences.Keys, Is.EquivalentTo(new[] { gggr.FullSequence }));
            }
            finally
            {
                Directory.Delete(outputFolder, true);
            }
        }

        /// <summary>
        /// The localized isomer is identified in the file, but its peptide fails the peptide-level q-value, so it is not
        /// reported. It then cannot stand for the feature, so the ambiguous match beside it is quantified there under the
        /// sites' name. (GGGR passes, so the list of peptides to report is not empty: an empty list means "report every
        /// peptide".)
        /// </summary>
        [Test]
        public static void AnAmbiguousMatchBesideAnUnreportedIsomerIsQuantifiedInItsFile()
        {
            string outputFolder = NewOutputFolder("LocalizationAmbiguousQuant_BesideUnreportedIsomer");
            string mzmlPath = Path.Combine(outputFolder, "file.mzML");
            var protein = Protein("ACMDMEKGGGR", "P1");
            var oxidizedAt3 = Peptide(protein, 1, 7, (3, OxidationOnM));
            var oxidizedAt5 = Peptide(protein, 1, 7, (5, OxidationOnM));
            var gggr = Peptide(protein, 8, 11);
            WriteMs1Fixture(mzmlPath, oxidizedAt3, gggr);
            var weakLocalized = Psm(mzmlPath, 1, oxidizedAt3);
            weakLocalized.PeptideFdrInfo = new FdrInfo { QValue = 0.5, QValueNotch = 0.5, PEP = 0.5, PEP_QValue = 0.5 };

            try
            {
                var results = Quantify(outputFolder, mzmlPath, new List<SpectralMatch>
                {
                    weakLocalized,
                    Psm(mzmlPath, 2, oxidizedAt3, oxidizedAt5),
                    Psm(mzmlPath, 3, gggr),
                }).FlashLfqResults;

                string name = string.Join("|", new[] { oxidizedAt3.FullSequence, oxidizedAt5.FullSequence }.OrderBy(s => s, System.StringComparer.Ordinal));
                var file = results.SpectraFiles.Single();
                Assert.That(results.PeptideModifiedSequences.Keys, Is.EquivalentTo(new[] { name, gggr.FullSequence }),
                    "the isomer fails the peptide-level q-value, so it is not reported");
                Assert.That(results.PeptideModifiedSequences[name].GetIntensity(file), Is.GreaterThan(0),
                    "the ambiguous match measures the feature in this file");
                Assert.That(results.PeptideModifiedSequences[name].GetDetectionType(file), Is.EqualTo(FlashLFQ.DetectionType.MSMS));
            }
            finally
            {
                Directory.Delete(outputFolder, true);
            }
        }

        /// <summary>
        /// SILAC rebuilds each match from its first candidate, which would give a localization-ambiguous match one site it
        /// never had. So with SILAC labels the ambiguous match stays out, and no row claims M3 or M5.
        /// </summary>
        [Test]
        public static void WithSilacLabelsALocalizationAmbiguousMatchIsLeftOut()
        {
            string outputFolder = NewOutputFolder("LocalizationAmbiguousQuant_Silac");
            string mzmlPath = Path.Combine(outputFolder, "file.mzML");
            var protein = Protein("ACMDMEKGGGR", "P1");
            var oxidizedAt3 = Peptide(protein, 1, 7, (3, OxidationOnM));
            var oxidizedAt5 = Peptide(protein, 1, 7, (5, OxidationOnM));
            var gggr = Peptide(protein, 8, 11);
            WriteMs1Fixture(mzmlPath, oxidizedAt3, gggr);

            try
            {
                var results = Quantify(outputFolder, mzmlPath, new List<SpectralMatch>
                    {
                        Psm(mzmlPath, 1, oxidizedAt3, oxidizedAt5),
                        Psm(mzmlPath, 2, gggr),
                    },
                    configure: search => search.SilacLabels = new List<SilacLabel> { new('K', 'a', "C{13}6", 6.020129) }).FlashLfqResults;

                Assert.That(results.PeptideModifiedSequences.Keys, Is.Not.Empty, "GGGR is quantified");
                Assert.That(results.PeptideModifiedSequences.Keys, Has.None.Contains("Oxidation"),
                    "no row carries a site the ambiguous match left open");
            }
            finally
            {
                Directory.Delete(outputFolder, true);
            }
        }

        /// <summary>
        /// Without parsimony each protein is its own group. A localization-ambiguous match whose candidates all come from
        /// P1 counts for P1; one whose candidates come from P1 and P2 (the same sequence in both) spans two groups and
        /// stays out, so the row names P1 alone.
        /// </summary>
        [Test]
        public static void WithoutParsimonyAnAmbiguousMatchCountsForItsOneProteinOnly()
        {
            string outputFolder = NewOutputFolder("LocalizationAmbiguousQuant_NoParsimony");
            string mzmlPath = Path.Combine(outputFolder, "file.mzML");
            var p1 = Protein("ACMDMEKGGGR", "P1");
            var p2 = Protein("ACMDMEKWWWR", "P2");
            WriteMs1Fixture(mzmlPath, Peptide(p1, 1, 7, (3, OxidationOnM)));

            try
            {
                var results = Quantify(outputFolder, mzmlPath, new List<SpectralMatch>
                    {
                        Psm(mzmlPath, 1, Peptide(p1, 1, 7, (3, OxidationOnM)), Peptide(p1, 1, 7, (5, OxidationOnM))),
                        Psm(mzmlPath, 2, Peptide(p1, 1, 7, (3, OxidationOnM)), Peptide(p1, 1, 7, (5, OxidationOnM)),
                            Peptide(p2, 1, 7, (3, OxidationOnM)), Peptide(p2, 1, 7, (5, OxidationOnM))),
                    },
                    configure: search => search.DoParsimony = false).FlashLfqResults;

                var peptide = results.PeptideModifiedSequences.Values.Single();
                Assert.That(peptide.Sequence.Split('|'), Has.Length.EqualTo(2));
                Assert.That(peptide.ProteinGroups.Select(g => g.ProteinGroupName), Is.EqualTo(new[] { "P1" }),
                    "the match spanning P1 and P2 is left out");
                Assert.That(peptide.GetIntensity(results.SpectraFiles.Single()), Is.GreaterThan(0));
            }
            finally
            {
                Directory.Delete(outputFolder, true);
            }
        }

        /// <summary>
        /// With match-between-runs on, decoys reach FlashLFQ (it needs them to score transfers), and a decoy's ambiguous
        /// match does too, under its sites' name, as a target's does. Only targets are reported, so the table holds the
        /// target's row alone.
        /// </summary>
        [Test]
        public static void WithMatchBetweenRunsADecoysAmbiguousMatchReachesFlashLfqButIsNotReported()
        {
            string outputFolder = NewOutputFolder("LocalizationAmbiguousQuant_Mbr");
            string first = Path.Combine(outputFolder, "first.mzML");
            string second = Path.Combine(outputFolder, "second.mzML");
            var target = Protein("ACMDMEKGGGR", "P1");
            var decoy = new Protein("WCMEMAKGGGR", "DECOY_P9", null, new List<System.Tuple<string, string>>(),
                new Dictionary<int, List<Modification>>(), isDecoy: true);
            var targetAt3 = Peptide(target, 1, 7, (3, OxidationOnM));
            var targetAt5 = Peptide(target, 1, 7, (5, OxidationOnM));
            var decoyAt3 = Peptide(decoy, 1, 7, (3, OxidationOnM));
            var decoyAt5 = Peptide(decoy, 1, 7, (5, OxidationOnM));
            WriteMs1Fixture(first, targetAt3, decoyAt3);
            WriteMs1Fixture(second, targetAt3, decoyAt3);

            try
            {
                // The target outscores the decoy, so its row passes the peptide-level q-value.
                var results = Quantify(outputFolder, new List<string> { first, second }, new List<SpectralMatch>
                    {
                        ScoredPsm(20, first, 1, targetAt3, targetAt5),
                        Psm(first, 2, decoyAt3, decoyAt5),
                    },
                    configure: search => search.MatchBetweenRuns = true).FlashLfqResults;

                string name = string.Join("|", new[] { targetAt3.FullSequence, targetAt5.FullSequence }.OrderBy(s => s, System.StringComparer.Ordinal));
                string decoyName = string.Join("|", new[] { decoyAt3.FullSequence, decoyAt5.FullSequence }.OrderBy(s => s, System.StringComparer.Ordinal));
                var files = results.SpectraFiles.ToDictionary(f => Path.GetFileName(f.FullFilePathWithExtension));
                Assert.That(results.PeptideModifiedSequences.Keys, Is.EquivalentTo(new[] { name }), "the decoy's row is not reported");
                Assert.That(results.PeptideModifiedSequences[name].GetIntensity(files["first.mzML"]), Is.GreaterThan(0));
                Assert.That(results.Peaks[files["first.mzML"]].SelectMany(peak => peak.Identifications).Select(id => id.ModifiedSequence),
                    Does.Contain(decoyName), "the decoy's ambiguous match reaches FlashLFQ, under its sites' name");
            }
            finally
            {
                Directory.Delete(outputFolder, true);
            }
        }

        /// <summary>
        /// A row of a localization-ambiguous peptide is held to the peptide-level q-value, as a localized peptide is
        /// (user, 2026-10-25). Each set of such matches counts as one peptide beside the full-sequence peptides, targets
        /// and decoys alike. Three decoy peptides score 20. ACMDMEK's set scores 25 and YMDMEK's 22, both above them, so
        /// their peptide-level q-values are 0 and both are reported: each set is its own peptide, not one entry for every
        /// ambiguous match. WMDMEK's set scores 5, below the decoys: its match passed the PSM-level q-value, but as a
        /// peptide it ranks under three decoys (q-value 1), so it is not reported, as a localized peptide there would not be.
        /// </summary>
        [Test]
        public static void AnAmbiguousRowIsReportedOnlyWhenItPassesThePeptideLevelQValue()
        {
            string outputFolder = NewOutputFolder("LocalizationAmbiguousQuant_PeptideQValue");
            string mzmlPath = Path.Combine(outputFolder, "file.mzML");
            var strongProtein = Protein("ACMDMEKGGGR", "P1");
            var secondProtein = Protein("YMDMEKR", "P3");
            var weakProtein = Protein("WMDMEKR", "P2");
            var decoy = new Protein("FFFFRLLLLKHHHHR", "DECOY_P9", null, new List<System.Tuple<string, string>>(),
                new Dictionary<int, List<Modification>>(), isDecoy: true);
            var strongAt3 = Peptide(strongProtein, 1, 7, (3, OxidationOnM));
            var strongAt5 = Peptide(strongProtein, 1, 7, (5, OxidationOnM));
            var secondAt2 = Peptide(secondProtein, 1, 6, (2, OxidationOnM));
            var secondAt4 = Peptide(secondProtein, 1, 6, (4, OxidationOnM));
            var weakAt2 = Peptide(weakProtein, 1, 6, (2, OxidationOnM));
            var weakAt4 = Peptide(weakProtein, 1, 6, (4, OxidationOnM));
            var gggr = Peptide(strongProtein, 8, 11);
            WriteMs1Fixture(mzmlPath, strongAt3, secondAt2, weakAt2, gggr);

            try
            {
                var weak = ScoredPsm(5, mzmlPath, 2, weakAt2, weakAt4);
                Assume.That(weak.PsmFdrInfo.QValue, Is.Zero, "the weak match passes at the PSM level");

                var results = Quantify(outputFolder, mzmlPath, new List<SpectralMatch>
                {
                    ScoredPsm(25, mzmlPath, 1, strongAt3, strongAt5),
                    weak,
                    ScoredPsm(30, mzmlPath, 3, gggr),
                    ScoredPsm(20, mzmlPath, 4, Peptide(decoy, 1, 5)),
                    ScoredPsm(20, mzmlPath, 5, Peptide(decoy, 6, 10)),
                    ScoredPsm(20, mzmlPath, 6, Peptide(decoy, 11, 15)),
                    ScoredPsm(22, mzmlPath, 7, secondAt2, secondAt4),
                }).FlashLfqResults;

                string strongName = string.Join("|", new[] { strongAt3.FullSequence, strongAt5.FullSequence }.OrderBy(s => s, System.StringComparer.Ordinal));
                string secondName = string.Join("|", new[] { secondAt2.FullSequence, secondAt4.FullSequence }.OrderBy(s => s, System.StringComparer.Ordinal));
                Assert.That(results.PeptideModifiedSequences.Keys, Is.EquivalentTo(new[] { strongName, secondName, gggr.FullSequence }),
                    "WMDMEK's set fails the peptide-level q-value");
                Assert.That(results.PeptideModifiedSequences[strongName].GetIntensity(results.SpectraFiles.Single()), Is.GreaterThan(0));
                Assert.That(weak.PeptideFdrInfo.QValue, Is.Zero, "the stored peptide-level values are left as they were");
            }
            finally
            {
                Directory.Delete(outputFolder, true);
            }
        }
    }
}
