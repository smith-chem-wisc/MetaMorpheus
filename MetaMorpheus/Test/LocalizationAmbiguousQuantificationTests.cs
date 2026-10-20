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
        private static SpectralMatch Psm(string filePath, int scanNumber, params PeptideWithSetModifications[] peptides)
        {
            var common = new CommonParameters(digestionParams: Digestion);
            SpectralMatch psm = new PeptideSpectralMatch(peptides[0], 0, 10, scanNumber, Scan(filePath, scanNumber), common, new List<MatchedFragmentIon>());
            foreach (var peptide in peptides.Skip(1))
            {
                psm.AddOrReplace(peptide, 10, 0, true, new List<MatchedFragmentIon>());
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
        private static PostSearchAnalysisParameters Quantify(string outputFolder, string mzmlPath, List<SpectralMatch> psms) =>
            Quantify(outputFolder, new List<string> { mzmlPath }, psms);

        private static PostSearchAnalysisParameters Quantify(string outputFolder, List<string> mzmlPaths, List<SpectralMatch> psms)
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
            foreach (string step in new[] { "ProteinAnalysis", "QuantificationAnalysis" })
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
        /// Which matches are ambiguous only in where their modifications sit. The same modification is judged by
        /// chemical formula, one per modification: an oxygen on M or on D qualifies (one formula, two residues), two
        /// oxidations against one dioxidation do not (one total formula, different modifications).
        /// </summary>
        [Test]
        public static void OnlyAMatchWhoseCandidatesCarryTheSameModificationsIsLocalizationAmbiguous()
        {
            var protein = Protein("ACMDMEKGGGR", "P1");
            var unambiguous = Psm("f", 1, Peptide(protein, 1, 7, (3, OxidationOnM)));
            var oxygenOnMOrD = Psm("f", 2, Peptide(protein, 1, 7, (3, OxidationOnM)), Peptide(protein, 1, 7, (4, HydroxylationOnD)));
            var twoOxidationsOrOneDioxidation = Psm("f", 3,
                Peptide(protein, 1, 7, (3, OxidationOnM), (5, OxidationOnM)), Peptide(protein, 1, 7, (3, DioxidationOnM)));
            var twoPeptides = Psm("f", 4, Peptide(protein, 1, 7), Peptide(protein, 8, 11));

            Assume.That(twoOxidationsOrOneDioxidation.BioPolymerWithSetModsMonoisotopicMass, Is.Not.Null, "one mass, one total formula");
            Assert.That(unambiguous.IsLocalizationAmbiguous, Is.False, "one full sequence: quantified under it, as before");
            Assert.That(oxygenOnMOrD.IsLocalizationAmbiguous, Is.True);
            Assert.That(oxygenOnMOrD.GetSiteFreeIdentity(), Is.EqualTo(unambiguous.GetSiteFreeIdentity()), "an isomer of the localized match");
            Assert.That(twoOxidationsOrOneDioxidation.IsLocalizationAmbiguous, Is.False);
            Assert.That(twoPeptides.IsLocalizationAmbiguous, Is.False, "two base sequences");
        }
    }
}
