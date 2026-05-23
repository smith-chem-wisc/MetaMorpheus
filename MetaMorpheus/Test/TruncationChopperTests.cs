using System.Collections.Generic;
using System.Linq;
using System.Text;
using Chemistry;
using EngineLayer;
using EngineLayer.Truncation;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Omics.Digestion;
using Omics.Fragmentation;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;

namespace Test
{
    /// <summary>
    /// Phase 2 tests: residue chopping (#9, #10) and Pass 3 standard scoring (#12, #13). Extends the
    /// Phase 1 synthetic style with modified parents (P2 phospho, P3 N-terminal acetyl) so PTM-aware
    /// chopping is exercised.
    /// </summary>
    [TestFixture]
    public class TruncationChopperTests
    {
        private const string AlphabetP1 = "ACDEFGHIKLNPQRSTVWY";
        private const string AlphabetP2 = "YWVTSRQPNLKIHGFEDCA"; // residue 24 is 'S' (phospho site)
        private const string AlphabetP3 = "ADEKRHSTNQGVLIFYWPC"; // residue 1 is 'A' (acetyl site)

        private CommonParameters _cp;
        private DotMassDiffAcceptor _exactAcceptor;
        private Modification _phospho;
        private Modification _acetyl;

        private PeptideWithSetModifications _p1; // 50 residues, no mods
        private PeptideWithSetModifications _p2; // 80 residues, phospho on residue 24
        private PeptideWithSetModifications _p3; // 30 residues, N-terminal acetyl on residue 1 ('A')

        [OneTimeSetUp]
        public void BuildFixture()
        {
            _cp = new CommonParameters(
                dissociationType: DissociationType.HCD,
                precursorMassTolerance: new PpmTolerance(10),
                productMassTolerance: new PpmTolerance(20));

            _exactAcceptor = new DotMassDiffAcceptor("exact", new[] { 0.0 }, new PpmTolerance(10));

            _phospho = MakePhospho();
            _acetyl = MakeAcetyl();

            _p1 = MakeProteoform(RepeatTo(AlphabetP1, 50), "P1", null);
            _p2 = MakeProteoform(RepeatTo(AlphabetP2, 80), "P2", new Dictionary<int, Modification> { { 25, _phospho } }); // residue 24 -> key 25
            _p3 = MakeProteoform(RepeatTo(AlphabetP3, 30), "P3", new Dictionary<int, Modification> { { 1, _acetyl } });  // N-terminal
        }

        [Test]
        public void ChopFromCTerm_RecoversS2_CorrectStartEnd()
        {
            double target = MakeProteoform(RepeatTo(AlphabetP1, 50).Substring(0, 40), "x", null).MonoisotopicMass;

            ChopResult result = ProteoformChopper.ChopUntilMassMatches(_p1, FragmentationTerminus.C, target, _exactAcceptor);

            Assert.That(result, Is.Not.Null);
            Assert.That(result.ResiduesChopped, Is.EqualTo(10));
            Assert.That(result.TruncatedForm.OneBasedStartResidueInProtein, Is.EqualTo(1));
            Assert.That(result.TruncatedForm.OneBasedEndResidueInProtein, Is.EqualTo(40));
        }

        [Test]
        public void ChopFromNTerm_RecoversS3_PhosphoRetained()
        {
            // P2[6..80] keeps the phospho (residue 24 is well inside the retained region).
            var expectedMods = new Dictionary<int, Modification> { { 20, _phospho } }; // residue 19 of the 75-mer
            double target = MakeProteoform(RepeatTo(AlphabetP2, 80).Substring(5), "x", expectedMods).MonoisotopicMass;

            ChopResult result = ProteoformChopper.ChopUntilMassMatches(_p2, FragmentationTerminus.N, target, _exactAcceptor);

            Assert.That(result, Is.Not.Null);
            Assert.That(result.ResiduesChopped, Is.EqualTo(5));
            Assert.That(result.TruncatedForm.OneBasedStartResidueInProtein, Is.EqualTo(6));
            Assert.That(result.TruncatedForm.OneBasedEndResidueInProtein, Is.EqualTo(80));
            var mod = result.TruncatedForm.AllModsOneIsNterminus.Single();
            Assert.That(mod.Key, Is.EqualTo(20));               // re-indexed onto the retained residue
            Assert.That(mod.Value.OriginalId, Is.EqualTo("phospho"));
        }

        [Test]
        public void ChopFromNTerm_RecoversS4_AcetylAndAlaRemovedTogether()
        {
            // P3[2..30] with the N-terminal acetyl-Ala removed as one step (no residual mod).
            double target = MakeProteoform(RepeatTo(AlphabetP3, 30).Substring(1), "x", null).MonoisotopicMass;

            ChopResult result = ProteoformChopper.ChopUntilMassMatches(_p3, FragmentationTerminus.N, target, _exactAcceptor);

            Assert.That(result, Is.Not.Null);
            Assert.That(result.ResiduesChopped, Is.EqualTo(1));
            Assert.That(result.TruncatedForm.OneBasedStartResidueInProtein, Is.EqualTo(2));
            Assert.That(result.TruncatedForm.OneBasedEndResidueInProtein, Is.EqualTo(30));
            Assert.That(result.TruncatedForm.AllModsOneIsNterminus, Is.Empty); // acetyl left with the chopped Ala
            Assert.That(result.TruncatedForm.MonoisotopicMass, Is.EqualTo(target).Within(0.01));
        }

        [Test]
        public void NoMatchScan_ReturnsNull()
        {
            // Half a dalton below the parent mass: no integer-residue chop can reach it.
            double unreachable = _p1.MonoisotopicMass - 0.5;

            ChopResult result = ProteoformChopper.ChopUntilMassMatches(_p1, FragmentationTerminus.C, unreachable, _exactAcceptor);

            Assert.That(result, Is.Null);
        }

        [Test]
        public void NotchOnePlus_AcceptedWithinTolerance()
        {
            const double c13Spacing = 1.0033548378;
            var notchAcceptor = new DotMassDiffAcceptor("1mm", new[] { 0.0, c13Spacing }, new PpmTolerance(10));

            double truncatedMass = MakeProteoform(RepeatTo(AlphabetP1, 50).Substring(0, 40), "x", null).MonoisotopicMass;
            double target = truncatedMass + c13Spacing; // M_obs = M_truncated_theo + 1.00335

            ChopResult result = ProteoformChopper.ChopUntilMassMatches(_p1, FragmentationTerminus.C, target, notchAcceptor);

            Assert.That(result, Is.Not.Null);
            Assert.That(result.ResiduesChopped, Is.EqualTo(10));
            Assert.That(result.Notch, Is.Not.EqualTo(0)); // matched at the +1 isotope notch, not notch 0
        }

        [Test]
        public void TwoParentsCollapseToSameTruncated()
        {
            var oxidation = MakeOxidation();
            var parent1 = MakeProteoform("PEPTIDEMK", "A1", null);                                            // -> PEPTIDE
            var parent2 = MakeProteoform("PEPTIDEMK", "A2", new Dictionary<int, Modification> { { 9, oxidation } }); // M(ox) at residue 8 -> PEPTIDE

            double peptideMass = MakeProteoform("PEPTIDE", "x", null).MonoisotopicMass;
            Ms2ScanWithSpecificMass scan = BuildScan(1, peptideMass, new List<double> { 500.0 }, _cp);

            TruncationPsm psm1 = TruncationPass3.ScoreTruncation(WinnerChoppingCTerm(parent1, "A1", scan), _cp, _exactAcceptor);
            TruncationPsm psm2 = TruncationPass3.ScoreTruncation(WinnerChoppingCTerm(parent2, "A2", scan), _cp, _exactAcceptor);

            Assert.That(psm1, Is.Not.Null);
            Assert.That(psm2, Is.Not.Null);
            Assert.That(psm1.TruncatedForm.BaseSequence, Is.EqualTo("PEPTIDE"));
            Assert.That(psm2.TruncatedForm.BaseSequence, Is.EqualTo("PEPTIDE"));

            List<SpectralMatch> collapsed = TruncationPass3.CollapseDuplicateTruncations(new[] { psm1, psm2 });

            Assert.That(collapsed.Count, Is.EqualTo(1));
            collapsed[0].ResolveAllAmbiguities();
            var accessions = collapsed[0].BestMatchingBioPolymersWithSetMods
                .Select(b => b.SpecificBioPolymer.Parent.Accession).Distinct().OrderBy(a => a).ToList();
            Assert.That(accessions, Is.EqualTo(new[] { "A1", "A2" })); // both parents preserved (#10)
        }

        [Test]
        public void Pass3StandardScore_UsesBothSeries()
        {
            var truncated = MakeProteoform(RepeatTo(AlphabetP1, 50).Substring(0, 40), "P1", null);
            Ms2ScanWithSpecificMass scan = BuildScan(1, truncated.MonoisotopicMass, SeriesMasses(truncated, FragmentationTerminus.Both, _cp), _cp);

            var bothSeriesProducts = new List<Product>();
            truncated.Fragment(_cp.DissociationType, FragmentationTerminus.Both, bothSeriesProducts, _cp.FragmentationParameters);
            double expectedScore = MetaMorpheusEngine.CalculatePeptideScore(scan.TheScan,
                MetaMorpheusEngine.MatchFragmentIons(scan, bothSeriesProducts, _cp));

            TruncationPsm psm = TruncationPass3.ScoreTruncation(WinnerChoppingCTerm(_p1, "P1", scan), _cp, _exactAcceptor);

            Assert.That(psm, Is.Not.Null);
            Assert.That(psm.Score, Is.EqualTo(expectedScore).Within(1e-6));
        }

        [Test]
        public void ChopFromCTerm_RetainsNTerminalMod()
        {
            // Chopping the C-terminus must keep the N-terminal acetyl.
            double target = MakeProteoform(RepeatTo(AlphabetP3, 30).Substring(0, 25), "x",
                new Dictionary<int, Modification> { { 1, _acetyl } }).MonoisotopicMass;

            ChopResult result = ProteoformChopper.ChopUntilMassMatches(_p3, FragmentationTerminus.C, target, _exactAcceptor);

            Assert.That(result, Is.Not.Null);
            Assert.That(result.ResiduesChopped, Is.EqualTo(5));
            Assert.That(result.TruncatedForm.OneBasedStartResidueInProtein, Is.EqualTo(1));
            Assert.That(result.TruncatedForm.OneBasedEndResidueInProtein, Is.EqualTo(25));
            Assert.That(result.TruncatedForm.AllModsOneIsNterminus.ContainsKey(1), Is.True);
            Assert.That(result.TruncatedForm.AllModsOneIsNterminus[1].OriginalId, Is.EqualTo("acetyl"));
        }

        [Test]
        public void ChopFromNTerm_RetainsCTerminalMod()
        {
            var cTermMod = MakeCTerminalMod();
            var parent = MakeProteoform("PEPTIDES", "C1", new Dictionary<int, Modification> { { 10, cTermMod } }); // C-term key = 8 + 2
            double target = MakeProteoform("PTIDES", "x", new Dictionary<int, Modification> { { 8, cTermMod } }).MonoisotopicMass; // 6-mer C-term key = 8

            ChopResult result = ProteoformChopper.ChopUntilMassMatches(parent, FragmentationTerminus.N, target, _exactAcceptor);

            Assert.That(result, Is.Not.Null);
            Assert.That(result.ResiduesChopped, Is.EqualTo(2));
            Assert.That(result.TruncatedForm.AllModsOneIsNterminus.ContainsKey(8), Is.True);
            Assert.That(result.TruncatedForm.AllModsOneIsNterminus[8].OriginalId, Is.EqualTo("cterm-amide"));
        }

        [Test]
        public void ScoreTruncation_ReturnsNull_WhenNoChopMatches()
        {
            Ms2ScanWithSpecificMass scan = BuildScan(1, _p1.MonoisotopicMass - 0.5, new List<double> { 500.0 }, _cp);

            TruncationPsm psm = TruncationPass3.ScoreTruncation(WinnerChoppingCTerm(_p1, "P1", scan), _cp, _exactAcceptor);

            Assert.That(psm, Is.Null);
        }

        [Test]
        public void ScoreTruncation_LabelsNTerminalTruncation()
        {
            // Winning C-terminal-ion series -> N-terminus chopped -> "N-terminal truncation".
            double target = MakeProteoform(RepeatTo(AlphabetP2, 80).Substring(5), "x",
                new Dictionary<int, Modification> { { 20, _phospho } }).MonoisotopicMass;
            Ms2ScanWithSpecificMass scan = BuildScan(1, target, new List<double> { 500.0 }, _cp);

            TruncationPsm psm = TruncationPass3.ScoreTruncation(WinnerChoppingNTerm(_p2, "P2", scan), _cp, _exactAcceptor);

            Assert.That(psm, Is.Not.Null);
            Assert.That(psm.TruncationProductType, Is.EqualTo(TruncationPass3.NTerminalTruncation));
            Assert.That(psm.TruncatedForm.OneBasedStartResidueInProtein, Is.EqualTo(6));
        }

        [Test]
        public void ChopFromNTerm_SingleInitiatorMet_LabeledMetExcision()
        {
            // Clean removal of the initiator Met (1 residue, parent starts at protein pos 1 with M) is NME,
            // not a 1-residue N-terminal truncation (#13).
            var parent = MakeProteoform("MDEKRHSTNQGVLIF", "MX", null);
            double target = MakeProteoform("DEKRHSTNQGVLIF", "x", null).MonoisotopicMass;

            ChopResult result = ProteoformChopper.ChopUntilMassMatches(parent, FragmentationTerminus.N, target, _exactAcceptor);

            Assert.That(result, Is.Not.Null);
            Assert.That(result.ResiduesChopped, Is.EqualTo(1));
            Assert.That(result.TruncatedForm.Description, Does.StartWith(TruncationPass3.NTerminalMetExcision));
            Assert.That(result.TruncatedForm.Description, Does.Not.StartWith(TruncationPass3.NTerminalTruncation));
        }

        [Test]
        public void ChopFromNTerm_SingleNonMet_IsNTerminalTruncation()
        {
            // One residue chopped from the N-terminus, but it is not the initiator Met -> real truncation.
            var parent = MakeProteoform("ADEKRHSTNQGVLIF", "AX", null);
            double target = MakeProteoform("DEKRHSTNQGVLIF", "x", null).MonoisotopicMass;

            ChopResult result = ProteoformChopper.ChopUntilMassMatches(parent, FragmentationTerminus.N, target, _exactAcceptor);

            Assert.That(result, Is.Not.Null);
            Assert.That(result.ResiduesChopped, Is.EqualTo(1));
            Assert.That(result.TruncatedForm.Description, Does.StartWith(TruncationPass3.NTerminalTruncation));
        }

        [Test]
        public void ChopFromNTerm_MetPlusMoreResidues_IsNTerminalTruncation()
        {
            // Initiator Met plus additional residues removed -> a real N-terminal truncation, not NME.
            var parent = MakeProteoform("MDEKRHSTNQGVLIF", "MX", null);
            double target = MakeProteoform("HSTNQGVLIF", "x", null).MonoisotopicMass; // residues 6-15 (chopped M,D,E,K,R)

            ChopResult result = ProteoformChopper.ChopUntilMassMatches(parent, FragmentationTerminus.N, target, _exactAcceptor);

            Assert.That(result, Is.Not.Null);
            Assert.That(result.ResiduesChopped, Is.EqualTo(5));
            Assert.That(result.TruncatedForm.Description, Does.StartWith(TruncationPass3.NTerminalTruncation));
        }

        [Test]
        public void ChopFromNTerm_MetExcisionWithAcetyl_LabeledMetExcisionPlusAcetyl()
        {
            // Parent M-A(acetyl)-...; chopping the initiator Met yields the Met-excised, N-acetylated form.
            var acetyl = MakeAcetyl();
            var parent = MakeProteoform("MADEKRHSTNQGVLIF", "MAcX",
                new Dictionary<int, Modification> { { 3, acetyl } });   // acetyl on parent residue 2 (A) -> key 3
            double target = MakeProteoform("ADEKRHSTNQGVLIF", "x",
                new Dictionary<int, Modification> { { 2, acetyl } }).MonoisotopicMass; // acetyl on residue 1 -> key 2

            ChopResult result = ProteoformChopper.ChopUntilMassMatches(parent, FragmentationTerminus.N, target, _exactAcceptor);

            Assert.That(result, Is.Not.Null);
            Assert.That(result.ResiduesChopped, Is.EqualTo(1));
            Assert.That(result.TruncatedForm.Description, Does.StartWith(TruncationPass3.NTerminalMetExcisionPlusAcetyl));
        }

        [Test]
        public void ChopInternal_RecoversInternalSpan_LabeledInternal()
        {
            // P1[6..40]: 5 residues chopped from the N-terminus AND 10 from the C-terminus -> internal fragment.
            double target = MakeProteoform(RepeatTo(AlphabetP1, 50).Substring(5, 35), "x", null).MonoisotopicMass;

            List<ChopResult> results = ProteoformChopper.ChopInternalCandidates(_p1, target, _exactAcceptor);
            ChopResult hit = results.FirstOrDefault(r => r.TruncatedForm.OneBasedStartResidueInProtein == 6
                                                       && r.TruncatedForm.OneBasedEndResidueInProtein == 40);

            Assert.That(hit, Is.Not.Null, "internal span 6-40 not recovered");
            Assert.That(hit.ResiduesChopped, Is.EqualTo(15));        // 5 from N + 10 from C
            Assert.That(hit.Notch, Is.EqualTo(0));
            Assert.That(hit.TruncatedForm.Description, Does.StartWith(TruncationPass3.InternalTruncation));
        }

        // ---------- helpers ----------

        // A winning N-terminal-ion series means Pass 3 chops the C-terminus.
        private static TruncationParentSelection WinnerChoppingCTerm(PeptideWithSetModifications parent, string accession, Ms2ScanWithSpecificMass scan) =>
            new()
            {
                ScanIndex = 0,
                Scan = scan,
                Outcome = TruncationScanOutcome.Winner,
                WinningParent = new TruncationParent(parent, accession, accession, false),
                WinningSeries = FragmentationTerminus.N
            };

        // A winning C-terminal-ion series means Pass 3 chops the N-terminus.
        private static TruncationParentSelection WinnerChoppingNTerm(PeptideWithSetModifications parent, string accession, Ms2ScanWithSpecificMass scan) =>
            new()
            {
                ScanIndex = 0,
                Scan = scan,
                Outcome = TruncationScanOutcome.Winner,
                WinningParent = new TruncationParent(parent, accession, accession, false),
                WinningSeries = FragmentationTerminus.C
            };

        private static Modification MakeCTerminalMod()
        {
            ModificationMotif.TryGetMotif("S", out var motif);
            return new Modification(_originalId: "cterm-amide", _modificationType: "testMod", _target: motif,
                _chemicalFormula: ChemicalFormula.ParseFormula("H1"), _locationRestriction: "C-terminal.");
        }

        private static Modification MakePhospho()
        {
            ModificationMotif.TryGetMotif("S", out var motif);
            return new Modification(_originalId: "phospho", _modificationType: "testMod", _target: motif,
                _chemicalFormula: ChemicalFormula.ParseFormula("H1O3P1"), _locationRestriction: "Anywhere.");
        }

        private static Modification MakeAcetyl()
        {
            ModificationMotif.TryGetMotif("A", out var motif);
            return new Modification(_originalId: "acetyl", _modificationType: "testMod", _target: motif,
                _chemicalFormula: ChemicalFormula.ParseFormula("C2H2O1"), _locationRestriction: "N-terminal.");
        }

        private static Modification MakeOxidation()
        {
            ModificationMotif.TryGetMotif("M", out var motif);
            return new Modification(_originalId: "oxidation", _modificationType: "testMod", _target: motif,
                _chemicalFormula: ChemicalFormula.ParseFormula("O1"), _locationRestriction: "Anywhere.");
        }

        private static string RepeatTo(string alphabet, int length)
        {
            var sb = new StringBuilder(length);
            while (sb.Length < length)
            {
                sb.Append(alphabet);
            }
            return sb.ToString().Substring(0, length);
        }

        private static PeptideWithSetModifications MakeProteoform(string sequence, string accession, Dictionary<int, Modification> mods)
        {
            var protein = new Protein(sequence, accession);
            var digestionParams = new DigestionParams(protease: "top-down", minPeptideLength: 1, maxPeptideLength: 100000);
            return new PeptideWithSetModifications(protein, digestionParams, 1, sequence.Length,
                CleavageSpecificity.Full, "top-down", 0, mods ?? new Dictionary<int, Modification>(), mods?.Count ?? 0);
        }

        private static List<double> SeriesMasses(PeptideWithSetModifications form, FragmentationTerminus terminus, CommonParameters cp)
        {
            var products = new List<Product>();
            form.Fragment(cp.DissociationType, terminus, products, cp.FragmentationParameters);
            return products.Where(p => !double.IsNaN(p.NeutralMass)).Select(p => p.NeutralMass).ToList();
        }

        private static Ms2ScanWithSpecificMass BuildScan(int scanNumber, double precursorMass, IReadOnlyList<double> fragmentNeutralMasses, CommonParameters cp)
        {
            const double intensity = 1000.0;
            var ordered = fragmentNeutralMasses.OrderBy(m => m).ToList();

            double[] mz = ordered.Select(m => m.ToMz(1)).ToArray();
            double[] intensities = ordered.Select(_ => intensity).ToArray();
            var spectrum = new MzSpectrum(mz, intensities, false);

            double tic = intensities.Sum();
            var msDataScan = new MsDataScan(
                massSpectrum: spectrum, oneBasedScanNumber: scanNumber, msnOrder: 2, isCentroid: true,
                polarity: Polarity.Positive, retentionTime: scanNumber, scanWindowRange: new MzRange(0, 1_000_000),
                scanFilter: "f", mzAnalyzer: MZAnalyzerType.Orbitrap, totalIonCurrent: tic <= 0 ? 1 : tic,
                injectionTime: 1.0, noiseData: null, nativeId: $"scan={scanNumber}");

            var envelopes = ordered
                .Select(m => new IsotopicEnvelope(new List<(double mz, double intensity)> { (m.ToMz(1), intensity) }, m, 1, intensity, 0))
                .ToArray();

            return new Ms2ScanWithSpecificMass(msDataScan, precursorMass.ToMz(1), 1, "synthetic", cp,
                neutralExperimentalFragments: envelopes);
        }
    }
}
