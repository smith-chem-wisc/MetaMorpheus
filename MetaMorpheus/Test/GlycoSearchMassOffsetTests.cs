using EngineLayer;
using EngineLayer.GlycoSearch;
using EngineLayer.Indexing;
using Nett;
using NUnit.Framework;
using Omics.Modifications;
using Proteomics;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using TaskLayer;
using UsefulProteomicsDatabases;

namespace Test
{
    /// <summary>
    /// Glyco search can search a precursor mass offset (a missed monoisotopic peak) the way the Classic and Modern searches do:
    /// each offset of the chosen MassDiffAcceptorType is tried as the observed precursor mass minus the offset, and the match
    /// keeps the offset's notch. A large glycopeptide's monoisotopic peak is small, so its precursor is often picked a 13C peak heavy.
    /// </summary>
    [TestFixture]
    public static class GlycoSearchMassOffsetTests
    {
        private static readonly double C13 = Chemistry.Constants.C13MinusC12;

        [Test]
        public static void PrecursorMassOffsetsOf_ReadsEachNotchAndItsOffset_ExactFirst()
        {
            var tolerance = new MzLibUtil.PpmTolerance(10);

            Assert.That(GlycoSearchEngine.PrecursorMassOffsetsOf(null), Is.EqualTo(new[] { (0, 0.0) }), "unset searches the observed precursor mass only");

            var exact = GlycoSearchEngine.PrecursorMassOffsetsOf(SearchTask.GetMassDiffAcceptor(tolerance, MassDiffAcceptorType.Exact, null));
            Assert.That(exact.Length, Is.EqualTo(1));
            Assert.That(exact[0].Notch, Is.EqualTo(0));
            Assert.That(exact[0].Offset, Is.EqualTo(0).Within(1e-9));

            var plusMinusThree = GlycoSearchEngine.PrecursorMassOffsetsOf(SearchTask.GetMassDiffAcceptor(tolerance, MassDiffAcceptorType.PlusOrMinusThreeMM, null));
            Assert.That(plusMinusThree.Select(o => o.Offset), Is.EqualTo(new[] { 0, C13, -C13, 2 * C13, -2 * C13, 3 * C13, -3 * C13 }).Within(1e-9),
                "exact first, then by size, the positive offset before the negative one");
            Assert.That(plusMinusThree.Select(o => o.Notch), Is.EqualTo(plusMinusThree.Select(o => (int)Math.Round(o.Offset * MassDiffAcceptor.NotchScalar))),
                "the notch is the one the Classic and Modern searches give the same offset");

            var oneMm = GlycoSearchEngine.PrecursorMassOffsetsOf(SearchTask.GetMassDiffAcceptor(tolerance, MassDiffAcceptorType.OneMM, null));
            Assert.That(oneMm.Select(o => o.Offset), Is.EqualTo(new[] { 0, 1.0029 }).Within(1e-9));
        }

        [Test]
        public static void PrecursorMassOffsetsOf_RejectsAnOpenAcceptor()
        {
            Assert.Throws<MetaMorpheusException>(() => GlycoSearchEngine.PrecursorMassOffsetsOf(new OpenSearchMode()));
        }

        [Test]
        public static void MassDiffAcceptorType_DefaultsToExact_AcceptsOnlyIsotopeOffsets_AndRoundTripsThroughToml()
        {
            var parameters = new GlycoSearchParameters();
            Assert.That(parameters.MassDiffAcceptorType, Is.EqualTo(MassDiffAcceptorType.Exact));

            foreach (MassDiffAcceptorType type in Enum.GetValues<MassDiffAcceptorType>())
            {
                if (GlycoSearchParameters.SupportedMassDiffAcceptorTypes.Contains(type))
                {
                    Assert.DoesNotThrow(() => parameters.MassDiffAcceptorType = type, type.ToString());
                }
                else
                {
                    Assert.Throws<MetaMorpheusException>(() => parameters.MassDiffAcceptorType = type, type.ToString());
                }
            }

            var task = new GlycoSearchTask();
            task._glycoSearchParameters.MassDiffAcceptorType = MassDiffAcceptorType.PlusOrMinusThreeMM;
            string toml = Path.Combine(TestContext.CurrentContext.TestDirectory, "GlycoMassOffsetRoundTrip.toml");
            try
            {
                Toml.WriteFile(task, toml, MetaMorpheusTask.tomlConfig);
                var read = Toml.ReadFile<GlycoSearchTask>(toml, MetaMorpheusTask.tomlConfig);
                Assert.That(read._glycoSearchParameters.MassDiffAcceptorType, Is.EqualTo(MassDiffAcceptorType.PlusOrMinusThreeMM));
            }
            finally
            {
                File.Delete(toml);
            }

            // a TOML written before the setting existed omits it and gets Exact
            var old = Toml.ReadFile<GlycoSearchTask>(Path.Combine(TestContext.CurrentContext.TestDirectory, "GlycoTestData", "NGlycanSearchTaskconfig.toml"), MetaMorpheusTask.tomlConfig);
            Assert.That(old._glycoSearchParameters.MassDiffAcceptorType, Is.EqualTo(MassDiffAcceptorType.Exact));
        }

        /// <summary>
        /// A real N-glycopeptide scan with its precursor moved one 13C peak heavy, as when the second isotope peak is picked as the
        /// monoisotopic one. Exact loses the glycopeptide; PlusOrMinusThreeMM finds it again at notch +1 with exactly the score it
        /// has at the true precursor mass, which also shows the glycan Y ions are built from the corrected precursor mass.
        /// </summary>
        [Test]
        public static void NGlycanSearch_PrecursorPickedOneIsotopeHeavy_IsFoundAtNotchOneWithTheSameScore()
        {
            var task = Toml.ReadFile<GlycoSearchTask>(Path.Combine(TestContext.CurrentContext.TestDirectory, "GlycoTestData", "NGlycanSearchTaskconfig.toml"), MetaMorpheusTask.tomlConfig);
            CommonParameters parameters = task.CommonParameters;
            var glyco = task._glycoSearchParameters;
            string spectra = Path.Combine(TestContext.CurrentContext.TestDirectory, "GlycoTestData", "yeast_glycan_25170.mgf");
            Ms2ScanWithSpecificMass[] scans = LoadScans(spectra, parameters);
            Ms2ScanWithSpecificMass[] heavyScans = MovePrecursors(scans, C13, parameters);

            var fixedMods = GlobalVariables.AllModsKnown.Where(m => parameters.ListOfModsFixed.Contains((m.ModificationType, m.IdWithMotif))).ToList();
            var variableMods = GlobalVariables.AllModsKnown.Where(m => parameters.ListOfModsVariable.Contains((m.ModificationType, m.IdWithMotif))).ToList();
            List<Protein> proteins = ProteinDbLoader.LoadProteinFasta(Path.Combine(TestContext.CurrentContext.TestDirectory, "GlycoTestData", "Q9C0Y4.fasta"),
                true, DecoyType.Reverse, false, out _);
            var index = BuildIndex(proteins, parameters, fixedMods, variableMods, DecoyType.Reverse);

            // the task file's own search settings (top N 10, oxonium filter on)
            GlycoSpectralMatch Best(Ms2ScanWithSpecificMass[] searched, MassDiffAcceptorType type) =>
                Search(searched, index, parameters, GlycoSearchType.NGlycanSearch, glyco.OGlycanDatabasefile, glyco.NGlycanDatabasefile, type,
                        glyco.GlycoSearchTopNum, glyco.MaximumOGlycanAllowed, glyco.OxoniumIonFilt)
                    .Where(g => g.NGlycan != null).OrderByDescending(g => g.Score).FirstOrDefault();

            GlycoSpectralMatch atTrueMass = Best(scans, MassDiffAcceptorType.Exact);
            Assert.That(atTrueMass, Is.Not.Null, "premise: the scan is identified as an N-glycopeptide at its true precursor mass");
            Assert.That(atTrueMass.Notch, Is.EqualTo(0));

            GlycoSpectralMatch heavyExact = Best(heavyScans, MassDiffAcceptorType.Exact);
            Assert.That(heavyExact?.FullSequence, Is.Not.EqualTo(atTrueMass.FullSequence), "premise: without an offset the heavy precursor loses the glycopeptide");

            GlycoSpectralMatch heavyWithOffsets = Best(heavyScans, MassDiffAcceptorType.PlusOrMinusThreeMM);
            Assert.That(heavyWithOffsets, Is.Not.Null);
            Assert.That(heavyWithOffsets.FullSequence, Is.EqualTo(atTrueMass.FullSequence));
            Assert.That(heavyWithOffsets.Notch, Is.EqualTo((int)Math.Round(C13 * MassDiffAcceptor.NotchScalar)), "found at the +1 isotope notch");
            Assert.That(heavyWithOffsets.Score, Is.EqualTo(atTrueMass.Score).Within(1e-9),
                "the corrected precursor gives the glycan Y ions their true masses, so the score is the one at the true precursor mass");
            Assert.That(NotchCell(heavyWithOffsets), Is.EqualTo("1.00335"), "the written Notch column carries the offset in Da");

            // offsets do not change a scan whose precursor was picked right
            GlycoSpectralMatch trueWithOffsets = Best(scans, MassDiffAcceptorType.PlusOrMinusThreeMM);
            Assert.That(trueWithOffsets.FullSequence, Is.EqualTo(atTrueMass.FullSequence));
            Assert.That(trueWithOffsets.Notch, Is.EqualTo(0));
            Assert.That(trueWithOffsets.Score, Is.EqualTo(atTrueMass.Score).Within(1e-9));
            Assert.That(NotchCell(trueWithOffsets), Is.EqualTo("0"));
        }

        [Test]
        public static void GlycoSearchTask_BuildsTheChosenAcceptorAroundTheFilesPrecursorTolerance()
        {
            var task = new GlycoSearchTask();
            task._glycoSearchParameters.MassDiffAcceptorType = MassDiffAcceptorType.TwoMM;
            var parameters = new CommonParameters(precursorMassTolerance: new MzLibUtil.PpmTolerance(7));

            MassDiffAcceptor acceptor = task.GetPrecursorMassDiffAcceptor(parameters);

            Assert.That(acceptor, Is.InstanceOf<DotMassDiffAcceptor>());
            Assert.That(acceptor.NumNotches, Is.EqualTo(3));
            Assert.That(acceptor.ToProseString(), Does.Contain("7"), "the file's precursor tolerance is the one used");
        }

        /// <summary>
        /// The O-glycan path (MatchGlycopeptide -> FindOGlycan -> CreateGsm) on a real O-glycopeptide scan with its precursor moved one
        /// 13C peak heavy: Exact loses it; PlusOrMinusThreeMM finds it at notch +1 with the score it has at the true precursor mass.
        /// </summary>
        [Test]
        public static void OGlycanSearch_PrecursorPickedOneIsotopeHeavy_IsFoundAtNotchOneWithTheSameScore()
        {
            var parameters = new CommonParameters(dissociationType: MassSpectrometry.DissociationType.EThcD, trimMsMsPeaks: false);
            string spectra = Path.Combine(TestContext.CurrentContext.TestDirectory, "GlycoTestData", "2019_09_16_StcEmix_35trig_EThcD25_rep1_4565.mgf");
            Ms2ScanWithSpecificMass[] scans = LoadScans(spectra, parameters);
            Ms2ScanWithSpecificMass[] heavyScans = MovePrecursors(scans, C13, parameters);

            // the peptide this scan is known to be (see GlycoSearchEngineTest), on its own protein
            var index = BuildIndex(new List<Protein> { new Protein("TTGSLEPSSGASGPQVSSVK", "P_O") }, parameters, new List<Modification>(), new List<Modification>(), DecoyType.None);

            GlycoSpectralMatch Best(Ms2ScanWithSpecificMass[] searched, MassDiffAcceptorType type) =>
                Search(searched, index, parameters, GlycoSearchType.OGlycanSearch, "OGlycan.gdb", "NGlycan.gdb", type)
                    .Where(g => g.BaseSequence == "TTGSLEPSSGASGPQVSSVK" && g.FullSequence.Contains("O-linked glycosylation"))
                    .OrderByDescending(g => g.Score).FirstOrDefault();

            GlycoSpectralMatch atTrueMass = Best(scans, MassDiffAcceptorType.Exact);
            Assert.That(atTrueMass, Is.Not.Null, "premise: the scan is identified as the O-glycopeptide at its true precursor mass");
            Assert.That(atTrueMass.Notch, Is.EqualTo(0));

            Assert.That(Best(heavyScans, MassDiffAcceptorType.Exact), Is.Null, "premise: without an offset the heavy precursor loses the O-glycopeptide");

            GlycoSpectralMatch heavyWithOffsets = Best(heavyScans, MassDiffAcceptorType.PlusOrMinusThreeMM);
            Assert.That(heavyWithOffsets, Is.Not.Null);
            Assert.That(heavyWithOffsets.FullSequence, Is.EqualTo(atTrueMass.FullSequence));
            Assert.That(heavyWithOffsets.Notch, Is.EqualTo((int)Math.Round(C13 * MassDiffAcceptor.NotchScalar)), "found at the +1 isotope notch");
            Assert.That(heavyWithOffsets.Score, Is.EqualTo(atTrueMass.Score).Within(1e-9));
            Assert.That(NotchCell(heavyWithOffsets), Is.EqualTo("1.00335"), "the written Notch column carries the offset in Da");
        }

        /// <summary>
        /// A non-glycopeptide (peptide-only match, FindSingle) whose precursor was picked one 13C peak light. The peptide is then heavier
        /// than the observed precursor, so the first pass must look above it by more than 1 Da (FirstPassMassMargin), the corrected mass is
        /// observed + 1.00335 (notch -1), and the match keeps that notch.
        /// </summary>
        [Test]
        public static void PeptideOnlyMatch_PrecursorPickedOneIsotopeLight_IsFoundAtNotchMinusOneWithTheSameScore()
        {
            var task = Toml.ReadFile<GlycoSearchTask>(Path.Combine(TestContext.CurrentContext.TestDirectory, "GlycoTestData", "GlycoSnip.toml"), MetaMorpheusTask.tomlConfig);
            CommonParameters parameters = task.CommonParameters.CloneWithNewTotalPartitions(1);
            string spectra = Path.Combine(TestContext.CurrentContext.TestDirectory, "GlycoTestData", "GlycoPepMix_snip.mzML");
            Ms2ScanWithSpecificMass[] scans = LoadScans(spectra, parameters);
            Ms2ScanWithSpecificMass[] lightScans = MovePrecursors(scans, -C13, parameters);

            var proteins = new List<Protein>();
            foreach (string fasta in new[] { Path.Combine("GlycoTestData", "GlycoProteinFASTA_7proteins.fasta"), Path.Combine("TestData", "hela_snip_for_unitTest.fasta") })
            {
                proteins.AddRange(ProteinDbLoader.LoadProteinFasta(Path.Combine(TestContext.CurrentContext.TestDirectory, fasta), true, DecoyType.Reverse, false, out _));
            }
            var fixedMods = GlobalVariables.AllModsKnown.Where(m => parameters.ListOfModsFixed.Contains((m.ModificationType, m.IdWithMotif))).ToList();
            var variableMods = GlobalVariables.AllModsKnown.Where(m => parameters.ListOfModsVariable.Contains((m.ModificationType, m.IdWithMotif))).ToList();
            var index = BuildIndex(proteins, parameters, fixedMods, variableMods, DecoyType.Reverse);

            List<GlycoSpectralMatch> PeptideOnly(Ms2ScanWithSpecificMass[] searched, MassDiffAcceptorType type) =>
                Search(searched, index, parameters, GlycoSearchType.NGlycanSearch, "OGlycan.gdb", "NGlycan.gdb", type)
                    .Where(g => g.NGlycan == null && g.LocalizationGraphs == null && !g.IsDecoy).ToList();

            // the best peptide-only match at the true precursor masses
            GlycoSpectralMatch atTrueMass = PeptideOnly(scans, MassDiffAcceptorType.Exact).OrderByDescending(g => g.Score).FirstOrDefault();
            Assert.That(atTrueMass, Is.Not.Null, "premise: the data hold a peptide-only match at the true precursor mass");
            Assert.That(atTrueMass.Notch, Is.EqualTo(0));

            bool SameMatch(GlycoSpectralMatch g) => g.ScanNumber == atTrueMass.ScanNumber && g.FullSequence == atTrueMass.FullSequence;
            Assert.That(PeptideOnly(lightScans, MassDiffAcceptorType.Exact).Any(SameMatch), Is.False,
                "premise: without an offset the light precursor loses the peptide");

            GlycoSpectralMatch lightWithOffsets = PeptideOnly(lightScans, MassDiffAcceptorType.PlusOrMinusThreeMM).FirstOrDefault(SameMatch);
            Assert.That(lightWithOffsets, Is.Not.Null, "the peptide, heavier than the light precursor, must reach the first pass and match at notch -1");
            Assert.That(lightWithOffsets.Notch, Is.EqualTo(-(int)Math.Round(C13 * MassDiffAcceptor.NotchScalar)), "found at the -1 isotope notch");
            Assert.That(lightWithOffsets.Score, Is.EqualTo(atTrueMass.Score).Within(1e-9));
            Assert.That(NotchCell(lightWithOffsets), Is.EqualTo("-1.00335"), "the written Notch column carries the offset in Da");
        }

        private static Ms2ScanWithSpecificMass[] LoadScans(string spectra, CommonParameters parameters)
        {
            var file = new MyFileManager(true).LoadFile(spectra, parameters);
            return MetaMorpheusTask.GetMs2Scans(file, spectra, parameters).OrderBy(s => s.PrecursorMass).ToArray();
        }

        /// <summary> The same scans with each precursor moved by <paramref name="shift"/> Da, as when the wrong isotope peak is picked. </summary>
        private static Ms2ScanWithSpecificMass[] MovePrecursors(Ms2ScanWithSpecificMass[] scans, double shift, CommonParameters parameters)
        {
            var moved = scans
                .Select(s => new Ms2ScanWithSpecificMass(s.TheScan, s.PrecursorMonoisotopicPeakMz + shift / s.PrecursorCharge, s.PrecursorCharge, s.FullFilePath, parameters, s.ExperimentalFragments)
                {
                    ChildScans = s.ChildScans
                })
                .OrderBy(s => s.PrecursorMass).ToArray();
            Assert.That(moved.Select(s => s.PrecursorMass).OrderBy(m => m), Is.EqualTo(scans.Select(s => s.PrecursorMass + shift).OrderBy(m => m)).Within(1e-6), "premise: moved by the shift");
            return moved;
        }

        private static IndexingResults BuildIndex(List<Protein> proteins, CommonParameters parameters, List<Modification> fixedMods, List<Modification> variableMods, DecoyType decoyType)
        {
            return (IndexingResults)new IndexingEngine(proteins, variableMods, fixedMods, null, null, null, 0, decoyType, parameters, null, 30000.0, false,
                new List<FileInfo>(), TargetContaminantAmbiguity.RemoveContaminant, new List<string>()).Run();
        }

        private static List<GlycoSpectralMatch> Search(Ms2ScanWithSpecificMass[] scans, IndexingResults index, CommonParameters parameters, GlycoSearchType type,
            string oGlycanDatabase, string nGlycanDatabase, MassDiffAcceptorType acceptorType, int topN = 50, int maxOGlycans = 4, bool oxoniumIonFilter = false)
        {
            var gsms = new List<GlycoSpectralMatch>[scans.Length];
            new GlycoSearchEngine(gsms, scans, index.PeptideIndex, index.FragmentIndex, null, 0, parameters, null, oGlycanDatabase, nGlycanDatabase, type,
                topN, maxOGlycans, oxoniumIonFilter, new List<string>())
            {
                PrecursorMassDiffAcceptor = SearchTask.GetMassDiffAcceptor(parameters.PrecursorMassTolerance, acceptorType, null)
            }.Run();
            return gsms.Where(g => g != null).SelectMany(g => g).ToList();
        }

        /// <summary> The Notch cell of the row GlycoSpectralMatch writes to the psmtsv files. </summary>
        private static string NotchCell(GlycoSpectralMatch gsm)
        {
            int column = Array.IndexOf(GlycoSpectralMatch.GetTabSepHeaderSingle().Split('\t'), "Notch");
            Assert.That(column, Is.GreaterThanOrEqualTo(0), "the header has a Notch column");
            return gsm.SingleToString().Split('\t')[column];
        }
    }
}
