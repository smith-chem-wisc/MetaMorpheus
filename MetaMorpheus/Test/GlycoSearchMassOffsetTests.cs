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
            string spectra = Path.Combine(TestContext.CurrentContext.TestDirectory, "GlycoTestData", "yeast_glycan_25170.mgf");
            var file = new MyFileManager(true).LoadFile(spectra, parameters);
            Ms2ScanWithSpecificMass[] scans = MetaMorpheusTask.GetMs2Scans(file, spectra, parameters).OrderBy(s => s.PrecursorMass).ToArray();

            // the same scans, with each precursor one 13C peak heavier
            Ms2ScanWithSpecificMass[] heavyScans = scans
                .Select(s => new Ms2ScanWithSpecificMass(s.TheScan, s.PrecursorMonoisotopicPeakMz + C13 / s.PrecursorCharge, s.PrecursorCharge, s.FullFilePath, parameters, s.ExperimentalFragments))
                .ToArray();
            Assert.That(heavyScans.Select(s => s.PrecursorMass), Is.EqualTo(scans.Select(s => s.PrecursorMass + C13)).Within(1e-6), "premise: one 13C heavier");

            var fixedMods = GlobalVariables.AllModsKnown.Where(m => parameters.ListOfModsFixed.Contains((m.ModificationType, m.IdWithMotif))).ToList();
            var variableMods = GlobalVariables.AllModsKnown.Where(m => parameters.ListOfModsVariable.Contains((m.ModificationType, m.IdWithMotif))).ToList();
            List<Protein> proteins = ProteinDbLoader.LoadProteinFasta(Path.Combine(TestContext.CurrentContext.TestDirectory, "GlycoTestData", "Q9C0Y4.fasta"),
                true, DecoyType.Reverse, false, out _);
            var index = (IndexingResults)new IndexingEngine(proteins, variableMods, fixedMods, null, null, null, 0, DecoyType.Reverse, parameters, null, 30000.0, false,
                new List<FileInfo>(), TargetContaminantAmbiguity.RemoveContaminant, new List<string>()).Run();

            GlycoSpectralMatch Best(Ms2ScanWithSpecificMass[] searched, MassDiffAcceptorType type)
            {
                var gsms = new List<GlycoSpectralMatch>[searched.Length];
                new GlycoSearchEngine(gsms, searched, index.PeptideIndex, index.FragmentIndex, null, 0, parameters, null,
                    task._glycoSearchParameters.OGlycanDatabasefile, task._glycoSearchParameters.NGlycanDatabasefile, GlycoSearchType.NGlycanSearch,
                    task._glycoSearchParameters.GlycoSearchTopNum, task._glycoSearchParameters.MaximumOGlycanAllowed, task._glycoSearchParameters.OxoniumIonFilt, new List<string>())
                {
                    PrecursorMassDiffAcceptor = SearchTask.GetMassDiffAcceptor(parameters.PrecursorMassTolerance, type, null)
                }.Run();
                return gsms.Where(g => g != null).SelectMany(g => g).Where(g => g.NGlycan != null).OrderByDescending(g => g.Score).FirstOrDefault();
            }

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

            // offsets do not change a scan whose precursor was picked right
            GlycoSpectralMatch trueWithOffsets = Best(scans, MassDiffAcceptorType.PlusOrMinusThreeMM);
            Assert.That(trueWithOffsets.FullSequence, Is.EqualTo(atTrueMass.FullSequence));
            Assert.That(trueWithOffsets.Notch, Is.EqualTo(0));
            Assert.That(trueWithOffsets.Score, Is.EqualTo(atTrueMass.Score).Within(1e-9));
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
    }
}
