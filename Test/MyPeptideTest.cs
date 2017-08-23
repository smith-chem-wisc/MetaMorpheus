using Chemistry;
using EngineLayer;
using EngineLayer.ClassicSearch;
using IO.MzML;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;

using Proteomics;
using System;
using System.Collections.Generic;
using System.Linq;
using TaskLayer;
using static Chemistry.PeriodicTable;

namespace Test
{
    [TestFixture]
    public class MyPeptideTest
    {
        #region Public Methods

        [Test]
        public static void TestGoodPeptide()
        {
            var prot = new Protein("MNNNKQQQQ", null);
            var protease = new Protease("Custom Protease", new List<string> { "K" }, new List<string>(), TerminusType.C, CleavageSpecificity.Full, null, null, null);

            var ye = prot.Digest(protease, 0, null, null, InitiatorMethionineBehavior.Retain, new List<ModificationWithMass>()).ToList();

            Assert.AreEqual(2, ye.Count);

            List<ModificationWithMass> variableModifications = new List<ModificationWithMass>();
            var pep1 = ye[0].GetPeptidesWithSetModifications(variableModifications, 4096, 3).First();
            Assert.IsTrue(pep1.MonoisotopicMass > 0);
            foreach (var huh in pep1.CompactPeptide(TerminusType.None).ProductMassesMightHaveDuplicatesAndNaNs(new List<ProductType> { ProductType.B, ProductType.Y }))
                Assert.IsTrue(huh > 0);
            var pep2 = ye[1].GetPeptidesWithSetModifications(variableModifications, 4096, 3).First();
            Assert.IsTrue(pep2.MonoisotopicMass > 0);
            foreach (var huh in pep2.CompactPeptide(TerminusType.None).ProductMassesMightHaveDuplicatesAndNaNs(new List<ProductType> { ProductType.B, ProductType.Y }))
                Assert.IsTrue(huh > 0);
        }

        [Test]
        public static void TestIdenticalPeaks()
        {
            IDictionary<int, List<Modification>> mods = new Dictionary<int, List<Modification>>();
            ModificationMotif.TryGetMotif("M", out ModificationMotif motif);
            mods.Add(1, new List<Modification> { new ModificationWithMass("Hehe", null, motif, TerminusLocalization.NProt, 18.010565, null, null, null, null) });
            var prot = new Protein("MMMM", null, null, mods);
            var ye = prot.Digest(GlobalTaskLevelSettings.ProteaseDictionary["trypsin"], 0, null, null, InitiatorMethionineBehavior.Retain, new List<ModificationWithMass>()).First();
            var thePep = ye.GetPeptidesWithSetModifications(new List<ModificationWithMass>(), 2, 1).Last();

            var massArray = thePep.CompactPeptide(TerminusType.None).ProductMassesMightHaveDuplicatesAndNaNs(new List<ProductType> { ProductType.B, ProductType.Y });
            Array.Sort(massArray);
            double[] intensities = new double[] { 1, 1, 1, 1 };
            double[] mz = new double[] { massArray[0].ToMz(1), massArray[2].ToMz(1), massArray[4].ToMz(1), 10000 };
            MzmlMzSpectrum massSpectrum = new MzmlMzSpectrum(mz, intensities, false);
            IMsDataScanWithPrecursor<IMzSpectrum<IMzPeak>> scan = new MzmlScanWithPrecursor(1, massSpectrum, 1, true, Polarity.Positive, 1, new MzRange(300, 2000), "", MZAnalyzerType.Unknown, massSpectrum.SumOfAllY, 0, null, null, 0, null, DissociationType.Unknown, 1, null, null);

            Psm[][] globalPsms = new Psm[1][];
            globalPsms[0] = new Psm[1];
            Ms2ScanWithSpecificMass[] arrayOfSortedMS2Scans = { new Ms2ScanWithSpecificMass(scan, new MzPeak(0, 0), 0, null) };
            CommonParameters CommonParameters = new CommonParameters
            {
                ProductMassTolerance = new PpmTolerance(5),
                Protease = GlobalTaskLevelSettings.ProteaseDictionary["trypsin"],
                MaxMissedCleavages = 0,
                MinPeptideLength = null,
                MaxPeptideLength = null,
                MaxModificationIsoforms = int.MaxValue,
                ConserveMemory = false,
                InitiatorMethionineBehavior = InitiatorMethionineBehavior.Retain,
                ScoreCutoff = 0
            };
            ClassicSearchEngine cse = new ClassicSearchEngine(globalPsms, arrayOfSortedMS2Scans, new List<ModificationWithMass>(), new List<ModificationWithMass>(), new List<Protein> { prot }, new List<ProductType> { ProductType.B, ProductType.Y }, new List<MassDiffAcceptor> { new OpenSearchMode() }, false, CommonParameters, null);

            cse.Run();

            Assert.Less(globalPsms[0][0].Score, 4);
            Assert.Greater(globalPsms[0][0].Score, 3);
        }

        [Test]
        public static void TestLastPeaks()
        {
            IDictionary<int, List<Modification>> mods = new Dictionary<int, List<Modification>>();
            ModificationMotif.TryGetMotif("M", out ModificationMotif motif);
            var prot = new Protein("MMMM", null, null, mods);
            var ye = prot.Digest(GlobalTaskLevelSettings.ProteaseDictionary["trypsin"], 0, null, null, InitiatorMethionineBehavior.Retain, new List<ModificationWithMass>()).First();
            var thePep = ye.GetPeptidesWithSetModifications(new List<ModificationWithMass>(), 2, 1).Last();

            var massArray = thePep.CompactPeptide(TerminusType.None).ProductMassesMightHaveDuplicatesAndNaNs(new List<ProductType> { ProductType.B, ProductType.Y });
            Array.Sort(massArray);
            double[] intensities = new double[] { 1, 1, 1 };
            double[] mz = new double[] { 1, 2, massArray[4].ToMz(1) };
            MzmlMzSpectrum massSpectrum = new MzmlMzSpectrum(mz, intensities, false);
            IMsDataScanWithPrecursor<IMzSpectrum<IMzPeak>> scan = new MzmlScanWithPrecursor(1, massSpectrum, 1, true, Polarity.Positive, 1, new MzRange(300, 2000), "", MZAnalyzerType.Unknown, massSpectrum.SumOfAllY, 0, null, null, 0, null, DissociationType.Unknown, 1, null, null);

            Psm[][] globalPsms = new Psm[1][];
            globalPsms[0] = new Psm[1];
            Ms2ScanWithSpecificMass[] arrayOfSortedMS2Scans = { new Ms2ScanWithSpecificMass(scan, new MzPeak(0, 0), 0, null) };
            CommonParameters CommonParameters = new CommonParameters
            {
                ProductMassTolerance = new PpmTolerance(5),
                Protease = GlobalTaskLevelSettings.ProteaseDictionary["trypsin"],
                MaxMissedCleavages = 0,
                MinPeptideLength = null,
                MaxPeptideLength = null,
                MaxModificationIsoforms = int.MaxValue,
                ConserveMemory = false,
                InitiatorMethionineBehavior = InitiatorMethionineBehavior.Retain,
                ScoreCutoff = 0
            };
            ClassicSearchEngine cse = new ClassicSearchEngine(globalPsms, arrayOfSortedMS2Scans, new List<ModificationWithMass>(), new List<ModificationWithMass>(), new List<Protein> { prot }, new List<ProductType> { ProductType.B, ProductType.Y }, new List<MassDiffAcceptor> { new OpenSearchMode() }, false, CommonParameters, null);

            cse.Run();
            Assert.Less(globalPsms[0][0].Score, 2);
            Assert.Greater(globalPsms[0][0].Score, 1);
        }

        [Test]
        public static void TestVeryCloseExperimentals()
        {
            IDictionary<int, List<Modification>> mods = new Dictionary<int, List<Modification>>();
            ModificationMotif.TryGetMotif("M", out ModificationMotif motif);
            var prot = new Protein("MMMM", null, null, mods);
            var ye = prot.Digest(GlobalTaskLevelSettings.ProteaseDictionary["trypsin"], 0, null, null, InitiatorMethionineBehavior.Retain, new List<ModificationWithMass>()).First();
            var thePep = ye.GetPeptidesWithSetModifications(new List<ModificationWithMass>(), 2, 1).Last();

            var massArray = thePep.CompactPeptide(TerminusType.None).ProductMassesMightHaveDuplicatesAndNaNs(new List<ProductType> { ProductType.B, ProductType.Y });
            Array.Sort(massArray);
            double[] intensities = new double[] { 1, 1, 1, 1 };
            double[] mz = new double[] { 1, 2, massArray[4].ToMz(1), massArray[4].ToMz(1) + 1e-9 };
            MzmlMzSpectrum massSpectrum = new MzmlMzSpectrum(mz, intensities, false);
            IMsDataScanWithPrecursor<IMzSpectrum<IMzPeak>> scan = new MzmlScanWithPrecursor(1, massSpectrum, 1, true, Polarity.Positive, 1, new MzRange(300, 2000), "", MZAnalyzerType.Unknown, massSpectrum.SumOfAllY, 0, null, null, 0, null, DissociationType.Unknown, 1, null, null);

            Psm[][] globalPsms = new Psm[1][];
            globalPsms[0] = new Psm[1];
            Ms2ScanWithSpecificMass[] arrayOfSortedMS2Scans = { new Ms2ScanWithSpecificMass(scan, new MzPeak(0, 0), 0, null) };
            CommonParameters CommonParameters = new CommonParameters
            {
                ProductMassTolerance = new PpmTolerance(5),
                Protease = GlobalTaskLevelSettings.ProteaseDictionary["trypsin"],
                MaxMissedCleavages = 0,
                MinPeptideLength = null,
                MaxPeptideLength = null,
                MaxModificationIsoforms = int.MaxValue,
                ConserveMemory = false,
                InitiatorMethionineBehavior = InitiatorMethionineBehavior.Retain,
                ScoreCutoff = 0
            };
            ClassicSearchEngine cse = new ClassicSearchEngine(globalPsms, arrayOfSortedMS2Scans, new List<ModificationWithMass>(), new List<ModificationWithMass>(), new List<Protein> { prot }, new List<ProductType> { ProductType.B, ProductType.Y }, new List<MassDiffAcceptor> { new OpenSearchMode() }, false, CommonParameters, null);

            cse.Run();
            Assert.Less(globalPsms[0][0].Score, 2);
            Assert.Greater(globalPsms[0][0].Score, 1);
        }

        [Test]
        public static void TestAllNaN()
        {
            IDictionary<int, List<Modification>> mods = new Dictionary<int, List<Modification>>();
            var prot = new Protein("XMMM", null, null, mods);
            var ye = prot.Digest(GlobalTaskLevelSettings.ProteaseDictionary["trypsin"], 0, null, null, InitiatorMethionineBehavior.Retain, new List<ModificationWithMass>()).First();
            var thePep = ye.GetPeptidesWithSetModifications(new List<ModificationWithMass>(), 2, 1).Last();

            var massArray = thePep.CompactPeptide(TerminusType.None).ProductMassesMightHaveDuplicatesAndNaNs(new List<ProductType> { ProductType.B });
            Array.Sort(massArray);
            double[] intensities = new double[] { 1, 1, 1, 1 };
            double[] mz = new double[] { 1, 2, 3, 4 };
            MzmlMzSpectrum massSpectrum = new MzmlMzSpectrum(mz, intensities, false);
            IMsDataScanWithPrecursor<IMzSpectrum<IMzPeak>> scan = new MzmlScanWithPrecursor(1, massSpectrum, 1, true, Polarity.Positive, 1, new MzRange(300, 2000), "", MZAnalyzerType.Unknown, massSpectrum.SumOfAllY, 0, null, null, 0, null, DissociationType.Unknown, 1, null, null);

            Psm[][] globalPsms = new Psm[1][];
            globalPsms[0] = new Psm[1];
            Ms2ScanWithSpecificMass[] arrayOfSortedMS2Scans = { new Ms2ScanWithSpecificMass(scan, new MzPeak(0, 0), 0, null) };
            CommonParameters CommonParameters = new CommonParameters
            {
                ProductMassTolerance = new PpmTolerance(5),
                Protease = GlobalTaskLevelSettings.ProteaseDictionary["trypsin"],
                MaxMissedCleavages = 0,
                MinPeptideLength = null,
                MaxPeptideLength = null,
                MaxModificationIsoforms = int.MaxValue,
                ConserveMemory = false,
                InitiatorMethionineBehavior = InitiatorMethionineBehavior.Retain,
                ScoreCutoff = 0
            };
            ClassicSearchEngine cse = new ClassicSearchEngine(globalPsms, arrayOfSortedMS2Scans, new List<ModificationWithMass>(), new List<ModificationWithMass>(), new List<Protein> { prot }, new List<ProductType> { ProductType.B, ProductType.Y }, new List<MassDiffAcceptor> { new OpenSearchMode() }, false, CommonParameters, null);

            cse.Run();
            Assert.IsNull(globalPsms[0][0]);
        }

        [Test]
        public static void TestNoCleavage()
        {
            List<ModificationWithMass> fixedModifications = new List<ModificationWithMass>();
            var prot = new Protein("MNNNKQQQQ", null, null, new Dictionary<int, List<Modification>>(), new List<ProteolysisProduct> { new ProteolysisProduct(5, 6, "lala") });
            var protease = new Protease("Custom Protease", null, null, TerminusType.None, CleavageSpecificity.None, null, null, null);

            var ye = prot.Digest(protease, int.MaxValue, null, null, InitiatorMethionineBehavior.Variable, fixedModifications).ToList();

            Assert.AreEqual(3, ye.Count);
        }

        [Test]
        public static void TestBadPeptide()
        {
            var prot = new Protein("MNNNKQQXQ", null);
            var protease = new Protease("Custom Protease", new List<string> { "K" }, new List<string>(), TerminusType.C, CleavageSpecificity.Full, null, null, null);

            var ye = prot.Digest(protease, 0, null, null, InitiatorMethionineBehavior.Retain, new List<ModificationWithMass>()).ToList();

            Assert.AreEqual(2, ye.Count);
            var pep1 = ye[0].GetPeptidesWithSetModifications(new List<ModificationWithMass>(), 4096, 3).First();
            Assert.IsTrue(pep1.MonoisotopicMass > 0);
            foreach (var huh in pep1.CompactPeptide(TerminusType.None).ProductMassesMightHaveDuplicatesAndNaNs(new List<ProductType> { ProductType.B, ProductType.Y }))
                Assert.IsTrue(huh > 0);

            var pep2 = ye[1].GetPeptidesWithSetModifications(new List<ModificationWithMass>(), 4096, 3).First();
            Assert.IsNaN(pep2.MonoisotopicMass);
            var cool = pep2.CompactPeptide(TerminusType.None).ProductMassesMightHaveDuplicatesAndNaNs(new List<ProductType> { ProductType.Y });
            Assert.IsTrue(cool[0] > 0);
            Assert.IsTrue(double.IsNaN(cool[1]));
            Assert.IsTrue(double.IsNaN(cool[2]));
            Assert.IsTrue(cool.Length == 3);
        }

        [Test]
        public static void TestPeptideWithSetModifications()
        {
            var prot = new Protein("M", null);
            var protease = new Protease("Custom Protease", new List<string> { "K" }, new List<string>(), TerminusType.C, CleavageSpecificity.Full, null, null, null);
            var ye = prot.Digest(protease, 0, null, null, InitiatorMethionineBehavior.Retain, new List<ModificationWithMass>()).First();
            List<ModificationWithMass> variableModifications = new List<ModificationWithMass>();
            ModificationMotif.TryGetMotif("M", out ModificationMotif motif);
            variableModifications.Add(new ModificationWithMassAndCf("ProtNmod", null, motif, TerminusLocalization.NProt, Chemistry.ChemicalFormula.ParseFormula("H"), GetElement(1).PrincipalIsotope.AtomicMass));
            variableModifications.Add(new ModificationWithMassAndCf("pepNmod", null, motif, TerminusLocalization.NPep, Chemistry.ChemicalFormula.ParseFormula("H"), GetElement(1).PrincipalIsotope.AtomicMass));
            variableModifications.Add(new ModificationWithMassAndCf("resMod", null, motif, TerminusLocalization.Any, Chemistry.ChemicalFormula.ParseFormula("H"), GetElement(1).PrincipalIsotope.AtomicMass));
            variableModifications.Add(new ModificationWithMassAndCf("PepCmod", null, motif, TerminusLocalization.PepC, Chemistry.ChemicalFormula.ParseFormula("H"), GetElement(1).PrincipalIsotope.AtomicMass));
            variableModifications.Add(new ModificationWithMassAndCf("ProtCmod", null, motif, TerminusLocalization.ProtC, Chemistry.ChemicalFormula.ParseFormula("H"), GetElement(1).PrincipalIsotope.AtomicMass));
            var ok = ye.GetPeptidesWithSetModifications(variableModifications, 4096, 5).ToList();
            Assert.AreEqual(8, ok.Count);

            Assert.AreEqual("[:ProtNmod]M[:resMod][:PepCmod]", ok.Last().Sequence);
            Assert.AreEqual("[H]M[H][H]", ok.Last().SequenceWithChemicalFormulas);
            Assert.AreEqual(5 * GetElement("H").PrincipalIsotope.AtomicMass + Residue.ResidueMonoisotopicMass['M'] + GetElement("O").PrincipalIsotope.AtomicMass, ok.Last().MonoisotopicMass, 1e-9);
        }

        [Test]
        public static void TestPeptideWithFixedModifications()
        {
            var prot = new Protein("M", null);
            var protease = new Protease("Custom Protease", new List<string> { "K" }, new List<string>(), TerminusType.C, CleavageSpecificity.Full, null, null, null);
            List<ModificationWithMass> fixedMods = new List<ModificationWithMass>();
            ModificationMotif.TryGetMotif("M", out ModificationMotif motif);
            fixedMods.Add(new ModificationWithMassAndCf("ProtNmod", null, motif, TerminusLocalization.NProt, Chemistry.ChemicalFormula.ParseFormula("H"), GetElement(1).PrincipalIsotope.AtomicMass));
            fixedMods.Add(new ModificationWithMassAndCf("PepNmod", null, motif, TerminusLocalization.NPep, Chemistry.ChemicalFormula.ParseFormula("H"), GetElement(1).PrincipalIsotope.AtomicMass));
            fixedMods.Add(new ModificationWithMassAndCf("resMod", null, motif, TerminusLocalization.Any, Chemistry.ChemicalFormula.ParseFormula("H"), GetElement(1).PrincipalIsotope.AtomicMass));
            fixedMods.Add(new ModificationWithMassAndCf("PepCmod", null, motif, TerminusLocalization.PepC, Chemistry.ChemicalFormula.ParseFormula("H"), GetElement(1).PrincipalIsotope.AtomicMass));
            fixedMods.Add(new ModificationWithMassAndCf("ProtCmod", null, motif, TerminusLocalization.ProtC, Chemistry.ChemicalFormula.ParseFormula("H"), GetElement(1).PrincipalIsotope.AtomicMass));

            var ye = prot.Digest(protease, 0, null, null, InitiatorMethionineBehavior.Retain, fixedMods).First();
            var ok = ye.GetPeptidesWithSetModifications(new List<ModificationWithMass>(), 4096, 5).ToList();

            Assert.AreEqual(1, ok.Count);

            Assert.AreEqual("[:PepNmod]M[:resMod][:ProtCmod]", ok.Last().Sequence);
            Assert.AreEqual("[H]M[H][H]", ok.Last().SequenceWithChemicalFormulas);
            Assert.AreEqual(5 * GetElement("H").PrincipalIsotope.AtomicMass + Residue.ResidueMonoisotopicMass['M'] + GetElement("O").PrincipalIsotope.AtomicMass, ok.Last().MonoisotopicMass, 1e-9);
        }

        [Test]
        public static void TestDigestIndices()
        {
            ModificationMotif.TryGetMotif("X", out ModificationMotif motif);
            Modification mod = new ModificationWithMass(null, null, motif, TerminusLocalization.Any, double.NaN);
            IDictionary<int, List<Modification>> modDict = new Dictionary<int, List<Modification>>
            {
                {2, new List<Modification> {mod } },
                {8, new List<Modification> {mod } }
            };
            var prot = new Protein("MNNNNKRRRRR", null, null, modDict);
            var protease = new Protease("Custom Protease", new List<string> { "K" }, new List<string>(), TerminusType.C, CleavageSpecificity.Full, null, null, null);
            var ye1 = prot.Digest(protease, 0, null, null, InitiatorMethionineBehavior.Retain, new List<ModificationWithMass>()).First();
            var ye2 = prot.Digest(protease, 0, null, null, InitiatorMethionineBehavior.Retain, new List<ModificationWithMass>()).Last();
            var ok1 = ye1.GetPeptidesWithSetModifications(new List<ModificationWithMass>(), 2, 1).Last();
            var ok2 = ye2.GetPeptidesWithSetModifications(new List<ModificationWithMass>(), 2, 1).Last();

            Assert.AreEqual(1, ok1.NumMods);
            Assert.IsTrue(ok1.allModsOneIsNterminus.ContainsKey(3));
            Assert.AreEqual(1, ok2.NumMods);
            Assert.IsTrue(ok2.allModsOneIsNterminus.ContainsKey(3));
        }

        [Test]
        public static void TestDigestDecoy()
        {
            ModificationMotif.TryGetMotif("Abcdefg", out ModificationMotif motif);
            Modification mod = new ModificationWithMass(null, null, motif, TerminusLocalization.Any, double.NaN);
            IDictionary<int, List<Modification>> modDict = new Dictionary<int, List<Modification>>
            {
                {2, new List<Modification> {mod } },
                {8, new List<Modification> {mod } }
            };
            var prot = new Protein("MNNNNKRRRRR", null, null, modDict, isDecoy: true);
            var protease = new Protease("Custom Protease", new List<string> { "K" }, new List<string>(), TerminusType.C, CleavageSpecificity.Full, null, null, null);
            var ye1 = prot.Digest(protease, 0, null, null, InitiatorMethionineBehavior.Retain, new List<ModificationWithMass>()).First();
            var ye2 = prot.Digest(protease, 0, null, null, InitiatorMethionineBehavior.Retain, new List<ModificationWithMass>()).Last();
            var ok1 = ye1.GetPeptidesWithSetModifications(new List<ModificationWithMass>(), 2, 1).Last();
            var ok2 = ye2.GetPeptidesWithSetModifications(new List<ModificationWithMass>(), 2, 1).Last();

            Assert.AreEqual(1, ok1.NumMods);
            Assert.IsTrue(ok1.allModsOneIsNterminus.ContainsKey(3));
            Assert.AreEqual(1, ok2.NumMods);
            Assert.IsTrue(ok2.allModsOneIsNterminus.ContainsKey(3));

            prot = new Protein("MNNNNKRRRRR", null, null, modDict);
            ye1 = prot.Digest(protease, 0, null, null, InitiatorMethionineBehavior.Retain, new List<ModificationWithMass>()).First();
            ye2 = prot.Digest(protease, 0, null, null, InitiatorMethionineBehavior.Retain, new List<ModificationWithMass>()).Last();
            ok1 = ye1.GetPeptidesWithSetModifications(new List<ModificationWithMass>(), 2, 1).Last();
            ok2 = ye2.GetPeptidesWithSetModifications(new List<ModificationWithMass>(), 2, 1).Last();

            Assert.AreEqual(0, ok1.NumMods);
            Assert.IsFalse(ok1.allModsOneIsNterminus.Any());
            Assert.AreEqual(0, ok2.NumMods);
            Assert.IsFalse(ok2.allModsOneIsNterminus.Any());
        }

        [Test]
        public static void TestGoodPeptideWithLength()
        {
            var prot = new Protein("MNNNKQQQQMNNNKQQQQ", null);
            var protease = new Protease("Custom Protease", new List<string> { "K" }, new List<string>(), TerminusType.C, CleavageSpecificity.Full, null, null, null);

            var ye = prot.Digest(protease, 0, null, null, InitiatorMethionineBehavior.Retain, new List<ModificationWithMass>()).ToList();
            var ye1 = prot.Digest(protease, 0, 5, null, InitiatorMethionineBehavior.Retain, new List<ModificationWithMass>()).ToList();
            var ye2 = prot.Digest(protease, 0, null, 5, InitiatorMethionineBehavior.Retain, new List<ModificationWithMass>()).ToList();
            var ye3 = prot.Digest(protease, 0, 5, 8, InitiatorMethionineBehavior.Retain, new List<ModificationWithMass>()).ToList();
            Assert.AreEqual(3, ye.Count);
            Assert.AreEqual(2, ye1.Count);
            Assert.AreEqual(2, ye2.Count);
            Assert.AreEqual(1, ye3.Count);
        }

        #endregion Public Methods
    }
}