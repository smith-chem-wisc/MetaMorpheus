using Chemistry;
using EngineLayer;
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
            foreach (var huh in pep1.CompactPeptide.ProductMassesMightHaveDuplicatesAndNaNs(new List<ProductType> { ProductType.B, ProductType.Y }))
                Assert.IsTrue(huh > 0);

            var pep2 = ye[1].GetPeptidesWithSetModifications(variableModifications, 4096, 3).First();
            Assert.IsTrue(pep2.MonoisotopicMass > 0);
            foreach (var huh in pep2.CompactPeptide.ProductMassesMightHaveDuplicatesAndNaNs(new List<ProductType> { ProductType.B, ProductType.Y }))
                Assert.IsTrue(huh > 0);
        }

        [Test]
        public static void TestIdenticalPeaks()
        {
            IDictionary<int, List<Modification>> mods = new Dictionary<int, List<Modification>>();
            ModificationMotif motif;
            ModificationMotif.TryGetMotif("M", out motif);
            mods.Add(1, new List<Modification> { new ModificationWithMass("Hehe", null, motif, ModificationSites.NProt, 18.010565, null, null, null, null) });
            var prot = new Protein("MMMM", null, null, mods);
            var ye = prot.Digest(GlobalTaskLevelSettings.ProteaseDictionary["trypsin"], 0, null, null, InitiatorMethionineBehavior.Retain, new List<ModificationWithMass>()).First();
            var thePep = ye.GetPeptidesWithSetModifications(new List<ModificationWithMass>(), 2, 1).Last();

            var massArray = thePep.CompactPeptide.ProductMassesMightHaveDuplicatesAndNaNs(new List<ProductType> { ProductType.B, ProductType.Y });
            double[] matchedIonMassesListPositiveIsMatch = new double[massArray.Length];
            Array.Sort(massArray);
            double[] intensities = new double[] { 1, 1, 1, 1 };
            double[] mz = new double[] { massArray[0].ToMz(1), massArray[2].ToMz(1), massArray[4].ToMz(1), 10000 };
            MzmlMzSpectrum massSpectrum = new MzmlMzSpectrum(mz, intensities, false);
            IMsDataScan<IMzSpectrum<IMzPeak>> scan = new MzmlScan(1, massSpectrum, 1, true, Polarity.Positive, 1, new MzRange(300, 2000), "", MZAnalyzerType.Unknown, massSpectrum.SumOfAllY, null);
            var score = PsmParent.MatchIons(scan, new PpmTolerance(5), massArray, matchedIonMassesListPositiveIsMatch);

            Assert.Less(score, 4);
            Assert.Greater(score, 3);
        }

        [Test]
        public static void TestLastPeaks()
        {
            IDictionary<int, List<Modification>> mods = new Dictionary<int, List<Modification>>();
            ModificationMotif motif;
            ModificationMotif.TryGetMotif("M", out motif);
            var prot = new Protein("MMMM", null, null, mods);
            var ye = prot.Digest(GlobalTaskLevelSettings.ProteaseDictionary["trypsin"], 0, null, null, InitiatorMethionineBehavior.Retain, new List<ModificationWithMass>()).First();
            var thePep = ye.GetPeptidesWithSetModifications(new List<ModificationWithMass>(), 2, 1).Last();

            var massArray = thePep.CompactPeptide.ProductMassesMightHaveDuplicatesAndNaNs(new List<ProductType> { ProductType.B, ProductType.Y });
            double[] matchedIonMassesListPositiveIsMatch = new double[massArray.Length];
            Array.Sort(massArray);
            double[] intensities = new double[] { 1, 1, 1 };
            double[] mz = new double[] { 1, 2, massArray[4].ToMz(1) };
            MzmlMzSpectrum massSpectrum = new MzmlMzSpectrum(mz, intensities, false);
            IMsDataScan<IMzSpectrum<IMzPeak>> scan = new MzmlScan(1, massSpectrum, 1, true, Polarity.Positive, 1, new MzRange(300, 2000), "", MZAnalyzerType.Unknown, massSpectrum.SumOfAllY, null);
            var score = PsmParent.MatchIons(scan, new PpmTolerance(5), massArray, matchedIonMassesListPositiveIsMatch);

            Assert.Less(score, 2);
            Assert.Greater(score, 1);
        }

        [Test]
        public static void TestVeryCloseExperimentals()
        {
            IDictionary<int, List<Modification>> mods = new Dictionary<int, List<Modification>>();
            ModificationMotif motif;
            ModificationMotif.TryGetMotif("M", out motif);
            var prot = new Protein("MMMM", null, null, mods);
            var ye = prot.Digest(GlobalTaskLevelSettings.ProteaseDictionary["trypsin"], 0, null, null, InitiatorMethionineBehavior.Retain, new List<ModificationWithMass>()).First();
            var thePep = ye.GetPeptidesWithSetModifications(new List<ModificationWithMass>(), 2, 1).Last();

            var massArray = thePep.CompactPeptide.ProductMassesMightHaveDuplicatesAndNaNs(new List<ProductType> { ProductType.B, ProductType.Y });
            double[] matchedIonMassesListPositiveIsMatch = new double[massArray.Length];
            Array.Sort(massArray);
            double[] intensities = new double[] { 1, 1, 1, 1 };
            double[] mz = new double[] { 1, 2, massArray[4].ToMz(1), massArray[4].ToMz(1) + 1e-9 };
            MzmlMzSpectrum massSpectrum = new MzmlMzSpectrum(mz, intensities, false);
            IMsDataScan<IMzSpectrum<IMzPeak>> scan = new MzmlScan(1, massSpectrum, 1, true, Polarity.Positive, 1, new MzRange(300, 2000), "", MZAnalyzerType.Unknown, massSpectrum.SumOfAllY, null);
            var score = PsmParent.MatchIons(scan, new PpmTolerance(5), massArray, matchedIonMassesListPositiveIsMatch);

            Assert.Less(score, 2);
            Assert.Greater(score, 1);
        }

        [Test]
        public static void TestAllNaN()
        {
            IDictionary<int, List<Modification>> mods = new Dictionary<int, List<Modification>>();
            var prot = new Protein("XMMM", null, null, mods);
            var ye = prot.Digest(GlobalTaskLevelSettings.ProteaseDictionary["trypsin"], 0, null, null, InitiatorMethionineBehavior.Retain, new List<ModificationWithMass>()).First();
            var thePep = ye.GetPeptidesWithSetModifications(new List<ModificationWithMass>(), 2, 1).Last();

            var massArray = thePep.CompactPeptide.ProductMassesMightHaveDuplicatesAndNaNs(new List<ProductType> { ProductType.B });
            double[] matchedIonMassesListPositiveIsMatch = new double[massArray.Length];
            Array.Sort(massArray);
            double[] intensities = new double[] { 1, 1, 1, 1 };
            double[] mz = new double[] { 1, 2, 3, 4 };
            MzmlMzSpectrum massSpectrum = new MzmlMzSpectrum(mz, intensities, false);
            IMsDataScan<IMzSpectrum<IMzPeak>> scan = new MzmlScan(1, massSpectrum, 1, true, Polarity.Positive, 1, new MzRange(300, 2000), "", MZAnalyzerType.Unknown, massSpectrum.SumOfAllY, null);
            var score = PsmParent.MatchIons(scan, new PpmTolerance(5), massArray, matchedIonMassesListPositiveIsMatch);

            Assert.AreEqual(0, score);
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
            foreach (var huh in pep1.CompactPeptide.ProductMassesMightHaveDuplicatesAndNaNs(new List<ProductType> { ProductType.B, ProductType.Y }))
                Assert.IsTrue(huh > 0);

            var pep2 = ye[1].GetPeptidesWithSetModifications(new List<ModificationWithMass>(), 4096, 3).First();
            Assert.IsNaN(pep2.MonoisotopicMass);
            var cool = pep2.CompactPeptide.ProductMassesMightHaveDuplicatesAndNaNs(new List<ProductType> { ProductType.Y });
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
            ModificationMotif motif;
            ModificationMotif.TryGetMotif("M", out motif);
            variableModifications.Add(new ModificationWithMassAndCf("ProtNmod", null, motif, ModificationSites.NProt, Chemistry.ChemicalFormula.ParseFormula("H"), GetElement(1).PrincipalIsotope.AtomicMass, null, new List<double> { 0 }, null, null));
            variableModifications.Add(new ModificationWithMassAndCf("pepNmod", null, motif, ModificationSites.NPep, Chemistry.ChemicalFormula.ParseFormula("H"), GetElement(1).PrincipalIsotope.AtomicMass, null, new List<double> { 0 }, null, null));
            variableModifications.Add(new ModificationWithMassAndCf("resMod", null, motif, ModificationSites.Any, Chemistry.ChemicalFormula.ParseFormula("H"), GetElement(1).PrincipalIsotope.AtomicMass, null, new List<double> { 0 }, null, null));
            variableModifications.Add(new ModificationWithMassAndCf("PepCmod", null, motif, ModificationSites.PepC, Chemistry.ChemicalFormula.ParseFormula("H"), GetElement(1).PrincipalIsotope.AtomicMass, null, new List<double> { 0 }, null, null));
            variableModifications.Add(new ModificationWithMassAndCf("ProtCmod", null, motif, ModificationSites.ProtC, Chemistry.ChemicalFormula.ParseFormula("H"), GetElement(1).PrincipalIsotope.AtomicMass, null, new List<double> { 0 }, null, null));
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
            ModificationMotif motif;
            ModificationMotif.TryGetMotif("M", out motif);
            fixedMods.Add(new ModificationWithMassAndCf("ProtNmod", null, motif, ModificationSites.NProt, Chemistry.ChemicalFormula.ParseFormula("H"), GetElement(1).PrincipalIsotope.AtomicMass, null, new List<double> { 0 }, null, null));
            fixedMods.Add(new ModificationWithMassAndCf("PepNmod", null, motif, ModificationSites.NPep, Chemistry.ChemicalFormula.ParseFormula("H"), GetElement(1).PrincipalIsotope.AtomicMass, null, new List<double> { 0 }, null, null));
            fixedMods.Add(new ModificationWithMassAndCf("resMod", null, motif, ModificationSites.Any, Chemistry.ChemicalFormula.ParseFormula("H"), GetElement(1).PrincipalIsotope.AtomicMass, null, new List<double> { 0 }, null, null));
            fixedMods.Add(new ModificationWithMassAndCf("PepCmod", null, motif, ModificationSites.PepC, Chemistry.ChemicalFormula.ParseFormula("H"), GetElement(1).PrincipalIsotope.AtomicMass, null, new List<double> { 0 }, null, null));
            fixedMods.Add(new ModificationWithMassAndCf("ProtCmod", null, motif, ModificationSites.ProtC, Chemistry.ChemicalFormula.ParseFormula("H"), GetElement(1).PrincipalIsotope.AtomicMass, null, new List<double> { 0 }, null, null));

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
            ModificationMotif motif;
            ModificationMotif.TryGetMotif("X", out motif);
            Modification mod = new ModificationWithMass(null, null, motif, ModificationSites.Any, double.NaN, null, new List<double> { double.NaN }, null, null);
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
            ModificationMotif motif;
            ModificationMotif.TryGetMotif("Abcdefg", out motif);
            Modification mod = new ModificationWithMass(null, null, motif, ModificationSites.Any, double.NaN, null, new List<double> { double.NaN }, null, null);
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