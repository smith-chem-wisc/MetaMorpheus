using EngineLayer;
using NUnit.Framework;

using Proteomics;
using System.Collections.Generic;
using System.Linq;
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
            var prot = new Protein("MNNNKQQQQ", null, new Dictionary<int, List<Modification>>(), new int?[0], new int?[0], new string[0], null, null, 0, false, false);
            var protease = new Protease("Custom Protease", new List<string> { "K" }, new List<string>(), TerminusType.C, CleavageSpecificity.Full, null, null, null);

            var ye = prot.Digest(protease, 0, InitiatorMethionineBehavior.Retain, new List<ModificationWithMass>()).ToList();

            Assert.AreEqual(2, ye.Count);

            List<ModificationWithMass> variableModifications = new List<ModificationWithMass>();
            var pep1 = ye[0].GetPeptidesWithSetModifications(variableModifications, 4096, 3).First();
            Assert.IsTrue(pep1.MonoisotopicMass > 0);
            foreach (var huh in pep1.FastSortedProductMasses(new List<ProductType> { ProductType.B, ProductType.Y }))
                Assert.IsTrue(huh > 0);

            var pep2 = ye[1].GetPeptidesWithSetModifications(variableModifications, 4096, 3).First();
            Assert.IsTrue(pep2.MonoisotopicMass > 0);
            foreach (var huh in pep2.FastSortedProductMasses(new List<ProductType> { ProductType.B, ProductType.Y }))
                Assert.IsTrue(huh > 0);
        }

        [Test]
        public static void TestNoCleavage()
        {
            List<ModificationWithMass> fixedModifications = new List<ModificationWithMass>();
            var prot = new Protein("MNNNKQQQQ", null, new Dictionary<int, List<Modification>>(), new int?[] { 5 }, new int?[] { 6 }, new string[] { "lala" }, null, null, 0, false, false);
            var protease = new Protease("Custom Protease", null, null, TerminusType.None, CleavageSpecificity.None, null, null, null);

            var ye = prot.Digest(protease, int.MaxValue, InitiatorMethionineBehavior.Variable, fixedModifications).ToList();

            Assert.AreEqual(3, ye.Count);
        }

        [Test]
        public static void TestBadPeptide()
        {
            var prot = new Protein("MNNNKQQXQ", null, new Dictionary<int, List<Modification>>(), new int?[0], new int?[0], new string[0], null, null, 0, false, false);
            var protease = new Protease("Custom Protease", new List<string> { "K" }, new List<string>(), TerminusType.C, CleavageSpecificity.Full, null, null, null);

            var ye = prot.Digest(protease, 0, InitiatorMethionineBehavior.Retain, new List<ModificationWithMass>()).ToList();

            Assert.AreEqual(2, ye.Count);
            var pep1 = ye[0].GetPeptidesWithSetModifications(new List<ModificationWithMass>(), 4096, 3).First();
            Assert.IsTrue(pep1.MonoisotopicMass > 0);
            foreach (var huh in pep1.FastSortedProductMasses(new List<ProductType> { ProductType.B, ProductType.Y }))
                Assert.IsTrue(huh > 0);

            var pep2 = ye[1].GetPeptidesWithSetModifications(new List<ModificationWithMass>(), 4096, 3).First();
            Assert.IsNaN(pep2.MonoisotopicMass);
            var cool = pep2.FastSortedProductMasses(new List<ProductType> { ProductType.Y });
            Assert.IsTrue(cool[0] > 0);
            Assert.IsNaN(cool[1]);
            Assert.IsNaN(cool[2]);
        }

        [Test]
        public static void TestPeptideWithSetModifications()
        {
            var prot = new Protein("M", null, new Dictionary<int, List<Modification>>(), new int?[0], new int?[0], new string[0], null, null, 0, false, false);
            var protease = new Protease("Custom Protease", new List<string> { "K" }, new List<string>(), TerminusType.C, CleavageSpecificity.Full, null, null, null);
            var ye = prot.Digest(protease, 0, InitiatorMethionineBehavior.Retain, new List<ModificationWithMass>()).First();
            List<ModificationWithMass> variableModifications = new List<ModificationWithMass>();
            ModificationMotif motif;
            ModificationMotif.TryGetMotif("M", out motif);
            variableModifications.Add(new ModificationWithMassAndCf("ProtNmod", null, motif, ModificationSites.NProt, Chemistry.ChemicalFormula.ParseFormula("H"), GetElement(1).PrincipalIsotope.AtomicMass, null, 0, null, null, null));
            variableModifications.Add(new ModificationWithMassAndCf("pepNmod", null, motif, ModificationSites.NPep, Chemistry.ChemicalFormula.ParseFormula("H"), GetElement(1).PrincipalIsotope.AtomicMass, null, 0, null, null, null));
            variableModifications.Add(new ModificationWithMassAndCf("resMod", null, motif, ModificationSites.Any, Chemistry.ChemicalFormula.ParseFormula("H"), GetElement(1).PrincipalIsotope.AtomicMass, null, 0, null, null, null));
            variableModifications.Add(new ModificationWithMassAndCf("PepCmod", null, motif, ModificationSites.PepC, Chemistry.ChemicalFormula.ParseFormula("H"), GetElement(1).PrincipalIsotope.AtomicMass, null, 0, null, null, null));
            variableModifications.Add(new ModificationWithMassAndCf("ProtCmod", null, motif, ModificationSites.ProtC, Chemistry.ChemicalFormula.ParseFormula("H"), GetElement(1).PrincipalIsotope.AtomicMass, null, 0, null, null, null));
            var ok = ye.GetPeptidesWithSetModifications(variableModifications, 4096, 5).ToList();
            Assert.AreEqual(8, ok.Count);

            Assert.AreEqual("[:ProtNmod]M[:resMod][:PepCmod]", ok.Last().Sequence);
            Assert.AreEqual("[H]M[H][H]", ok.Last().SequenceWithChemicalFormulas);
            Assert.AreEqual(5 * GetElement("H").PrincipalIsotope.AtomicMass + Residue.ResidueMonoisotopicMass['M'] + GetElement("O").PrincipalIsotope.AtomicMass, ok.Last().MonoisotopicMass, 1e-9);
        }

        [Test]
        public static void TestPeptideWithFixedModifications()
        {
            var prot = new Protein("M", null, new Dictionary<int, List<Modification>>(), new int?[0], new int?[0], new string[0], null, null, 0, false, false);
            var protease = new Protease("Custom Protease", new List<string> { "K" }, new List<string>(), TerminusType.C, CleavageSpecificity.Full, null, null, null);
            List<ModificationWithMass> fixedMods = new List<ModificationWithMass>();
            ModificationMotif motif;
            ModificationMotif.TryGetMotif("M", out motif);
            fixedMods.Add(new ModificationWithMassAndCf("ProtNmod", null, motif, ModificationSites.NProt, Chemistry.ChemicalFormula.ParseFormula("H"), GetElement(1).PrincipalIsotope.AtomicMass, null, 0, null, null, null));
            fixedMods.Add(new ModificationWithMassAndCf("PepNmod", null, motif, ModificationSites.NPep, Chemistry.ChemicalFormula.ParseFormula("H"), GetElement(1).PrincipalIsotope.AtomicMass, null, 0, null, null, null));
            fixedMods.Add(new ModificationWithMassAndCf("resMod", null, motif, ModificationSites.Any, Chemistry.ChemicalFormula.ParseFormula("H"), GetElement(1).PrincipalIsotope.AtomicMass, null, 0, null, null, null));
            fixedMods.Add(new ModificationWithMassAndCf("PepCmod", null, motif, ModificationSites.PepC, Chemistry.ChemicalFormula.ParseFormula("H"), GetElement(1).PrincipalIsotope.AtomicMass, null, 0, null, null, null));
            fixedMods.Add(new ModificationWithMassAndCf("ProtCmod", null, motif, ModificationSites.ProtC, Chemistry.ChemicalFormula.ParseFormula("H"), GetElement(1).PrincipalIsotope.AtomicMass, null, 0, null, null, null));

            var ye = prot.Digest(protease, 0, InitiatorMethionineBehavior.Retain, fixedMods).First();
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
            Modification mod = new ModificationWithMass(null, null, motif, ModificationSites.Any, double.NaN, null, double.NaN, null, null, null);
            IDictionary<int, List<Modification>> modDict = new Dictionary<int, List<Modification>>
            {
                {2, new List<Modification> {mod } },
                {8, new List<Modification> {mod } }
            };
            var prot = new Protein("MNNNNKRRRRR", null, modDict, new int?[0], new int?[0], new string[0], null, null, 0, false, false);
            var protease = new Protease("Custom Protease", new List<string> { "K" }, new List<string>(), TerminusType.C, CleavageSpecificity.Full, null, null, null);
            var ye1 = prot.Digest(protease, 0, InitiatorMethionineBehavior.Retain, new List<ModificationWithMass>()).First();
            var ye2 = prot.Digest(protease, 0, InitiatorMethionineBehavior.Retain, new List<ModificationWithMass>()).Last();
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
            Modification mod = new ModificationWithMass(null, null, motif, ModificationSites.Any, double.NaN, null, double.NaN, null, null, null);
            IDictionary<int, List<Modification>> modDict = new Dictionary<int, List<Modification>>
            {
                {2, new List<Modification> {mod } },
                {8, new List<Modification> {mod } }
            };
            var prot = new Protein("MNNNNKRRRRR", null, modDict, new int?[0], new int?[0], new string[0], null, null, 0, true, false);
            var protease = new Protease("Custom Protease", new List<string> { "K" }, new List<string>(), TerminusType.C, CleavageSpecificity.Full, null, null, null);
            var ye1 = prot.Digest(protease, 0, InitiatorMethionineBehavior.Retain, new List<ModificationWithMass>()).First();
            var ye2 = prot.Digest(protease, 0, InitiatorMethionineBehavior.Retain, new List<ModificationWithMass>()).Last();
            var ok1 = ye1.GetPeptidesWithSetModifications(new List<ModificationWithMass>(), 2, 1).Last();
            var ok2 = ye2.GetPeptidesWithSetModifications(new List<ModificationWithMass>(), 2, 1).Last();

            Assert.AreEqual(1, ok1.NumMods);
            Assert.IsTrue(ok1.allModsOneIsNterminus.ContainsKey(3));
            Assert.AreEqual(1, ok2.NumMods);
            Assert.IsTrue(ok2.allModsOneIsNterminus.ContainsKey(3));


             prot = new Protein("MNNNNKRRRRR", null, modDict, new int?[0], new int?[0], new string[0], null, null, 0, false, false);
             ye1 = prot.Digest(protease, 0, InitiatorMethionineBehavior.Retain, new List<ModificationWithMass>()).First();
             ye2 = prot.Digest(protease, 0, InitiatorMethionineBehavior.Retain, new List<ModificationWithMass>()).Last();
             ok1 = ye1.GetPeptidesWithSetModifications(new List<ModificationWithMass>(), 2, 1).Last();
             ok2 = ye2.GetPeptidesWithSetModifications(new List<ModificationWithMass>(), 2, 1).Last();

            Assert.AreEqual(0, ok1.NumMods);
            Assert.IsFalse(ok1.allModsOneIsNterminus.Any());
            Assert.AreEqual(0, ok2.NumMods);
            Assert.IsFalse(ok2.allModsOneIsNterminus.Any());


        }
        #endregion Public Methods

    }
}