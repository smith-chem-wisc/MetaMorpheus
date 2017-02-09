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
            var prot = new Protein("MNNNKQQQQ", null, new Dictionary<int, List<ModificationWithMass>>(), new int[0], new int[0], new string[0], null, null, 0, false, false);
            var protease = new Protease("Custom Protease", new List<string> { "K" }, new List<string>(), TerminusType.C, CleavageSpecificity.Full, null, null, null);

            var ye = prot.Digest(protease, 0, InitiatorMethionineBehavior.Retain, new List<ModificationWithMass>()).ToList();

            Assert.AreEqual(2, ye.Count);
            
            List<ModificationWithMass> variableModifications = new List<ModificationWithMass>();
            var pep1 = ye[0].GetPeptideWithSetModifications(variableModifications, 4096, 3).First();
            Assert.IsTrue(pep1.MonoisotopicMass > 0);
            foreach (var huh in pep1.FastSortedProductMasses(new List<ProductType> { ProductType.B, ProductType.Y }))
                Assert.IsTrue(huh > 0);

            var pep2 = ye[1].GetPeptideWithSetModifications(variableModifications, 4096, 3).First();
            Assert.IsTrue(pep2.MonoisotopicMass > 0);
            foreach (var huh in pep2.FastSortedProductMasses(new List<ProductType> { ProductType.B, ProductType.Y }))
                Assert.IsTrue(huh > 0);
        }

        [Test]
        public static void TestNoCleavage()
        {
            List<ModificationWithMass> fixedModifications = new List<ModificationWithMass>();
            var prot = new Protein("MNNNKQQQQ", null, new Dictionary<int, List<ModificationWithMass>>(), new int[] { 5 }, new int[] { 6 }, new string[] { "lala" }, null, null, 0, false, false);
            var protease = new Protease("Custom Protease", null, null, TerminusType.None, CleavageSpecificity.None, null, null, null);

            var ye = prot.Digest(protease, int.MaxValue, InitiatorMethionineBehavior.Variable, fixedModifications).ToList();

            Assert.AreEqual(3, ye.Count);
        }

        [Test]
        public static void TestBadPeptide()
        {
            var prot = new Protein("MNNNKQQXQ", null, new Dictionary<int, List<ModificationWithMass>>(), new int[0], new int[0], new string[0], null, null, 0, false, false);
            var protease = new Protease("Custom Protease", new List<string> { "K" }, new List<string>(), TerminusType.C, CleavageSpecificity.Full, null, null, null);

            var ye = prot.Digest(protease, 0, InitiatorMethionineBehavior.Retain, new List<ModificationWithMass>()).ToList();

            Assert.AreEqual(2, ye.Count);
            var pep1 = ye[0].GetPeptideWithSetModifications(new List<ModificationWithMass>(), 4096, 3).First();
            Assert.IsTrue(pep1.MonoisotopicMass > 0);
            foreach (var huh in pep1.FastSortedProductMasses(new List<ProductType> { ProductType.B, ProductType.Y }))
                Assert.IsTrue(huh > 0);

            var pep2 = ye[1].GetPeptideWithSetModifications(new List<ModificationWithMass>(), 4096, 3).First();
            Assert.IsNaN(pep2.MonoisotopicMass);
            var cool = pep2.FastSortedProductMasses(new List<ProductType> { ProductType.Y });
            Assert.IsTrue(cool[0] > 0);
            Assert.IsNaN(cool[1]);
            Assert.IsNaN(cool[2]);
        }

        [Test]
        public static void TestPeptideWithSetModifications()
        {
            var prot = new Protein("M", null, new Dictionary<int, List<ModificationWithMass>>(), new int[0], new int[0], new string[0], null, null, 0, false, false);
            var protease = new Protease("Custom Protease", new List<string> { "K" }, new List<string>(), TerminusType.C, CleavageSpecificity.Full, null, null, null);
            var ye = prot.Digest(protease, 0, InitiatorMethionineBehavior.Retain, new List<ModificationWithMass>()).First();
            List<ModificationWithMass> variableModifications = new List<ModificationWithMass>();
            variableModifications.Add(new ModificationWithMass("ProtNmod", ModificationType.ProteinNTerminus, 'M', null, '\0', GetElement(1).PrincipalIsotope.AtomicMass, double.NaN, double.NaN, new Chemistry.ChemicalFormula("H")));
            variableModifications.Add(new ModificationWithMass("PepNmod", ModificationType.PeptideNTerminus, 'M', null, '\0', GetElement(1).PrincipalIsotope.AtomicMass, double.NaN, double.NaN, new Chemistry.ChemicalFormula("H")));
            variableModifications.Add(new ModificationWithMass("resMod", ModificationType.AminoAcidResidue, 'M', null, '\0', GetElement(1).PrincipalIsotope.AtomicMass, double.NaN, double.NaN, new Chemistry.ChemicalFormula("H")));
            variableModifications.Add(new ModificationWithMass("PepCmod", ModificationType.PeptideCTerminus, 'M', null, '\0', GetElement(1).PrincipalIsotope.AtomicMass, double.NaN, double.NaN, new Chemistry.ChemicalFormula("H")));
            variableModifications.Add(new ModificationWithMass("ProtCmod", ModificationType.ProteinCTerminus, 'M', null, '\0', GetElement(1).PrincipalIsotope.AtomicMass, double.NaN, double.NaN, new Chemistry.ChemicalFormula("H")));
            var ok = ye.GetPeptideWithSetModifications(variableModifications, 4096, 5).ToList();

            Assert.AreEqual(8, ok.Count);

            Assert.AreEqual("[:ProtNmod]M[:resMod][:PepCmod]", ok.Last().Sequence);
            Assert.AreEqual("[H]M[H][H]", ok.Last().SequenceWithChemicalFormulas);
            Assert.AreEqual(5 * GetElement("H").PrincipalIsotope.AtomicMass + Residue.ResidueMonoisotopicMass['M'] + GetElement("O").PrincipalIsotope.AtomicMass, ok.Last().MonoisotopicMass, 1e-9);
        }

        [Test]
        public static void TestPeptideWithFixedModifications()
        {
            var prot = new Protein("M", null, new Dictionary<int, List<ModificationWithMass>>(), new int[0], new int[0], new string[0], null, null, 0, false, false);
            var protease = new Protease("Custom Protease", new List<string> { "K" }, new List<string>(), TerminusType.C, CleavageSpecificity.Full, null, null, null);
            List<ModificationWithMass> fixedMods = new List<ModificationWithMass>();
            fixedMods.Add(new ModificationWithMass("ProtNmod", ModificationType.ProteinNTerminus, 'M', null, '\0', GetElement(1).PrincipalIsotope.AtomicMass, double.NaN, double.NaN, new Chemistry.ChemicalFormula("H")));
            fixedMods.Add(new ModificationWithMass("PepNmod", ModificationType.PeptideNTerminus, 'M', null, '\0', GetElement(1).PrincipalIsotope.AtomicMass, double.NaN, double.NaN, new Chemistry.ChemicalFormula("H")));
            fixedMods.Add(new ModificationWithMass("resMod", ModificationType.AminoAcidResidue, 'M', null, '\0', GetElement(1).PrincipalIsotope.AtomicMass, double.NaN, double.NaN, new Chemistry.ChemicalFormula("H")));
            fixedMods.Add(new ModificationWithMass("PepCmod", ModificationType.PeptideCTerminus, 'M', null, '\0', GetElement(1).PrincipalIsotope.AtomicMass, double.NaN, double.NaN, new Chemistry.ChemicalFormula("H")));
            fixedMods.Add(new ModificationWithMass("ProtCmod", ModificationType.ProteinCTerminus, 'M', null, '\0', GetElement(1).PrincipalIsotope.AtomicMass, double.NaN, double.NaN, new Chemistry.ChemicalFormula("H")));

            var ye = prot.Digest(protease, 0, InitiatorMethionineBehavior.Retain, fixedMods).First();
            var ok = ye.GetPeptideWithSetModifications(new List<ModificationWithMass>(), 4096, 5).ToList();

            Assert.AreEqual(1, ok.Count);

            Assert.AreEqual("[:PepNmod]M[:resMod][:ProtCmod]", ok.Last().Sequence);
            Assert.AreEqual("[H]M[H][H]", ok.Last().SequenceWithChemicalFormulas);
            Assert.AreEqual(5 * GetElement("H").PrincipalIsotope.AtomicMass + Residue.ResidueMonoisotopicMass['M'] + GetElement("O").PrincipalIsotope.AtomicMass, ok.Last().MonoisotopicMass, 1e-9);
        }

        #endregion Public Methods

    }
}