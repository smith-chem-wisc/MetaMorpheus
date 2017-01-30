using NUnit.Framework;
using OldInternalLogic;
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
            var prot = new Protein("MNNNKQQQQ", null, new Dictionary<int, List<MorpheusModification>>(), new int[0], new int[0], new string[0], null, null, 0, false, false);
            var protease = new Protease("Custom Protease", new List<string> { "K" }, new List<string>(), OldLogicTerminus.C, CleavageSpecificity.Full, null, null, null);

            var ye = prot.Digest(protease, 0, InitiatorMethionineBehavior.Retain).ToList();

            Assert.AreEqual(2, ye.Count);

            List<MorpheusModification> variableModifications = new List<MorpheusModification>();
            var pep1 = ye[0].GetPeptideWithSetModifications(variableModifications, 4096, 3, new List<MorpheusModification>()).First();
            Assert.IsTrue(pep1.MonoisotopicMass > 0);
            foreach (var huh in pep1.FastSortedProductMasses(new List<ProductType> { ProductType.B, ProductType.Y }))
                Assert.IsTrue(huh > 0);

            var pep2 = ye[1].GetPeptideWithSetModifications(variableModifications, 4096, 3, new List<MorpheusModification>()).First();
            Assert.IsTrue(pep2.MonoisotopicMass > 0);
            foreach (var huh in pep2.FastSortedProductMasses(new List<ProductType> { ProductType.B, ProductType.Y }))
                Assert.IsTrue(huh > 0);
        }

        [Test]
        public static void TestBadPeptide()
        {
            var prot = new Protein("MNNNKQQXQ", null, new Dictionary<int, List<MorpheusModification>>(), new int[0], new int[0], new string[0], null, null, 0, false, false);
            var protease = new Protease("Custom Protease", new List<string> { "K" }, new List<string>(), OldLogicTerminus.C, CleavageSpecificity.Full, null, null, null);

            var ye = prot.Digest(protease, 0, InitiatorMethionineBehavior.Retain).ToList();

            Assert.AreEqual(2, ye.Count);
            var pep1 = ye[0].GetPeptideWithSetModifications(new List<MorpheusModification>(), 4096, 3, new List<MorpheusModification>()).First();
            Assert.IsTrue(pep1.MonoisotopicMass > 0);
            foreach (var huh in pep1.FastSortedProductMasses(new List<ProductType> { ProductType.B, ProductType.Y }))
                Assert.IsTrue(huh > 0);

            var pep2 = ye[1].GetPeptideWithSetModifications(new List<MorpheusModification>(), 4096, 3, new List<MorpheusModification>()).First();
            Assert.IsNaN(pep2.MonoisotopicMass);
            var cool = pep2.FastSortedProductMasses(new List<ProductType> { ProductType.Y });
            Assert.IsTrue(cool[0] > 0);
            Assert.IsNaN(cool[1]);
            Assert.IsNaN(cool[2]);
        }

        [Test]
        public static void TestPeptideWithSetModifications()
        {
            var prot = new Protein("M", null, new Dictionary<int, List<MorpheusModification>>(), new int[0], new int[0], new string[0], null, null, 0, false, false);
            var protease = new Protease("Custom Protease", new List<string> { "K" }, new List<string>(), OldLogicTerminus.C, CleavageSpecificity.Full, null, null, null);
            var ye = prot.Digest(protease, 0, InitiatorMethionineBehavior.Retain).First();
            List<MorpheusModification> variableModifications = new List<MorpheusModification>();
            variableModifications.Add(new MorpheusModification("ProtNmod", ModificationType.ProteinNTerminus, 'M', null, '\0', GetElement(1).PrincipalIsotope.AtomicMass, double.NaN, double.NaN, new Chemistry.ChemicalFormula("H")));
            variableModifications.Add(new MorpheusModification("PepNmod", ModificationType.PeptideNTerminus, 'M', null, '\0', GetElement(1).PrincipalIsotope.AtomicMass, double.NaN, double.NaN, new Chemistry.ChemicalFormula("H")));
            variableModifications.Add(new MorpheusModification("resMod", ModificationType.AminoAcidResidue, 'M', null, '\0', GetElement(1).PrincipalIsotope.AtomicMass, double.NaN, double.NaN, new Chemistry.ChemicalFormula("H")));
            variableModifications.Add(new MorpheusModification("PepCmod", ModificationType.PeptideCTerminus, 'M', null, '\0', GetElement(1).PrincipalIsotope.AtomicMass, double.NaN, double.NaN, new Chemistry.ChemicalFormula("H")));
            variableModifications.Add(new MorpheusModification("ProtCmod", ModificationType.ProteinCTerminus, 'M', null, '\0', GetElement(1).PrincipalIsotope.AtomicMass, double.NaN, double.NaN, new Chemistry.ChemicalFormula("H")));
            var ok = ye.GetPeptideWithSetModifications(variableModifications, 4096, 5, new List<MorpheusModification>()).ToList();

            Assert.AreEqual(8, ok.Count);

            Assert.AreEqual("[:ProtNmod]M[:resMod][:PepCmod]", ok.Last().Sequence);
            Assert.AreEqual("[H]M[H][H]", ok.Last().SequenceWithChemicalFormulas);
            Assert.AreEqual(5 * GetElement("H").PrincipalIsotope.AtomicMass + Residue.ResidueMonoisotopicMass['M'] + GetElement("O").PrincipalIsotope.AtomicMass, ok.Last().MonoisotopicMass, 1e-9);
        }

        [Test]
        public static void TestPeptideWithFixedModifications()
        {
            var prot = new Protein("M", null, new Dictionary<int, List<MorpheusModification>>(), new int[0], new int[0], new string[0], null, null, 0, false, false);
            var protease = new Protease("Custom Protease", new List<string> { "K" }, new List<string>(), OldLogicTerminus.C, CleavageSpecificity.Full, null, null, null);
            var ye = prot.Digest(protease, 0, InitiatorMethionineBehavior.Retain).First();
            List<MorpheusModification> fixedMods = new List<MorpheusModification>();
            fixedMods.Add(new MorpheusModification("ProtNmod", ModificationType.ProteinNTerminus, 'M', null, '\0', GetElement(1).PrincipalIsotope.AtomicMass, double.NaN, double.NaN, new Chemistry.ChemicalFormula("H")));
            fixedMods.Add(new MorpheusModification("PepNmod", ModificationType.PeptideNTerminus, 'M', null, '\0', GetElement(1).PrincipalIsotope.AtomicMass, double.NaN, double.NaN, new Chemistry.ChemicalFormula("H")));
            fixedMods.Add(new MorpheusModification("resMod", ModificationType.AminoAcidResidue, 'M', null, '\0', GetElement(1).PrincipalIsotope.AtomicMass, double.NaN, double.NaN, new Chemistry.ChemicalFormula("H")));
            fixedMods.Add(new MorpheusModification("PepCmod", ModificationType.PeptideCTerminus, 'M', null, '\0', GetElement(1).PrincipalIsotope.AtomicMass, double.NaN, double.NaN, new Chemistry.ChemicalFormula("H")));
            fixedMods.Add(new MorpheusModification("ProtCmod", ModificationType.ProteinCTerminus, 'M', null, '\0', GetElement(1).PrincipalIsotope.AtomicMass, double.NaN, double.NaN, new Chemistry.ChemicalFormula("H")));

            var ok = ye.GetPeptideWithSetModifications(new List<MorpheusModification>(), 4096, 5, fixedMods).ToList();

            Assert.AreEqual(1, ok.Count);

            Assert.AreEqual("[:ProtNmod]M[:resMod][:PepCmod]", ok.Last().Sequence);
            Assert.AreEqual("[H]M[H][H]", ok.Last().SequenceWithChemicalFormulas);
            Assert.AreEqual(5 * GetElement("H").PrincipalIsotope.AtomicMass + Residue.ResidueMonoisotopicMass['M'] + GetElement("O").PrincipalIsotope.AtomicMass, ok.Last().MonoisotopicMass, 1e-9);
        }

        #endregion Public Methods

    }
}