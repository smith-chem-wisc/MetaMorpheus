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
            var prot = new Protein("MNNNKQQQQ", null, new Dictionary<int, List<MetaMorpheusModification>>(), new int[0], new int[0], new string[0], null, null, 0, false, false);
            var protease = new Protease("Custom Protease", new List<string> { "K" }, new List<string>(), TerminusType.C, CleavageSpecificity.Full, null, null, null);

            var ye = prot.Digest(protease, 0, InitiatorMethionineBehavior.Retain).ToList();

            Assert.AreEqual(2, ye.Count);

            List<MetaMorpheusModification> variableModifications = new List<MetaMorpheusModification>();
            var pep1 = ye[0].GetPeptideWithSetModifications(variableModifications, 4096, 3, new List<MetaMorpheusModification>()).First();
            Assert.IsTrue(pep1.MonoisotopicMass > 0);
            foreach (var huh in pep1.FastSortedProductMasses(new List<ProductType> { ProductType.B, ProductType.Y }))
                Assert.IsTrue(huh > 0);

            var pep2 = ye[1].GetPeptideWithSetModifications(variableModifications, 4096, 3, new List<MetaMorpheusModification>()).First();
            Assert.IsTrue(pep2.MonoisotopicMass > 0);
            foreach (var huh in pep2.FastSortedProductMasses(new List<ProductType> { ProductType.B, ProductType.Y }))
                Assert.IsTrue(huh > 0);
        }

        [Test]
        public static void TestBadPeptide()
        {
            var prot = new Protein("MNNNKQQXQ", null, new Dictionary<int, List<MetaMorpheusModification>>(), new int[0], new int[0], new string[0], null, null, 0, false, false);
            var protease = new Protease("Custom Protease", new List<string> { "K" }, new List<string>(), TerminusType.C, CleavageSpecificity.Full, null, null, null);

            var ye = prot.Digest(protease, 0, InitiatorMethionineBehavior.Retain).ToList();

            Assert.AreEqual(2, ye.Count);
            var pep1 = ye[0].GetPeptideWithSetModifications(new List<MetaMorpheusModification>(), 4096, 3, new List<MetaMorpheusModification>()).First();
            Assert.IsTrue(pep1.MonoisotopicMass > 0);
            foreach (var huh in pep1.FastSortedProductMasses(new List<ProductType> { ProductType.B, ProductType.Y }))
                Assert.IsTrue(huh > 0);

            var pep2 = ye[1].GetPeptideWithSetModifications(new List<MetaMorpheusModification>(), 4096, 3, new List<MetaMorpheusModification>()).First();
            Assert.IsNaN(pep2.MonoisotopicMass);
            var cool = pep2.FastSortedProductMasses(new List<ProductType> { ProductType.Y });
            Assert.IsTrue(cool[0] > 0);
            Assert.IsNaN(cool[1]);
            Assert.IsNaN(cool[2]);
        }

        [Test]
        public static void TestPeptideWithSetModifications()
        {
            var prot = new Protein("M", null, new Dictionary<int, List<MetaMorpheusModification>>(), new int[0], new int[0], new string[0], null, null, 0, false, false);
            var protease = new Protease("Custom Protease", new List<string> { "K" }, new List<string>(), TerminusType.C, CleavageSpecificity.Full, null, null, null);
            var ye = prot.Digest(protease, 0, InitiatorMethionineBehavior.Retain).First();
            List<MetaMorpheusModification> variableModifications = new List<MetaMorpheusModification>();
            variableModifications.Add(new MetaMorpheusModification("ProtNmod", ModificationType.ProteinNTerminus, 'M', null, '\0', GetElement(1).PrincipalIsotope.AtomicMass, double.NaN, double.NaN, new Chemistry.ChemicalFormula("H")));
            variableModifications.Add(new MetaMorpheusModification("PepNmod", ModificationType.PeptideNTerminus, 'M', null, '\0', GetElement(1).PrincipalIsotope.AtomicMass, double.NaN, double.NaN, new Chemistry.ChemicalFormula("H")));
            variableModifications.Add(new MetaMorpheusModification("resMod", ModificationType.AminoAcidResidue, 'M', null, '\0', GetElement(1).PrincipalIsotope.AtomicMass, double.NaN, double.NaN, new Chemistry.ChemicalFormula("H")));
            variableModifications.Add(new MetaMorpheusModification("PepCmod", ModificationType.PeptideCTerminus, 'M', null, '\0', GetElement(1).PrincipalIsotope.AtomicMass, double.NaN, double.NaN, new Chemistry.ChemicalFormula("H")));
            variableModifications.Add(new MetaMorpheusModification("ProtCmod", ModificationType.ProteinCTerminus, 'M', null, '\0', GetElement(1).PrincipalIsotope.AtomicMass, double.NaN, double.NaN, new Chemistry.ChemicalFormula("H")));
            var ok = ye.GetPeptideWithSetModifications(variableModifications, 4096, 5, new List<MetaMorpheusModification>()).ToList();

            Assert.AreEqual(8, ok.Count);

            Assert.AreEqual("[:ProtNmod]M[:resMod][:PepCmod]", ok.Last().Sequence);
            Assert.AreEqual("[H]M[H][H]", ok.Last().SequenceWithChemicalFormulas);
            Assert.AreEqual(5 * GetElement("H").PrincipalIsotope.AtomicMass + Residue.ResidueMonoisotopicMass['M'] + GetElement("O").PrincipalIsotope.AtomicMass, ok.Last().MonoisotopicMass, 1e-9);
        }

        [Test]
        public static void TestPeptideWithFixedModifications()
        {
            var prot = new Protein("M", null, new Dictionary<int, List<MetaMorpheusModification>>(), new int[0], new int[0], new string[0], null, null, 0, false, false);
            var protease = new Protease("Custom Protease", new List<string> { "K" }, new List<string>(), TerminusType.C, CleavageSpecificity.Full, null, null, null);
            var ye = prot.Digest(protease, 0, InitiatorMethionineBehavior.Retain).First();
            List<MetaMorpheusModification> fixedMods = new List<MetaMorpheusModification>();
            fixedMods.Add(new MetaMorpheusModification("ProtNmod", ModificationType.ProteinNTerminus, 'M', null, '\0', GetElement(1).PrincipalIsotope.AtomicMass, double.NaN, double.NaN, new Chemistry.ChemicalFormula("H")));
            fixedMods.Add(new MetaMorpheusModification("PepNmod", ModificationType.PeptideNTerminus, 'M', null, '\0', GetElement(1).PrincipalIsotope.AtomicMass, double.NaN, double.NaN, new Chemistry.ChemicalFormula("H")));
            fixedMods.Add(new MetaMorpheusModification("resMod", ModificationType.AminoAcidResidue, 'M', null, '\0', GetElement(1).PrincipalIsotope.AtomicMass, double.NaN, double.NaN, new Chemistry.ChemicalFormula("H")));
            fixedMods.Add(new MetaMorpheusModification("PepCmod", ModificationType.PeptideCTerminus, 'M', null, '\0', GetElement(1).PrincipalIsotope.AtomicMass, double.NaN, double.NaN, new Chemistry.ChemicalFormula("H")));
            fixedMods.Add(new MetaMorpheusModification("ProtCmod", ModificationType.ProteinCTerminus, 'M', null, '\0', GetElement(1).PrincipalIsotope.AtomicMass, double.NaN, double.NaN, new Chemistry.ChemicalFormula("H")));

            var ok = ye.GetPeptideWithSetModifications(new List<MetaMorpheusModification>(), 4096, 5, fixedMods).ToList();

            Assert.AreEqual(1, ok.Count);

            Assert.AreEqual("[:ProtNmod]M[:resMod][:PepCmod]", ok.Last().Sequence);
            Assert.AreEqual("[H]M[H][H]", ok.Last().SequenceWithChemicalFormulas);
            Assert.AreEqual(5 * GetElement("H").PrincipalIsotope.AtomicMass + Residue.ResidueMonoisotopicMass['M'] + GetElement("O").PrincipalIsotope.AtomicMass, ok.Last().MonoisotopicMass, 1e-9);
        }

        #endregion Public Methods

    }
}