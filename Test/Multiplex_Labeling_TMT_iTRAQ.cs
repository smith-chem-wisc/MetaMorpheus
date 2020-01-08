using Chemistry;
using EngineLayer;
using MassSpectrometry;
using NUnit.Framework;
using Proteomics;
using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System.Collections.Generic;
using System.Linq;

namespace Test
{
    [TestFixture]
    internal static class Multiplex_Labeling_TMT_iTRAQ
    {
        [Test]
        [TestCase("C8 N1 H16", 126.128274520)]
        [TestCase("C8 H16 N{15}1", 127.125309415)]
        [TestCase("C7 H16 C{13}1 N{15}1", 128.128664250)]
        [TestCase("C6 H16 C{13}2 N{15}1", 129.132019085)]
        [TestCase("C5 H16 C{13}3 N{15}1", 130.135373920)]
        [TestCase("C3 N1 H16 C{13}5", 131.145048695)]
        [TestCase("C7 N1 H16 C{13}1", 127.131629355)]
        [TestCase("C6 N1 H16 C{13}2", 128.134984190)]
        [TestCase("C5 N1 H16 C{13}3", 129.138339025)]
        [TestCase("C4 N1 H16 C{13}4", 130.141693860)]
        [TestCase("C4 H16 C{13}4 N{15}1", 131.138728755)]
        public static void TestChemicalFormulaWithIsotopesTMT(string formula, double mass)
        {
            ChemicalFormula cf = ChemicalFormula.ParseFormula(formula);
            Assert.AreEqual(mass, ClassExtensions.RoundedDouble(cf.MonoisotopicMass));
        }

        [Test]
        [TestCase("PEPTIDE", 1029.5302)]
        [TestCase("PEPTIDEK", 1386.7881)]
        public static void TestPeptideLabelledWithTMT(string peptide, double totalMass)
        {
            List<Modification> gptmdModifications = new List<Modification>();
            gptmdModifications.AddRange(GlobalVariables.AllModsKnown);
            List<Modification> tmt10Mods = gptmdModifications.Where(m => m.ModificationType == "Multiplex Label" && m.IdWithMotif.Contains("TMT10")).ToList();

            Protein P = new Protein(peptide, "", "", null, null, null, null, null, false, false, null, null, null, null);
            CommonParameters CommonParameters = new CommonParameters(digestionParams: new DigestionParams(minPeptideLength: 1));
            var p = P.Digest(CommonParameters.DigestionParams, tmt10Mods, new List<Modification>()).First();
            var f = p.Fragment(DissociationType.HCD, FragmentationTerminus.Both);

            List<double> productMasses = f.Select(m => m.NeutralMass.ToMz(1)).ToList();
            productMasses.Distinct();
            productMasses.Sort();

            Assert.AreEqual(totalMass, ClassExtensions.RoundedDouble(p.MonoisotopicMass.ToMz(1), 4));
        }

        [Test]
        [TestCase("PEPTIDE", 944.4712)]
        [TestCase("PEPTIDEK", 1216.6702)]
        public static void TestPeptideLabelledWith_iTRAQ_4plex(string peptide, double totalMass)
        {
            List<Modification> gptmdModifications = new List<Modification>();
            gptmdModifications.AddRange(GlobalVariables.AllModsKnown);
            List<Modification> itraq4plex = gptmdModifications.Where(m => m.ModificationType == "Multiplex Label" && m.IdWithMotif.Contains("iTRAQ-4plex")).ToList();

            Protein P = new Protein(peptide, "", "", null, null, null, null, null, false, false, null, null, null, null);
            CommonParameters CommonParameters = new CommonParameters(digestionParams: new DigestionParams(minPeptideLength: 1));
            var p = P.Digest(CommonParameters.DigestionParams, itraq4plex, new List<Modification>()).First();
            var f = p.Fragment(DissociationType.HCD, FragmentationTerminus.Both);

            List<double> productMasses = f.Select(m => m.NeutralMass.ToMz(1)).ToList();
            productMasses.Distinct();
            productMasses.Sort();

            Assert.AreEqual(totalMass, ClassExtensions.RoundedDouble(p.MonoisotopicMass.ToMz(1), 4));
        }

        [Test]
        [TestCase("PEPTIDE", 945.4975)]
        [TestCase("PEPTIDEK", 1218.7227)]
        public static void TestPeptideLabelledWith_DiLeu_4plex(string peptide, double totalMass)
        {
            List<Modification> gptmdModifications = new List<Modification>();
            gptmdModifications.AddRange(GlobalVariables.AllModsKnown);
            List<Modification> itraq4plex = gptmdModifications.Where(m => m.ModificationType == "Multiplex Label" && m.IdWithMotif.Contains("DiLeu-4plex")).ToList();

            Protein P = new Protein(peptide, "", "", null, null, null, null, null, false, false, null, null, null, null);
            CommonParameters CommonParameters = new CommonParameters(digestionParams: new DigestionParams(minPeptideLength: 1));
            var p = P.Digest(CommonParameters.DigestionParams, itraq4plex, new List<Modification>()).First();
            var f = p.Fragment(DissociationType.HCD, FragmentationTerminus.Both);

            List<double> productMasses = f.Select(m => m.NeutralMass.ToMz(1)).ToList();
            productMasses.Distinct();
            productMasses.Sort();

            Assert.AreEqual(totalMass, ClassExtensions.RoundedDouble(p.MonoisotopicMass.ToMz(1), 4));
        }

        [Test]
        [TestCase("PEPTIDE", 1104.5694)]
        [TestCase("PEPTIDEK", 1536.8666)]
        public static void TestPeptideLabelledWith_iTRAQ_8plex(string peptide, double totalMass)
        {
            List<Modification> gptmdModifications = new List<Modification>();
            gptmdModifications.AddRange(GlobalVariables.AllModsKnown);
            List<Modification> itraq8plex = gptmdModifications.Where(m => m.ModificationType == "Multiplex Label" && m.IdWithMotif.Contains("iTRAQ-8plex")).ToList();

            Protein P = new Protein(peptide, "", "", null, null, null, null, null, false, false, null, null, null, null);
            CommonParameters CommonParameters = new CommonParameters(digestionParams: new DigestionParams(minPeptideLength: 1));
            var p = P.Digest(CommonParameters.DigestionParams, itraq8plex, new List<Modification>()).First();
            var f = p.Fragment(DissociationType.HCD, FragmentationTerminus.Both);

            List<double> productMasses = f.Select(m => m.NeutralMass.ToMz(1)).ToList();
            productMasses.Distinct();
            productMasses.Sort();

            Assert.AreEqual(totalMass, ClassExtensions.RoundedDouble(p.MonoisotopicMass.ToMz(1), 4));
        }

        [Test]
        [TestCase("C5 N2 H12 C{13}1", 114.110679698, true)]
        //[TestCase("C4 N2 H12 C{13}2", 115.114034533, true)] this is old style (ABI or Thermo, don't know). no longer used
        [TestCase("C5 C{13}1 N1 N{15}1 H12", 115.107714592, true)]
        [TestCase("C4 N1 H12 C{13}2 N{15}1", 116.111069427, true)]
        [TestCase("C3 N1 H12 C{13}3 N{15}1", 117.114424262, true)]
        [TestCase("C5 N2 H12 C{13}2 O{18}1", 144.105917679, false)]
        [TestCase("C4 N1O1 H12 C{13}3 N{15}1", 144.102062415, false)]
        [TestCase("C7 N3O3 H24 C{13}7 N{15}1", 304.205359390, false)]//for tmt 8-plex 113, 114, 116, 117
        [TestCase("C8 N2O3 H24 C{13}6 N{15}2", 304.199039449, false)]//for tmt 8-plex 115, 118, 119, 121
        public static void TestChemicalFormulaWithIsotopes_iTRAQ(string formula, double mass, bool mz)
        {
            ChemicalFormula cf = ChemicalFormula.ParseFormula(formula);
            if (mz)
            {
                Assert.AreEqual(mass, ClassExtensions.RoundedDouble(cf.MonoisotopicMass.ToMz(1)));
            }
            else
            {
                Assert.AreEqual(mass, ClassExtensions.RoundedDouble(cf.MonoisotopicMass));
            }
        }

        [Test]
        [TestCase("C7 N{15}1 H15", 115.124760849, true)]
        [TestCase("C7 N1 H13 H{2}2", 116.140279447, true)]
        [TestCase("C7 N{15}1 H13 H{2}2", 117.137314341, true)]
        [TestCase("C7 N1 H11 H{2}4", 118.152832938, true)]
        [TestCase("C7 N{15}1 H15 C{13}1 O{18}1", 145.119998830, false)]
        [TestCase("C8 N1 H13 H{2}2 O{18}1", 145.132162593, false)]
        [TestCase("C7 N{15}1 H13 H{2}2 C{13}1 O1", 145.128307329, false)]
        [TestCase("C8 N1 H11 H{2}4 O1", 145.140471091, false)]
        public static void TestChemicalFormulaWithIsotopes_DiLeu4plex(string formula, double mass, bool mz)
        {
            ChemicalFormula cf = ChemicalFormula.ParseFormula(formula);
            if (mz)
            {
                Assert.AreEqual(mass, ClassExtensions.RoundedDouble(cf.MonoisotopicMass.ToMz(1)));
            }
            else
            {
                Assert.AreEqual(mass, ClassExtensions.RoundedDouble(cf.MonoisotopicMass));
            }
        }
        [Test]
        [TestCase("C7 H15 N{15}1 ", 115.12476, true)]
        [TestCase("C6 C{13}1 H15 N{15}1 ", 116.12812, true)]
        [TestCase("C5 C{13}2 H15 N1 ", 116.13444, true)]
        [TestCase("C7 H13 H{2}2 N1 ", 116.14028, true)]
        [TestCase("C6 C{13}1 H15 N1 ", 115.13108, true)]
        [TestCase("C5 C{13}2 H15 N{15}1 ", 117.13147, true)]
        [TestCase("C7 H13 H{2}2 N{15}1 ", 117.13731, true)]
        [TestCase("C6 C{13}1 H13 H{2}2 N1 ", 117.14363, true)]
        [TestCase("C4 C{13}3 H15 N{15}1 ", 118.13483, true)]
        [TestCase("C6 C{13}1 H13 H{2}2 N{15}1 ", 118.14067, true)]
        [TestCase("C5 C{13}2 H13 H{2}2 N1 ", 118.14699, true)]
        [TestCase("C7 H11 H{2}4 N1 ", 118.15283, true)]

        [TestCase("C7 C{13}1 H15 N{15}1 O{18}1 ", 145.119998830, false)]
        [TestCase("C6 C{13}2 H15 N1 O{18}1 ", 145.126318771, false)]
        [TestCase("C7 C{13}1 H15 N{15}1 O{18}1 ", 145.119998830, false)]
        [TestCase("C6 C{13}2 H15 N1 O{18}1 ", 145.126318771, false)]
        [TestCase("C8 H13 H{2}2 N1 O{18}1 ", 145.132162593, false)]
        [TestCase("C5 C{13}3 H15 N{15}1 O1 ", 145.122463507, false)]
        [TestCase("C7 C{13}1 H13 H{2}2 N{15}1 O1 ", 145.128307329, false)]
        [TestCase("C6 C{13}2 H13 H{2}2 N1 O1 ", 145.134627269, false)]
        [TestCase("C5 C{13}3 H15 N{15}1 O1 ", 145.122463507, false)]
        [TestCase("C7 C{13}1 H13 H{2}2 N{15}1 O1 ", 145.128307329, false)]
        [TestCase("C6 C{13}2 H13 H{2}2 N1 O1 ", 145.134627269, false)]
        [TestCase("C8 H11 H{2}4 N1 O1 ", 145.140471091, false)]

        public static void TestChemicalFormulaWithIsotopes_DiLeu12plex(string formula, double mass, bool mz)
        {
            ChemicalFormula cf = ChemicalFormula.ParseFormula(formula);
            if (mz)
            {
                Assert.AreEqual(mass, ClassExtensions.RoundedDouble(cf.MonoisotopicMass.ToMz(1),5));
            }
            else
            {
                Assert.AreEqual(mass, ClassExtensions.RoundedDouble(cf.MonoisotopicMass));
            }
        }
    }
}