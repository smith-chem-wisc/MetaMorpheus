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
        [TestCase("PEPTIDE", 1)]
        //[TestCase("PEPTIDEK", 1)]
        //[TestCase("PK", 1)]
        //[TestCase("KP", 1)]
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

            Assert.AreEqual(1029.5302, ClassExtensions.RoundedDouble(p.MonoisotopicMass.ToMz(1),4));

            Assert.AreEqual(11, f.Where(z => z.ProductType == ProductType.D).ToList().Count());
            Assert.AreEqual(11, f.Where(z => z.ProductType == ProductType.M).ToList().Count());
        }

        [Test]
        [TestCase("C5 N2 H13 C{13}1", 114.111228263)]
        [TestCase("C4 N2 H13 C{13}2", 115.114583098)]
        [TestCase("C4 N1 H13 C{13}2 N{15}1", 116.111617992)]
        [TestCase("C3 N1 H13 C{13}3 N{15}1", 117.114972828)]
        [TestCase("C5 N2 H13 C{13}2 O{18}1", 145.113742711)]
        [TestCase("C4 N1O1 H13 C{13}3 N{15}1", 145.109887447)]
        [TestCase("C7 N3O3 H25 C{13}7 N{15}1", 305.213184422)]
        [TestCase("C8 N2O3 H25 C{13}6 N{15}2", 305.206864482)]
        public static void TestChemicalFormulaWithIsotopes_iTRAQ(string formula, double mass)
        {
            ChemicalFormula cf = ChemicalFormula.ParseFormula(formula);
            Assert.AreEqual(mass, ClassExtensions.RoundedDouble(cf.MonoisotopicMass));
        }
    }
}