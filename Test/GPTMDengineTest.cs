using Chemistry;
using EngineLayer;
using EngineLayer.ClassicSearch;
using EngineLayer.Gptmd;
using IO.MzML;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Proteomics;
using System;
using System.Collections.Generic;
using System.Linq;

namespace Test
{
    [TestFixture]
    public class GptmdEngineTest
    {

        #region Public Methods

        [Test]
        public static void TestGptmdEngine()
        {
            List<NewPsmWithFdr> allResultingIdentifications = null;
            ModificationMotif motifN;
            ModificationMotif.TryGetMotif("N", out motifN);
            var gptmdModifications = new List<ModificationWithMass> { new ModificationWithMass("21", null, motifN, ModificationSites.Any, 21.981943, null, new List<double> { 0 }, new List<double> { 21.981943 }, null) };
            IEnumerable<Tuple<double, double>> combos = new List<Tuple<double, double>>();
            Tolerance precursorMassTolerance = new Tolerance(ToleranceUnit.PPM, 10);

            allResultingIdentifications = new List<NewPsmWithFdr>();
            var engine = new GptmdEngine(allResultingIdentifications, gptmdModifications, combos, precursorMassTolerance);
            var res = (GptmdResults)engine.Run();
            Assert.AreEqual(0, res.Mods.Count);

            //PsmParent newPsm = new TestParentSpectrumMatch(588.22520189093 + 21.981943);
            Ms2ScanWithSpecificMass scan = new Ms2ScanWithSpecificMass(new MzmlScanWithPrecursor(0, new MzmlMzSpectrum(new double[] { 1 }, new double[] { 1 }, false), 1, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, double.NaN, null, null, double.NaN, null, DissociationType.AnyActivationType, 0, null, null), new MzPeak((588.22520189093 + 21.981943).ToMz(1), 1), 1, null);

            var parentProtein = new Protein("NNNNN", "accession", null, new Dictionary<int, List<Modification>>(), null, null, null, null, null, false, false, null);
            var protease = new Protease("Custom Protease", new List<string> { "K" }, new List<string>(), TerminusType.C, CleavageSpecificity.Full, null, null, null);

            var modPep = parentProtein.Digest(protease, 0, null, null, InitiatorMethionineBehavior.Variable, new List<ModificationWithMass>()).First();
            //var twoBasedVariableAndLocalizeableModificationss = new Dictionary<int, MorpheusModification>();
            List<ModificationWithMass> variableModifications = new List<ModificationWithMass>();
            var peptidesWithSetModifications = new List<PeptideWithSetModifications> { modPep.GetPeptidesWithSetModifications(variableModifications, 4096, 3).First() };
            PsmParent newPsm = new PsmClassic(peptidesWithSetModifications.First(), 0, 0, 0, scan);

            Dictionary<ModificationWithMass, ushort> modsDictionary = new Dictionary<ModificationWithMass, ushort>();
            Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>> matching = new Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>>
            {
                {newPsm.GetCompactPeptide(modsDictionary), new HashSet<PeptideWithSetModifications>{ peptidesWithSetModifications.First() } }
            };
            List<ProductType> lp = new List<ProductType> { ProductType.B, ProductType.Y };
            Tolerance fragmentTolerance = new Tolerance(ToleranceUnit.Absolute, 0.01);
            newPsm.ComputeProteinLevelInfo(matching, fragmentTolerance, scan, lp, modsDictionary);

            NewPsmWithFdr thePsmwithfdr = new NewPsmWithFdr(newPsm);
            thePsmwithfdr.SetValues(1, 0, 0, 1, 0, 0);
            allResultingIdentifications.Add(thePsmwithfdr);

            engine = new GptmdEngine(allResultingIdentifications, gptmdModifications, combos, precursorMassTolerance);
            res = (GptmdResults)engine.Run();
            Assert.AreEqual(1, res.Mods.Count);
            Assert.AreEqual(5, res.Mods["accession"].Count);
        }

        [Test]
        public static void TestCombos()
        {
            List<NewPsmWithFdr> allIdentifications = null;
            ModificationMotif motifN;
            ModificationMotif.TryGetMotif("N", out motifN);
            ModificationMotif motifP;
            ModificationMotif.TryGetMotif("P", out motifP);
            var gptmdModifications = new List<ModificationWithMass> { new ModificationWithMass("21", null, motifN, ModificationSites.Any, 21.981943,null, new List<double> { 0 }, new List<double> { 21.981943 },  null),
                                                                      new ModificationWithMass("16", null, motifP, ModificationSites.Any, 15.994915,null, new List<double> { 0 }, new List<double> { 15.994915 },  null) };
            IEnumerable<Tuple<double, double>> combos = new List<Tuple<double, double>> { new Tuple<double, double>(21.981943, 15.994915) };
            Tolerance precursorMassTolerance = new Tolerance(ToleranceUnit.PPM, 10);
            var protease = new Protease("Custom Protease", new List<string> { "K" }, new List<string>(), TerminusType.C, CleavageSpecificity.Full, null, null, null);

            IMsDataScanWithPrecursor<IMzSpectrum<IMzPeak>> dfd = new MzmlScanWithPrecursor(0, new MzmlMzSpectrum(new double[] { 1 }, new double[] { 1 }, false), 1, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, double.NaN, null, null, double.NaN, null, DissociationType.AnyActivationType, 0, null, null);
            Ms2ScanWithSpecificMass scan = new Ms2ScanWithSpecificMass(dfd, new MzPeak((651.297638557 + 21.981943 + 15.994915).ToMz(1), 1), 1, null);

            var parentProtein = new Protein("NNNPPP", "accession", null, new Dictionary<int, List<Modification>>(), null, null, null, null, null, false, false, null);
            var modPep = parentProtein.Digest(protease, 0, null, null, InitiatorMethionineBehavior.Variable, new List<ModificationWithMass>()).First();

            List<ModificationWithMass> variableModifications = new List<ModificationWithMass>();
            var peptidesWithSetModifications = new List<PeptideWithSetModifications> { modPep.GetPeptidesWithSetModifications(variableModifications, 4096, 3).First() };
            PsmParent match = new PsmClassic(peptidesWithSetModifications.First(), 0, 0, 0, scan);
            PsmParent newPsm = new PsmClassic(peptidesWithSetModifications.First(), 0, 0, 0, scan);
            Dictionary<ModificationWithMass, ushort> modsDictionary = new Dictionary<ModificationWithMass, ushort>();
            Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>> matching = new Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>>
            {
                {newPsm.GetCompactPeptide(modsDictionary), new HashSet<PeptideWithSetModifications>{ peptidesWithSetModifications.First() } }
            };
            List<ProductType> lp = new List<ProductType> { ProductType.B, ProductType.Y };
            Tolerance fragmentTolerance = new Tolerance(ToleranceUnit.Absolute, 0.01);
            match.ComputeProteinLevelInfo(matching, fragmentTolerance, scan, lp, modsDictionary);

            NewPsmWithFdr thePsmwithfdr = new NewPsmWithFdr(match);
            thePsmwithfdr.SetValues(1, 0, 0, 1, 0, 0);
            allIdentifications = new List<NewPsmWithFdr> { thePsmwithfdr };

            var engine = new GptmdEngine(allIdentifications, gptmdModifications, combos, precursorMassTolerance);
            var res = (GptmdResults)engine.Run();
            Assert.AreEqual(1, res.Mods.Count);
            Assert.AreEqual(6, res.Mods["accession"].Count);
            Assert.AreEqual(3, res.Mods["accession"].Where(b => b.Item2.id.Equals("21")).Count());
            Assert.AreEqual(3, res.Mods["accession"].Where(b => b.Item2.id.Equals("16")).Count());
        }

        #endregion Public Methods

    }
}