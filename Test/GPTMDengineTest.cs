using Chemistry;
using EngineLayer;
using EngineLayer.Gptmd;
using IO.MzML;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Proteomics;
using System;
using System.Collections.Generic;
using System.Linq;
using TaskLayer;

namespace Test
{
    [TestFixture]
    public static class GptmdEngineTest
    {
        #region Public Methods

        [Test]
        public static void TestGptmdEngine()
        {
            List<Psm> allResultingIdentifications = null;
            ModificationMotif.TryGetMotif("N", out ModificationMotif motifN);
            var gptmdModifications = new List<ModificationWithMass> { new ModificationWithMass("21", "mt", motifN, TerminusLocalization.Any, 21.981943) };
            IEnumerable<Tuple<double, double>> combos = new List<Tuple<double, double>>();
            Tolerance precursorMassTolerance = new PpmTolerance(10);

            allResultingIdentifications = new List<Psm>();
            var engine = new GptmdEngine(allResultingIdentifications, gptmdModifications, combos, precursorMassTolerance, new List<string>());
            var res = (GptmdResults)engine.Run();
            Assert.AreEqual(0, res.Mods.Count);

            //PsmParent newPsm = new TestParentSpectrumMatch(588.22520189093 + 21.981943);
            Ms2ScanWithSpecificMass scan = new Ms2ScanWithSpecificMass(new MzmlScanWithPrecursor(0, new MzmlMzSpectrum(new double[] { 1 }, new double[] { 1 }, false), 1, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, double.NaN, null, null, double.NaN, null, DissociationType.AnyActivationType, 0, null, null, "scan=1"), (588.22520189093 + 21.981943).ToMz(1), 1, null);

            var parentProtein = new Protein("NNNNN", "accession");
            var protease = new Protease("Custom Protease", new List<string> { "K" }, new List<string>(), TerminusType.C, CleavageSpecificity.Full, null, null, null);

            DigestionParams digestionParams = new DigestionParams();
            List<ModificationWithMass> variableModifications = new List<ModificationWithMass>();
            var modPep = parentProtein.Digest(digestionParams, new List<ModificationWithMass>(), variableModifications).First();

            var peptidesWithSetModifications = new List<PeptideWithSetModifications> { modPep };
            Psm newPsm = new Psm(peptidesWithSetModifications.First().CompactPeptide(TerminusType.None), 0, 0, 0, scan);

            Dictionary<ModificationWithMass, ushort> modsDictionary = new Dictionary<ModificationWithMass, ushort>();
            Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>> matching = new Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>>
            {
                {peptidesWithSetModifications.First().CompactPeptide(TerminusType.None), new HashSet<PeptideWithSetModifications>{ peptidesWithSetModifications.First() } }
            };
            List<ProductType> lp = new List<ProductType> { ProductType.B, ProductType.Y };
            Tolerance fragmentTolerance = new AbsoluteTolerance(0.01);
            newPsm.MatchToProteinLinkedPeptides(matching);

            newPsm.SetFdrValues(1, 0, 0, 1, 0, 0, 0, 0, 0, false);
            allResultingIdentifications.Add(newPsm);

            engine = new GptmdEngine(allResultingIdentifications, gptmdModifications, combos, precursorMassTolerance, new List<string>());
            res = (GptmdResults)engine.Run();
            Assert.AreEqual(1, res.Mods.Count);
            Assert.AreEqual(5, res.Mods["accession"].Count);
        }

        [Test]
        public static void TestCombos()
        {
            List<Psm> allIdentifications = null;
            ModificationMotif.TryGetMotif("N", out ModificationMotif motifN);
            ModificationMotif.TryGetMotif("P", out ModificationMotif motifP);
            var gptmdModifications = new List<ModificationWithMass> { new ModificationWithMass("21", "mt", motifN, TerminusLocalization.Any, 21.981943,null),
                                                                      new ModificationWithMass("16",  "mt", motifP, TerminusLocalization.Any, 15.994915,null) };
            IEnumerable<Tuple<double, double>> combos = new List<Tuple<double, double>> { new Tuple<double, double>(21.981943, 15.994915) };
            Tolerance precursorMassTolerance = new PpmTolerance(10);
            var protease = new Protease("Custom Protease", new List<string> { "K" }, new List<string>(), TerminusType.C, CleavageSpecificity.Full, null, null, null);

            IMsDataScanWithPrecursor<IMzSpectrum<IMzPeak>> dfd = new MzmlScanWithPrecursor(0, new MzmlMzSpectrum(new double[] { 1 }, new double[] { 1 }, false), 1, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, double.NaN, null, null, double.NaN, null, DissociationType.AnyActivationType, 0, null, null, "scan=1");
            Ms2ScanWithSpecificMass scan = new Ms2ScanWithSpecificMass(dfd, (651.297638557 + 21.981943 + 15.994915).ToMz(1), 1, null);

            var parentProtein = new Protein("NNNPPP", "accession");
            DigestionParams digestionParams = new DigestionParams();
            List<ModificationWithMass> variableModifications = new List<ModificationWithMass>();
            var modPep = parentProtein.Digest(digestionParams, new List<ModificationWithMass>(), variableModifications).First();

            var peptidesWithSetModifications = new List<PeptideWithSetModifications> { modPep };
            Psm match = new Psm(peptidesWithSetModifications.First().CompactPeptide(TerminusType.None), 0, 0, 0, scan);
            Psm newPsm = new Psm(peptidesWithSetModifications.First().CompactPeptide(TerminusType.None), 0, 0, 0, scan);
            Dictionary<ModificationWithMass, ushort> modsDictionary = new Dictionary<ModificationWithMass, ushort>();
            Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>> matching = new Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>>
            {
                {peptidesWithSetModifications.First().CompactPeptide(TerminusType.None), new HashSet<PeptideWithSetModifications>{ peptidesWithSetModifications.First() } }
            };
            List<ProductType> lp = new List<ProductType> { ProductType.B, ProductType.Y };

            Tolerance fragmentTolerance = new AbsoluteTolerance(0.01);

            match.MatchToProteinLinkedPeptides(matching);

            match.SetFdrValues(1, 0, 0, 1, 0, 0, 0, 0, 0, false);
            allIdentifications = new List<Psm> { match };

            var engine = new GptmdEngine(allIdentifications, gptmdModifications, combos, precursorMassTolerance, new List<string>());
            var res = (GptmdResults)engine.Run();
            Assert.AreEqual(1, res.Mods.Count);
            Assert.AreEqual(6, res.Mods["accession"].Count);
            Assert.AreEqual(3, res.Mods["accession"].Where(b => b.Item2.id.Equals("21")).Count());
            Assert.AreEqual(3, res.Mods["accession"].Where(b => b.Item2.id.Equals("16")).Count());
        }

        #endregion Public Methods
    }
}