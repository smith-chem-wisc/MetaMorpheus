using Chemistry;
using EngineLayer;
using EngineLayer.Gptmd;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Linq;

namespace Test
{
    [TestFixture]
    public static class GptmdEngineTest
    {
        [Test]
        public static void TestGptmdEngine()
        {
            List<PeptideSpectralMatch> allResultingIdentifications = null;
            ModificationMotif.TryGetMotif("N", out ModificationMotif motifN);
            var gptmdModifications = new List<ModificationWithMass> { new ModificationWithMass("21", "mt", motifN, TerminusLocalization.Any, 21.981943) };
            IEnumerable<Tuple<double, double>> combos = new List<Tuple<double, double>>();
            Tolerance precursorMassTolerance = new PpmTolerance(10);

            allResultingIdentifications = new List<PeptideSpectralMatch>();
            var engine = new GptmdEngine(allResultingIdentifications, gptmdModifications, combos, new Dictionary<string, Tolerance> { { "filepath", precursorMassTolerance } }, new CommonParameters(), new List<string>());
            var res = (GptmdResults)engine.Run();
            Assert.AreEqual(0, res.Mods.Count);

            //PsmParent newPsm = new TestParentSpectrumMatch(588.22520189093 + 21.981943);
            Ms2ScanWithSpecificMass scan = new Ms2ScanWithSpecificMass(new MsDataScan(new MzSpectrum(new double[] { 1 }, new double[] { 1 }, false), 0, 1, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, null, null, "scan=1", double.NaN, null, null, double.NaN, null, DissociationType.AnyActivationType, 0, null), (588.22520189093 + 21.981943).ToMz(1), 1, "filepath");

            var parentProtein = new Protein("NNNNN", "accession");

            DigestionParams digestionParams = new DigestionParams(minPeptideLength: 5);
            List<ModificationWithMass> variableModifications = new List<ModificationWithMass>();
            var modPep = parentProtein.Digest(digestionParams, new List<ModificationWithMass>(), variableModifications).First();

            var peptidesWithSetModifications = new List<PeptideWithSetModifications> { modPep };
            PeptideSpectralMatch newPsm = new PeptideSpectralMatch(peptidesWithSetModifications.First().CompactPeptide(TerminusType.None), 0, 0, 0, scan, digestionParams);

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

            engine = new GptmdEngine(allResultingIdentifications, gptmdModifications, combos, new Dictionary<string, Tolerance> { { "filepath", precursorMassTolerance } }, new CommonParameters(), new List<string>());
            res = (GptmdResults)engine.Run();
            Assert.AreEqual(1, res.Mods.Count);
            Assert.AreEqual(5, res.Mods["accession"].Count);
        }

        [Test]
        public static void TestCombos()
        {
            List<PeptideSpectralMatch> allIdentifications = null;
            ModificationMotif.TryGetMotif("N", out ModificationMotif motifN);
            ModificationMotif.TryGetMotif("P", out ModificationMotif motifP);
            var gptmdModifications = new List<ModificationWithMass> { new ModificationWithMass("21", "mt", motifN, TerminusLocalization.Any, 21.981943,null),
                                                                      new ModificationWithMass("16",  "mt", motifP, TerminusLocalization.Any, 15.994915,null) };
            IEnumerable<Tuple<double, double>> combos = new List<Tuple<double, double>> { new Tuple<double, double>(21.981943, 15.994915) };
            Tolerance precursorMassTolerance = new PpmTolerance(10);

            MsDataScan dfd = new MsDataScan(new MzSpectrum(new double[] { 1 }, new double[] { 1 }, false), 0, 1, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, null, null, "scan=1", double.NaN, null, null, double.NaN, null, DissociationType.AnyActivationType, 0, null);
            Ms2ScanWithSpecificMass scan = new Ms2ScanWithSpecificMass(dfd, (651.297638557 + 21.981943 + 15.994915).ToMz(1), 1, "filepath");

            var parentProtein = new Protein("NNNPPP", "accession");
            DigestionParams digestionParams = new DigestionParams(minPeptideLength: 5);
            List<ModificationWithMass> variableModifications = new List<ModificationWithMass>();
            var modPep = parentProtein.Digest(digestionParams, new List<ModificationWithMass>(), variableModifications).First();

            var peptidesWithSetModifications = new List<PeptideWithSetModifications> { modPep };
            PeptideSpectralMatch match = new PeptideSpectralMatch(peptidesWithSetModifications.First().CompactPeptide(TerminusType.None), 0, 0, 0, scan, digestionParams);
            PeptideSpectralMatch newPsm = new PeptideSpectralMatch(peptidesWithSetModifications.First().CompactPeptide(TerminusType.None), 0, 0, 0, scan, digestionParams);
            Dictionary<ModificationWithMass, ushort> modsDictionary = new Dictionary<ModificationWithMass, ushort>();
            Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>> matching = new Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>>
            {
                {peptidesWithSetModifications.First().CompactPeptide(TerminusType.None), new HashSet<PeptideWithSetModifications>{ peptidesWithSetModifications.First() } }
            };

            List<ProductType> lp = new List<ProductType> { ProductType.B, ProductType.Y };

            Tolerance fragmentTolerance = new AbsoluteTolerance(0.01);

            match.MatchToProteinLinkedPeptides(matching);

            match.SetFdrValues(1, 0, 0, 1, 0, 0, 0, 0, 0, false);
            allIdentifications = new List<PeptideSpectralMatch> { match };

            var engine = new GptmdEngine(allIdentifications, gptmdModifications, combos, new Dictionary<string, Tolerance> { { "filepath", precursorMassTolerance } }, new CommonParameters(), new List<string>());
            var res = (GptmdResults)engine.Run();
            Assert.AreEqual(1, res.Mods.Count);
            Assert.AreEqual(6, res.Mods["accession"].Count);
            Assert.AreEqual(3, res.Mods["accession"].Where(b => b.Item2.id.Equals("21")).Count());
            Assert.AreEqual(3, res.Mods["accession"].Where(b => b.Item2.id.Equals("16")).Count());
        }
    }
}