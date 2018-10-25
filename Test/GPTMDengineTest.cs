using Chemistry;
using EngineLayer;
using EngineLayer.Gptmd;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Proteomics;
using Proteomics.Fragmentation;
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
            var gptmdModifications = new List<Modification> { new Modification(_originalId: "21", _modificationType: "mt", _target: motifN, _locationRestriction: "Anywhere.", _monoisotopicMass: 21.981943) };
            IEnumerable<Tuple<double, double>> combos = new List<Tuple<double, double>>();
            Tolerance precursorMassTolerance = new PpmTolerance(10);

            allResultingIdentifications = new List<PeptideSpectralMatch>();
            var engine = new GptmdEngine(allResultingIdentifications, gptmdModifications, combos, new Dictionary<string, Tolerance> { { "filepath", precursorMassTolerance } }, new CommonParameters(), new List<string>());
            var res = (GptmdResults)engine.Run();
            Assert.AreEqual(0, res.Mods.Count);

            //PsmParent newPsm = new TestParentSpectrumMatch(588.22520189093 + 21.981943);
            Ms2ScanWithSpecificMass scan = new Ms2ScanWithSpecificMass(new MsDataScan(new MzSpectrum(new double[] { 1 }, new double[] { 1 }, false), 0, 1, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, null, null, "scan=1", double.NaN, null, null, double.NaN, null, DissociationType.AnyActivationType, 0, null), (588.22520189093 + 21.981943).ToMz(1), 1, "filepath", new CommonParameters());

            var parentProtein = new Protein("NNNNN", "accession");

            DigestionParams digestionParams = new DigestionParams(minPeptideLength: 5);
            List<Modification> variableModifications = new List<Modification>();
            var modPep = parentProtein.Digest(digestionParams, new List<Modification>(), variableModifications).First();

            var peptidesWithSetModifications = new List<PeptideWithSetModifications> { modPep };
            PeptideSpectralMatch newPsm = new PeptideSpectralMatch(peptidesWithSetModifications.First(), 0, 0, 0, scan, digestionParams, new List<MatchedFragmentIon>());
            
            Tolerance fragmentTolerance = new AbsoluteTolerance(0.01);

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
            var gptmdModifications = new List<Modification> { new Modification(_originalId: "21", _modificationType: "mt", _target: motifN, _locationRestriction: "Anywhere.", _monoisotopicMass: 21.981943),
                                                                      new Modification(_originalId: "16",  _modificationType: "mt", _target: motifP, _locationRestriction: "Anywhere.", _monoisotopicMass: 15.994915) };
            IEnumerable<Tuple<double, double>> combos = new List<Tuple<double, double>> { new Tuple<double, double>(21.981943, 15.994915) };
            Tolerance precursorMassTolerance = new PpmTolerance(10);

            MsDataScan dfd = new MsDataScan(new MzSpectrum(new double[] { 1 }, new double[] { 1 }, false), 0, 1, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, null, null, "scan=1", double.NaN, null, null, double.NaN, null, DissociationType.AnyActivationType, 0, null);
            Ms2ScanWithSpecificMass scan = new Ms2ScanWithSpecificMass(dfd, (651.297638557 + 21.981943 + 15.994915).ToMz(1), 1, "filepath", new CommonParameters());

            var parentProtein = new Protein("NNNPPP", "accession");
            DigestionParams digestionParams = new DigestionParams(minPeptideLength: 5);
            List<Modification> variableModifications = new List<Modification>();
            var modPep = parentProtein.Digest(digestionParams, new List<Modification>(), variableModifications).First();

            var peptidesWithSetModifications = new List<PeptideWithSetModifications> { modPep };
            PeptideSpectralMatch match = new PeptideSpectralMatch(peptidesWithSetModifications.First(), 0, 0, 0, scan, digestionParams, new List<MatchedFragmentIon>());
            PeptideSpectralMatch newPsm = new PeptideSpectralMatch(peptidesWithSetModifications.First(), 0, 0, 0, scan, digestionParams, new List<MatchedFragmentIon>());
            
            Tolerance fragmentTolerance = new AbsoluteTolerance(0.01);

            match.SetFdrValues(1, 0, 0, 1, 0, 0, 0, 0, 0, false);
            allIdentifications = new List<PeptideSpectralMatch> { match };

            var engine = new GptmdEngine(allIdentifications, gptmdModifications, combos, new Dictionary<string, Tolerance> { { "filepath", precursorMassTolerance } }, new CommonParameters(), new List<string>());
            var res = (GptmdResults)engine.Run();
            Assert.AreEqual(1, res.Mods.Count);
            Assert.AreEqual(6, res.Mods["accession"].Count);
            Assert.AreEqual(3, res.Mods["accession"].Where(b => b.Item2.OriginalId.Equals("21")).Count());
            Assert.AreEqual(3, res.Mods["accession"].Where(b => b.Item2.OriginalId.Equals("16")).Count());
        }

        [Test]
        [TestCase("P", "PETID", "junk", 1, 5, 1, false)]
        [TestCase("P", "PETID", "Unassigned.", 1, 5, 1, false)]
        [TestCase("P", "PETID", "Anywhere.", 1, 5, 1, true)]
        [TestCase("P", "PETID", "N-terminal.", 1, 5, 1, true)]
        [TestCase("P", "PETID", "Peptide N-terminal.", 1, 5, 1, true)]
        [TestCase("P", "PETID", "C-terminal.", 1, 5, 1, false)]
        [TestCase("P", "PETID", "Peptide C-terminal.", 1, 5, 1, false)]
        [TestCase("E", "PETID", "Anywhere.", 2, 5, 2, true)]
        [TestCase("E", "PETID", "N-terminal.", 2, 5, 2, true)]
        [TestCase("E", "PETID", "Peptide N-terminal.", 2, 5, 2, false)]
        [TestCase("E", "PETID", "C-terminal.", 2, 5, 2, false)]
        [TestCase("E", "PETID", "Peptide C-terminal.", 2, 5, 2, false)]
        [TestCase("D", "PETID", "Anywhere.", 5, 5, 5, true)]
        [TestCase("D", "PETID", "N-terminal.", 5, 5, 5, false)]
        [TestCase("D", "PETID", "Peptide N-terminal.", 5, 5, 5, false)]
        [TestCase("D", "PETID", "C-terminal.", 5, 5, 5, true)]
        [TestCase("D", "PETID", "Peptide C-terminal.", 5, 5, 5, true)]
        public static void Test_GptmdEngineModFits(string targetAminoAcid, string proteinSequence, string locationRestriction, int peptideOneBasedIndex, int peptideLength, int proteinOneBasedIndex, bool result)
        {
            ModificationMotif.TryGetMotif(targetAminoAcid, out ModificationMotif motif);
            Modification attemptToLocalize = new Modification(null, null, null, null, _target: motif, _locationRestriction: locationRestriction, _chemicalFormula: null, _monoisotopicMass: 1, _databaseReference: null, _taxonomicRange: null, _keywords: null, _neutralLosses: null, _diagnosticIons: null, _fileOrigin: null);
            Dictionary<int, List<Modification>> oneBasedModifications = new Dictionary<int, List<Modification>>();
            oneBasedModifications.Add(proteinOneBasedIndex, new List<Modification>() { attemptToLocalize });
            Protein protein = new Protein(proteinSequence, null, null, null, oneBasedModifications, null, null, null, false, false, null, null, null, "");

            Assert.AreEqual(result, GptmdEngine.ModFits(attemptToLocalize, protein, peptideOneBasedIndex, peptideLength, proteinOneBasedIndex));
        }
    }
}