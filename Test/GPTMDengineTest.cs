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
        [TestCase("NNNNN", "accession", @"not applied", 5)]
        [TestCase("NNNNN", "accession", @"1\t50000000\t.\tA\tG\t.\tPASS\tANN=\tGT:AD:DP\t1/1:30,30:30", 4)]
        public static void TestGptmdEngine(string proteinSequence, string accession, string sequenceVariantDescription, int numModifiedResidues)
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

            var parentProtein = new Protein(proteinSequence, accession, sequenceVariations: new List<SequenceVariation> { new SequenceVariation(1, "N", "A", sequenceVariantDescription) });
            var variantProteins = parentProtein.GetVariantProteins();

            DigestionParams digestionParams = new DigestionParams(minPeptideLength: 5);
            List<Modification> variableModifications = new List<Modification>();
            var modPep = variantProteins.SelectMany(p => p.Digest(digestionParams, new List<Modification>(), variableModifications)).First();

            //PsmParent newPsm = new TestParentSpectrumMatch(588.22520189093 + 21.981943);
            Ms2ScanWithSpecificMass scan = new Ms2ScanWithSpecificMass(new MsDataScan(new MzSpectrum(new double[] { 1 }, new double[] { 1 }, false), 0, 1, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, null, null, "scan=1", double.NaN, null, null, double.NaN, null, DissociationType.AnyActivationType, 0, null), (new Proteomics.AminoAcidPolymer.Peptide(modPep.BaseSequence).MonoisotopicMass + 21.981943).ToMz(1), 1, "filepath");

            var peptidesWithSetModifications = new List<PeptideWithSetModifications> { modPep };
            PeptideSpectralMatch newPsm = new PeptideSpectralMatch(peptidesWithSetModifications.First(), 0, 0, 0, scan, digestionParams, new List<MatchedFragmentIon>());
            
            Tolerance fragmentTolerance = new AbsoluteTolerance(0.01);

            newPsm.SetFdrValues(1, 0, 0, 1, 0, 0, 0, 0, 0, false);
            allResultingIdentifications.Add(newPsm);

            engine = new GptmdEngine(allResultingIdentifications, gptmdModifications, combos, new Dictionary<string, Tolerance> { { "filepath", precursorMassTolerance } }, new CommonParameters(), new List<string>());
            res = (GptmdResults)engine.Run();
            Assert.AreEqual(1, res.Mods.Count);
            Assert.AreEqual(numModifiedResidues, res.Mods["accession"].Count);
        }

        [Test]
        [TestCase("NNNPPP", "accession", @"not applied", 6, 3, 3)]
        [TestCase("NNNPPP", "accession", @"1\t50000000\t.\tA\tG\t.\tPASS\tANN=\tGT:AD:DP\t1/1:30,30:30", 5, 2, 3)]
        public static void TestCombos(string proteinSequence, string accession, string sequenceVariantDescription, int numModifiedResidues, int numModifiedResiduesN, int numModifiedResiduesP)
        {
            List<PeptideSpectralMatch> allIdentifications = null;
            ModificationMotif.TryGetMotif("N", out ModificationMotif motifN);
            ModificationMotif.TryGetMotif("P", out ModificationMotif motifP);
            var gptmdModifications = new List<Modification> { new Modification(_originalId: "21", _modificationType: "mt", _target: motifN, _locationRestriction: "Anywhere.", _monoisotopicMass: 21.981943),
                                                                      new Modification(_originalId: "16",  _modificationType: "mt", _target: motifP, _locationRestriction: "Anywhere.", _monoisotopicMass: 15.994915) };
            IEnumerable<Tuple<double, double>> combos = new List<Tuple<double, double>> { new Tuple<double, double>(21.981943, 15.994915) };
            Tolerance precursorMassTolerance = new PpmTolerance(10);

            var parentProtein = new Protein(proteinSequence, accession, sequenceVariations: new List<SequenceVariation> { new SequenceVariation(1, "N", "A", sequenceVariantDescription) });
            var variantProteins = parentProtein.GetVariantProteins();

            DigestionParams digestionParams = new DigestionParams(minPeptideLength: 5);
            List<Modification> variableModifications = new List<Modification>();
            var modPep = variantProteins.SelectMany(p => p.Digest(digestionParams, new List<Modification>(), variableModifications)).First();

            MsDataScan dfd = new MsDataScan(new MzSpectrum(new double[] { 1 }, new double[] { 1 }, false), 0, 1, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, null, null, "scan=1", double.NaN, null, null, double.NaN, null, DissociationType.AnyActivationType, 0, null);
            Ms2ScanWithSpecificMass scan = new Ms2ScanWithSpecificMass(dfd, (new Proteomics.AminoAcidPolymer.Peptide(modPep.BaseSequence).MonoisotopicMass + 21.981943 + 15.994915).ToMz(1), 1, "filepath");

            var peptidesWithSetModifications = new List<PeptideWithSetModifications> { modPep };
            PeptideSpectralMatch match = new PeptideSpectralMatch(peptidesWithSetModifications.First(), 0, 0, 0, scan, digestionParams, new List<MatchedFragmentIon>());
            PeptideSpectralMatch newPsm = new PeptideSpectralMatch(peptidesWithSetModifications.First(), 0, 0, 0, scan, digestionParams, new List<MatchedFragmentIon>());
            
            Tolerance fragmentTolerance = new AbsoluteTolerance(0.01);

            match.SetFdrValues(1, 0, 0, 1, 0, 0, 0, 0, 0, false);
            allIdentifications = new List<PeptideSpectralMatch> { match };

            var engine = new GptmdEngine(allIdentifications, gptmdModifications, combos, new Dictionary<string, Tolerance> { { "filepath", precursorMassTolerance } }, new CommonParameters(), new List<string>());
            var res = (GptmdResults)engine.Run();
            Assert.AreEqual(1, res.Mods.Count);
            Assert.AreEqual(numModifiedResidues, res.Mods["accession"].Count);
            Assert.AreEqual(numModifiedResiduesN, res.Mods["accession"].Where(b => b.Item2.OriginalId.Equals("21")).Count());
            Assert.AreEqual(numModifiedResiduesP, res.Mods["accession"].Where(b => b.Item2.OriginalId.Equals("16")).Count());
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