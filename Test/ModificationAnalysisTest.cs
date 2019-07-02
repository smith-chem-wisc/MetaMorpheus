using EngineLayer;
using EngineLayer.FdrAnalysis;
using EngineLayer.ModificationAnalysis;
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
    public static class ModificationAnalysisTest
    {
        [Test]
        public static void TestModificationAnalysis()
        {
            Ms2ScanWithSpecificMass scan = new Ms2ScanWithSpecificMass(new MsDataScan(new MzSpectrum(new double[,] { }), 0, 0, true, Polarity.Positive,
                0, new MzLibUtil.MzRange(0, 0), "", MZAnalyzerType.FTICR, 0, null, null, ""), 0, 0, "", new CommonParameters());

            ModificationMotif.TryGetMotif("N", out ModificationMotif motif1);
            Modification mod1 = new Modification(_originalId: "mod1", _modificationType: "myModType", _target: motif1, _locationRestriction: "Anywhere.", _monoisotopicMass: 10);

            ModificationMotif.TryGetMotif("L", out ModificationMotif motif2);
            Modification mod2 = new Modification(_originalId: "mod2", _modificationType: "myModType", _target: motif2, _locationRestriction: "Anywhere.", _monoisotopicMass: 10);

            IDictionary<int, List<Modification>> oneBasedModifications = new Dictionary<int, List<Modification>>
            {
                {2, new List<Modification>{ mod1 }},
                {5, new List<Modification>{ mod2 }},
                {7, new List<Modification>{ mod1 }},
            };
            Protein protein1 = new Protein("MNLDLDNDL", "prot1", oneBasedModifications: oneBasedModifications);

            Dictionary<int, Modification> allModsOneIsNterminus1 = new Dictionary<int, Modification>
            {
                {2, mod1},
            };
            PeptideWithSetModifications pwsm1 = new PeptideWithSetModifications(protein1, new DigestionParams(), 2, 9, CleavageSpecificity.Unknown, null, 0, allModsOneIsNterminus1, 0);

            Dictionary<int, Modification> allModsOneIsNterminus2 = new Dictionary<int, Modification>
            {
                {2, mod1},
                {7, mod1},
            };
            PeptideWithSetModifications pwsm2 = new PeptideWithSetModifications(protein1, new DigestionParams(), 2, 9, CleavageSpecificity.Unknown, null, 0, allModsOneIsNterminus2, 0);

            Dictionary<int, Modification> allModsOneIsNterminus3 = new Dictionary<int, Modification>
            {
                {7, mod1},
            };
            PeptideWithSetModifications pwsm3 = new PeptideWithSetModifications(protein1, new DigestionParams(), 2, 9, CleavageSpecificity.Unknown, null, 0, allModsOneIsNterminus3, 0);

            Dictionary<int, Modification> allModsOneIsNterminus4 = new Dictionary<int, Modification>
            {
                {8, mod1},
            };
            PeptideWithSetModifications pwsm4 = new PeptideWithSetModifications(protein1, new DigestionParams(), 1, 9, CleavageSpecificity.Unknown, null, 0, allModsOneIsNterminus4, 0);

            CommonParameters CommonParameters = new CommonParameters(
                digestionParams: new DigestionParams(
                    maxMissedCleavages: 0,
                    minPeptideLength: 1,
                    maxModificationIsoforms: int.MaxValue),
                scoreCutoff: 1);

            var newPsms = new List<PeptideSpectralMatch>
            {
                new PeptideSpectralMatch(pwsm1, 0, 10, 0, scan, CommonParameters.DigestionParams, new List<MatchedFragmentIon>()),
                new PeptideSpectralMatch(pwsm1, 0, 10, 0, scan, CommonParameters.DigestionParams, new List<MatchedFragmentIon>()),
                new PeptideSpectralMatch(pwsm2, 0, 10, 0, scan, CommonParameters.DigestionParams, new List<MatchedFragmentIon>()),
                new PeptideSpectralMatch(pwsm3, 0, 10, 0, scan, CommonParameters.DigestionParams, new List<MatchedFragmentIon>()),
                new PeptideSpectralMatch(pwsm4, 0, 10, 0, scan, CommonParameters.DigestionParams, new List<MatchedFragmentIon>()),
            };

            foreach (var psm in newPsms)
            {
                psm.ResolveAllAmbiguities();
            }

            MassDiffAcceptor searchMode = new SinglePpmAroundZeroSearchMode(5);
            List<Protein> proteinList = new List<Protein> { protein1 };

            FdrAnalysisEngine fdrAnalysisEngine = new FdrAnalysisEngine(newPsms, searchMode.NumNotches, CommonParameters, new List<string>());
            fdrAnalysisEngine.Run();
            ModificationAnalysisEngine modificationAnalysisEngine = new ModificationAnalysisEngine(newPsms, new CommonParameters(), new List<string>());
            var res = (ModificationAnalysisResults)modificationAnalysisEngine.Run();

            Assert.AreEqual(2, res.CountOfEachModSeenOnProteins.Count());
            Assert.AreEqual(2, res.CountOfEachModSeenOnProteins[mod1.IdWithMotif]);
            Assert.AreEqual(1, res.CountOfEachModSeenOnProteins[mod2.IdWithMotif]);

            Assert.AreEqual(1, res.CountOfModsSeenAndLocalized.Count());
            Assert.AreEqual(2, res.CountOfModsSeenAndLocalized[mod1.IdWithMotif]);

            Assert.AreEqual(0, res.CountOfAmbiguousButLocalizedModsSeen.Count());

            Assert.AreEqual(0, res.CountOfUnlocalizedMods.Count());

            Assert.AreEqual(0, res.CountOfUnlocalizedFormulas.Count());
        }

        [Test]
        public static void TestModificationAnalysisWithNonLocalizedPtms()
        {
            Ms2ScanWithSpecificMass scan = new Ms2ScanWithSpecificMass(new MsDataScan(new MzSpectrum(new double[,] { }), 0, 0, true, Polarity.Positive,
                0, new MzLibUtil.MzRange(0, 0), "", MZAnalyzerType.FTICR, 0, null, null, ""), 0, 0, "", new CommonParameters());

            ModificationMotif.TryGetMotif("N", out ModificationMotif motif1);
            Modification mod1 = new Modification(_originalId: "mod1", _modificationType: "mt", _target: motif1, _locationRestriction: "Anywhere.", _monoisotopicMass: 10, _neutralLosses: new Dictionary<DissociationType, List<double>> { { MassSpectrometry.DissociationType.AnyActivationType, new List<double> { 10 } } });

            IDictionary<int, List<Modification>> oneBasedModifications = new Dictionary<int, List<Modification>>
            {
                {2, new List<Modification>{ mod1 }},
                {7, new List<Modification>{ mod1 }},
            };
            Protein protein1 = new Protein("MNLDLDNDL", "prot1", oneBasedModifications: oneBasedModifications);

            Dictionary<int, Modification> allModsOneIsNterminus1 = new Dictionary<int, Modification>
            {
                {2, mod1},
            };
            PeptideWithSetModifications pwsm1 = new PeptideWithSetModifications(protein1, new DigestionParams(), 2, 9, CleavageSpecificity.Unknown, null, 0, allModsOneIsNterminus1, 0);

            Dictionary<int, Modification> allModsOneIsNterminus3 = new Dictionary<int, Modification>
            {
                {7, mod1},
            };
            PeptideWithSetModifications pwsm2 = new PeptideWithSetModifications(protein1, new DigestionParams(), 2, 9, CleavageSpecificity.Unknown, null, 0, allModsOneIsNterminus3, 0);

            CommonParameters CommonParameters = new CommonParameters(digestionParams: new DigestionParams(maxMissedCleavages: 0, minPeptideLength: 1), scoreCutoff: 1);

            PeptideSpectralMatch myPsm = new PeptideSpectralMatch(pwsm1, 0, 10, 0, scan, new DigestionParams(), new List<MatchedFragmentIon>());
            myPsm.AddOrReplace(pwsm2, 10, 0, true, new List<MatchedFragmentIon>(),0);
            
            myPsm.ResolveAllAmbiguities();

            MassDiffAcceptor searchMode = new SinglePpmAroundZeroSearchMode(5);
            List<Protein> proteinList = new List<Protein> { protein1 };

            FdrAnalysisEngine fdrAnalysisEngine = new FdrAnalysisEngine(new List<PeptideSpectralMatch> { myPsm }, searchMode.NumNotches, CommonParameters, new List<string>());
            fdrAnalysisEngine.Run();
            ModificationAnalysisEngine modificationAnalysisEngine = new ModificationAnalysisEngine(new List<PeptideSpectralMatch> { myPsm }, new CommonParameters(), new List<string>());
            var res = (ModificationAnalysisResults)modificationAnalysisEngine.Run();

            Assert.AreEqual(1, res.CountOfEachModSeenOnProteins.Count());
            Assert.AreEqual(2, res.CountOfEachModSeenOnProteins[mod1.IdWithMotif]);
            Assert.AreEqual(0, res.CountOfModsSeenAndLocalized.Count());
            Assert.AreEqual(0, res.CountOfAmbiguousButLocalizedModsSeen.Count);
            Assert.AreEqual(1, res.CountOfUnlocalizedMods[mod1.IdWithMotif]); // Saw it, but not sure where!
            Assert.AreEqual(0, res.CountOfUnlocalizedFormulas.Count());
        }
    }
}