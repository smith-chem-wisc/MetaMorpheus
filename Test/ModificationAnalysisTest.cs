using EngineLayer;
using EngineLayer.FdrAnalysis;
using EngineLayer.ModificationAnalysis;
using MassSpectrometry;
using NUnit.Framework;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System.Collections.Generic;

namespace Test
{
    [TestFixture]
    public static class ModificationAnalysisTest
    {
        [Test]
        public static void TestModificationAnalysis()
        {
            IScan scan = new ThisTestScan();

            ModificationMotif.TryGetMotif("N", out ModificationMotif motif1);
            ModificationWithMass mod1 = new ModificationWithMass("mod1", null, motif1, TerminusLocalization.Any, 10);

            ModificationMotif.TryGetMotif("L", out ModificationMotif motif2);
            ModificationWithMass mod2 = new ModificationWithMass("mod2", null, motif2, TerminusLocalization.Any, 10);

            IDictionary<int, List<Modification>> oneBasedModifications = new Dictionary<int, List<Modification>>
            {
                {2, new List<Modification>{ mod1 }},
                {5, new List<Modification>{ mod2 }},
                {7, new List<Modification>{ mod1 }},
            };
            Protein protein1 = new Protein("MNLDLDNDL", "prot1", oneBasedModifications: oneBasedModifications);

            Dictionary<int, ModificationWithMass> allModsOneIsNterminus1 = new Dictionary<int, ModificationWithMass>
            {
                {2, mod1},
            };
            PeptideWithSetModifications pwsm1 = new PeptideWithSetModifications(0, protein1, 2, 9, allModsOneIsNterminus1);
            CompactPeptideBase pep1 = new CompactPeptide(pwsm1, TerminusType.None);

            Dictionary<int, ModificationWithMass> allModsOneIsNterminus2 = new Dictionary<int, ModificationWithMass>
            {
                {2, mod1},
                {7, mod1},
            };
            PeptideWithSetModifications pwsm2 = new PeptideWithSetModifications(0, protein1, 2, 9, allModsOneIsNterminus2);
            CompactPeptideBase pep2 = new CompactPeptide(pwsm2, TerminusType.None);

            Dictionary<int, ModificationWithMass> allModsOneIsNterminus3 = new Dictionary<int, ModificationWithMass>
            {
                {7, mod1},
            };
            PeptideWithSetModifications pwsm3 = new PeptideWithSetModifications(0, protein1, 2, 9, allModsOneIsNterminus3);
            CompactPeptideBase pep3 = new CompactPeptide(pwsm3, TerminusType.None);

            Dictionary<int, ModificationWithMass> allModsOneIsNterminus4 = new Dictionary<int, ModificationWithMass>
            {
                {8, mod1},
            };
            PeptideWithSetModifications pwsm4 = new PeptideWithSetModifications(0, protein1, 1, 9, allModsOneIsNterminus4);
            CompactPeptideBase pep4 = new CompactPeptide(pwsm4, TerminusType.None);

            CommonParameters CommonParameters = new CommonParameters(
                digestionParams: new DigestionParams(
                    maxMissedCleavages: 0,
                    minPeptideLength: 1,
                    maxModificationIsoforms: int.MaxValue),
                scoreCutoff: 1);

            var newPsms = new List<PeptideSpectralMatch>
            {
                new PeptideSpectralMatch(pep1, 0,10,0,scan,CommonParameters.DigestionParams),
                new PeptideSpectralMatch(pep1, 0,10,0,scan, CommonParameters.DigestionParams),
                new PeptideSpectralMatch(pep2, 0,10,0,scan, CommonParameters.DigestionParams),
                new PeptideSpectralMatch(pep3, 0,10,0,scan, CommonParameters.DigestionParams),
                new PeptideSpectralMatch(pep4, 0,10,0,scan,CommonParameters.DigestionParams),
            };

            MassDiffAcceptor searchMode = new SinglePpmAroundZeroSearchMode(5);
            List<Protein> proteinList = new List<Protein> { protein1 };

            SequencesToActualProteinPeptidesEngine sequencesToActualProteinPeptidesEngine =
                new SequencesToActualProteinPeptidesEngine(newPsms, proteinList, new List<ModificationWithMass>(), new List<ModificationWithMass>(),
                    new List<ProductType> { ProductType.B, ProductType.Y }, new List<DigestionParams> { CommonParameters.DigestionParams }, CommonParameters.ReportAllAmbiguity,
                    CommonParameters, new List<string>());
            var nice = (SequencesToActualProteinPeptidesEngineResults)sequencesToActualProteinPeptidesEngine.Run();
            foreach (var psm in newPsms)
            {
                psm.MatchToProteinLinkedPeptides(nice.CompactPeptideToProteinPeptideMatching);
            }
            FdrAnalysisEngine fdrAnalysisEngine = new FdrAnalysisEngine(newPsms, searchMode.NumNotches, CommonParameters, new List<string>());
            fdrAnalysisEngine.Run();
            ModificationAnalysisEngine modificationAnalysisEngine = new ModificationAnalysisEngine(newPsms, new CommonParameters(), new List<string>());
            var res = (ModificationAnalysisResults)modificationAnalysisEngine.Run();

            Assert.AreEqual(2, res.AllModsOnProteins.Count);
            Assert.AreEqual(2, res.AllModsOnProteins[mod1.id]);
            Assert.AreEqual(1, res.AllModsOnProteins[mod2.id]);

            Assert.AreEqual(1, res.ModsSeenAndLocalized.Count);
            Assert.AreEqual(2, res.ModsSeenAndLocalized[mod1.id]);

            Assert.AreEqual(0, res.AmbiguousButLocalizedModsSeen.Count);

            Assert.AreEqual(0, res.UnlocalizedMods.Count);

            Assert.AreEqual(0, res.UnlocalizedFormulas.Count);
        }

        [Test]
        public static void TestModificationAnalysisWithNonLocalizedPtms()
        {
            IScan scan = new ThisTestScan();

            ModificationMotif.TryGetMotif("N", out ModificationMotif motif1);
            ModificationWithMass mod1 = new ModificationWithMass("mod1", "mt", motif1, TerminusLocalization.Any, 10, neutralLosses: new List<double> { 10 });

            IDictionary<int, List<Modification>> oneBasedModifications = new Dictionary<int, List<Modification>>
            {
                {2, new List<Modification>{ mod1 }},
                {7, new List<Modification>{ mod1 }},
            };
            Protein protein1 = new Protein("MNLDLDNDL", "prot1", oneBasedModifications: oneBasedModifications);

            Dictionary<int, ModificationWithMass> allModsOneIsNterminus1 = new Dictionary<int, ModificationWithMass>
            {
                {2, mod1},
            };
            PeptideWithSetModifications pwsm1 = new PeptideWithSetModifications(0, protein1, 2, 9, allModsOneIsNterminus1);
            CompactPeptideBase pep1 = new CompactPeptide(pwsm1, TerminusType.None);

            Dictionary<int, ModificationWithMass> allModsOneIsNterminus3 = new Dictionary<int, ModificationWithMass>
            {
                {7, mod1},
            };
            PeptideWithSetModifications pwsm3 = new PeptideWithSetModifications(0, protein1, 2, 9, allModsOneIsNterminus3);
            CompactPeptideBase pep3 = new CompactPeptide(pwsm3, TerminusType.None);

            CommonParameters CommonParameters = new CommonParameters(digestionParams: new DigestionParams(maxMissedCleavages: 0, minPeptideLength: 1), scoreCutoff: 1);

            var newPsms = new List<PeptideSpectralMatch>
            {
                new PeptideSpectralMatch(pep1, 0,10,0,scan, CommonParameters.DigestionParams),
                new PeptideSpectralMatch(pep3, 0,10,0,scan, CommonParameters.DigestionParams),
            };

            MassDiffAcceptor searchMode = new SinglePpmAroundZeroSearchMode(5);
            List<Protein> proteinList = new List<Protein> { protein1 };

            SequencesToActualProteinPeptidesEngine sequencesToActualProteinPeptidesEngine = new SequencesToActualProteinPeptidesEngine(newPsms, proteinList, new List<ModificationWithMass>(),
                new List<ModificationWithMass>(), new List<ProductType> { ProductType.B, ProductType.Y }, new List<DigestionParams> { CommonParameters.DigestionParams },
                CommonParameters.ReportAllAmbiguity, CommonParameters, new List<string>());

            var nice = (SequencesToActualProteinPeptidesEngineResults)sequencesToActualProteinPeptidesEngine.Run();
            foreach (var psm in newPsms)
            {
                psm.MatchToProteinLinkedPeptides(nice.CompactPeptideToProteinPeptideMatching);
            }

            Assert.AreEqual(2, nice.CompactPeptideToProteinPeptideMatching[pep1].Count);

            FdrAnalysisEngine fdrAnalysisEngine = new FdrAnalysisEngine(newPsms, searchMode.NumNotches, CommonParameters, new List<string>());
            fdrAnalysisEngine.Run();
            ModificationAnalysisEngine modificationAnalysisEngine = new ModificationAnalysisEngine(newPsms, new CommonParameters(), new List<string>());
            var res = (ModificationAnalysisResults)modificationAnalysisEngine.Run();

            Assert.AreEqual(1, res.AllModsOnProteins.Count);
            Assert.AreEqual(2, res.AllModsOnProteins[mod1.id]);

            Assert.AreEqual(0, res.ModsSeenAndLocalized.Count);

            Assert.AreEqual(0, res.AmbiguousButLocalizedModsSeen.Count);

            Assert.AreEqual(1, res.UnlocalizedMods[mod1.id]); // Saw it, but not sure where!

            Assert.AreEqual(0, res.UnlocalizedFormulas.Count);
        }
    }

    internal class ThisTestScan : IScan
    {
        public string FullFilePath => null;
        public int OneBasedScanNumber => 0;
        public int? OneBasedPrecursorScanNumber => 0;
        public double RetentionTime => 0;
        public int NumPeaks => 0;
        public double TotalIonCurrent => 0;
        public int PrecursorCharge => 0;
        public double PrecursorMonoisotopicPeakMz => 0;
        public double PrecursorMass => 0;
    }
}