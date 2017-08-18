﻿using EngineLayer;
using EngineLayer.Analysis;
using MassSpectrometry;
using NUnit.Framework;
using Proteomics;
using System.Collections.Generic;
using System.Linq;
using TaskLayer;
using System;

namespace Test
{
    [TestFixture]
    public static class ModificationAnalysisTest
    {
        #region Public Methods

        [Test]
        public static void TestModificationAnalysis()
        {
            List<Psm>[] newPsms = new List<Psm>[1];
            IScan scan = new ThisTestScan();

            ModificationMotif.TryGetMotif("N", out ModificationMotif motif1);
            ModificationWithMass mod1 = new ModificationWithMass("mod1", null, motif1, TerminusLocalization.Any, 10, null, null, null, "testType");

            ModificationMotif.TryGetMotif("L", out ModificationMotif motif2);
            ModificationWithMass mod2 = new ModificationWithMass("mod2", null, motif2, TerminusLocalization.Any, 10, null, null, null, "testType");

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
            PeptideWithSetModifications pwsm1 = new PeptideWithSetModifications(allModsOneIsNterminus1, 0, protein1, 2, 9, null);
            CompactPeptideBase pep1 = new CompactPeptide(pwsm1, TerminusType.None);

            Dictionary<int, ModificationWithMass> allModsOneIsNterminus2 = new Dictionary<int, ModificationWithMass>
            {
                {2, mod1},
                {7, mod1},
            };
            PeptideWithSetModifications pwsm2 = new PeptideWithSetModifications(allModsOneIsNterminus2, 0, protein1, 2, 9, null);
            CompactPeptideBase pep2 = new CompactPeptide(pwsm2, TerminusType.None);

            Dictionary<int, ModificationWithMass> allModsOneIsNterminus3 = new Dictionary<int, ModificationWithMass>
            {
                {7, mod1},
            };
            PeptideWithSetModifications pwsm3 = new PeptideWithSetModifications(allModsOneIsNterminus3, 0, protein1, 2, 9, null);
            CompactPeptideBase pep3 = new CompactPeptide(pwsm3, TerminusType.None);

            Dictionary<int, ModificationWithMass> allModsOneIsNterminus4 = new Dictionary<int, ModificationWithMass>
            {
                {8, mod1},
            };
            PeptideWithSetModifications pwsm4 = new PeptideWithSetModifications(allModsOneIsNterminus4, 0, protein1, 1, 9, null);
            CompactPeptideBase pep4 = new CompactPeptide(pwsm4, TerminusType.None);

            newPsms[0] = new List<Psm>
            {
                new Psm(pep1, 0,10,0,scan),
                new Psm(pep1, 0,10,0,scan),
                new Psm(pep2, 0,10,0,scan),
                new Psm(pep3, 0,10,0,scan),
                new Psm(pep4, 0,10,0,scan),
            };

            List<MassDiffAcceptor> searchModes = new List<MassDiffAcceptor> { new SinglePpmAroundZeroSearchMode(5) };
            List<Protein> proteinList = new List<Protein> { protein1 };
            Protease protease = GlobalTaskLevelSettings.ProteaseDictionary["trypsin"];
            SequencesToActualProteinPeptidesEngine sequencesToActualProteinPeptidesEngine = new SequencesToActualProteinPeptidesEngine(newPsms, proteinList, searchModes, protease, 0, null, null, InitiatorMethionineBehavior.Variable, new List<ModificationWithMass>(), new List<ModificationWithMass>(), int.MaxValue, new List<string>(), TerminusType.None);
            var nice = (SequencesToActualProteinPeptidesEngineResults)sequencesToActualProteinPeptidesEngine.Run();
            foreach (var psm in newPsms[0])
                psm.MatchToProteinLinkedPeptides(nice.CompactPeptideToProteinPeptideMatching);
            FdrAnalysisEngine fdrAnalysisEngine = new FdrAnalysisEngine(newPsms, searchModes, new List<string>());
            fdrAnalysisEngine.Run();
            ModificationAnalysisEngine modificationAnalysisEngine = new ModificationAnalysisEngine(newPsms, searchModes.Count, new List<string>());
            var res = (ModificationAnalysisResults)modificationAnalysisEngine.Run();

            Assert.AreEqual(1, res.ModsSeenAndLocalized[0].Count());
            Assert.AreEqual(2, res.ModsSeenAndLocalized[0][mod1.id]);

            Assert.AreEqual(0, res.NonLocalizedModsSeen[0].Count()); // Weird

            Assert.AreEqual(2, res.AllModsOnProteins[0].Count());
            Assert.AreEqual(2, res.AllModsOnProteins[0][mod1.id]);
            Assert.AreEqual(1, res.AllModsOnProteins[0][mod2.id]);
        }

        [Test]
        public static void TestModificationAnalysisWithNonLocalizedPtms()
        {
            List<Psm>[] newPsms = new List<Psm>[1];
            IScan scan = new ThisTestScan();

            ModificationMotif.TryGetMotif("N", out ModificationMotif motif1);
            ModificationWithMass mod1 = new ModificationWithMass("mod1", null, motif1, TerminusLocalization.Any, 10, null, new List<double> { 10 }, null, "testType");

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
            PeptideWithSetModifications pwsm1 = new PeptideWithSetModifications(allModsOneIsNterminus1, 0, protein1, 2, 9, null);
            CompactPeptideBase pep1 = new CompactPeptide(pwsm1, TerminusType.None);

            Dictionary<int, ModificationWithMass> allModsOneIsNterminus3 = new Dictionary<int, ModificationWithMass>
            {
                {7, mod1},
            };
            PeptideWithSetModifications pwsm3 = new PeptideWithSetModifications(allModsOneIsNterminus3, 0, protein1, 2, 9, null);
            CompactPeptideBase pep3 = new CompactPeptide(pwsm3, TerminusType.None);

            newPsms[0] = new List<Psm>
            {
                new Psm(pep1, 0,10,0,scan),
                new Psm(pep3, 0,10,0,scan),
            };

            List<MassDiffAcceptor> searchModes = new List<MassDiffAcceptor> { new SinglePpmAroundZeroSearchMode(5) };
            List<Protein> proteinList = new List<Protein> { protein1 };
            Protease protease = GlobalTaskLevelSettings.ProteaseDictionary["trypsin"];
            SequencesToActualProteinPeptidesEngine sequencesToActualProteinPeptidesEngine = new SequencesToActualProteinPeptidesEngine(newPsms, proteinList, searchModes, protease, 0, null, null, InitiatorMethionineBehavior.Variable, new List<ModificationWithMass>(), new List<ModificationWithMass>(), int.MaxValue, new List<string>(), TerminusType.None);
            var nice = (SequencesToActualProteinPeptidesEngineResults)sequencesToActualProteinPeptidesEngine.Run();
            foreach (var psm in newPsms[0])
                psm.MatchToProteinLinkedPeptides(nice.CompactPeptideToProteinPeptideMatching);

            Assert.AreEqual(2, nice.CompactPeptideToProteinPeptideMatching[pep1].Count);

            FdrAnalysisEngine fdrAnalysisEngine = new FdrAnalysisEngine(newPsms, searchModes, new List<string>());
            fdrAnalysisEngine.Run();
            ModificationAnalysisEngine modificationAnalysisEngine = new ModificationAnalysisEngine(newPsms, searchModes.Count, new List<string>());
            var res = (ModificationAnalysisResults)modificationAnalysisEngine.Run();

            Assert.AreEqual(0, res.ModsSeenAndLocalized[0].Count()); // Not localized

            Assert.AreEqual(1, res.NonLocalizedModsSeen[0][mod1.id]); // Saw it, but not sure where!

            Assert.AreEqual(1, res.AllModsOnProteins[0].Count());
            Assert.AreEqual(2, res.AllModsOnProteins[0][mod1.id]);
        }

        #endregion Public Methods

        #region Private Classes

        private class ThisTestScan : IScan
        {
            #region Public Properties

            public string FullFilePath => null;

            public int OneBasedScanNumber => 0;

            public int OneBasedPrecursorScanNumber => 0;

            public double RetentionTime => 0;

            public int NumPeaks => 0;

            public double TotalIonCurrent => 0;

            public int PrecursorCharge => 0;

            public IMzPeak PrecursorMonoisotopicPeak => null;

            public double PrecursorMass => 0;

            #endregion Public Properties
        }

        #endregion Private Classes
    }
}