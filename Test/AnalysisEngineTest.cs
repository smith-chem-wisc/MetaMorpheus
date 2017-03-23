using EngineLayer;
using EngineLayer.Analysis;
using EngineLayer.ModernSearch;
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
    public class AnalysisEngineTests
    {

        #region Public Methods

        [Test]
        public static void TestAnalysisEngineTests()
        {
            List<ModificationWithMass> localizeableModifications = null;
            List<ModificationWithMass> variableModifications = new List<ModificationWithMass>();
            List<ModificationWithMass> fixedModifications = new List<ModificationWithMass>();

            PsmParent[][] newPsms = new PsmParent[1][];

            Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>> compactPeptideToProteinPeptideMatching = new Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>>();

            var proteinList = new List<Protein> { new Protein("MNNNKQQQ", "accession", null, new Dictionary<int, List<Modification>>(), new int?[0], new int?[0], new string[0], null, null, false, false, null) };

            var protease = new Protease("Custom Protease", new List<string> { "K" }, new List<string>(), TerminusType.C, CleavageSpecificity.Full, null, null, null);

            PeptideWithPossibleModifications modPep = proteinList.First().Digest(protease, 0, InitiatorMethionineBehavior.Variable, fixedModifications).Last();
            HashSet<PeptideWithSetModifications> value = new HashSet<PeptideWithSetModifications> { modPep.GetPeptidesWithSetModifications(variableModifications, 4096, 3).First() };
            CompactPeptide compactPeptide1 = new CompactPeptide(value.First(), variableModifications, localizeableModifications, fixedModifications);

            Assert.AreEqual("QQQ", value.First().BaseSequence);
            PeptideWithPossibleModifications modPep2 = proteinList.First().Digest(protease, 0, InitiatorMethionineBehavior.Variable, fixedModifications).First();
            HashSet<PeptideWithSetModifications> value2 = new HashSet<PeptideWithSetModifications> { modPep2.GetPeptidesWithSetModifications(variableModifications, 4096, 3).First() };
            CompactPeptide compactPeptide2 = new CompactPeptide(value2.First(), variableModifications, localizeableModifications, fixedModifications);

            Assert.AreEqual("MNNNK", value2.First().BaseSequence);

            PeptideWithPossibleModifications modPep3 = proteinList.First().Digest(protease, 0, InitiatorMethionineBehavior.Variable, fixedModifications).ToList()[1];
            HashSet<PeptideWithSetModifications> value3 = new HashSet<PeptideWithSetModifications> { modPep3.GetPeptidesWithSetModifications(variableModifications, 4096, 3).First() };
            CompactPeptide compactPeptide3 = new CompactPeptide(value3.First(), variableModifications, localizeableModifications, fixedModifications);
            Assert.AreEqual("NNNK", value3.First().BaseSequence);

            newPsms[0] = new PsmParent[] { new PsmModern(compactPeptide1, null, 1, 1, 1, 2, 2, 1,1, 1, 1, 3,0),
                                           new PsmModern(compactPeptide2, null, 2,2,2+132.040,3,3,2,2,2,2,2,0),
                                           new PsmModern(compactPeptide3, null, 3,3,3,4,3, 3, 3,3,3,3,0) };

            compactPeptideToProteinPeptideMatching.Add(compactPeptide1, value);
            compactPeptideToProteinPeptideMatching.Add(compactPeptide2, value2);
            compactPeptideToProteinPeptideMatching.Add(compactPeptide3, value3);

            Action<BinTreeStructure, string> action1 = (BinTreeStructure l, string s) => {; };
            Tolerance fragmentTolerance = new Tolerance(ToleranceUnit.PPM, 10);
            IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> myMsDataFile = new TestDataFile(new List<PeptideWithSetModifications> { value.First(), value2.First(), value3.First() });

            var searchModes = new List<SearchMode> { new SinglePpmAroundZeroSearchMode(5) };
            Action<List<ProteinGroup>, string> action3 = null;
            Action<List<NewPsmWithFdr>, string> action2 = (List<NewPsmWithFdr> l, string s) => {; };
            bool doParsimony = false;
            bool noOneHitWonders = false;
            bool modPepsAreUnique = false;
            bool quant = false;
            double quantRtTol = 0;
            double quantPpmTol = 0;
            AnalysisEngine engine = new AnalysisEngine(newPsms, compactPeptideToProteinPeptideMatching, proteinList, variableModifications, fixedModifications, localizeableModifications, protease, searchModes, myMsDataFile, fragmentTolerance, action1, action2, action3, doParsimony, noOneHitWonders, modPepsAreUnique, 2, 4096, true, new List<ProductType> { ProductType.B, ProductType.Y }, 0.003, InitiatorMethionineBehavior.Variable, new List<string> { "ff" }, quant, quantRtTol, quantPpmTol);

            engine.Run();
        }

        #endregion Public Methods

    }
}