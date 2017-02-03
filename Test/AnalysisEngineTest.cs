using EngineLayer;
using EngineLayer.Analysis;
using EngineLayer.ModernSearch;
using MassSpectrometry;
using NUnit.Framework;

using Spectra;
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
            List<MetaMorpheusModification> localizeableModifications = null;
            List<MetaMorpheusModification> variableModifications = new List<MetaMorpheusModification>();
            List<MetaMorpheusModification> fixedModifications = new List<MetaMorpheusModification>();

            PsmParent[][] newPsms = new PsmParent[1][];

            Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>> compactPeptideToProteinPeptideMatching = new Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>>();

            var proteinList = new List<Protein> { new Protein("MNNNKQQQ", "accession", new Dictionary<int, List<MetaMorpheusModification>>(), new int[0], new int[0], new string[0], null, null, 0, false, false) };

            PeptideWithPossibleModifications modPep = new PeptideWithPossibleModifications(6, 8, proteinList.First(), 0, "ya");
            HashSet<PeptideWithSetModifications> value = new HashSet<PeptideWithSetModifications> { modPep.GetPeptideWithSetModifications(variableModifications, 4096, 3, new List<MetaMorpheusModification>()).First() };
            CompactPeptide compactPeptide1 = new CompactPeptide(value.First(), variableModifications, localizeableModifications, fixedModifications);

            Assert.AreEqual("QQQ", value.First().BaseSequence);
            PeptideWithPossibleModifications modPep2 = new PeptideWithPossibleModifications(1, 5, proteinList.First(), 0, "ya");
            HashSet<PeptideWithSetModifications> value2 = new HashSet<PeptideWithSetModifications> { modPep2.GetPeptideWithSetModifications(variableModifications, 4096, 3, new List<MetaMorpheusModification>()).First() };
            CompactPeptide compactPeptide2 = new CompactPeptide(value2.First(), variableModifications, localizeableModifications, fixedModifications);

            Assert.AreEqual("MNNNK", value2.First().BaseSequence);

            PeptideWithPossibleModifications modPep3 = new PeptideWithPossibleModifications(2, 5, proteinList.First(), 0, "ya");
            HashSet<PeptideWithSetModifications> value3 = new HashSet<PeptideWithSetModifications> { modPep3.GetPeptideWithSetModifications(variableModifications, 4096, 3, new List<MetaMorpheusModification>()).First() };
            CompactPeptide compactPeptide3 = new CompactPeptide(value3.First(), variableModifications, localizeableModifications, fixedModifications);
            Assert.AreEqual("NNNK", value3.First().BaseSequence);

            newPsms[0] = new PsmParent[] { new PsmModern(compactPeptide1, null,1, 1, 1, 2, 1, 1, 1, 1, 3,0),
                                                     new PsmModern(compactPeptide2, null, 2,2,2+132.040,3,2,2,2,2,2,0),
                                                     new PsmModern(compactPeptide3, null, 3,3,3,4,3,3,3,3,3,0) };

            compactPeptideToProteinPeptideMatching.Add(compactPeptide1, value);
            compactPeptideToProteinPeptideMatching.Add(compactPeptide2, value2);
            compactPeptideToProteinPeptideMatching.Add(compactPeptide3, value3);

            Action<BinTreeStructure, string> action1 = (BinTreeStructure l, string s) => {; };
            Tolerance fragmentTolerance = new Tolerance(ToleranceUnit.PPM, 10);
            IMsDataFile<IMzSpectrum<MzPeak>> myMsDataFile = new TestDataFile(new List<PeptideWithSetModifications> { value.First(), value2.First(), value3.First() });

            var protease = new Protease("Custom Protease", new List<string> { "K" }, new List<string>(), TerminusType.C, CleavageSpecificity.Full, null, null, null);

            var searchModes = new List<SearchMode> { new SinglePpmAroundZeroSearchMode(5) };
            Action<List<ProteinGroup>, string> action3 = null;
            Action<List<NewPsmWithFdr>, string> action2 = (List<NewPsmWithFdr> l, string s) => {; };
            bool doParsimony = false;
            AnalysisEngine engine = new AnalysisEngine(newPsms, compactPeptideToProteinPeptideMatching, proteinList, variableModifications, fixedModifications, localizeableModifications, protease, searchModes, myMsDataFile, fragmentTolerance, action1, action2, action3, doParsimony, 2, 4096, true, new List<ProductType> { ProductType.B, ProductType.Y }, 0.003);

            engine.Run();
        }

        #endregion Public Methods

    }
}