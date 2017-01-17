using InternalLogicEngineLayer;
using MassSpectrometry;
using NUnit.Framework;
using OldInternalLogic;
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
        public void TestAnalysisEngineTests()
        {
            List<MorpheusModification> localizeableModifications = null;
            List<MorpheusModification> variableModifications = new List<MorpheusModification>();

            ParentSpectrumMatch[][] newPsms = new ParentSpectrumMatch[1][];
            Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>> compactPeptideToProteinPeptideMatching = new Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>>();

            var proteinList = new List<Protein> { new Protein("MNNNKQQQ", "accession", null, new Dictionary<int, List<MorpheusModification>>(), new int[0], new int[0], new string[0], null, null, 0, false) };

            Dictionary<int, MorpheusModification> twoBasedVariableAndLocalizeableModificationss = new Dictionary<int, MorpheusModification>(); ;
            PeptideWithPossibleModifications modPep = new PeptideWithPossibleModifications(6, 8, proteinList.First(), 0, "ya");

            HashSet<PeptideWithSetModifications> value = new HashSet<PeptideWithSetModifications> { new PeptideWithSetModifications(modPep, twoBasedVariableAndLocalizeableModificationss) };

            CompactPeptide key = new CompactPeptide(value.First(), variableModifications, localizeableModifications);

            compactPeptideToProteinPeptideMatching.Add(key, value);

            List<MorpheusModification> fixedModifications = new List<MorpheusModification>();
            Action<BinTreeStructure, string> action1 = null;
            Tolerance fragmentTolerance = null;
            IMsDataFile<IMzSpectrum<MzPeak>> myMsDataFile = null;

            var protease = new Protease("Custom Protease", new List<string> { "K" }, new List<string>(), Terminus.C, CleavageSpecificity.Full, null, null, null);

            var searchModes = new List<SearchMode> { new SinglePpmAroundZeroSearchMode("", 5) };
            Action<List<ProteinGroup>, string> action3 = null;
            Action<List<NewPsmWithFDR>, string> action2 = null;
            bool doParsimony = false;
            AnalysisEngine engine = new AnalysisEngine(newPsms, compactPeptideToProteinPeptideMatching, proteinList, variableModifications, fixedModifications, localizeableModifications, protease, searchModes, myMsDataFile, fragmentTolerance, action1, action2, action3, doParsimony);

            var res = (AnalysisResults)engine.Run();
        }

        #endregion Public Methods

    }
}