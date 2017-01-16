using InternalLogicEngineLayer;
using MassSpectrometry;
using NUnit.Framework;
using OldInternalLogic;
using Spectra;
using System;
using System.Collections.Generic;

namespace Test
{
    [TestFixture]
    public class AnalysisEngineTests
    {

        #region Public Methods

        [Test]
        public void TestAnalysisEngineTests()
        {
            ParentSpectrumMatch[][] newPsms = null;
            Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>> compactPeptideToProteinPeptideMatching = null;
            List<Protein> proteinList = null;
            List<MorpheusModification> variableModifications = null;
            List<MorpheusModification> fixedModifications = null;
            Action<BinTreeStructure, string> action1 = null;
            Tolerance fragmentTolerance = null;
            IMsDataFile<IMzSpectrum<MzPeak>> myMsDataFile = null;
            Protease protease = null;
            List<MorpheusModification> localizeableModifications = null;
            List<SearchMode> searchModes = null;
            Action<List<ProteinGroup>, string> action3 = null;
            Action<List<NewPsmWithFDR>, string> action2 = null;
            bool doParsimony = false;
            AnalysisEngine engine = new AnalysisEngine(newPsms, compactPeptideToProteinPeptideMatching, proteinList, variableModifications, fixedModifications, localizeableModifications, protease, searchModes, myMsDataFile, fragmentTolerance, action1, action2, action3, doParsimony);
        }

        #endregion Public Methods

    }
}