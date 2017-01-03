using FragmentGeneration;
using IndexSearchAndAnalyze;
using MassSpectrometry;
using MetaMorpheus;
using NUnit.Framework;
using Proteomics;
using Spectra;
using System.Collections.Generic;
using System.IO;

namespace Test
{
    [TestFixture]
    public class AnalysisEngineTest
    {
        [Test]
        public void TestAnalysis()
        {

            List<NewPsm>[] newPsms = null;
            Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>> compactPeptideToProteinPeptideMatching = null;
            List<Protein> proteinList = null;
            List<MorpheusModification> variableModifications = null;
            List<MorpheusModification> fixedModifications = null;
            List<MorpheusModification> localizeableModifications = null;
            IMsDataFile<IMzSpectrum<MzPeak>> iMsDataFile = null;
            Protease protease = null;
            List<SearchMode> searchModes = null;
            double fragmentTolerance = 0;
            UsefulProteomicsDatabases.Generated.unimod unimodDeserialized = null;
            Dictionary<int, ChemicalFormulaModification> uniprotDeseralized = null;

            AnalysisParams analysisParams = new AnalysisParams(newPsms, compactPeptideToProteinPeptideMatching, proteinList, variableModifications, fixedModifications, localizeableModifications, protease, searchModes, iMsDataFile, fragmentTolerance, unimodDeserialized, uniprotDeseralized, (BinTreeStructure myTreeStructure, string s) => { }, (List<NewPsmWithFDR> h, string s) => { });
            AnalysisEngine analysisEngine = new AnalysisEngine(analysisParams);

            Assert.That(() => analysisEngine.Run(), Throws.TypeOf<ValidationException>()
                    .With.Property("Message").EqualTo("newPsms is null"));


        }
    }
}