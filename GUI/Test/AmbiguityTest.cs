using EngineLayer;
using EngineLayer.ClassicSearch;
using MzLibUtil;
using NUnit.Framework;
using Proteomics;
using System.Collections.Generic;
using System.Linq;
using TaskLayer;

namespace Test
{
    [TestFixture]
    internal static class AmbiguityTest
    {
        #region Public Methods

        [Test]
        public static void TestResolveAmbiguities()
        {
            CommonParameters CommonParameters = new CommonParameters
            {
                DigestionParams = new DigestionParams
                {
                    Protease = new Protease("Custom Protease", new List<string> { "K" }, new List<string>(), TerminusType.C, CleavageSpecificity.Full, null, null, null),
                    MinPeptideLength = null,
                },
                ConserveMemory = false,
                ScoreCutoff = 1,
                ReportAllAmbiguity = false
            };

            var myMsDataFile = new TestDataFile();
            var variableModifications = new List<ModificationWithMass>();
            var fixedModifications = new List<ModificationWithMass>();
            var proteinList = new List<Protein> { new Protein("MNNKNKNKQQQ", "Prot1", null, null, null, null, null, true), new Protein("MNNNKQQQ", "Prot2") };

            var searchModes = new SinglePpmAroundZeroSearchMode(5);

            bool DoPrecursorDeconvolution = true;
            bool UseProvidedPrecursorInfo = true;
            double DeconvolutionIntensityRatio = 4;
            int DeconvolutionMaxAssumedChargeState = 10;
            Tolerance DeconvolutionMassTolerance = new PpmTolerance(5);

            var listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(myMsDataFile, null, DoPrecursorDeconvolution, UseProvidedPrecursorInfo, DeconvolutionIntensityRatio, DeconvolutionMaxAssumedChargeState, DeconvolutionMassTolerance).OrderBy(b => b.PrecursorMass).ToArray();

            Psm[] allPsmsArrayt = new Psm[listOfSortedms2Scans.Length];
            Psm[] allPsmsArrayf = new Psm[listOfSortedms2Scans.Length];
            new ClassicSearchEngine(allPsmsArrayt, listOfSortedms2Scans, variableModifications, fixedModifications, proteinList, new List<ProductType> { ProductType.B, ProductType.Y }, searchModes, false, CommonParameters, new List<string>()).Run();
            new ClassicSearchEngine(allPsmsArrayf, listOfSortedms2Scans, variableModifications, fixedModifications, proteinList, new List<ProductType> { ProductType.B, ProductType.Y }, searchModes, false, CommonParameters, new List<string>()).Run();

            var haht = (SequencesToActualProteinPeptidesEngineResults)new SequencesToActualProteinPeptidesEngine(new List<Psm> { allPsmsArrayt[0] }, proteinList, fixedModifications, variableModifications, new List<ProductType> { ProductType.B, ProductType.Y }, new List<DigestionParams> { CommonParameters.DigestionParams }, true, new List<string>()).Run();
            var hahf = (SequencesToActualProteinPeptidesEngineResults)new SequencesToActualProteinPeptidesEngine(new List<Psm> { allPsmsArrayf[0] }, proteinList, fixedModifications, variableModifications, new List<ProductType> { ProductType.B, ProductType.Y }, new List<DigestionParams> { CommonParameters.DigestionParams }, CommonParameters.ReportAllAmbiguity, new List<string>()).Run();

            foreach (var huh in allPsmsArrayt)
                if (huh != null && huh.MostProbableProteinInfo == null)
                    huh.MatchToProteinLinkedPeptides(haht.CompactPeptideToProteinPeptideMatching);

            foreach (var huh in allPsmsArrayf)
                if (huh != null && huh.MostProbableProteinInfo == null)
                    huh.MatchToProteinLinkedPeptides(hahf.CompactPeptideToProteinPeptideMatching);

            Assert.AreEqual("QQQ", allPsmsArrayt[0].BaseSequence);
            Assert.AreEqual("QQQ", allPsmsArrayf[0].BaseSequence);
            Assert.IsTrue(allPsmsArrayt[0].ProteinLength == null);
            Assert.IsTrue(allPsmsArrayf[0].ProteinLength == 8);
            Assert.IsTrue(allPsmsArrayt[0].OneBasedStartResidueInProtein == null);
            Assert.IsTrue(allPsmsArrayf[0].OneBasedStartResidueInProtein == 6);
        }

        #endregion Public Methods
    }
}