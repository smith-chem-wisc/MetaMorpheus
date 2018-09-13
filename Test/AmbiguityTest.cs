using EngineLayer;
using EngineLayer.ClassicSearch;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Linq;
using TaskLayer;

namespace Test
{
    [TestFixture]
    internal static class AmbiguityTest
    {
        [Test]
        public static void TestResolveAmbiguities()
        {
            Protease protease = new Protease("Custom Protease4", new List<Tuple<string, TerminusType>> { new Tuple<string, TerminusType>("K", TerminusType.C) }, new List<Tuple<string, TerminusType>>(), CleavageSpecificity.Full, null, null, null);
            ProteaseDictionary.Dictionary.Add(protease.Name, protease);
            CommonParameters CommonParameters = new CommonParameters(
                digestionParams: new DigestionParams(
                    protease: protease.Name,
                    minPeptideLength: 1),
                scoreCutoff: 1,
                reportAllAmbiguity: false);

            var myMsDataFile = new TestDataFile();
            var variableModifications = new List<ModificationWithMass>();
            var fixedModifications = new List<ModificationWithMass>();
            var proteinList = new List<Protein> { new Protein("MNNKNKNKQQQ", "Prot1", "organism", null, null, null, null, null, true), new Protein("MNNNKQQQ", "Prot2") };

            var searchModes = new SinglePpmAroundZeroSearchMode(5);

            bool DoPrecursorDeconvolution = true;
            bool UseProvidedPrecursorInfo = true;
            double DeconvolutionIntensityRatio = 4;
            int DeconvolutionMaxAssumedChargeState = 10;
            Tolerance DeconvolutionMassTolerance = new PpmTolerance(5);

            var listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(myMsDataFile, null, DoPrecursorDeconvolution, UseProvidedPrecursorInfo, DeconvolutionIntensityRatio, DeconvolutionMaxAssumedChargeState, DeconvolutionMassTolerance).OrderBy(b => b.PrecursorMass).ToArray();

            PeptideSpectralMatch[] allPsmsArrayt = new PeptideSpectralMatch[listOfSortedms2Scans.Length];
            PeptideSpectralMatch[] allPsmsArrayf = new PeptideSpectralMatch[listOfSortedms2Scans.Length];
            new ClassicSearchEngine(allPsmsArrayt, null, listOfSortedms2Scans, variableModifications, fixedModifications, proteinList, new List<ProductType> { ProductType.B, ProductType.Y }, searchModes, CommonParameters, new List<string>()).Run();
            new ClassicSearchEngine(allPsmsArrayf, null, listOfSortedms2Scans, variableModifications, fixedModifications, proteinList, new List<ProductType> { ProductType.B, ProductType.Y }, searchModes, CommonParameters, new List<string>()).Run();

            var haht = (SequencesToActualProteinPeptidesEngineResults)new SequencesToActualProteinPeptidesEngine(new List<PeptideSpectralMatch>
            { allPsmsArrayt[0] }, proteinList, fixedModifications, variableModifications, new List<ProductType>
            { ProductType.B, ProductType.Y }, new List<DigestionParams> { CommonParameters.DigestionParams }, true, CommonParameters, new List<string>()).Run();
            var hahf = (SequencesToActualProteinPeptidesEngineResults)new SequencesToActualProteinPeptidesEngine(new List<PeptideSpectralMatch>
            { allPsmsArrayf[0] }, proteinList, fixedModifications, variableModifications, new List<ProductType>
            { ProductType.B, ProductType.Y }, new List<DigestionParams> { CommonParameters.DigestionParams }, CommonParameters.ReportAllAmbiguity, CommonParameters, new List<string>()).Run();

            foreach (var huh in allPsmsArrayt)
            {
                if (huh != null)
                {
                    huh.MatchToProteinLinkedPeptides(haht.CompactPeptideToProteinPeptideMatching);
                }
            }

            foreach (var huh in allPsmsArrayf)
            {
                if (huh != null)
                {
                    huh.MatchToProteinLinkedPeptides(hahf.CompactPeptideToProteinPeptideMatching);
                }
            }

            Assert.AreEqual("QQQ", allPsmsArrayt[0].BaseSequence);
            Assert.AreEqual("QQQ", allPsmsArrayf[0].BaseSequence);
            Assert.IsTrue(allPsmsArrayt[0].ProteinLength == null);
            Assert.IsTrue(allPsmsArrayf[0].ProteinLength == 8);
            Assert.IsTrue(allPsmsArrayt[0].OneBasedStartResidueInProtein == null);
            Assert.IsTrue(allPsmsArrayf[0].OneBasedStartResidueInProtein == 6);
        }
    }
}