using EngineLayer;
using EngineLayer.ClassicSearch;
using MzLibUtil;
using NUnit.Framework;
using Proteomics;
using Proteomics.Fragmentation;
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
            Protease protease = new Protease("Custom Protease4", CleavageSpecificity.Full, null, null, new List<DigestionMotif> { new DigestionMotif("K", null, 1, "") });
            ProteaseDictionary.Dictionary.Add(protease.Name, protease);
            CommonParameters CommonParameters_t = new CommonParameters(
                dissociationType: MassSpectrometry.DissociationType.HCD,
                digestionParams: new DigestionParams(
                    protease: protease.Name,
                    minPeptideLength: 1),
                scoreCutoff: 1,
                reportAllAmbiguity: true);

            CommonParameters CommonParameters_f = new CommonParameters(
                dissociationType: MassSpectrometry.DissociationType.HCD,
                digestionParams: new DigestionParams(
                    protease: protease.Name,
                    minPeptideLength: 1),
                scoreCutoff: 1,
                reportAllAmbiguity: false);

            var myMsDataFile = new TestDataFile();
            var variableModifications = new List<Modification>();
            var fixedModifications = new List<Modification>();
            var proteinList = new List<Protein> { new Protein("MNNKNKNKQQQ", "Prot1"), new Protein("MNNNKQQQ", "Prot2") };
            var searchModes = new SinglePpmAroundZeroSearchMode(5);

            Tolerance DeconvolutionMassTolerance = new PpmTolerance(5);

            var listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(myMsDataFile, null, new CommonParameters()).OrderBy(b => b.PrecursorMass).ToArray();
            
            PeptideSpectralMatch[] allPsmsArray_withAmbiguity = new PeptideSpectralMatch[listOfSortedms2Scans.Length];
            
            PeptideSpectralMatch[] allPsmsArray_withOutAmbiguity = new PeptideSpectralMatch[listOfSortedms2Scans.Length];

            new ClassicSearchEngine(allPsmsArray_withAmbiguity, listOfSortedms2Scans, variableModifications, fixedModifications, proteinList, searchModes, CommonParameters_t, new List<string>()).Run(); //report all ambiguity TRUE
            new ClassicSearchEngine(allPsmsArray_withOutAmbiguity, listOfSortedms2Scans, variableModifications, fixedModifications, proteinList, searchModes, CommonParameters_f, new List<string>()).Run(); //report all ambiguity FALSE

            Assert.AreEqual("QQQ", allPsmsArray_withAmbiguity[0].BaseSequence);
            Assert.AreEqual("QQQ", allPsmsArray_withOutAmbiguity[0].BaseSequence);
            Assert.IsTrue(allPsmsArray_withAmbiguity[0].ProteinLength == null);
            Assert.IsTrue(allPsmsArray_withOutAmbiguity[0].ProteinLength != null);
            Assert.IsTrue(allPsmsArray_withAmbiguity[0].OneBasedStartResidueInProtein == null);
            Assert.IsTrue(allPsmsArray_withOutAmbiguity[0].OneBasedStartResidueInProtein != null);
        }
    }
}