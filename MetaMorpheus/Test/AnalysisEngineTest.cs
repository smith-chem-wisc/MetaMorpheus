using EngineLayer;
using EngineLayer.FdrAnalysis;
using EngineLayer.HistogramAnalysis;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Proteomics;
using Omics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Linq;
using EngineLayer.DatabaseLoading;
using Omics.Digestion;
using Omics.Modifications;
using TaskLayer;
using Omics;

namespace Test
{
    [TestFixture]
    public static class AnalysisEngineTests
    {
        [Test]
        public static void TestAnalysisEngineTests()
        {
            List<DigestionMotif> motifs = new List<DigestionMotif> { new DigestionMotif("K", null, 1, null) };
            Protease protease = new Protease("Custom Protease5", CleavageSpecificity.Full, null, null, motifs);
            ProteaseDictionary.Dictionary.Add(protease.Name, protease);
            CommonParameters CommonParameters = new CommonParameters(
                digestionParams: new DigestionParams(
                    protease: protease.Name,
                    maxMissedCleavages: 0,
                    minPeptideLength: 1,
                    maxModificationIsoforms: 1042),
                scoreCutoff: 1,
                productMassTolerance: new PpmTolerance(10));
            var fsp = new List<(string fileName, CommonParameters fileSpecificParameters)>();
            fsp.Add(("", CommonParameters));

            List<Modification> localizeableModifications = new List<Modification>();
            List<Modification> variableModifications = new List<Modification>();
            List<Modification> fixedModifications = new List<Modification>();

            Dictionary<Modification, ushort> modsDictionary = new Dictionary<Modification, ushort>();
            foreach (var mod in fixedModifications)
                modsDictionary.Add(mod, 0);
            int i = 1;
            foreach (var mod in variableModifications)
            {
                modsDictionary.Add(mod, (ushort)i);
                i++;
            }
            foreach (var mod in localizeableModifications)
            {
                modsDictionary.Add(mod, (ushort)i);
                i++;
            }

            var proteinList = new List<Protein> { new Protein("MNNNKQQQ", "accession") };
            var modPep = proteinList.First().Digest(CommonParameters.DigestionParams, fixedModifications, variableModifications).Last();
            HashSet<IBioPolymerWithSetMods> value1 = new HashSet<IBioPolymerWithSetMods> { modPep };
            var compactPeptide1 = value1.First();

            Assert.That(value1.First().BaseSequence, Is.EqualTo("QQQ"));
            var modPep2 = proteinList.First().Digest(CommonParameters.DigestionParams, fixedModifications, variableModifications).First();
            HashSet<IBioPolymerWithSetMods> value2 = new HashSet<IBioPolymerWithSetMods> { modPep2 };
            var compactPeptide2 = value2.First();

            Assert.That(value2.First().BaseSequence, Is.EqualTo("MNNNK"));

            var modPep3 = proteinList.First().Digest(CommonParameters.DigestionParams, fixedModifications, variableModifications).ToList()[1];
            HashSet<IBioPolymerWithSetMods> value3 = new HashSet<IBioPolymerWithSetMods> { modPep3 };
            var compactPeptide3 = value3.First();
            Assert.That(value3.First().BaseSequence, Is.EqualTo("NNNK"));


            Ms2ScanWithSpecificMass scanA = new Ms2ScanWithSpecificMass(new MsDataScan(new MzSpectrum(new double[] { 1 }, new double[] { 1 }, false), 2, 1, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, null, null, "scan=1", double.NaN, null, null, double.NaN, null, DissociationType.AnyActivationType, 1, null), 1, 1, null, new CommonParameters());
            Ms2ScanWithSpecificMass scanB = new Ms2ScanWithSpecificMass(new MsDataScan(new MzSpectrum(new double[] { 1 }, new double[] { 1 }, false), 3, 1, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, null, null, "scan=2", double.NaN, null, null, double.NaN, null, DissociationType.AnyActivationType, 1, null), 2 + 132.040, 1, null, new CommonParameters());
            Ms2ScanWithSpecificMass scanC = new Ms2ScanWithSpecificMass(new MsDataScan(new MzSpectrum(new double[] { 1 }, new double[] { 1 }, false), 4, 1, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, null, null, "scan=3", double.NaN, null, null, double.NaN, null, DissociationType.AnyActivationType, 1, null), 3, 1, null, new CommonParameters());

            SpectralMatch matchA = new PeptideSpectralMatch(compactPeptide1, 0, 10, 0, scanA, CommonParameters, new List<MatchedFragmentIon>());
            SpectralMatch matchB = new PeptideSpectralMatch(compactPeptide2, 0, 10, 0, scanB, CommonParameters, new List<MatchedFragmentIon>());
            SpectralMatch matchC = new PeptideSpectralMatch(compactPeptide3, 0, 10, 0, scanC, CommonParameters, new List<MatchedFragmentIon>());

            var newPsms = new List<SpectralMatch> { matchA, matchB, matchC };

            MsDataFile myMsDataFile = new TestDataFile(new List<IBioPolymerWithSetMods> { value1.First(), value2.First(), value3.First() });

            var searchMode = new SinglePpmAroundZeroSearchMode(5);
            Action<List<SpectralMatch>, string, List<string>> action2 = (List<SpectralMatch> l, string s, List<string> sdf) => {; };

            Tolerance DeconvolutionMassTolerance = new PpmTolerance(5);

            var arrayOfMs2ScansSortedByMass = MetaMorpheusTask.GetMs2Scans(myMsDataFile, null, new CommonParameters()).OrderBy(b => b.PrecursorMass).ToArray();

            Action<BinTreeStructure, string> action1 = (BinTreeStructure l, string s) =>
            {
                Assert.That(l.FinalBins.Count, Is.EqualTo(1));
            };

            FdrAnalysisEngine engine = new FdrAnalysisEngine(newPsms, searchMode.NumNotches, CommonParameters, fsp, new List<string> { "ff" });

            engine.Run();
        }

        [Test]
        public static void TestRunSpecificPostSearchAnalysis()
        {
            //code coverage unit test for an unused abstract method in post search analysis
            var task = new PostSearchAnalysisTask();
            task.RunTask(TestContext.CurrentContext.TestDirectory, new List<DbForTask>(), new List<string>(), "");
        }
    }
}