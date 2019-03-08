using EngineLayer;
using EngineLayer.FdrAnalysis;
using EngineLayer.HistogramAnalysis;
using MassSpectrometry;
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
            HashSet<PeptideWithSetModifications> value1 = new HashSet<PeptideWithSetModifications> { modPep };
            PeptideWithSetModifications compactPeptide1 = value1.First();

            Assert.AreEqual("QQQ", value1.First().BaseSequence);
            var modPep2 = proteinList.First().Digest(CommonParameters.DigestionParams, fixedModifications, variableModifications).First();
            HashSet<PeptideWithSetModifications> value2 = new HashSet<PeptideWithSetModifications> { modPep2 };
            PeptideWithSetModifications compactPeptide2 = value2.First();

            Assert.AreEqual("MNNNK", value2.First().BaseSequence);

            var modPep3 = proteinList.First().Digest(CommonParameters.DigestionParams, fixedModifications, variableModifications).ToList()[1];
            HashSet<PeptideWithSetModifications> value3 = new HashSet<PeptideWithSetModifications> { modPep3 };
            PeptideWithSetModifications compactPeptide3 = value3.First();
            Assert.AreEqual("NNNK", value3.First().BaseSequence);


            Ms2ScanWithSpecificMass scanA = new Ms2ScanWithSpecificMass(new MsDataScan(new MzSpectrum(new double[] { 1 }, new double[] { 1 }, false), 2, 1, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, null, null, "scan=1", double.NaN, null, null, double.NaN, null, DissociationType.AnyActivationType, 1, null), 1, 1, null, new CommonParameters());
            Ms2ScanWithSpecificMass scanB = new Ms2ScanWithSpecificMass(new MsDataScan(new MzSpectrum(new double[] { 1 }, new double[] { 1 }, false), 3, 1, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, null, null, "scan=2", double.NaN, null, null, double.NaN, null, DissociationType.AnyActivationType, 1, null), 2 + 132.040, 1, null, new CommonParameters());
            Ms2ScanWithSpecificMass scanC = new Ms2ScanWithSpecificMass(new MsDataScan(new MzSpectrum(new double[] { 1 }, new double[] { 1 }, false), 4, 1, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, null, null, "scan=3", double.NaN, null, null, double.NaN, null, DissociationType.AnyActivationType, 1, null), 3, 1, null, new CommonParameters());

            PeptideSpectralMatch matchA = new PeptideSpectralMatch(compactPeptide1, 0, 0, 0, scanA, CommonParameters.DigestionParams, new List<MatchedFragmentIon>());
            PeptideSpectralMatch matchB = new PeptideSpectralMatch(compactPeptide2, 0, 0, 0, scanB, CommonParameters.DigestionParams, new List<MatchedFragmentIon>());
            PeptideSpectralMatch matchC = new PeptideSpectralMatch(compactPeptide3, 0, 0, 0, scanC, CommonParameters.DigestionParams, new List<MatchedFragmentIon>());

            var newPsms = new List<PeptideSpectralMatch> { matchA, matchB, matchC };

            MsDataFile myMsDataFile = new TestDataFile(new List<PeptideWithSetModifications> { value1.First(), value2.First(), value3.First() });

            var searchMode = new SinglePpmAroundZeroSearchMode(5);
            Action<List<PeptideSpectralMatch>, string, List<string>> action2 = (List<PeptideSpectralMatch> l, string s, List<string> sdf) => {; };

            Tolerance DeconvolutionMassTolerance = new PpmTolerance(5);

            var arrayOfMs2ScansSortedByMass = MetaMorpheusTask.GetMs2Scans(myMsDataFile, null, new CommonParameters()).OrderBy(b => b.PrecursorMass).ToArray();

            Action<BinTreeStructure, string> action1 = (BinTreeStructure l, string s) =>
            {
                Assert.AreEqual(1, l.FinalBins.Count);
            };

            FdrAnalysisEngine engine = new FdrAnalysisEngine(newPsms, searchMode.NumNotches, CommonParameters, new List<string> { "ff" });

            engine.Run();
        }
    }
}