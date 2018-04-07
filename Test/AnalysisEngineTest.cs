using EngineLayer;
using EngineLayer.FdrAnalysis;
using EngineLayer.HistogramAnalysis;
using IO.MzML;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Proteomics;
using System;
using System.Collections.Generic;
using System.Linq;
using TaskLayer;

namespace Test
{
    [TestFixture]
    public static class AnalysisEngineTests
    {
        #region Public Methods

        [Test]
        public static void TestAnalysisEngineTests()
        {
            CommonParameters CommonParameters = new CommonParameters
            {
                DigestionParams = new DigestionParams
                {
                    Protease = new Protease("Custom Protease", new List<string> { "K" }, new List<string>(), TerminusType.C, CleavageSpecificity.Full, null, null, null),
                    MinPeptideLength = 1,
                    MaxMissedCleavages = 0,
                    MaxModificationIsoforms = 1042,
                },
                ConserveMemory = false,
                ScoreCutoff = 1,
                ProductMassTolerance = new PpmTolerance(10),
            };

            List<ModificationWithMass> localizeableModifications = new List<ModificationWithMass>();
            List<ModificationWithMass> variableModifications = new List<ModificationWithMass>();
            List<ModificationWithMass> fixedModifications = new List<ModificationWithMass>();

            Dictionary<ModificationWithMass, ushort> modsDictionary = new Dictionary<ModificationWithMass, ushort>();
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
            CompactPeptide compactPeptide1 = new CompactPeptide(value1.First(), TerminusType.None);

            Assert.AreEqual("QQQ", value1.First().BaseSequence);
            var modPep2 = proteinList.First().Digest(CommonParameters.DigestionParams, fixedModifications, variableModifications).First();
            HashSet<PeptideWithSetModifications> value2 = new HashSet<PeptideWithSetModifications> { modPep2 };
            CompactPeptide compactPeptide2 = new CompactPeptide(value2.First(), TerminusType.None);

            Assert.AreEqual("MNNNK", value2.First().BaseSequence);

            var modPep3 = proteinList.First().Digest(CommonParameters.DigestionParams, fixedModifications, variableModifications).ToList()[1];
            HashSet<PeptideWithSetModifications> value3 = new HashSet<PeptideWithSetModifications> { modPep3 };
            CompactPeptide compactPeptide3 = new CompactPeptide(value3.First(), TerminusType.None);
            Assert.AreEqual("NNNK", value3.First().BaseSequence);

            //newPsms[0] = new List<PsmParent>[] { new List<PsmParent>{ new PsmModern(compactPeptide1, null, 1,  1, 2, 2, 1,1, 1, 1, 3,0) },
            //                                     new List<PsmParent>{  new PsmModern(compactPeptide2, null, 2,2+132.040,3,3,2,2,2,2,2,0) },
            //                                     new List<PsmParent>{ new PsmModern(compactPeptide3, null, 3, 3, 4, 3, 3, 3, 3, 3, 3, 0)} };

            Ms2ScanWithSpecificMass scanA = new Ms2ScanWithSpecificMass(new MzmlScanWithPrecursor(2, new MzmlMzSpectrum(new double[] { 1 }, new double[] { 1 }, false), 1, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, double.NaN, null, null, double.NaN, null, DissociationType.AnyActivationType, 1, null, null, "scan=1"), 1, 1, null);
            Ms2ScanWithSpecificMass scanB = new Ms2ScanWithSpecificMass(new MzmlScanWithPrecursor(3, new MzmlMzSpectrum(new double[] { 1 }, new double[] { 1 }, false), 1, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, double.NaN, null, null, double.NaN, null, DissociationType.AnyActivationType, 1, null, null, "scan=2"), 2 + 132.040, 1, null);
            Ms2ScanWithSpecificMass scanC = new Ms2ScanWithSpecificMass(new MzmlScanWithPrecursor(4, new MzmlMzSpectrum(new double[] { 1 }, new double[] { 1 }, false), 1, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, double.NaN, null, null, double.NaN, null, DissociationType.AnyActivationType, 1, null, null, "scan=3"), 3, 1, null);

            PeptideSpectralMatch matchA = new PeptideSpectralMatch(compactPeptide1, 0, 0, 0, scanA);
            PeptideSpectralMatch matchB = new PeptideSpectralMatch(compactPeptide2, 0, 0, 0, scanB);
            PeptideSpectralMatch matchC = new PeptideSpectralMatch(compactPeptide3, 0, 0, 0, scanC);

            var newPsms = new List<PeptideSpectralMatch> { matchA, matchB, matchC };

            IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> myMsDataFile = new TestDataFile(new List<PeptideWithSetModifications> { value1.First(), value2.First(), value3.First() });

            var searchMode = new SinglePpmAroundZeroSearchMode(5);
            Action<List<PeptideSpectralMatch>, string, List<string>> action2 = (List<PeptideSpectralMatch> l, string s, List<string> sdf) => {; };

            bool DoPrecursorDeconvolution = true;
            bool UseProvidedPrecursorInfo = true;
            double DeconvolutionIntensityRatio = 4;
            int DeconvolutionMaxAssumedChargeState = 10;
            Tolerance DeconvolutionMassTolerance = new PpmTolerance(5);

            var arrayOfMs2ScansSortedByMass = MetaMorpheusTask.GetMs2Scans(myMsDataFile, null, DoPrecursorDeconvolution, UseProvidedPrecursorInfo, DeconvolutionIntensityRatio, DeconvolutionMaxAssumedChargeState, DeconvolutionMassTolerance).OrderBy(b => b.PrecursorMass).ToArray();

            Action<BinTreeStructure, string> action1 = (BinTreeStructure l, string s) =>
            {
                Assert.AreEqual(1, l.FinalBins.Count);
            };

            SequencesToActualProteinPeptidesEngine sequencesToActualProteinPeptidesEngine = new SequencesToActualProteinPeptidesEngine(newPsms, proteinList, fixedModifications, variableModifications, new List<ProductType> { ProductType.B, ProductType.Y }, new List<IDigestionParams> { CommonParameters.DigestionParams }, CommonParameters.ReportAllAmbiguity, new List<string>());

            var res = (SequencesToActualProteinPeptidesEngineResults)sequencesToActualProteinPeptidesEngine.Run();
            var compactPeptideToProteinPeptideMatching = res.CompactPeptideToProteinPeptideMatching;

            foreach (var huh in newPsms)
                if (huh != null)
                    huh.MatchToProteinLinkedPeptides(compactPeptideToProteinPeptideMatching);

            FdrAnalysisEngine engine = new FdrAnalysisEngine(newPsms, searchMode.NumNotches, CommonParameters, new List<string> { "ff" });

            engine.Run();
        }

        #endregion Public Methods
    }
}