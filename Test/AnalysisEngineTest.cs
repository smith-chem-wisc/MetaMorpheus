using EngineLayer;
using EngineLayer.Analysis;
using EngineLayer.ClassicSearch;
using IO.MzML;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Proteomics;
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

            List<PsmParent>[] newPsms = new List<PsmParent>[1];

            var proteinList = new List<Protein> { new Protein("MNNNKQQQ", "accession") };

            var protease = new Protease("Custom Protease", new List<string> { "K" }, new List<string>(), TerminusType.C, CleavageSpecificity.Full, null, null, null);

            PeptideWithPossibleModifications modPep = proteinList.First().Digest(protease, 0, null, null, InitiatorMethionineBehavior.Variable, fixedModifications).Last();
            HashSet<PeptideWithSetModifications> value1 = new HashSet<PeptideWithSetModifications> { modPep.GetPeptidesWithSetModifications(variableModifications, 4096, 3).First() };
            CompactPeptide compactPeptide1 = new CompactPeptide(value1.First(), modsDictionary);

            Assert.AreEqual("QQQ", value1.First().BaseSequence);
            PeptideWithPossibleModifications modPep2 = proteinList.First().Digest(protease, 0, null, null, InitiatorMethionineBehavior.Variable, fixedModifications).First();
            HashSet<PeptideWithSetModifications> value2 = new HashSet<PeptideWithSetModifications> { modPep2.GetPeptidesWithSetModifications(variableModifications, 4096, 3).First() };
            CompactPeptide compactPeptide2 = new CompactPeptide(value2.First(), modsDictionary);

            Assert.AreEqual("MNNNK", value2.First().BaseSequence);

            PeptideWithPossibleModifications modPep3 = proteinList.First().Digest(protease, 0, null, null, InitiatorMethionineBehavior.Variable, fixedModifications).ToList()[1];
            HashSet<PeptideWithSetModifications> value3 = new HashSet<PeptideWithSetModifications> { modPep3.GetPeptidesWithSetModifications(variableModifications, 4096, 3).First() };
            CompactPeptide compactPeptide3 = new CompactPeptide(value3.First(), modsDictionary);
            Assert.AreEqual("NNNK", value3.First().BaseSequence);

            //newPsms[0] = new List<PsmParent>[] { new List<PsmParent>{ new PsmModern(compactPeptide1, null, 1,  1, 2, 2, 1,1, 1, 1, 3,0) },
            //                                     new List<PsmParent>{  new PsmModern(compactPeptide2, null, 2,2+132.040,3,3,2,2,2,2,2,0) },
            //                                     new List<PsmParent>{ new PsmModern(compactPeptide3, null, 3, 3, 4, 3, 3, 3, 3, 3, 3, 0)} };

            IMzPeak peakA = new MzPeak(1, 1);
            IMzPeak peakB = new MzPeak(2 + 132.040, 1);
            IMzPeak peakC = new MzPeak(3, 1);

            Ms2ScanWithSpecificMass scanA = new Ms2ScanWithSpecificMass(new MzmlScanWithPrecursor(2, new MzmlMzSpectrum(new double[] { 1 }, new double[] { 1 }, false), 1, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, double.NaN, null, null, double.NaN, null, DissociationType.AnyActivationType, 1, null, null), peakA, 1, null);
            Ms2ScanWithSpecificMass scanB = new Ms2ScanWithSpecificMass(new MzmlScanWithPrecursor(3, new MzmlMzSpectrum(new double[] { 1 }, new double[] { 1 }, false), 1, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, double.NaN, null, null, double.NaN, null, DissociationType.AnyActivationType, 1, null, null), peakB, 1, null);
            Ms2ScanWithSpecificMass scanC = new Ms2ScanWithSpecificMass(new MzmlScanWithPrecursor(4, new MzmlMzSpectrum(new double[] { 1 }, new double[] { 1 }, false), 1, true, Polarity.Positive, double.NaN, null, null, MZAnalyzerType.Orbitrap, double.NaN, double.NaN, null, null, double.NaN, null, DissociationType.AnyActivationType, 1, null, null), peakC, 1, null);

            PsmParent matchA = new PsmClassic(value1.First(), 0, 0, 0, scanA);
            PsmParent matchB = new PsmClassic(value2.First(), 0, 0, 0, scanB);
            PsmParent matchC = new PsmClassic(value3.First(), 0, 0, 0, scanC);

            newPsms[0] = new List<PsmParent> { matchA, matchB, matchC };

            Tolerance fragmentTolerance = new Tolerance(ToleranceUnit.PPM, 10);
            IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> myMsDataFile = new TestDataFile(new List<PeptideWithSetModifications> { value1.First(), value2.First(), value3.First() });

            var searchModes = new List<MassDiffAcceptor> { new SinglePpmAroundZeroSearchMode(5) };
            Action<List<PsmParent>, string, List<string>> action2 = (List<PsmParent> l, string s, List<string> sdf) => {; };

            bool useProvidedPrecursorInfo = true;
            bool findAllPrecursors = true;
            var intensityRatio = 4;
            var arrayOfMs2ScansSortedByMass = MetaMorpheusEngine.GetMs2Scans(myMsDataFile, findAllPrecursors, useProvidedPrecursorInfo, intensityRatio, null).OrderBy(b => b.PrecursorMass).ToArray();

            Action<BinTreeStructure, string> action1 = (BinTreeStructure l, string s) =>
            {
                Assert.AreEqual(1, l.FinalBins.Count);
            };

            SequencesToActualProteinPeptidesEngine sequencesToActualProteinPeptidesEngine = new SequencesToActualProteinPeptidesEngine(newPsms, modsDictionary, proteinList, searchModes, protease, 2, null, null, InitiatorMethionineBehavior.Variable, fixedModifications, variableModifications, 1024, null);
            var res = (SequencesToActualProteinPeptidesEngineResults)sequencesToActualProteinPeptidesEngine.Run();
            var compactPeptideToProteinPeptideMatching = res.CompactPeptideToProteinPeptideMatching;

            FdrAnalysisEngine engine = new FdrAnalysisEngine(newPsms, searchModes, new List<string> { "ff" });

            engine.Run();
        }

        #endregion Public Methods

    }
}