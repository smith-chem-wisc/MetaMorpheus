﻿using Chemistry;
using EngineLayer;
using EngineLayer.ClassicSearch;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Proteomics;
using System.Collections.Generic;
using System.Linq;
using TaskLayer;

namespace Test
{
    [TestFixture]
    public static class CoIsolationTests
    {
        #region Public Methods

        [Test]
        public static void TestCoIsolation()
        {
            Protease protease = new Protease("CustProtease", new List<string> { "K" }, new List<string>(), TerminusType.C, CleavageSpecificity.Full, null, null, null);
            GlobalVariables.ProteaseDictionary.Add(protease.Name, protease);
            CommonParameters CommonParameters = new CommonParameters
            {
                DigestionParams = new DigestionParams(protease.Name, MinPeptideLength: 1),
                ConserveMemory = false,
                ScoreCutoff = 1,
                DeconvolutionIntensityRatio = 50
            };

            var variableModifications = new List<ModificationWithMass>();
            var fixedModifications = new List<ModificationWithMass>();
            var proteinList = new List<Protein> { new Protein("MNNNKNDNK", null) };

            var searchModes = new SinglePpmAroundZeroSearchMode(5);

            Proteomics.Peptide pep1 = new Proteomics.Peptide("NNNK");
            Proteomics.Peptide pep2 = new Proteomics.Peptide("NDNK");

            var dist1 = IsotopicDistribution.GetDistribution(pep1.GetChemicalFormula(), 0.1, 0.01);

            var dist2 = IsotopicDistribution.GetDistribution(pep2.GetChemicalFormula(), 0.1, 0.01);

            MsDataScan[] Scans = new MsDataScan[2];
            double[] ms1intensities = new double[] { 0.8, 0.8, 0.2, 0.02, 0.2, 0.02 };
            double[] ms1mzs = dist1.Masses.Concat(dist2.Masses).OrderBy(b => b).Select(b => b.ToMz(1)).ToArray();

            double selectedIonMz = ms1mzs[1];

            MzSpectrum MS1 = new MzSpectrum(ms1mzs, ms1intensities, false);

            Scans[0] = new MsDataScan(MS1, 1, 1, false, Polarity.Positive, 1.0, new MzRange(300, 2000), "first spectrum", MZAnalyzerType.Unknown, MS1.SumOfAllY, null, null, "scan=1");

            double[] ms2intensities = new double[] { 1, 1, 1, 1, 1 };
            double[] ms2mzs = new double[] { 146.106.ToMz(1), 228.086.ToMz(1), 229.07.ToMz(1), 260.148.ToMz(1), 342.129.ToMz(1) };
            MzSpectrum MS2 = new MzSpectrum(ms2mzs, ms2intensities, false);
            double isolationMZ = selectedIonMz;
            Scans[1] = new MsDataScan(MS2, 2, 2, false, Polarity.Positive, 2.0, new MzRange(100, 1500), "second spectrum", MZAnalyzerType.Unknown, MS2.SumOfAllY, null, null, "scan=2", selectedIonMz, null, null, isolationMZ, 2.5, DissociationType.HCD, 1, null);

            var myMsDataFile = new MsDataFile(Scans, null);

            bool DoPrecursorDeconvolution = true;
            bool UseProvidedPrecursorInfo = true;
            double DeconvolutionIntensityRatio = 50;
            int DeconvolutionMaxAssumedChargeState = 10;
            Tolerance DeconvolutionMassTolerance = new PpmTolerance(5);

            var listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(myMsDataFile, null, DoPrecursorDeconvolution, UseProvidedPrecursorInfo, DeconvolutionIntensityRatio, DeconvolutionMaxAssumedChargeState, DeconvolutionMassTolerance).OrderBy(b => b.PrecursorMass).ToArray();

            PeptideSpectralMatch[] allPsmsArray = new PeptideSpectralMatch[listOfSortedms2Scans.Length];

            List<ProductType> lp = new List<ProductType> { ProductType.B, ProductType.Y };
            new ClassicSearchEngine(allPsmsArray, listOfSortedms2Scans, variableModifications, fixedModifications, proteinList, lp, searchModes, false, CommonParameters, new List<string>()).Run();

            // Two matches for this single scan! Corresponding to two co-isolated masses
            Assert.AreEqual(2, allPsmsArray.Length);

            Assert.IsTrue(allPsmsArray[0].Score > 1);
            Assert.AreEqual(2, allPsmsArray[0].ScanNumber);

            var ojdfkj = (SequencesToActualProteinPeptidesEngineResults)new SequencesToActualProteinPeptidesEngine(new List<PeptideSpectralMatch>
            { allPsmsArray[0], allPsmsArray[1] }, proteinList, fixedModifications, variableModifications, lp, new List<DigestionParams> { CommonParameters.DigestionParams }, CommonParameters.ReportAllAmbiguity, new List<string>()).Run();

            foreach (var huh in allPsmsArray)
            {
                if (huh != null)
                {
                    huh.MatchToProteinLinkedPeptides(ojdfkj.CompactPeptideToProteinPeptideMatching);
                }
            }

            Assert.AreEqual("NNNK", allPsmsArray[0].BaseSequence);
            Assert.AreEqual("NDNK", allPsmsArray[1].BaseSequence);
        }

        #endregion Public Methods
    }
}