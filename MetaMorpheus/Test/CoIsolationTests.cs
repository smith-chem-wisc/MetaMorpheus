using Chemistry;
using EngineLayer;
using EngineLayer.ClassicSearch;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Proteomics;
using Omics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Linq;
using Omics.Digestion;
using Omics.Modifications;
using Readers;
using TaskLayer;

namespace Test
{
    [TestFixture]
    public static class CoIsolationTests
    {
        [Test]
        public static void TestCoIsolation()
        {
            List<DigestionMotif> motifs = new List<DigestionMotif> { new DigestionMotif("K", null, 1, null) };
            Protease protease = new Protease("CustProtease", CleavageSpecificity.Full, null, null, motifs);
            ProteaseDictionary.Dictionary.Add(protease.Name, protease);
            CommonParameters CommonParameters = new CommonParameters(scoreCutoff: 1, deconvolutionIntensityRatio: 50, digestionParams: new DigestionParams(protease.Name, minPeptideLength: 1));
            var fsp = new List<(string fileName, CommonParameters fileSpecificParameters)>();
            fsp.Add(("", CommonParameters));

            var variableModifications = new List<Modification>();
            var fixedModifications = new List<Modification>();
            var proteinList = new List<Protein> { new Protein("MNNNKNDNK", null) };

            var searchModes = new SinglePpmAroundZeroSearchMode(5);

            Proteomics.AminoAcidPolymer.Peptide pep1 = new Proteomics.AminoAcidPolymer.Peptide("NNNK");
            Proteomics.AminoAcidPolymer.Peptide pep2 = new Proteomics.AminoAcidPolymer.Peptide("NDNK");

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

            var myMsDataFile = new GenericMsDataFile(Scans, null);
            
            var listOfSortedms2Scans = MetaMorpheusTask.GetMs2Scans(myMsDataFile, null, new CommonParameters(deconvolutionIntensityRatio: 50)).OrderBy(b => b.PrecursorMass).ToArray();
            
            SpectralMatch[] allPsmsArray = new PeptideSpectralMatch[listOfSortedms2Scans.Length]; ;

            bool writeSpectralLibrary = false;
            new ClassicSearchEngine(allPsmsArray, listOfSortedms2Scans, variableModifications, fixedModifications, null, null, null, 
                proteinList, searchModes, CommonParameters, fsp, null, new List<string>(), writeSpectralLibrary).Run();

            // Two matches for this single scan! Corresponding to two co-isolated masses
            Assert.That(allPsmsArray.Length, Is.EqualTo(2));

            Assert.That(allPsmsArray[0].Score > 1);
            Assert.That(allPsmsArray[0].ScanNumber, Is.EqualTo(2));

            Assert.That(allPsmsArray[0].BaseSequence, Is.EqualTo("NNNK"));
            Assert.That(allPsmsArray[1].BaseSequence, Is.EqualTo("NDNK"));
        }

        // Chimeric low-CID precursors share one spectrum but get one wrapper each, and XCorr rewrites that
        // spectrum once. Every wrapper has to re-read the count, not just the one that did the rewriting,
        // or Num Experimental Peaks reports a pre-XCorr value for the rest.
        [Test]
        public static void GetMs2ScansRefreshesPeakCountForEveryChimericPrecursor()
        {
            Proteomics.AminoAcidPolymer.Peptide pep1 = new Proteomics.AminoAcidPolymer.Peptide("NNNK");
            Proteomics.AminoAcidPolymer.Peptide pep2 = new Proteomics.AminoAcidPolymer.Peptide("NDNK");
            var dist1 = IsotopicDistribution.GetDistribution(pep1.GetChemicalFormula(), 0.1, 0.01);
            var dist2 = IsotopicDistribution.GetDistribution(pep2.GetChemicalFormula(), 0.1, 0.01);

            double[] ms1intensities = new double[] { 0.8, 0.8, 0.2, 0.02, 0.2, 0.02 };
            double[] ms1mzs = dist1.Masses.Concat(dist2.Masses).OrderBy(b => b).Select(b => b.ToMz(1)).ToArray();
            double selectedIonMz = ms1mzs[1];
            MzSpectrum ms1 = new MzSpectrum(ms1mzs, ms1intensities, false);

            // 228.486 shares a nominal mass bin with 228.086, so XCorr pre-processing collapses the two
            // and the peak count actually moves. Without that the refresh has nothing to catch.
            double[] ms2intensities = new double[] { 1, 1, 1, 1, 1, 1 };
            double[] ms2mzs = new double[] { 146.106.ToMz(1), 228.086.ToMz(1), 228.486.ToMz(1), 229.07.ToMz(1), 260.148.ToMz(1), 342.129.ToMz(1) };
            MzSpectrum ms2 = new MzSpectrum(ms2mzs, ms2intensities, false);
            int peakCountBeforeXcorr = ms2.Size;

            MsDataScan[] scans = new MsDataScan[2];
            scans[0] = new MsDataScan(ms1, 1, 1, false, Polarity.Positive, 1.0, new MzRange(300, 2000), "first spectrum",
                MZAnalyzerType.Unknown, ms1.SumOfAllY, null, null, "scan=1");
            scans[1] = new MsDataScan(ms2, 2, 2, false, Polarity.Positive, 2.0, new MzRange(100, 1500), "second spectrum",
                MZAnalyzerType.Unknown, ms2.SumOfAllY, null, null, "scan=2", selectedIonMz, null, null, selectedIonMz, 2.5,
                DissociationType.LowCID, 1, null);

            // Single-threaded on purpose: the two wrappers have to be visited in sequence, so the second
            // one sees a spectrum that is already XcorrProcessed. Let them land in separate partitions
            // and both would take the pre-processing branch, hiding the regression.
            var commonParameters = new CommonParameters(dissociationType: DissociationType.LowCID,
                deconvolutionIntensityRatio: 50, maxThreadsToUsePerFile: 1);
            var wrappers = MetaMorpheusTask.GetMs2Scans(new GenericMsDataFile(scans, null), null, commonParameters).ToList();

            // Without two wrappers over one rewritten spectrum there is no regression to catch.
            Assert.That(wrappers.Count, Is.EqualTo(2));
            Assert.That(wrappers.Select(w => w.TheScan).Distinct().Count(), Is.EqualTo(1));
            Assert.That(wrappers[0].TheScan.MassSpectrum.XcorrProcessed, Is.True);
            Assert.That(wrappers[0].TheScan.MassSpectrum.Size, Is.Not.EqualTo(peakCountBeforeXcorr));

            foreach (var wrapper in wrappers)
            {
                Assert.That(wrapper.ScanMetadata.NumPeaks, Is.EqualTo(wrapper.TheScan.MassSpectrum.Size));
            }
        }
    }
}