using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Text.RegularExpressions;
using System.Threading.Tasks;
using System.Windows.Controls;
using Chemistry;
using Easy.Common.Extensions;
using EngineLayer;
using EngineLayer.PsmTsv;
using FlashLFQ;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Proteomics;
using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using TaskLayer;
using EngineLayer.SpectralRecovery;

namespace Test
{
    [TestFixture]
    internal class GenericReaderTest
    {
        public string SpectrumFilePath;
        public Dictionary<string, PsmFromTsv> MaxQuantPsms;
        public Dictionary<string, List<ChromatographicPeak>> MbrPeaks;
        public List<SpectraFileInfo> SpectraFiles;

        [Test]
        public static void TestHeLaSlice()
        {
            List<string> rawSlices = new List<string> {
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", @"SpectralRecoveryTest\K13_02ng_1min_frac1.mzML"),
                Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", @"SpectralRecoveryTest\K13_20ng_1min_frac1.mzML") };
            string experimentalDesignFilepath = Path.Combine(TestContext.CurrentContext.TestDirectory,
                "TestData", @"SpectralRecoveryTest\ExperimentalDesign.tsv");

            List<SpectraFileInfo> spectraFiles = ExperimentalDesign.ReadExperimentalDesign(
                experimentalDesignFilepath, rawSlices, out var errors);

            Assert.That(spectraFiles.Count, Is.EqualTo(2));

            string maxQuantEvidencePath = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", @"SpectralRecoveryTest\K13_slice_evidence.txt");
            Dictionary<string, List<ChromatographicPeak>> mbrPeaks = PsmGenericReader.ReadInMbrPeaks(
                maxQuantEvidencePath, silent: false, spectraFiles);

            string maxQuantMsmsPath = Path.Combine(TestContext.CurrentContext.TestDirectory, "TestData", @"SpectralRecoveryTest\K13_slice_msms.txt");
            Dictionary<string, PsmFromTsv> maxQuantPsms = PsmGenericReader.GetDonorPsms(
                maxQuantMsmsPath, spectraFiles, mbrPeaks, ignoreArtifactIons: false);

            Assert.That(mbrPeaks.Count == 22);
            Assert.That(maxQuantPsms.Count == 14);

            PpmTolerance testTolerance = new PpmTolerance(5);
            // Check that the PeptideMonoMass (derived from msms.txt mass columns) and the pwsm mass (derived from converted full sequence)
            // matches for every psm
            Assert.IsFalse(maxQuantPsms.Values.Any(p =>
                !testTolerance.Within(double.Parse(p.PeptideMonoMass), p.PeptideWithSetModifications.MonoisotopicMass)));

            string outputFolder = Path.Combine(TestContext.CurrentContext.TestDirectory, @"TestSpectralRecoveryOutput");
            Directory.CreateDirectory(outputFolder);
            // Writing a spectral library
            var spectraForSpectraLibrary = new List<LibrarySpectrum>();
            foreach (var psm in maxQuantPsms.Values)
            {
                var standardSpectrum = new LibrarySpectrum(psm.FullSequence, psm.PrecursorMz, psm.PrecursorCharge, psm.MatchedIons, psm.RetentionTime ?? 0);
                spectraForSpectraLibrary.Add(standardSpectrum);
            }
            string spectrumFilePath =Path.Combine(outputFolder, "SpectralLibrary.msp");
            if (File.Exists(spectrumFilePath)) File.Delete(spectrumFilePath);
            using (StreamWriter output = new StreamWriter(spectrumFilePath))
            {
                foreach (var librarySpectrum in spectraForSpectraLibrary)
                {
                    output.WriteLine(librarySpectrum.ToString());
                }
            }

            // Tolerances taken from MaxQuant defaults
            CommonParameters commonParams = new CommonParameters(dissociationType: DissociationType.Autodetect,
                productMassTolerance: new PpmTolerance(20), deconvolutionMassTolerance: new PpmTolerance(7));

            var mbrAnalysisResults = SpectralRecoveryRunner.RunSpectralRecoveryFromMaxQuant(
                spectraFiles,
                mbrPeaks,
                maxQuantPsms,
                spectrumFilePath,
                outputFolder,
                commonParams);

            mbrAnalysisResults.WritePeakQuantificationResultsToTsv(outputFolder, "PeakQuant_NoArtifact");
        }

        [Test]
        public static void TestMaxQuantEvidenceReaderHeLa()
        {
            string experimentalDesignFilepath = @"D:\HelaSingleCellQCmzML\rawFiles\ExperimentalDesign.tsv";
            List<string> rawFilePathList = new List<string>
            {
                @"D:\HelaSingleCellQCmzML\rawFiles\Ex_Auto_J3_30umTB_2ngQC_60m_2.raw",
                @"D:\HelaSingleCellQCmzML\rawFiles\Ex_Auto_K13_30umTA_02ngQC_60m_1.raw",
                @"D:\HelaSingleCellQCmzML\rawFiles\Ex_Auto_K13_30umTA_02ngQC_60m_2.raw",
                @"D:\HelaSingleCellQCmzML\rawFiles\Ex_Auto_K13_30umTA_2ngQC_60m_1.raw",
                @"D:\HelaSingleCellQCmzML\rawFiles\Ex_Auto_K13_30umTA_2ngQC_60m_2.raw",
                @"D:\HelaSingleCellQCmzML\rawFiles\Ex_Auto_W17_30umTA_02ngQC_60m_3.raw",
                @"D:\HelaSingleCellQCmzML\rawFiles\Ex_Auto_W17_30umTA_02ngQC_60m_4.raw",
                @"D:\HelaSingleCellQCmzML\rawFiles\Ex_Auto_W17_30umTB_2ngQC_60m_1.raw",
                @"D:\HelaSingleCellQCmzML\rawFiles\Ex_Auto_W17_30umTB_2ngQC_60m_2.raw",
                @"D:\HelaSingleCellQCmzML\rawFiles\Ex_Auto_J3_30umTB_02ngQC_60m_1.raw",
                @"D:\HelaSingleCellQCmzML\rawFiles\Ex_Auto_J3_30umTB_02ngQC_60m_2.raw",
                @"D:\HelaSingleCellQCmzML\rawFiles\Ex_Auto_J3_30umTB_2ngQC_60m_1.raw"
            };
            List<SpectraFileInfo> spectraFiles = ExperimentalDesign.ReadExperimentalDesign(
                experimentalDesignFilepath, rawFilePathList, out var errors);

            Assert.That(spectraFiles.Count, Is.EqualTo(12));

            string maxQuantEvidencePath = @"D:\HelaSingleCellQCmzML\rawFiles\combined\txt\evidence.txt";
            Dictionary<string, List<ChromatographicPeak>> mbrPeaks = PsmGenericReader.ReadInMbrPeaks(
                maxQuantEvidencePath, silent: false, spectraFiles);

            string maxQuantMsmsPath = @"D:\HelaSingleCellQCmzML\rawFiles\combined\txt\msms.txt";
            Dictionary<string, PsmFromTsv> maxQuantPsms = PsmGenericReader.GetDonorPsms(
                maxQuantMsmsPath, spectraFiles, mbrPeaks, ignoreArtifactIons: false);

            // Every MBR run ID should have a corresponding PSM in the msms.txt file
            Assert.That(mbrPeaks.Count == maxQuantPsms.Count);

            PpmTolerance testTolerance = new PpmTolerance(5);
            // Check that the PeptideMonoMass (derived from msms.txt mass columns) and the pwsm mass (derived from converted full sequence)
            // matches for every psm
            Assert.IsFalse(maxQuantPsms.Values.Any(p => 
                !testTolerance.Within(double.Parse(p.PeptideMonoMass), p.PeptideWithSetModifications.MonoisotopicMass)));

            // Writing a spectral library
            var spectraForSpectraLibrary = new List<LibrarySpectrum>();
            foreach (var psm in maxQuantPsms.Values)
            {
                var standardSpectrum = new LibrarySpectrum(psm.FullSequence, psm.PrecursorMz, psm.PrecursorCharge, psm.MatchedIons, psm.RetentionTime ?? 0);
                spectraForSpectraLibrary.Add(standardSpectrum);
            }
            string spectrumFilePath = @"D:\HelaSingleCellQCmzML\rawFiles\combined\txt\spectralLibrary.msp";
            if(File.Exists(spectrumFilePath)) File.Delete(spectrumFilePath);
            using (StreamWriter output = new StreamWriter(spectrumFilePath))
            {
                foreach (var librarySpectrum in spectraForSpectraLibrary)
                {
                    output.WriteLine(librarySpectrum.ToString());
                }
            }

            // Tolerances taken from MaxQuant defaults
            CommonParameters commonParams = new CommonParameters(dissociationType: DissociationType.Autodetect, 
                productMassTolerance: new PpmTolerance(20), deconvolutionMassTolerance: new PpmTolerance(7));

            string resultsPath = @"D:\SingleCellDataSets\MSV000087524\SpectralRecoveryMaxQuant";
            Directory.CreateDirectory(resultsPath);

            var mbrAnalysisResults = SpectralRecoveryRunner.RunSpectralRecoveryFromMaxQuant(
                spectraFiles,
                mbrPeaks,
                maxQuantPsms,
                spectrumFilePath,
                resultsPath,
                commonParams);

            mbrAnalysisResults.WritePeakQuantificationResultsToTsv(resultsPath, "PeakQuant_NoArtifact");
        }


        [Test]
        public static void TestMaxQuantEvidenceReader955()
        {
            string folderPath = @"D:\SingleCellDataSets\PXD031955";
            string experimentalDesignFilepath = @"D:\SingleCellDataSets\PXD031955\ExperimentalDesign_PXD031955.tsv";
            List<string> rawFilePathList = new List<string>
            {
                @"D:\SingleCellDataSets\PXD031955\2016-04-11_SC05_hip_cultured_neu0_125_ugul.raw",
                @"D:\SingleCellDataSets\PXD031955\2016-04-11_SC06_hip_cultured_neu0_25_ugul.raw",
                @"D:\SingleCellDataSets\PXD031955\2016-04-11_SC07_hip_cultured_neu0_25_ugul.raw",
                @"D:\SingleCellDataSets\PXD031955\2016-04-01_SC03_hip_cultured_neu0_5_ugul.raw",
                @"D:\SingleCellDataSets\PXD031955\2016-04-04_SC02_hip_cultured_neu0_5_ugul.raw",
                @"D:\SingleCellDataSets\PXD031955\2016-04-11_SC04_hip_cultured_neu0_125_ugul.raw",
            };
            List<SpectraFileInfo> spectraFiles = ExperimentalDesign.ReadExperimentalDesign(
                experimentalDesignFilepath, rawFilePathList, out var errors);

            string maxQuantEvidencePath = @"D:\SingleCellDataSets\PXD031955\combined\txt\evidence.txt";
            Dictionary<string, List<ChromatographicPeak>> mbrPeaks = PsmGenericReader.ReadInMbrPeaks(
                maxQuantEvidencePath, silent: false, spectraFiles);

            string maxQuantMsmsPath = @"D:\SingleCellDataSets\PXD031955\combined\txt\msms.txt";
            Dictionary<string, PsmFromTsv> maxQuantPsms = PsmGenericReader.GetDonorPsms(
                maxQuantMsmsPath, spectraFiles, mbrPeaks, ignoreArtifactIons: true);

            // Every MBR run ID should have a corresponding PSM in the msms.txt file
            Assert.That(mbrPeaks.Count == maxQuantPsms.Count);

            PpmTolerance testTolerance = new PpmTolerance(5);
            // Check that the PeptideMonoMass (derived from msms.txt mass columns) and the pwsm mass (derived from converted full sequence)
            // matches for every psm
            Assert.IsFalse(maxQuantPsms.Values.Any(p =>
                !testTolerance.Within(double.Parse(p.PeptideMonoMass), p.PeptideWithSetModifications.MonoisotopicMass)));

            // Writing a spectral library
            var spectraForSpectraLibrary = new List<LibrarySpectrum>();
            foreach (var psm in maxQuantPsms.Values)
            {
                var standardSpectrum = new LibrarySpectrum(psm.FullSequence, psm.PrecursorMz, psm.PrecursorCharge, psm.MatchedIons, psm.RetentionTime ?? 0);
                spectraForSpectraLibrary.Add(standardSpectrum);
            }
            string spectrumFilePath = @"D:\SingleCellDataSets\PXD031955\combined\txt\spectralLibrary.msp";
            if (File.Exists(spectrumFilePath)) File.Delete(spectrumFilePath);
            using (StreamWriter output = new StreamWriter(spectrumFilePath))
            {
                foreach (var librarySpectrum in spectraForSpectraLibrary)
                {
                    output.WriteLine(librarySpectrum.ToString());
                }
            }

            // Tolerances taken from MaxQuant defaults
            CommonParameters commonParams = new CommonParameters(dissociationType: DissociationType.Autodetect,
                productMassTolerance: new PpmTolerance(20), deconvolutionMassTolerance: new PpmTolerance(7));

            string resultsPath = Path.Combine(folderPath, "SpectralRecoveryMaxQuant");
            Directory.CreateDirectory(resultsPath);

            var mbrAnalysisResults = SpectralRecoveryRunner.RunSpectralRecoveryFromMaxQuant(
                spectraFiles,
                mbrPeaks,
                maxQuantPsms,
                spectrumFilePath,
                resultsPath + @"\",
                commonParams);

            mbrAnalysisResults.WritePeakQuantificationResultsToTsv(resultsPath + @"\", "PeakQuant_NoArtifact");
        }


        [Test]
        public static void TestMaxQuantEvidenceReader24017()
        {
            string folderPath = @"D:\SingleCellDataSets\PXD024017";
            string experimentalDesignFilepath = @"D:\SingleCellDataSets\PXD024017\ExperimentalDesign_PXD024017_raw.tsv";
            List<string> rawFilePathList = new List<string>
            {
                @"D:\SingleCellDataSets\PXD024017\20201207_Exploris_RSLC9_PepSep_PHF-3um_HeLaPi_500pg_03.raw",
                @"D:\SingleCellDataSets\PXD024017\20201207_Exploris_RSLC9_PepSep_PHF-3um_HeLaPi_500pg_02.raw",
                @"D:\SingleCellDataSets\PXD024017\20201207_Exploris_RSLC9_PepSep_PHF-3um_HeLaPi_500pg_01.raw",
                @"D:\SingleCellDataSets\PXD024017\20201207_Exploris_RSLC9_PepSep_PHF-3um_HeLaPi_2c5ng_04.raw",
                @"D:\SingleCellDataSets\PXD024017\20201207_Exploris_RSLC9_PepSep_PHF-3um_HeLaPi_2c5ng_03.raw",
                @"D:\SingleCellDataSets\PXD024017\20201207_Exploris_RSLC9_PepSep_PHF-3um_HeLaPi_2c5ng_02.raw",
                @"D:\SingleCellDataSets\PXD024017\20201207_Exploris_RSLC9_PepSep_PHF-3um_HeLaPi_2c5ng_01.raw",
                @"D:\SingleCellDataSets\PXD024017\20201207_Exploris_RSLC9_PepSep_PHF-3um_HeLaPi_250pg_04.raw",
                @"D:\SingleCellDataSets\PXD024017\20201207_Exploris_RSLC9_PepSep_PHF-3um_HeLaPi_250pg_03.raw",
                @"D:\SingleCellDataSets\PXD024017\20201207_Exploris_RSLC9_PepSep_PHF-3um_HeLaPi_250pg_02.raw",
                @"D:\SingleCellDataSets\PXD024017\20201207_Exploris_RSLC9_PepSep_PHF-3um_HeLaPi_250g_01.raw",
                @"D:\SingleCellDataSets\PXD024017\20201207_Exploris_RSLC9_PepSep_PHF-3um_HeLaPi_1ng_04.raw",
                @"D:\SingleCellDataSets\PXD024017\20201207_Exploris_RSLC9_PepSep_PHF-3um_HeLaPi_1ng_03.raw",
                @"D:\SingleCellDataSets\PXD024017\20201207_Exploris_RSLC9_PepSep_PHF-3um_HeLaPi_1ng_01.raw",
                @"D:\SingleCellDataSets\PXD024017\20201207_Exploris_RSLC9_PepSep_PHF-3um_HeLaPi_1ng_02.raw",
                @"D:\SingleCellDataSets\PXD024017\20201207_Exploris_RSLC9_PepSep_PHF-3um_HeLaPi_5ng_04.raw",
                @"D:\SingleCellDataSets\PXD024017\20201207_Exploris_RSLC9_PepSep_PHF-3um_HeLaPi_5ng_03.raw",
                @"D:\SingleCellDataSets\PXD024017\20201207_Exploris_RSLC9_PepSep_PHF-3um_HeLaPi_5ng_02.raw",
                @"D:\SingleCellDataSets\PXD024017\20201207_Exploris_RSLC9_PepSep_PHF-3um_HeLaPi_5ng_01.raw",
                @"D:\SingleCellDataSets\PXD024017\20201207_Exploris_RSLC9_PepSep_PHF-3um_HeLaPi_500pg_04.raw"
            };
            List<SpectraFileInfo> spectraFiles = ExperimentalDesign.ReadExperimentalDesign(
                experimentalDesignFilepath, rawFilePathList, out var errors);

            string maxQuantEvidencePath = Path.Join(folderPath, @"combined\txt\evidence.txt");
            Dictionary<string, List<ChromatographicPeak>> mbrPeaks = PsmGenericReader.ReadInMbrPeaks(
                maxQuantEvidencePath, silent: false, spectraFiles);

            string maxQuantMsmsPath = Path.Join(folderPath, @"combined\txt\msms.txt");
            Dictionary<string, PsmFromTsv> maxQuantPsms = PsmGenericReader.GetDonorPsms(
                maxQuantMsmsPath, spectraFiles, mbrPeaks, ignoreArtifactIons: true);

            // Every MBR run ID should have a corresponding PSM in the msms.txt file
            Assert.AreEqual(mbrPeaks.Count, maxQuantPsms.Count, 1); // There is one donor where the matched fragments are missing

            PpmTolerance testTolerance = new PpmTolerance(5);
            // Check that the PeptideMonoMass (derived from msms.txt mass columns) and the pwsm mass (derived from converted full sequence)
            // matches for every psm
            Assert.IsFalse(maxQuantPsms.Values.Any(p =>
                !testTolerance.Within(double.Parse(p.PeptideMonoMass), p.PeptideWithSetModifications.MonoisotopicMass)));

            // Writing a spectral library
            var spectraForSpectraLibrary = new List<LibrarySpectrum>();
            foreach (var psm in maxQuantPsms.Values)
            {
                var standardSpectrum = new LibrarySpectrum(psm.FullSequence, psm.PrecursorMz, psm.PrecursorCharge, psm.MatchedIons, psm.RetentionTime ?? 0);
                spectraForSpectraLibrary.Add(standardSpectrum);
            }
            string spectrumFilePath = Path.Join(folderPath, @"combined\txt\spectralLibrary.msp");
            if (File.Exists(spectrumFilePath)) File.Delete(spectrumFilePath);
            using (StreamWriter output = new StreamWriter(spectrumFilePath))
            {
                foreach (var librarySpectrum in spectraForSpectraLibrary)
                {
                    output.WriteLine(librarySpectrum.ToString());
                }
            }

            // Tolerances taken from MaxQuant defaults
            CommonParameters commonParams = new CommonParameters(dissociationType: DissociationType.Autodetect,
                productMassTolerance: new PpmTolerance(20), deconvolutionMassTolerance: new PpmTolerance(7));

            string resultsPath = Path.Combine(folderPath, "SpectralRecoveryMaxQuant");
            Directory.CreateDirectory(resultsPath);

            var mbrAnalysisResults = SpectralRecoveryRunner.RunSpectralRecoveryFromMaxQuant(
                spectraFiles,
                mbrPeaks,
                maxQuantPsms,
                spectrumFilePath,
                resultsPath + @"\",
                commonParams);

            mbrAnalysisResults.WritePeakQuantificationResultsToTsv(resultsPath, "PeakQuant_NoArtifact");
        }

        [Test]
        public static void TestMaxQuantEvidenceReader19515()
        {
            string folderPath = @"D:\SingleCellDataSets\PXD019515";
            string experimentalDesignFilepath = @"D:\SingleCellDataSets\PXD019515\ExperimentalDesign_PXD019515.tsv";
            List<string> rawFilePathList = new List<string>
            {
                @"D:\SingleCellDataSets\PXD019515\FAIMS_2CV_OTIT_HCD_300ITMS2_SingleMotorNeuron_2.raw",
                @"D:\SingleCellDataSets\PXD019515\FAIMS_2CV_OTIT_HCD_300ITMS2_SingleMotorNeuron_3.raw",
                @"D:\SingleCellDataSets\PXD019515\FAIMS_2CV_OTIT_HCD_300ITMS2_Single_HeLa_1.raw",
                @"D:\SingleCellDataSets\PXD019515\FAIMS_2CV_OTIT_HCD_300ITMS2_Single_HeLa_2.raw",
                @"D:\SingleCellDataSets\PXD019515\FAIMS_2CV_OTIT_HCD_300ITMS2_Single_HeLa_3.raw",
                @"D:\SingleCellDataSets\PXD019515\FAIMS_2CV_OTIT_HCD_300ITMS2_SingleInterNeuron_1.raw",
                @"D:\SingleCellDataSets\PXD019515\FAIMS_2CV_OTIT_HCD_300ITMS2_SingleInterNeuron_2.raw",
                @"D:\SingleCellDataSets\PXD019515\FAIMS_2CV_OTIT_HCD_300ITMS2_SingleInterNeuron_3.raw",
                @"D:\SingleCellDataSets\PXD019515\FAIMS_2CV_OTIT_HCD_300ITMS2_SingleMotorNeuron_1.raw"
            };
            List<SpectraFileInfo> spectraFiles = ExperimentalDesign.ReadExperimentalDesign(
                experimentalDesignFilepath, rawFilePathList, out var errors);

            string maxQuantEvidencePath = Path.Join(folderPath, @"combined\txt\evidence.txt");
            Dictionary<string, List<ChromatographicPeak>> mbrPeaks = PsmGenericReader.ReadInMbrPeaks(
                maxQuantEvidencePath, silent: false, spectraFiles);

            string maxQuantMsmsPath = Path.Join(folderPath, @"combined\txt\msms.txt");
            Dictionary<string, PsmFromTsv> maxQuantPsms = PsmGenericReader.GetDonorPsms(
                maxQuantMsmsPath, spectraFiles, mbrPeaks, ignoreArtifactIons: true);

            // Every MBR run ID should have a corresponding PSM in the msms.txt file
            Assert.AreEqual(mbrPeaks.Count, maxQuantPsms.Count, 1); // There is one donor where the matched fragments are missing

            PpmTolerance testTolerance = new PpmTolerance(5);
            // Check that the PeptideMonoMass (derived from msms.txt mass columns) and the pwsm mass (derived from converted full sequence)
            // matches for every psm
            Assert.IsFalse(maxQuantPsms.Values.Any(p =>
                !testTolerance.Within(double.Parse(p.PeptideMonoMass), p.PeptideWithSetModifications.MonoisotopicMass)));

            // Writing a spectral library
            var spectraForSpectraLibrary = new List<LibrarySpectrum>();
            foreach (var psm in maxQuantPsms.Values)
            {
                var standardSpectrum = new LibrarySpectrum(psm.FullSequence, psm.PrecursorMz, psm.PrecursorCharge, psm.MatchedIons, psm.RetentionTime ?? 0);
                spectraForSpectraLibrary.Add(standardSpectrum);
            }
            string spectrumFilePath = Path.Join(folderPath, @"combined\txt\spectralLibrary.msp");
            if (File.Exists(spectrumFilePath)) File.Delete(spectrumFilePath);
            using (StreamWriter output = new StreamWriter(spectrumFilePath))
            {
                foreach (var librarySpectrum in spectraForSpectraLibrary)
                {
                    output.WriteLine(librarySpectrum.ToString());
                }
            }

            // Tolerances taken from MaxQuant defaults
            CommonParameters commonParams = new CommonParameters(dissociationType: DissociationType.Autodetect,
                productMassTolerance: new PpmTolerance(20), deconvolutionMassTolerance: new PpmTolerance(7));

            string resultsPath = Path.Combine(folderPath, "SpectralRecoveryMaxQuant");
            Directory.CreateDirectory(resultsPath);

            var mbrAnalysisResults = SpectralRecoveryRunner.RunSpectralRecoveryFromMaxQuant(
                spectraFiles,
                mbrPeaks,
                maxQuantPsms,
                spectrumFilePath,
                resultsPath,
                commonParams);

            mbrAnalysisResults.WritePeakQuantificationResultsToTsv(resultsPath, "PeakQuant_NoArtifact");
        }

        [Test]
        public static void TestMaxQuantEvidenceReader85937_HeLa()
        {
            string folderPath = @"D:\SingleCellDataSets\MSV000085937\HeLa_SingleCellProteomics";
            string experimentalDesignFilepath = @"D:\SingleCellDataSets\MSV000085937\HeLa_SingleCellProteomics\ExperimentalDesign_MSV000085937_HeLa_SCP.tsv";
            List<string> rawFilePathList = new List<string>
            {
                @"D:\SingleCellDataSets\MSV000085937\HeLa_SingleCellProteomics\NanoPOTS_FAIMS_HeLa_Lib_Chip1_A4.raw",
                @"D:\SingleCellDataSets\MSV000085937\HeLa_SingleCellProteomics\NanoPOTS_FAIMS_HeLa_Lib_Chip1_A5.raw",
                @"D:\SingleCellDataSets\MSV000085937\HeLa_SingleCellProteomics\NanoPOTS_FAIMS_HeLa_SC_Chip1_A11.raw",
                @"D:\SingleCellDataSets\MSV000085937\HeLa_SingleCellProteomics\NanoPOTS_FAIMS_HeLa_SC_Chip1_B1.raw",
                @"D:\SingleCellDataSets\MSV000085937\HeLa_SingleCellProteomics\NanoPOTS_FAIMS_HeLa_SC_Chip1_B3.raw",
                @"D:\SingleCellDataSets\MSV000085937\HeLa_SingleCellProteomics\NanoPOTS_FAIMS_HeLa_SC_Chip1_B4.raw",
                @"D:\SingleCellDataSets\MSV000085937\HeLa_SingleCellProteomics\NanoPOTS_FAIMS_HeLa_SC_Chip1_B5.raw",
                @"D:\SingleCellDataSets\MSV000085937\HeLa_SingleCellProteomics\NanoPOTS_FAIMS_HeLa_SC_Chip1_B6.raw",
                @"D:\SingleCellDataSets\MSV000085937\HeLa_SingleCellProteomics\NanoPOTS_FAIMS_HeLa_SC_Chip1_B7.raw",
                @"D:\SingleCellDataSets\MSV000085937\HeLa_SingleCellProteomics\NanoPOTS_FAIMS_HeLa_SC_Chip1_B8.raw",
                @"D:\SingleCellDataSets\MSV000085937\HeLa_SingleCellProteomics\NanoPOTS_FAIMS_HeLa_SC_Chip1_B10.raw",
                @"D:\SingleCellDataSets\MSV000085937\HeLa_SingleCellProteomics\NanoPOTS_FAIMS_HeLa_SC_Chip1_B11.raw",
                @"D:\SingleCellDataSets\MSV000085937\HeLa_SingleCellProteomics\NanoPOTS_FAIMS_HeLa_Lib_Chip1_A1.raw",
                @"D:\SingleCellDataSets\MSV000085937\HeLa_SingleCellProteomics\NanoPOTS_FAIMS_HeLa_Lib_Chip1_A2.raw",
                @"D:\SingleCellDataSets\MSV000085937\HeLa_SingleCellProteomics\NanoPOTS_FAIMS_HeLa_Lib_Chip1_A3.raw"
            };
            List<SpectraFileInfo> spectraFiles = ExperimentalDesign.ReadExperimentalDesign(
                experimentalDesignFilepath, rawFilePathList, out var errors);

            string maxQuantEvidencePath = Path.Join(folderPath, @"combined\txt\evidence.txt");
            Dictionary<string, List<ChromatographicPeak>> mbrPeaks = PsmGenericReader.ReadInMbrPeaks(
                maxQuantEvidencePath, silent: false, spectraFiles);

            string maxQuantMsmsPath = Path.Join(folderPath, @"combined\txt\msms.txt");
            Dictionary<string, PsmFromTsv> maxQuantPsms = PsmGenericReader.GetDonorPsms(
                maxQuantMsmsPath, spectraFiles, mbrPeaks, ignoreArtifactIons: true);

            // Every MBR run ID should have a corresponding PSM in the msms.txt file
            Assert.AreEqual(mbrPeaks.Count, maxQuantPsms.Count, 1); // There is one donor where the matched fragments are missing

            PpmTolerance testTolerance = new PpmTolerance(5);
            // Check that the PeptideMonoMass (derived from msms.txt mass columns) and the pwsm mass (derived from converted full sequence)
            // matches for every psm
            Assert.IsFalse(maxQuantPsms.Values.Any(p =>
                !testTolerance.Within(double.Parse(p.PeptideMonoMass), p.PeptideWithSetModifications.MonoisotopicMass)));

            // Writing a spectral library
            var spectraForSpectraLibrary = new List<LibrarySpectrum>();
            foreach (var psm in maxQuantPsms.Values)
            {
                var standardSpectrum = new LibrarySpectrum(psm.FullSequence, psm.PrecursorMz, psm.PrecursorCharge, psm.MatchedIons, psm.RetentionTime ?? 0);
                spectraForSpectraLibrary.Add(standardSpectrum);
            }
            string spectrumFilePath = Path.Join(folderPath, @"combined\txt\spectralLibrary.msp");
            if (File.Exists(spectrumFilePath)) File.Delete(spectrumFilePath);
            using (StreamWriter output = new StreamWriter(spectrumFilePath))
            {
                foreach (var librarySpectrum in spectraForSpectraLibrary)
                {
                    output.WriteLine(librarySpectrum.ToString());
                }
            }

            // Tolerances taken from MaxQuant defaults
            CommonParameters commonParams = new CommonParameters(dissociationType: DissociationType.Autodetect,
                productMassTolerance: new PpmTolerance(20), deconvolutionMassTolerance: new PpmTolerance(7));

            //Dictionary<SpectraFileInfo, string> calibFileDictionary = new Dictionary<SpectraFileInfo, string>();
            //foreach(SpectraFileInfo spectraFile in spectraFiles)
            //{
            //    string calibFilePath = Path.Join(@"D:\SingleCellDataSets\MSV000085937\HeLa_SingleCellProteomics\Calibration\Task1-CalibrateTask",
            //        spectraFile.FilenameWithoutExtension + "-calib.mzML");
            //    calibFileDictionary.Add(spectraFile, calibFilePath);
            //}

            string resultsPath = Path.Combine(folderPath, "SpectralRecoveryMaxQuant");
            Directory.CreateDirectory(resultsPath);

            var mbrAnalysisResults = SpectralRecoveryRunner.RunSpectralRecoveryFromMaxQuant(
                spectraFiles,
                mbrPeaks,
                maxQuantPsms,
                spectrumFilePath,
                resultsPath + @"\",
                commonParams);

            mbrAnalysisResults.WritePeakQuantificationResultsToTsv(resultsPath, "PeakQuant_NoArtifact");
        }

        [Test]
        public static void TestMaxQuantEvidenceReader85937_Lung()
        {
            string folderPath = @"D:\SingleCellDataSets\MSV000085937\NonDepleted_Lung_SingleCellProteomics";
            string experimentalDesignFilepath = @"D:\SingleCellDataSets\MSV000085937\NonDepleted_Lung_SingleCellProteomics\ExperimentalDesign_MSV000085937_NonDepletedLung_SCP.tsv";
            List<string> rawFilePathList = new List<string>
            {
                @"D:\SingleCellDataSets\MSV000085937\NonDepleted_Lung_SingleCellProteomics\NanoPOTS_D022_Primary_SC_Non_Resting_Chip1_A7.raw",
                @"D:\SingleCellDataSets\MSV000085937\NonDepleted_Lung_SingleCellProteomics\NanoPOTS_D022_Primary_SC_Non_Resting_Chip1_A8.raw",
                @"D:\SingleCellDataSets\MSV000085937\NonDepleted_Lung_SingleCellProteomics\NanoPOTS_D022_Primary_SC_Non_Resting_Chip1_A10.raw",
                @"D:\SingleCellDataSets\MSV000085937\NonDepleted_Lung_SingleCellProteomics\NanoPOTS_D022_Primary_SC_Non_Resting_Chip1_A11.raw",
                @"D:\SingleCellDataSets\MSV000085937\NonDepleted_Lung_SingleCellProteomics\NanoPOTS_D022_Primary_SC_Non_Resting_Chip1_B02.raw",
                @"D:\SingleCellDataSets\MSV000085937\NonDepleted_Lung_SingleCellProteomics\NanoPOTS_D022_Primary_SC_Resting_Chip1_C1.raw",
                @"D:\SingleCellDataSets\MSV000085937\NonDepleted_Lung_SingleCellProteomics\NanoPOTS_D022_Primary_SC_Resting_Chip1_C3.raw",
                @"D:\SingleCellDataSets\MSV000085937\NonDepleted_Lung_SingleCellProteomics\NanoPOTS_D022_Primary_SC_Resting_Chip1_C4.raw",
                @"D:\SingleCellDataSets\MSV000085937\NonDepleted_Lung_SingleCellProteomics\NanoPOTS_D022_Primary_SC_Resting_Chip1_C5.raw",
                @"D:\SingleCellDataSets\MSV000085937\NonDepleted_Lung_SingleCellProteomics\NanoPOTS_D022_Primary_SC_Resting_Chip1_C6.raw",
                @"D:\SingleCellDataSets\MSV000085937\NonDepleted_Lung_SingleCellProteomics\NanoPOTS_D022_Primary_SC_Resting_Chip1_C7.raw",
                @"D:\SingleCellDataSets\MSV000085937\NonDepleted_Lung_SingleCellProteomics\NanoPOTS_D022_Primary_SC_Resting_Chip1_C8.raw",
                @"D:\SingleCellDataSets\MSV000085937\NonDepleted_Lung_SingleCellProteomics\NanoPOTS_D022_Primary_SC_Resting_Chip1_C9.raw",
                @"D:\SingleCellDataSets\MSV000085937\NonDepleted_Lung_SingleCellProteomics\NanoPOTS_D022_Primary_SC_Non_Resting_Chip1_A1.raw",
                @"D:\SingleCellDataSets\MSV000085937\NonDepleted_Lung_SingleCellProteomics\NanoPOTS_D022_Primary_SC_Non_Resting_Chip1_A2.raw",
                @"D:\SingleCellDataSets\MSV000085937\NonDepleted_Lung_SingleCellProteomics\NanoPOTS_D022_Primary_SC_Non_Resting_Chip1_A3.raw",
                @"D:\SingleCellDataSets\MSV000085937\NonDepleted_Lung_SingleCellProteomics\NanoPOTS_D022_Primary_SC_Non_Resting_Chip1_A4.raw",
                @"D:\SingleCellDataSets\MSV000085937\NonDepleted_Lung_SingleCellProteomics\NanoPOTS_D022_Primary_SC_Non_Resting_Chip1_A5.raw",
                @"D:\SingleCellDataSets\MSV000085937\NonDepleted_Lung_SingleCellProteomics\NanoPOTS_D022_Primary_SC_Non_Resting_Chip1_A6.raw"
            };
            List<SpectraFileInfo> spectraFiles = ExperimentalDesign.ReadExperimentalDesign(
                experimentalDesignFilepath, rawFilePathList, out var errors);

            Assert.That(errors.Count == 0);

            string maxQuantEvidencePath = Path.Join(folderPath, @"combined\txt\evidence.txt");
            Dictionary<string, List<ChromatographicPeak>> mbrPeaks = PsmGenericReader.ReadInMbrPeaks(
                maxQuantEvidencePath, silent: false, spectraFiles);

            string maxQuantMsmsPath = Path.Join(folderPath, @"combined\txt\msms.txt");
            Dictionary<string, PsmFromTsv> maxQuantPsms = PsmGenericReader.GetDonorPsms(
                maxQuantMsmsPath, spectraFiles, mbrPeaks, ignoreArtifactIons: true);

            // Every MBR run ID should have a corresponding PSM in the msms.txt file
            Assert.AreEqual(mbrPeaks.Count, maxQuantPsms.Count, 1); // There is one donor where the matched fragments are missing

            PpmTolerance testTolerance = new PpmTolerance(5);
            // Check that the PeptideMonoMass (derived from msms.txt mass columns) and the pwsm mass (derived from converted full sequence)
            // matches for every psm
            Assert.IsFalse(maxQuantPsms.Values.Any(p =>
                !testTolerance.Within(double.Parse(p.PeptideMonoMass), p.PeptideWithSetModifications.MonoisotopicMass)));

            // Writing a spectral library
            var spectraForSpectraLibrary = new List<LibrarySpectrum>();
            foreach (var psm in maxQuantPsms.Values)
            {
                var standardSpectrum = new LibrarySpectrum(psm.FullSequence, psm.PrecursorMz, psm.PrecursorCharge, psm.MatchedIons, psm.RetentionTime ?? 0);
                spectraForSpectraLibrary.Add(standardSpectrum);
            }
            string spectrumFilePath = Path.Join(folderPath, @"combined\txt\spectralLibrary.msp");
            if (File.Exists(spectrumFilePath)) File.Delete(spectrumFilePath);
            using (StreamWriter output = new StreamWriter(spectrumFilePath))
            {
                foreach (var librarySpectrum in spectraForSpectraLibrary)
                {
                    output.WriteLine(librarySpectrum.ToString());
                }
            }

            // Tolerances taken from MaxQuant defaults
            CommonParameters commonParams = new CommonParameters(dissociationType: DissociationType.Autodetect,
                productMassTolerance: new PpmTolerance(20), deconvolutionMassTolerance: new PpmTolerance(7));

            string resultsPath = Path.Combine(folderPath, "SpectralRecoveryMaxQuant");
            Directory.CreateDirectory(resultsPath);

            var mbrAnalysisResults = SpectralRecoveryRunner.RunSpectralRecoveryFromMaxQuant(
                spectraFiles,
                mbrPeaks,
                maxQuantPsms,
                spectrumFilePath,
                resultsPath + @"\",
                commonParams);

            mbrAnalysisResults.WritePeakQuantificationResultsToTsv(resultsPath, "PeakQuant_NoArtifact");
        }

        [Test]
        public static void ParseModificationsTest()
        {
            string fullSeq = "_DEM(Oxidation (M))LESICGVT(Phospho (T))PLEK_";
            string fullSeq2 = "_(Acetyl (Protein N-term))AAAAAAAAAAGAAGGR_";
            string fullSeq3 = "_CCCC_";
            string fullSeq4 = "_(Acetyl (Protein N-term))M(Oxidation (M))M(Oxidation (M))IGLPGAGK_";
            string fullSeq5 = "_(Acetyl (Protein N-term))CNTPTYCDLGK_";

            string mmFullSeq = PsmGenericReader.ConvertMaxQuantFullSequence(fullSeq, out var allKnownMods, out int numFixedMods);
            Assert.That(mmFullSeq.Equals("DEM[Common Variable:Oxidation on M]LESIC[Common Fixed:Carbamidomethyl on C]GVT[Unimod:Phospho on T]PLEK"));
            Assert.That(allKnownMods.Count == 3 & numFixedMods == 1);

            string mmFullSeq2 = PsmGenericReader.ConvertMaxQuantFullSequence(fullSeq2, out allKnownMods, out numFixedMods);
            Assert.That(mmFullSeq2.Equals("[Unimod:Acetyl on X]AAAAAAAAAAGAAGGR"));
            Assert.That(allKnownMods.Count == 1 & numFixedMods == 0);

            string mmFullSeq3 = PsmGenericReader.ConvertMaxQuantFullSequence(fullSeq3, out allKnownMods, out numFixedMods);
            Assert.That(mmFullSeq3.Equals("C[Common Fixed:Carbamidomethyl on C]C[Common Fixed:Carbamidomethyl on C]C[Common Fixed:Carbamidomethyl on C]C[Common Fixed:Carbamidomethyl on C]"));
            Assert.That(allKnownMods.Count == 1 & numFixedMods == 4);

            string mmFullSeq4 = PsmGenericReader.ConvertMaxQuantFullSequence(fullSeq4, out allKnownMods, out numFixedMods);
            Assert.That(mmFullSeq4.Equals("[Unimod:Acetyl on X]M[Common Variable:Oxidation on M]M[Common Variable:Oxidation on M]IGLPGAGK"));
            Assert.That(allKnownMods.Count == 2 & numFixedMods == 0);

            string mmFullSeq5 = PsmGenericReader.ConvertMaxQuantFullSequence(fullSeq5, out allKnownMods, out numFixedMods);
            Assert.That(mmFullSeq5.Equals("[Unimod:Acetyl on C]C[Common Fixed:Carbamidomethyl on C]NTPTYC[Common Fixed:Carbamidomethyl on C]DLGK"));
            Assert.That(allKnownMods.Count == 2 & numFixedMods == 2);

        }

        [Test]
        public static void ReadMaxQuantFragmentsTest()
        {
            // Test that neutral mass losses are read correctly
            string matches = "y1;y3;y4;y5;y6;y7;y8;y9;y10;y11;y12;y13;y5-NH3;b2;b3;b4;b5;b6;b7(10+);b8(2+)";
            string masses =
                "175.1192169189453;289.1623840332031;360.2002868652344;431.2371520996094;488.2598507332258;559.2946166992188;630.3342644496789;701.3716342218012;772.4094891178714;843.4456672801743;914.4825152186543;985.5216128341629;414.2112731933594;185.0926055908203;256.1300048828125;327.16729736328125;398.20452947726307;469.2414310041778;54.9350987;306.16195";
            string intensities =
                "1280.7003173828125;1053.8074951171875;1268.592041015625;1887.4310302734375;6766.16357421875;5207.5390625;9180.8564453125;12792.216796875;15985.1328125;13708.62109375;10652.765625;4033.37939453125;895.0408935546875;7100.369140625;12120.0615234375;14993.857421875;14382.29296875;9648.4521484375;5916.3203125;2529.598876953125";

            List<string> fragmentsAndMass = PsmFromTsv.BuildFragmentStringsFromMaxQuant(matches, masses);
            List<string> fragmentsAndIntensities = PsmFromTsv.BuildFragmentStringsFromMaxQuant(matches, intensities);

            Assert.That(fragmentsAndMass.Count == 20 & fragmentsAndIntensities.Count == 20);

            List<MatchedFragmentIon> ions = PsmFromTsv.ReadFragmentIonsFromList(fragmentsAndMass,
                fragmentsAndIntensities,
                peptideBaseSequence: "AAAAAAAAAAGAAGGR");

            // Test that the neutral mass losses was read correctly
            Assert.That(ions.
                Where(i => i.NeutralTheoreticalProduct.Annotation.Equals("y5-17.03")).
                Select(i => i.NeutralTheoreticalProduct.NeutralMass + i.NeutralTheoreticalProduct.NeutralLoss + Constants.ProtonMass).
                First(), Is.EqualTo(431.24).Within(0.01));

            // Test that fragment charge states are read correctly
            Assert.That(ions.
                Where(i => i.NeutralTheoreticalProduct.Annotation.Equals("b7")).
                Select(i => i.Charge).
                First() == 10);
            Assert.That(ions.
                Where(i => i.NeutralTheoreticalProduct.Annotation.Equals("b8")).
                Select(i => i.Charge).
                First() == 2);
            Assert.That(ions.Count(i => i.Charge == 1) == 18);

        }
    }
}
