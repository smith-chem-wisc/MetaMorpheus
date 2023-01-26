using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Text.RegularExpressions;
using System.Threading.Tasks;
using EngineLayer;
using EngineLayer.PsmTsv;
using FlashLFQ;
using NUnit.Framework;

namespace Test
{
    [TestFixture]
    internal class GenericReaderTest
    {

        [Test]
        public static void TestMaxQuantEvidenceReader()
        {
            string experimentalDesignFilepath = @"D:\HelaSingleCellQCmzML\ExperimentalDesign_rawFiles.tsv";
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

            string maxQuantEvidencePath = @"D:\HelaSingleCellQCmzML\rawFiles\combined\txt\evidence.txt";
            List<ChromatographicPeak> mbrPeaks = PsmGenericReader.ReadInMbrPeaks(
                maxQuantEvidencePath, silent: false, spectraFiles);

            List<ChromatographicPeak> duplicatePeaks = new();
            List<Identification> maxQuantIds = PsmGenericReader.ReadPsms(
                maxQuantEvidencePath, silent: false, spectraFiles, duplicatePeaks);

            Assert.That(mbrPeaks == duplicatePeaks);
        }

        [Test]
        public static void ParseModificationsTest()
        {
            string fullSeq = "_AM(Oxidation (M))LESICGVT(Phospho (T))PLEK_";
            string fullSeq2 = "_(Acetyl (Protein N-term))AAAAAAAAAAGAAGGR_";

            //string x = PsmGenericReader.ConvertMaxQuantFullSequence(fullSeq, out var allKnownMods, out int numFixedMods);
            string y = PsmGenericReader.ConvertMaxQuantFullSequence(fullSeq2, out var allKnownMods, out var numFixedMods);

            var testMod = GlobalVariables.AllModsKnown.Where(m => m.OriginalId.Contains("Phospho"));
            //var resultsDictionary = PsmGenericReader.ParseMaxQuantFullSeq(fullSeq);
            foreach (var mod in testMod)
            {
                int placeholder = 0;
            }


        }
    }
}
