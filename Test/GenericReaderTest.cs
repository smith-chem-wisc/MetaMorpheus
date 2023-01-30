using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Text.RegularExpressions;
using System.Threading.Tasks;
using Chemistry;
using Easy.Common.Extensions;
using EngineLayer;
using EngineLayer.PsmTsv;
using FlashLFQ;
using NUnit.Framework;
using Proteomics.Fragmentation;

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
            Dictionary<string, List<ChromatographicPeak>> mbrPeaks = PsmGenericReader.ReadInMbrPeaks(
                maxQuantEvidencePath, silent: false, spectraFiles);

            string maxQuantMsmsPath = @"D:\HelaSingleCellQCmzML\rawFiles\combined\txt\msms.txt";
            List<PsmFromTsv> maxQuantIds = PsmGenericReader.GetDonorPsms(
                maxQuantMsmsPath, spectraFiles, mbrPeaks);

            //Assert.That(mbrPeaks == duplicatePeaks);
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
            Assert.That(ions.Where(i => i.Charge == 1).Count() == 18);

        }
    }
}
