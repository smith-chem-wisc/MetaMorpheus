using System;
using System.Collections.Generic;
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
    internal class GenericReaderTest
    {

        [Test]
        public static void TestMaxQuantEvidenceReader()
        {
            string experimentalDesignFilepath = @"C/";
            List<string> rawFilePathList = new();
            List<SpectraFileInfo> spectraFiles = ExperimentalDesign.ReadExperimentalDesign(
                experimentalDesignFilepath, rawFilePathList, out var errors);

            string maxQuantEvidencePath = @"";
            List<ChromatographicPeak> mbrPeaks = new();
            List<Identification> maxQuantIds = PsmGenericReader.ReadPsms(
                maxQuantEvidencePath, silent: false, spectraFiles, mbrPeaks);

        }

        [Test]
        public static void ParseModificationsTest()
        {
            string fullSeq = "_AM(Oxidation (M))LESICGVT(Phospho (T))PLEK_";
            string fullSeq2 = "_(Acetyl (Protein N-term))AAAAAAAAAAGAAGGR_";

            string x = PsmGenericReader.MaxQuantToPWSM(fullSeq, null);
            string y = PsmGenericReader.MaxQuantToPWSM(fullSeq2, null);

            var testMod = GlobalVariables.AllModsKnown.Where(m => m.OriginalId.Contains("Phospho"));
            //var resultsDictionary = PsmGenericReader.ParseMaxQuantFullSeq(fullSeq);
            foreach (var mod in testMod)
            {
                int placeholder = 0;
            }


        }
    }
}
