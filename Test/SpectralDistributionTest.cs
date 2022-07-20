using EngineLayer;
using EngineLayer.ClassicSearch;
using MassSpectrometry;
using NUnit.Framework;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using TaskLayer;
using EngineLayer.MbrAnalysis;

namespace Test
{
    internal class SpectralDistributionTest
    {
        [Test]
        public static void TestSequenceConverter()
        {
            string pepStringNoMod = "LVEALCAEHQINLIK";
            pepStringNoMod = new string(TrueNegativeDistribution.ConvertSequence(pepStringNoMod));
            Assert.That(pepStringNoMod.Equals("LVEALCAEHQINLIK"));

            string pepStringOneMod = "LVEALC[Common Fixed:Carbamidomethyl on C]AEHQINLIK";
            pepStringOneMod = new string(TrueNegativeDistribution.ConvertSequence(pepStringOneMod));
            Assert.That(pepStringOneMod.Equals("LVEAL0AEHQINLIK"));

            string pepStringTwoMod = "PNM[Common Variable:Oxidation on M]VTPGHAC[Common Fixed:Carbamidomethyl on C]TQK";
            pepStringTwoMod = new string(TrueNegativeDistribution.ConvertSequence(pepStringTwoMod));
            Assert.That(pepStringTwoMod.Equals("PN0VTPGHA1TQK"));

        }

        [Test]
        public static void TestPercentHomology()
        {
            string pepA = "LVEALCAEHQINLIK";
            string pepB = "LVEALC[Common Fixed:Carbamidomethyl on C]AEHQINLIK";
            string pepC = "PNMVTP[Common Variable:Oxidation on M]GHAC[Common Fixed:Carbamidomethyl on C]TQLIK";

            double homology = TrueNegativeDistribution.GetPercentHomology(pepA, pepB);
            Assert.That(homology > 0.93 && homology < 0.95);

            // This shouldn't be true, is caused by local instantiation of dictionary, needs to be distribtion specific
            homology = TrueNegativeDistribution.GetPercentHomology(pepC, pepB);
            Assert.That(Math.Abs(homology-0.26666) < 0.001);

        }
    }
}
