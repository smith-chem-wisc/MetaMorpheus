using EngineLayer.DIA;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Linq;

namespace Test.DIATests
{
    public class TestPfPairAndPfGroup
    {
        [Test]
        public static void TestFragmentRank()
        {
            //Create a list of fake peaks as the precursor XIC
            var precursorPeaks = new List<IIndexedPeak>();
            double[] intensityMultipliers = { 1, 2, 3, 4, 3, 2, 1 };
            for (int i = 0; i < intensityMultipliers.Length; i++)
            {
                precursorPeaks.Add(new IndexedMassSpectralPeak(intensity: 1e5 * intensityMultipliers[i], retentionTime: 1 + i / 10, zeroBasedScanIndex: i, mz: 500.0));
            }
            var precursorXic = new ExtractedIonChromatogram(precursorPeaks);

            //Create a list of peak list as fragment XICs
            int numberOfFragments = 5;
            var fragmentXics = new List<ExtractedIonChromatogram>();
            for (int i = 0; i < numberOfFragments; i++)
            {
                var fragmentPeaks = new List<IIndexedPeak>();
                for (int j = 0; j < intensityMultipliers.Length; j++)
                {
                    fragmentPeaks.Add(new IndexedMassSpectralPeak(intensity: 1e5 * intensityMultipliers[j] * (i + 1), retentionTime: 1 + j / 10, zeroBasedScanIndex: j, mz: 500.0 + i));
                }
                fragmentXics.Add(new ExtractedIonChromatogram(fragmentPeaks));
            }

            //For each fragment XIC, try to pair with the precursor XIC
        }
    }
}
