using EngineLayer.DIA;
using MassSpectrometry;
using MathNet.Numerics.Statistics;
using MzLibUtil;
using NUnit.Framework;
using NUnit.Framework.Interfaces;
using System;
using System.Collections.Generic;
using System.Linq;
using static Org.BouncyCastle.Asn1.Cmp.Challenge;

namespace Test.DIATests
{
    public class TestPfPairAndPfGroup
    {
        [Test]
        public static void TestFragmentAndPrecursorRank()
        {
            //This test creates a set of precursor and fragment XICs with random intensities and group them.
            //The resulting precursor-fragment pairs will have random correlations depending on their generated intensities. However, the rankings should be done as expected regardless of the correlation values.
            var random = new Random(42);

            // Create a set of precursor XICs with random intensities
            int numberOfPrecursors = 10;
            var precursorXics = new List<ExtractedIonChromatogram>();
            for (int i = 0; i < numberOfPrecursors; i++)
            {
                var peakList = new List<IIndexedPeak>();
                for (int j = 0; j < 10; j++)
                {
                    peakList.Add(new IndexedMassSpectralPeak(500.0 + i, random.NextDouble() * 1e5,  j, 1 + j / 10.0));
                }
                precursorXics.Add(new ExtractedIonChromatogram(peakList));
            }

            // Create a set of fragment XICs with random intensities
            int numberOfFragments = 20;
            var fragmentXics = new List<ExtractedIonChromatogram>();
            for (int i = 0; i < numberOfFragments; i++)
            {
                var peakList = new List<IIndexedPeak>();
                for (int j = 0; j < 10; j++)
                {
                    peakList.Add(new IndexedMassSpectralPeak(300.0 + i, random.NextDouble() * 1e5, j, 1.01 + j / 10.0));
                }
                fragmentXics.Add(new ExtractedIonChromatogram(peakList));
            }

            //Try to group fragment Xics with each precursor Xic. A set of reasonable threshold values are used for grouping.
            var pfGroups = new List<PrecursorFragmentsGroup>();
            foreach (var precursorXic in precursorXics)
            {
                var pfGroup = XicGroupingEngine.GroupFragmentsForOnePrecursor(precursorXic, fragmentXics, 0.3f, 0.3, 0, 1);
                if (pfGroup != null)
                {
                    pfGroups.Add(pfGroup);
                    //In this precursor-fragments grouping method, all pfPairs should have a valid Correlation value higher than the threshold
                    Assert.That(pfGroup.PFpairs.All(pf => pf.Correlation.HasValue));
                }
            }

            //Test fragment ranking. Fragment ranking is performed within each pfGroup.
            foreach (var pfGroup in pfGroups)
            {
                pfGroup.SetFragmentRankForPfPairs();
                //All PFpairs should be ordered with decreasing correlation values.
                for (int i = 0; i < pfGroup.PFpairs.Count - 1; i++)
                {
                    Assert.That(pfGroup.PFpairs[i].Correlation, Is.GreaterThanOrEqualTo(pfGroup.PFpairs[i + 1].Correlation));
                    Assert.That(pfGroup.PFpairs[i].FragmentRank, Is.EqualTo(i + 1));
                }
            }

            //Test precursor ranking. Precursor ranking needs to be performed across all groups, based on each fragment Xic. 
            PrecursorFragmentPair.SetPrecursorRankForPfPairs(pfGroups.SelectMany(g => g.PFpairs));
            //pfPairs are grouped by fragment Xics to rank precursor Xics, so precursor ranks are tested by each fragment Xic
            foreach(var fragment in fragmentXics)
            {
                //Get all pfPairs associated with this fragment Xic. Skip if no pfPairs are found.
                var fragmentPairs = pfGroups.SelectMany(g => g.PFpairs.Where(pf => pf.FragmentXic == fragment)).ToList();
                if (fragmentPairs.IsNullOrEmpty()) continue;
                //precursor ranking is performed across all pfGroups so any pfPair should have a valid PrecursorRank value
                Assert.That(fragmentPairs.All(pf => pf.PrecursorRank.HasValue));
                //The largest PrecursorRank of these pfPairs should be equal to the total number of precursor XICs associated with this fragment Xic
                Assert.That(fragmentPairs.Max(pf => pf.PrecursorRank), Is.EqualTo(fragmentPairs.Count()));
                //If we rank these pfPairs based on their PrecursorRank, they should be ordered with decreasing Correlation values.
                fragmentPairs.Sort((a, b) => a.PrecursorRank.Value.CompareTo(b.PrecursorRank.Value));
                for (int i = 0; i < fragmentPairs.Count - 1; i++)
                {
                    Assert.That(fragmentPairs[i].Correlation, Is.GreaterThanOrEqualTo(fragmentPairs[i + 1].Correlation));
                }
            }
        }
    }
}
