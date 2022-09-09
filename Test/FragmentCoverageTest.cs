using Chemistry;
using EngineLayer;
using EngineLayer.ClassicSearch;
using EngineLayer.FdrAnalysis;
using EngineLayer.Indexing;
using EngineLayer.ModernSearch;
using MassSpectrometry;
using MzLibUtil;
using NUnit.Framework;
using Proteomics;
using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Drawing;
using System.IO;
using System.Linq;
using TaskLayer;
using UsefulProteomicsDatabases;

namespace Test
{
    [TestFixture]
    public static class FragmentCoverageTest
    {

        [Test]
        public static void TestFragmentCoverage()
        {
            var testPSMs = PsmTsvReader.ReadTsv(@"TestData\SequenceCoverageTestPSM.psmtsv", out var warnings);

            var peptideLength = testPSMs[0].BaseSeq.Length;

            var peptideLengthMinus1 = peptideLength - 1;

            var startAndEnd = testPSMs[0].StartAndEndResiduesInProtein;

            List<int> bIonPositions = new();
            List<int> yIonPositions = new();

            var temp = testPSMs[0].MatchedIons.Where(p => p.Annotation.Contains('b') || p.Annotation.Contains('y')).ToList();
            foreach (var ion in testPSMs[0].MatchedIons)
            {
                
                if (ion.NeutralTheoreticalProduct.ProductType == ProductType.b)
                {
                    bIonPositions.Add(ion.NeutralTheoreticalProduct.AminoAcidPosition);
                }

                if (ion.NeutralTheoreticalProduct.ProductType == ProductType.y)
                {
                    yIonPositions.Add(ion.NeutralTheoreticalProduct.AminoAcidPosition);
                }
            }

            var bIonPositionsUnique = bIonPositions.Distinct().ToList();
            var yIonPositionsUnique = yIonPositions.Distinct().ToList();

            bIonPositionsUnique.Sort();
            yIonPositionsUnique.Sort();

            List<int> Covered = new();

            for (int i = 0; i < bIonPositionsUnique.Count - 1; i++)
            {

                if (bIonPositionsUnique[i+1] - bIonPositionsUnique[i] == 1)
                {
                    Covered.Add(bIonPositionsUnique[i+1]);
                }

                if (bIonPositionsUnique[i] == peptideLengthMinus1)
                {
                    Covered.Add(bIonPositionsUnique[i]);
                }

                var yValNecessary = peptideLengthMinus1 - bIonPositionsUnique[i];

                if (yIonPositionsUnique.Contains(yValNecessary))
                {
                    Covered.Add(bIonPositionsUnique[i] + 1);
                }
            }

            for (int i = 0; i < yIonPositionsUnique.Count - 1; i++)
            {

                if (yIonPositionsUnique[i + 1] - yIonPositionsUnique[i] == 1)
                {
                    Covered.Add(peptideLength - yIonPositionsUnique[i]);
                }

                if (yIonPositionsUnique[i] == peptideLengthMinus1)
                {
                    Covered.Add(1);
                }
            }

            var CoveredUnique = Covered.Distinct().ToList();
            CoveredUnique.Sort();

            bool[] coveredBool = new bool[peptideLength];
            for (int i = 0; i < coveredBool.Length; i++)
            {
                if (CoveredUnique.Contains(i + 1))
                {
                    coveredBool[i] = true;
                }
                else
                {
                    coveredBool[i] = false;
                }
            }

            char[] aminoAcids = testPSMs[0].BaseSeq.ToCharArray();

            for (int i = 0; i < aminoAcids.Length; i++)
            {
                if (!coveredBool[i])
                {
                    aminoAcids[i] = char.ToLower(aminoAcids[i]);
                }
            }

            var fragmentCoverageMap = string.Join("", aminoAcids);
            Console.WriteLine(fragmentCoverageMap);
        }
    }
}