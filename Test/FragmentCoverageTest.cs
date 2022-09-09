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
using Newtonsoft.Json.Linq;
using TaskLayer;
using UsefulProteomicsDatabases;
using System.Text.RegularExpressions;

namespace Test
{
    [TestFixture]
    public static class FragmentCoverageTest
    {
        private static List<string> FragmentCoveragePeptides(List<PsmFromTsv> PSMs)
        {
            List<string> CoveredAminoAcids = new();

            foreach (var psm in PSMs)
            {
                var peptideLength = psm.BaseSeq.Length;

                var peptideLengthMinus1 = peptideLength - 1;

                var startAndEnd = psm.StartAndEndResiduesInProtein;

                List<int> bIonPositions = new();
                List<int> yIonPositions = new();

                var temp = psm.MatchedIons.Where(p => p.NeutralTheoreticalProduct.ProductType == ProductType.b ||p.NeutralTheoreticalProduct.ProductType == ProductType.y).ToList();
                foreach (var ion in psm.MatchedIons)
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

                for (int i = 0; i < bIonPositionsUnique.Count; i++)
                {
                    if (i + 1 < bIonPositionsUnique.Count)
                    {
                        if (bIonPositionsUnique[i+1] - bIonPositionsUnique[i] == 1)
                        {
                            Covered.Add(bIonPositionsUnique[i+1]);
                        }
                    }

                    var yValNecessary = peptideLengthMinus1 - bIonPositionsUnique[i];

                    if (yIonPositionsUnique.Contains(yValNecessary))
                    {
                        Covered.Add(bIonPositionsUnique[i] + 1);
                    }
                }

                if (bIonPositionsUnique[^1] == peptideLengthMinus1)
                {
                    Covered.Add(peptideLength);
                }

                for (int i = 0; i < yIonPositionsUnique.Count; i++)
                {
                    if (i + 1 < yIonPositionsUnique.Count)
                    {
                        if (yIonPositionsUnique[i + 1] - yIonPositionsUnique[i] == 1)
                        {
                            Covered.Add(peptideLength - yIonPositionsUnique[i]);
                        }
                    }

                }

                if (yIonPositionsUnique[^1] == peptideLengthMinus1)
                {
                    Covered.Add(1);
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

                char[] aminoAcids = psm.BaseSeq.ToCharArray();

                for (int i = 0; i < aminoAcids.Length; i++)
                {
                    if (!coveredBool[i])
                    {
                        aminoAcids[i] = char.ToLower(aminoAcids[i]);
                    }
                }

                var fragmentCoverageMap = string.Join("", aminoAcids);
                CoveredAminoAcids.Add(fragmentCoverageMap);
                Console.WriteLine(fragmentCoverageMap);
                var startResidue = Int32.Parse(Regex.Match(psm.StartAndEndResiduesInProtein, @"\d+").Value);
                Console.WriteLine("Start value " + startResidue);
                List<int> CoveredUniqueAminoAcidNumberInProt = new();
                for (int i = 0; i < CoveredUnique.Count; i++)
                {
                    CoveredUniqueAminoAcidNumberInProt.Add(CoveredUnique[i] + startResidue - 1);
                }
                Console.WriteLine(string.Join(", ", CoveredUniqueAminoAcidNumberInProt));
            }
            return CoveredAminoAcids;
        }

        [Test]
        public static void TestFragmentCoverage()
        {
            var testPSMs = PsmTsvReader.ReadTsv(@"TestData\oglyco.psmtsv", out var warnings);

            FragmentCoveragePeptides(testPSMs);

            //for (int j = 0; j < testPSMs.Count; j++)
            //{
            //    var peptideLength = testPSMs[j].BaseSeq.Length;

            //    var peptideLengthMinus1 = peptideLength - 1;

            //    var startAndEnd = testPSMs[j].StartAndEndResiduesInProtein;

            //    List<int> bIonPositions = new();
            //    List<int> yIonPositions = new();

            //    var temp = testPSMs[j].MatchedIons.Where(p => p.NeutralTheoreticalProduct.ProductType == ProductType.b ||p.NeutralTheoreticalProduct.ProductType == ProductType.y).ToList();
            //    foreach (var ion in testPSMs[j].MatchedIons)
            //    {

            //        if (ion.NeutralTheoreticalProduct.ProductType == ProductType.b)
            //        {
            //            bIonPositions.Add(ion.NeutralTheoreticalProduct.AminoAcidPosition);
            //        }

            //        if (ion.NeutralTheoreticalProduct.ProductType == ProductType.y)
            //        {
            //            yIonPositions.Add(ion.NeutralTheoreticalProduct.AminoAcidPosition);
            //        }
            //    }

            //    var bIonPositionsUnique = bIonPositions.Distinct().ToList();
            //    var yIonPositionsUnique = yIonPositions.Distinct().ToList();

            //    bIonPositionsUnique.Sort();
            //    yIonPositionsUnique.Sort();

            //    List<int> Covered = new();

            //    for (int i = 0; i < bIonPositionsUnique.Count - 1; i++)
            //    {

            //        if (bIonPositionsUnique[i+1] - bIonPositionsUnique[i] == 1)
            //        {
            //            Covered.Add(bIonPositionsUnique[i+1]);
            //        }

            //        if (bIonPositionsUnique[i] == peptideLengthMinus1)
            //        {
            //            Covered.Add(bIonPositionsUnique[i]);
            //        }

            //        var yValNecessary = peptideLengthMinus1 - bIonPositionsUnique[i];

            //        if (yIonPositionsUnique.Contains(yValNecessary))
            //        {
            //            Covered.Add(bIonPositionsUnique[i] + 1);
            //        }
            //    }

            //    for (int i = 0; i < yIonPositionsUnique.Count - 1; i++)
            //    {

            //        if (yIonPositionsUnique[i + 1] - yIonPositionsUnique[i] == 1)
            //        {
            //            Covered.Add(peptideLength - yIonPositionsUnique[i]);
            //        }

            //        if (yIonPositionsUnique[i] == peptideLengthMinus1)
            //        {
            //            Covered.Add(1);
            //        }
            //    }

            //    var CoveredUnique = Covered.Distinct().ToList();
            //    CoveredUnique.Sort();

            //    bool[] coveredBool = new bool[peptideLength];
            //    for (int i = 0; i < coveredBool.Length; i++)
            //    {
            //        if (CoveredUnique.Contains(i + 1))
            //        {
            //            coveredBool[i] = true;
            //        }
            //        else
            //        {
            //            coveredBool[i] = false;
            //        }
            //    }

            //    char[] aminoAcids = testPSMs[j].BaseSeq.ToCharArray();

            //    for (int i = 0; i < aminoAcids.Length; i++)
            //    {
            //        if (!coveredBool[i])
            //        {
            //            aminoAcids[i] = char.ToLower(aminoAcids[i]);
            //        }
            //    }

            //    var fragmentCoverageMap = string.Join("", aminoAcids);
            //    Console.WriteLine(fragmentCoverageMap);
            //}
        }

    }
}