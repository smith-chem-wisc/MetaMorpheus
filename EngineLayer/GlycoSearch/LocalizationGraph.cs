using System;
using System.Collections.Generic;
using System.Text;
using System.Linq;
using EngineLayer.ModernSearch;
using Proteomics.Fragmentation;

namespace EngineLayer.GlycoSearch
{
    public class LocalizationGraph
    {
        public AdjNode[][] array;
        public LocalizationGraph(int posCount, int glycanCount)
        {
            array = new AdjNode[posCount][];
            for (int i = 0; i < posCount; i++)
            {
                array[i] = new AdjNode[glycanCount];
            }
        }

        public void Localization(int[] modPos, GlycanBox glycanBox, GlycanBox[] boxes, HashSet<int> allPeaks, List<Product> products, int peptideLength)
        {
            for (int i = 0; i < modPos.Length; i++)
            {
                int maxLength = i + 1;
                int minlength = glycanBox.GlycanIds.Length - (modPos.Length - 1 - i);

                for (int j = 0; j < boxes.Length; j++)
                {
                    if (boxes[j].NumberOfGlycans <= maxLength && boxes[j].NumberOfGlycans >= minlength)
                    {
                        AdjNode adjNode = new AdjNode(i, j, modPos[i], boxes[j]);
                        var cost = LocalizationGraph.CalculateCost(allPeaks, products, peptideLength, modPos, i, glycanBox, boxes[j], 1000);
                        if (i == 0)
                        {
                            //Get cost                             
                            adjNode.maxCost = cost;

                        }
                        else
                        {
                            double maxCost = 0;
                            for (int prej = 0; prej <= j; prej++)
                            {
                                if (boxes[j].NumberOfGlycans <= boxes[prej].NumberOfGlycans + 1 &&
                                    (boxes[prej].GlycanIds.Length == 0 || !boxes[prej].GlycanIds.Except(boxes[j].GlycanIds).Any()) &&
                                    array[i - 1][prej] != null)
                                {
                                    var tempCost = cost + array[i - 1][prej].maxCost;
                                    if (tempCost > maxCost)
                                    {
                                        adjNode.Sources.Clear();
                                        adjNode.Costs.Clear();

                                        adjNode.Sources.Add(prej);
                                        adjNode.Costs.Add(tempCost);
                                        maxCost = tempCost;
                                    }
                                    else if (tempCost == maxCost)
                                    {
                                        adjNode.Sources.Add(prej);
                                        adjNode.Costs.Add(tempCost);
                                    }

                                }
                            }
                            //if (adjNode.Costs.Any())
                            {
                                adjNode.maxCost = adjNode.Costs.Max();
                            }
                        }

                        array[i][j] = adjNode;
                    }
                }

            }
        }

        public static double CalculateCost(HashSet<int> allPeaksForLocalization, List<Product> products, int peptideLength, int[] modPos, int modInd, GlycanBox OGlycanBox, GlycanBox box, int FragmentBinsPerDalton)
        {
            if (modInd == modPos.Length - 1)
            {
                return 0;
            }

            var fragmentHash = GlycoPeptides.GetFragmentHash(products, peptideLength, modPos, modInd, OGlycanBox, box, FragmentBinsPerDalton);

            int currentLocalizationScore = allPeaksForLocalization.Intersect(fragmentHash).Count();

            return (double)currentLocalizationScore;
        }
    }

    public class AdjNode
    {
        public AdjNode(int i, int j, int modPos, GlycanBox box)
        {
            PosX = i;
            PosY = j;
            ModPos = modPos;
            glycanBox = box;
        }

        public int PosX { get; }
        public int PosY { get; }

        public int ModPos { get; }
        public GlycanBox glycanBox { get; }

        //sources are represented by index.
        public List<int> Sources = new List<int>();

        public List<double> Costs = new List<double>();

        public double maxCost { get; set; }

    }

}
