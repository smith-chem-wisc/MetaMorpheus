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
