using System;
using System.Collections.Generic;
using System.Text;
using System.Linq;
using EngineLayer.ModernSearch;

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
        public static void CalculateCost(HashSet<int> allPeaksForLocalization, GlycanBox glycanBox, List<int> modPos, int FragmentBinsPerDalton)
        {
            //var fragmentHash = GlycoPeptides.GetFragmentHash(products, glycanBoxId_localization[i], glycanBox, FragmentBinsPerDalton);

            //int currentLocalizationScore = allPeaksForLocalization.Intersect(fragmentHash).Count();

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
