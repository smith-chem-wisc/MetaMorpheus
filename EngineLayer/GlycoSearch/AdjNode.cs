using System;
using System.Collections.Generic;
using System.Text;

namespace EngineLayer.GlycoSearch
{
    public class AdjNode
    {
        public AdjNode(int i, int j, int modPos, ModBox box)
        {
            PosX = i;
            PosY = j;
            ModPos = modPos;
            ModBox = box;
        }

        public int PosX { get; }
        public int PosY { get; }

        public int ModPos { get; }
        public ModBox ModBox { get; }

        //sources are represented by index. Only track ones with highest cummulative cost
        public List<int> CummulativeSources = new List<int>();

        public double maxCost { get; set; }

        public double CurrentCost { get; set; }

        public List<int> AllSources { get; set; } = new List<int>();
    }
}
