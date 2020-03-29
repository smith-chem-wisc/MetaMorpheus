using System;
using System.Collections.Generic;
using System.Text;

namespace EngineLayer.GlycoSearch
{
    public class Route
    { 

        public int ModBoxId { get; set; }

        //Tuple<int, int, double> mod pos, glycan id, local peak exist
        //For the local peak exist, the idea is that, in the localization graph matrix, if the node is detected as a mod, then the node score and the previous node has a current score >0.
        public List<Tuple<int, int, bool>> Mods { get; private set; } = new List<Tuple<int, int, bool>>();

        public double Score { get; set; }

        public double ReversePScore { get; set; }

        public void AddPos(int pos, int id, bool localPeakExist)
        {
            Mods.Add(new Tuple<int, int, bool>(pos, id, localPeakExist));
        }

    }
}
