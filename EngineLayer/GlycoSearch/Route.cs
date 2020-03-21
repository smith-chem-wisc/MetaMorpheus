using System;
using System.Collections.Generic;
using System.Text;

namespace EngineLayer.GlycoSearch
{
    public class Route
    {
        public int ModBoxId { get; set; }
        //public List<int> ModPos { get; set; } = new List<int>();

        //public List<int> ModId { get; set; } = new List<int>();

        public List<Tuple<int, int>> Mods { get; set; } = new List<Tuple<int, int>>();

        public double Score { get; set; }

        public double ReversePScore { get; set; }

        public void AddPos(int pos, int id)
        {
            //ModPos.Add(pos);
            //ModId.Add(id);

            Mods.Add(new Tuple<int, int>(pos, id));
        }

    }
}
