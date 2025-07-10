using System;
using System.Collections.Generic;
using System.Text;
using EngineLayer.ModSearch;
using Omics;

namespace EngineLayer.GlycoSearch
{
    /// <summary>
    /// Presents a possible hypothesis in the glycan localization graph, which contains a set of modification site pairs and their associated scores.
    /// </summary>
    public class Route
    { 

        public int ModBoxId { get; set; }

        //Tuple<int, int, double> mod pos, glycan id, local peak exist
        //For the local peak exist, the idea is that, in the localization graph matrix, if the node is detected as a mod, then the node score and the previous node has a current score >0.
        public List<ModSitePair> ModSitePairs { get; private set; } = new List<ModSitePair>();

        public double Score { get; set; }

        public double ReversePScore { get; set; }

        public void AddPos(int pos, int id, bool localPeakExist)
        {
            ModSitePairs.Add(new ModSitePair(pos, id, localPeakExist));
        }

    }
}
