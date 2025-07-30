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

        /// <summary>
        /// The all modification site pairs in this route.
        /// </summary>
        public List<ModSitePair> ModSitePairs { get; private set; } = new List<ModSitePair>();

        public double Score { get; set; }

        public double ReversePScore { get; set; }

        public void AddPos(int pos, int id, bool localPeakExist)
        {
            ModSitePairs.Add(new ModSitePair(pos, id, localPeakExist));
        }

        /// <summary>
        /// Set the reverse p score to each pair in the route.
        /// </summary>
        public void SetRPScoreToPair()
        {
            foreach (var pair in ModSitePairs)
            {
                pair.RouteScore = ReversePScore;
            }
        }

    }
}
