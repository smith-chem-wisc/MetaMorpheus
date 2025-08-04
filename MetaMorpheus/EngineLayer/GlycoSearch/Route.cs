using EngineLayer.ModSearch;
using Omics;
using Omics.Modifications;
using System;
using System.Collections.Generic;
using System.Data.Entity.Core.Mapping;
using System.Data.Entity.ModelConfiguration.Conventions;
using System.Text;

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
        /// Check if this route contains any N-linked glycans. And out the list of N-linked glycans if exists.
        /// </summary>
        /// <param name="nGlycans"></param>
        /// <returns></returns>
        public bool ContainNGlycan(out List<Modification> nGlycans)
        {
            nGlycans = new List<Modification>();
            foreach (var modSitePair in ModSitePairs)
            {
                if (ModBox.GlobalModifications[modSitePair.ModId] is Glycan glycan &&
                    glycan.ModificationType == "N-linked glycosylation" &&
                    glycan.Ions != null)
                {
                    nGlycans.Add(glycan);
                }
            }
            return nGlycans.Count > 0;
        }

    }
}
