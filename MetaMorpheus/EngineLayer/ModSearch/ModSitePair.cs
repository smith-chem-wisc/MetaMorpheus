using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace EngineLayer.ModSearch
{
    /// <summary>
    /// Represent a hypothesis of the pair of modification and site. Include the modification Id, motif index and the existece of the MS2 spectrum.
    /// </summary>
    public class ModSitePair(int siteIndex, int modId, bool hasMs2Spectrum)
    {
        /// <summary>
        /// The index of the site in the peptide sequence. The N-terminus is 1, the first aminoAcid is 1, and so on.
        /// </summary>
        public int SiteIndex { get; } = siteIndex;
        /// <summary>
        /// The ID of the modification, which corresponds to its index in the GlobalOGlycans List that defined in GlycanBox class.
        /// </summary>
        public int ModId { get; } = modId;
        /// <summary>
        /// Indicates whether the MS2 spectrum exists for this modification at the site.
        /// The idea is that, in the localization graph matrix, if the node is detected as a mod, then the node score and the previous node has a current score >0.
        /// </summary>
        public bool HasMs2Spectrum { get; } = hasMs2Spectrum;
        /// <summary>
        /// Indicates whether the modification is confidently localized at the site. If this pair occurs within all hypotheses (all paths), then it is confidently localized.
        /// </summary>
        public bool Confident { get; set; } = false;

        /// <summary>
        /// Probability of the modification being localized at the site.
        /// </summary>
        public double Probability { get; set; } = 0.0;

        /// <summary>
        /// The ReversePScore of the route that this pair belongs to.
        /// </summary>
        public double RouteRPScore { get; set; }

        public override bool Equals(object obj)
        {
            if (obj is ModSitePair other)
            {
                return SiteIndex == other.SiteIndex && ModId == other.ModId;
            }
            return false;
        }

        public override int GetHashCode()
        {
            return HashCode.Combine(SiteIndex, ModId);
        }

    }
}
