using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace EngineLayer.ModSearch
{
    /// <summary>
    /// Represent a hypothesis of the pair of modification and site. Include the modification Id, motif index and existece of the MS2 spectrum.
    /// </summary>
    public class ModSitePair(int siteIndex, int modId, bool hasMs2Spectrum)
    {
        public int SiteIndex { get; } = siteIndex; // The index of the site in the peptide sequence
        public int ModId { get; } = modId; // The ID of the modification, which corresponds to its index in the modBox
        public bool HasMs2Spectrum { get; } = hasMs2Spectrum; // Indicates whether the MS2 spectrum exists for this modification at the site
        public bool Confident { get; set; } = false; // Indicates whether the modification is confidently localized at the site, if this pair occurs within of all hypothesis, then it is confident localized.

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
