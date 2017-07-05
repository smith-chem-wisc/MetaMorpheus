using Proteomics;
using System.Collections.Generic;

namespace EngineLayer.ModernSearch
{
    public class PsmModern : PsmParent
    {

        #region Private Fields

        private readonly CompactPeptide compactPeptide;

        #endregion Private Fields

        #region Public Constructors

        public PsmModern(CompactPeptide theBestPeptide, int notch, double score, int scanIndex, Ms2ScanWithSpecificMass scan) : base(notch, score, scanIndex, scan, theBestPeptide.MonoisotopicMassIncludingFixedMods)
        {
            compactPeptide = theBestPeptide;
        }

        #endregion Public Constructors

        #region Public Methods

        public override CompactPeptide GetCompactPeptide(Dictionary<ModificationWithMass, ushort> modsDictionary)
        {
            return compactPeptide;
        }

        #endregion Public Methods

    }
}