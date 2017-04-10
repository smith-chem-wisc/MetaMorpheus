using Proteomics;
using System.Collections.Generic;

namespace EngineLayer.ModernSearch
{
    public class PsmModern : PsmParent
    {

        #region Private Fields

        private CompactPeptide compactPeptide;

        #endregion Private Fields

        #region Public Constructors

        public PsmModern(CompactPeptide theBestPeptide, string fileName, double scanRetentionTime, double scanPrecursorIntensity, double scanPrecursorMass, int scanNumber, int scanPrecursorNumber, int scanPrecursorCharge, int scanExperimentalPeaks, double totalIonCurrent, double scanPrecursorMZ, double score, int notch) : base(fileName, scanRetentionTime, scanPrecursorIntensity, scanPrecursorMass, scanNumber, scanPrecursorNumber, scanPrecursorCharge, scanExperimentalPeaks, totalIonCurrent, scanPrecursorMZ, score, notch)
        {
            compactPeptide = theBestPeptide;
        }

        #endregion Public Constructors

        #region Public Methods

        public override CompactPeptide GetCompactPeptide(List<ModificationWithMass> variableModifications, List<ModificationWithMass> localizeableModifications, List<ModificationWithMass> fixedModifications)
        {
            return compactPeptide;
        }

        #endregion Public Methods

    }
}