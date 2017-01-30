using System.Collections.Generic;

namespace EngineLayer.ModernSearch
{
    public class ModernSpectrumMatch : ParentSpectrumMatch
    {

        #region Private Fields

        private CompactPeptide compactPeptide;

        #endregion Private Fields

        #region Public Constructors

        public ModernSpectrumMatch(CompactPeptide theBestPeptide, string fileName, double scanRetentionTime, double scanPrecursorIntensity, double scanPrecursorMass, int scanNumber, int scanPrecursorCharge, int scanExperimentalPeaks, double totalIonCurrent, double scanPrecursorMZ, double score) : base(fileName, scanRetentionTime, scanPrecursorIntensity, scanPrecursorMass, scanNumber, scanPrecursorCharge, scanExperimentalPeaks, totalIonCurrent, scanPrecursorMZ, score)
        {
            compactPeptide = theBestPeptide;
        }

        #endregion Public Constructors

        #region Public Methods

        public override CompactPeptide GetCompactPeptide(List<MetaMorpheusModification> variableModifications, List<MetaMorpheusModification> localizeableModifications, List<MetaMorpheusModification> fixedModifications)
        {
            return compactPeptide;
        }

        #endregion Public Methods

    }
}