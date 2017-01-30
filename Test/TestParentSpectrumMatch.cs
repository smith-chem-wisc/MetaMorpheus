using EngineLayer;

using System;
using System.Collections.Generic;

namespace Test
{
    internal class TestParentSpectrumMatch : PsmParent
    {

        #region Public Constructors

        public TestParentSpectrumMatch(int scanNumber, int scanPrecursorCharge) : base(null, double.NaN, double.NaN, double.NaN, scanNumber, scanPrecursorCharge, 0, double.NaN, double.NaN, double.NaN)
        {
        }

        #endregion Public Constructors

        #region Public Methods

        public override CompactPeptide GetCompactPeptide(List<MetaMorpheusModification> variableModifications, List<MetaMorpheusModification> localizeableModifications, List<MetaMorpheusModification> fixedModifications)
        {
            throw new NotImplementedException();
        }

        #endregion Public Methods

    }
}