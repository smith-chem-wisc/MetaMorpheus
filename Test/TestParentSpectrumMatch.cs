using InternalLogicEngineLayer;
using OldInternalLogic;
using System;
using System.Collections.Generic;

namespace Test
{
    internal class TestParentSpectrumMatch : ParentSpectrumMatch
    {

        #region Public Constructors

        public TestParentSpectrumMatch(int scanNumber, int scanPrecursorCharge) : base(null, double.NaN, double.NaN, double.NaN, scanNumber, scanPrecursorCharge, 0, double.NaN, double.NaN, double.NaN)
        {
        }

        #endregion Public Constructors

        #region Public Methods

        public override CompactPeptide GetCompactPeptide(List<MorpheusModification> variableModifications, List<MorpheusModification> localizeableModifications, List<MorpheusModification> fixedModifications)
        {
            throw new NotImplementedException();
        }

        #endregion Public Methods

    }
}