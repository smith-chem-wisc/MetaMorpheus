using InternalLogicEngineLayer;
using OldInternalLogic;
using System;
using System.Collections.Generic;

namespace Test
{
    internal class TestParentSpectrumMatch : ParentSpectrumMatch
    {

        #region Public Constructors

        public TestParentSpectrumMatch(int scanNumber, int scanPrecursorCharge)
        {
            this.scanNumber = scanNumber;
            this.scanPrecursorCharge = scanPrecursorCharge;
        }

        #endregion Public Constructors

        #region Public Methods

        public override CompactPeptide GetCompactPeptide(List<MorpheusModification> variableModifications, List<MorpheusModification> localizeableModifications)
        {
            throw new NotImplementedException();
        }

        #endregion Public Methods

    }
}