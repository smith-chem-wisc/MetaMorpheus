using EngineLayer;
using Proteomics;
using System;
using System.Collections.Generic;

namespace Test
{
    internal class TestParentSpectrumMatch : PsmParent
    {

        #region Public Constructors

        public TestParentSpectrumMatch(int scanNumber, int scanPrecursorCharge, int precursorScanNumber) : base(null, double.NaN, double.NaN, double.NaN, scanNumber, precursorScanNumber, scanPrecursorCharge, 0, double.NaN, double.NaN, double.NaN, 1)
        {
        }

        #endregion Public Constructors

        #region Public Methods

        public override CompactPeptide GetCompactPeptide(Dictionary<ModificationWithMass, ushort> modsDictionary)
        {
            throw new NotImplementedException();
        }

        #endregion Public Methods

    }
}