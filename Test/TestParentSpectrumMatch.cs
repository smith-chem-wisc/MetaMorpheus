using EngineLayer;
using Proteomics;
using System;
using System.Collections.Generic;

namespace Test
{
    internal class TestParentSpectrumMatch : PsmParent
    {

        #region Public Constructors

        public TestParentSpectrumMatch(Ms2ScanWithSpecificMass scan) : base(scan, double.NaN, 1)
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