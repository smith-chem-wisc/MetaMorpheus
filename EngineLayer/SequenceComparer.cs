using System.Collections.Generic;

namespace EngineLayer
{
    internal class SequenceComparer : IEqualityComparer<PSMwithProteinHashSet>
    {

        #region Public Methods

        public bool Equals(PSMwithProteinHashSet x, PSMwithProteinHashSet y)
        {
            return x.FullSequence.Equals(y.FullSequence);
        }

        public int GetHashCode(PSMwithProteinHashSet obj)
        {
            return obj.FullSequence.GetHashCode();
        }

        #endregion Public Methods

    }
}