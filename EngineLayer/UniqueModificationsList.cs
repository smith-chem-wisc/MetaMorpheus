using System;
using System.Collections.Generic;

namespace EngineLayer
{
    public class UniqueModificationsCollection : List<MorpheusModification>
    {

        #region Internal Methods

        internal new void Add(MorpheusModification mod)
        {
            foreach (MorpheusModification modHere in this)
            {
                if (Math.Abs(modHere.PrecursorMassShift - mod.PrecursorMassShift) < 0.001)
                {
                    return;
                }
            }
            base.Add(mod);
        }

        #endregion Internal Methods

    }
}