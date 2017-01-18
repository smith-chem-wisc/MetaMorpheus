using System;
using System.Collections.Generic;

namespace OldInternalLogic
{
    public class UniqueModificationsCollection : List<MorpheusModification>
    {
        internal new void Add(MorpheusModification mod)
        {
            foreach (MorpheusModification modHere in this)
            {
                if (Math.Abs(modHere.MonoisotopicMassShift - mod.MonoisotopicMassShift) < 0.001)
                {
                    return;
                }
            }
            base.Add(mod);
        }
    }
}