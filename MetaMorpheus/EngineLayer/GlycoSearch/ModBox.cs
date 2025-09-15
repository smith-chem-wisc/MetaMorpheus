using System;
using Easy.Common.Extensions;

namespace EngineLayer
{
    public class ModBox //The superclass of GlycanBox
    {
        //One peptide can have several modifications. The combined modifications are grouped as a modification box. Used for localization. 
        //ModBox -- a defined combination of modifications will be considered to modify on one peptide. The box means the combined group of modification. 
        public ModBox(int[] ids)
        {
            ModIds = ids;
            NumberOfMods = ids.Length;
            TargetDecoy = true;
        }

        public ModBox(int[] oGlycanIds, int nGlycanId)
        {
            if (nGlycanId != 0) // if there is an N-glycan
            {
                int[] newArray = new int[oGlycanIds.Length + 1];
                Array.Copy(oGlycanIds, newArray, oGlycanIds.Length);
                newArray[oGlycanIds.Length] = nGlycanId;
                ModIds = newArray;
            }
            else // if there is no N-glycan
            {
                ModIds = oGlycanIds;
            }
            NumberOfMods = nGlycanId == 0 ? oGlycanIds.Length : oGlycanIds.Length + 1;
            TargetDecoy = true;
        }

        public int[] ModIds { get;  }

        public int NumberOfMods { get; }

        public double Mass { get; set; }

        public double DecoyMass { get; set; }

        public bool TargetDecoy { get; set; }

    }
}
