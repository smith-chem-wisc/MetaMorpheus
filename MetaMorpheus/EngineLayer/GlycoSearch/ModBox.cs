using MzLibUtil;

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

        // Simplified constructor logic for combining oGlycanIds and nGlycanId
            public ModBox(int[] oGlycanIds, int nGlycanId)
            {
                if (nGlycanId != 0)
                {
                    if (oGlycanIds.IsNullOrEmpty())
                    {
                        ModIds = [nGlycanId];
                    }
                    else
                    {
                        ModIds = [.. oGlycanIds, nGlycanId];
                    }
                }
                else
                {
                    ModIds = oGlycanIds;
                }
            NumberOfMods = ModIds.Length;
            TargetDecoy = true;
        }

        public int[] ModIds { get;  }

        public int NumberOfMods { get; }

        public double Mass { get; set; }

        public double DecoyMass { get; set; }

        public bool TargetDecoy { get; set; }

    }
}
