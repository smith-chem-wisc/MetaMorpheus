using System;
using System.Linq;
using System.Collections.Generic;

namespace EngineLayer
{
    public class ModBox
    {
        //One peptide can have several modifications. The combined modifications are grouped as a modification box. Used for localization. 
        //ModBox -- a defined combination of modifications will be considered to modify on one peptide. The box means the combined group of modification. 
        public ModBox(int[] ids)
        {
            ModIds = ids;
            NumberOfMods = ids.Length;
            TargetDecoy = true;
        }

        public int[] ModIds { get;  }

        public int NumberOfMods { get; }

        public double Mass { get; set; }

        public double DecoyMass { get; set; }

        public bool TargetDecoy { get; set; }

        public static IEnumerable<ModBox> BuildChildModBoxes(int[] modIds)
        {
            yield return new ModBox(new int[0]);
            HashSet<string> seen = new HashSet<string>();
            var maxNum = modIds.Length;
            for (int i = 1; i <= maxNum; i++)
            {
                foreach (var idCombine in Glycan.GetKCombs(Enumerable.Range(0, maxNum), i))
                {
                    List<int> ids = new List<int>();
                    foreach (var id in idCombine)
                    {
                        ids.Add(modIds[id]);
                    }

                    if (!seen.Contains(string.Join(",", ids.Select(p => p.ToString()))))
                    {
                        seen.Add(string.Join(",", ids.Select(p => p.ToString())));

                        ModBox modBox = new ModBox(ids.ToArray());

                        yield return modBox;
                    }

                }
            }
        }

    }
}
