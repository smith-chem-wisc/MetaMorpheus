using System.Collections.Generic;
using System.Linq;
using Omics.Modifications;

namespace EngineLayer
{
    public class ModBox //The superclass of GlycanBox
    {
        //One peptide can have several modifications. The combined modifications are grouped as a modification box. Used for localization. 
        //ModBox -- a defined combination of modifications will be considered to modify on one peptide. The box means the combined group of modification. 
        public ModBox(int[] ids, bool targetDecoy = true)
        {
            ModIds = ids;
            NumberOfMods = ids.Length;
            TargetDecoy = targetDecoy;
        }

        public int[] ModIds { get;  }

        public int NumberOfMods { get; }

        public double Mass { get; set; }

        public double DecoyMass { get; set; }

        public bool TargetDecoy { get; set; }

        public static Modification[] GlobalModifications { get; set; }

        public ModBox[] ChildModBoxes { get; set; }


        public static IEnumerable<ModBox> BuildModBoxes(int maxNum,  bool buildDecoy = false)
        {
            for (int i = 1; i <= maxNum; i++)
            {
                foreach (var idCombine in Glycan.GetKCombsWithRept(Enumerable.Range(0, GlobalModifications.Length), i))
                {
                    ModBox modBox = new ModBox(idCombine.ToArray());
                    modBox.TargetDecoy = true;
                    modBox.ChildModBoxes = BuildChildModBoxes(modBox.NumberOfMods, modBox.ModIds, modBox.TargetDecoy).ToArray();

                    yield return modBox;

                    if (buildDecoy)
                    {
                        ModBox modBox_decoy = new ModBox(idCombine.ToArray(), false); // decoy glycanBox
                        modBox_decoy.TargetDecoy = false;
                        modBox_decoy.ChildModBoxes = BuildChildModBoxes(modBox.NumberOfMods, modBox.ModIds, modBox.TargetDecoy).ToArray();
                        yield return modBox_decoy;
                    }
                }
            }
        }

        public static IEnumerable<ModBox> BuildChildModBoxes(int maxNum, int[] modIds, bool targetDecoy = true)
        {
            yield return new GlycanBox(new int[0], targetDecoy);
            HashSet<string> seen = new HashSet<string>();
            for (int i = 1; i <= maxNum; i++)
            {
                foreach (var idCombine in Glycan.GetKCombs(Enumerable.Range(0, maxNum), i)) //get all combinations of glycans on the peptide, ex. we have three glycosite and three glycan maybe on that (A,B,C) 
                {                                                                           //the combination of glycans on the peptide can be (A),(A+B),(A+C),(B+C),(A+B+C) totally six 
                    List<int> ids = new List<int>();
                    foreach (var id in idCombine)
                    {
                        ids.Add(modIds[id]);
                    }

                    if (!seen.Contains(string.Join(",", ids.Select(p => p.ToString()))))
                    {
                        seen.Add(string.Join(",", ids.Select(p => p.ToString())));

                        ModBox modBox = new ModBox(ids.ToArray(), targetDecoy);

                        yield return modBox;
                    }

                }
            }
        }


    }
}
