using Omics.Modifications;
using System;
using System.Collections.Generic;
using System.Linq;

namespace EngineLayer
{
    public class ModBox //The superclass of GlycanBox
    {
        //One peptide can have several modifications. The combined modifications are grouped as a modification box. Used for localization. 
        //ModBox -- a defined combination of modifications will be considered to modify on one peptide. The box means the combined group of modification. 
        public ModBox(int[] ids, bool targetDecoy = true)
        {
            // ModBox setting
            ModIds = ids;
            NumberOfMods = ids.Length;
            TargetDecoy = targetDecoy;

            // Glycan setting
            byte[] kind = new byte[] { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
            foreach (var id in ModIds.Where(id=> GlobalModifications[id] is Glycan)) // iterate through all glycan in the box to generate the glycan kind
            {
                for (int i = 0; i < kind.Length; i++)
                {
                    kind[i] += GlobalOGlycans[id].Kind[i]; //kind is the sum of all glycan Kind in the Box.
                }
            }
            Kind = kind;
        }

        public int[] ModIds { get;  }
        public int NumberOfMods { get; }
        public double Mass
        {
            get
            {
                if (TargetDecoy)
                {
                    return Glycan.GetMass(Kind) / 1E5;
                }
                // if decoy, we randomly shift the mass of the glycan by a random value from SugarShift.
                Random random = new Random();
                int shiftInd = random.Next(SugarShift.Length);
                return (Glycan.GetMass(Kind) + SugarShift[shiftInd]) / 1E5;
            }
        }
        public double DecoyMass { get; set; }
        public bool TargetDecoy { get; set; }
        public static Modification[] GlobalModifications { get; set; }
        public ModBox[] ChildBoxes { get; set; }

        // Glycans properties
        public byte[] Kind { get; private set; }
        /// <summary>
        /// The global list of O-glycans loaded from the glycan database file.
        /// </summary>
        public static Glycan[] GlobalOGlycans
        {
            get
            {
                return GlobalModifications.Where(p => p is Glycan).Cast<Glycan>().ToArray();
            }
        }

        //TO DO: Decoy O-glycan can be created, but the results need to be reasoned.
        //public static int[] SugarShift = new int[]{ -16205282, -20307937, -29109542, -14605791, -30709033, -15005282, -36513219, -40615874, 16205282, 20307937, 29109542, 14605791, 30709033, 15005282, 36513219, 40615874 };
        protected readonly static int[] SugarShift = new int[] //still unclear about the shift...
        {
            7103710, 10300920, 11502690, 12904260, 14706840, 5702150, 13705890, 12809500, 11308410, 13104050,
            11404290, 9705280, 12805860, 15610110, 8703200, 10104770, 9906840, 18607930, 16306330,
            -7103710, -10300920, -11502690, -12904260, -14706840, -5702150, -13705890, -12809500, -11308410, -13104050,
            -11404290, -9705280, -12805860, -15610110, -8703200, -10104770, -9906840, -18607930, -16306330,

        };


        /// <summary>
        /// Use the Mod from database to create all possible combination Mods set into ModBox. 
        /// </summary>
        /// <param name="maxNum"> The maxNum is maximum glycans allowed on one peptides </param>
        /// <returns> The glycanBox collection, glycanBox[]</returns>
        public static IEnumerable<ModBox> BuildModBoxes(int maxNum,  bool buildDecoy = false)
        {
            for (int i = 1; i <= maxNum; i++)
            {
                foreach (var idCombine in Glycan.GetKCombsWithRept(Enumerable.Range(0, GlobalModifications.Length), i))
                {
                    ModBox modBox = new ModBox(idCombine.ToArray());
                    modBox.TargetDecoy = true;
                    modBox.ChildBoxes = BuildChildBoxes(modBox.NumberOfMods, modBox.ModIds, modBox.TargetDecoy).ToArray();

                    yield return modBox;

                    if (buildDecoy)
                    {
                        ModBox modBox_decoy = new ModBox(idCombine.ToArray(), false); // decoy glycanBox
                        modBox_decoy.TargetDecoy = false;
                        modBox_decoy.ChildBoxes = BuildChildBoxes(modBox.NumberOfMods, modBox.ModIds, modBox.TargetDecoy).ToArray();
                        yield return modBox_decoy;
                    }
                }
            }
        }


        /// <summary>
        /// Generate all possible child/fragment box of the specific glycanBox. The childBoxes is uesd for LocalizationGraph.
        /// </summary>
        /// <param name="maxNum"></param>
        /// <param name="glycanIds"> The glycanBox, ex. [0,0,1] means glycan0 + glycan0 + glycan1 </param>
        /// <param name="targetDecoy"></param>
        /// <returns> The ChildBox collection, ChildBox[] </returns>
        public static IEnumerable<ModBox> BuildChildBoxes(int maxNum, int[] modIds, bool targetDecoy = true)
        {
            yield return new ModBox(new int[0], targetDecoy);
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
