using System.Collections.Generic;
using System.IO;
using System.Linq;
using Chemistry;
using System;
using Proteomics;
using MassSpectrometry;
using Omics.Modifications;

namespace EngineLayer
{

    /// <summary>
    /// A defined combination of glycans to modify on one peptide. Ex. if we have 3 glycans on one peptide (g1,g2,g3), the GlycanBoxMass is the sum of the three glycans.(glycanBox: [g1,g2,g3])
    /// </summary>
    public class GlycanBox:ModBox
    {

            //TO DO: Decoy O-glycan can be created, but the results need to be reasoned.
            //public static int[] SugarShift = new int[]{ -16205282, -20307937, -29109542, -14605791, -30709033, -15005282, -36513219, -40615874, 16205282, 20307937, 29109542, 14605791, 30709033, 15005282, 36513219, 40615874 };
            private readonly static int[] SugarShift = new int[] //still unclear about the shift...
            {
                7103710, 10300920, 11502690, 12904260, 14706840, 5702150, 13705890, 12809500, 11308410, 13104050,
                11404290, 9705280, 12805860, 15610110, 8703200, 10104770, 9906840, 18607930, 16306330,
                -7103710, -10300920, -11502690, -12904260, -14706840, -5702150, -13705890, -12809500, -11308410, -13104050,
                -11404290, -9705280, -12805860, -15610110, -8703200, -10104770, -9906840, -18607930, -16306330,

            };

        /// <summary>
        /// Use the glycan from the database to create all possible combination glycans set into GlycanBox.
        /// </summary>
        /// <param name="maxNum"></param>
        /// <param name="buildDecoy"></param>
        /// <returns></returns>
        public static IEnumerable<GlycanBox> BuildModBoxes(int maxNum, bool buildDecoy = false)
        {
            for (int i = 1; i <= maxNum; i++)
            {
                foreach (var idCombine in Glycan.GetKCombsWithRept(Enumerable.Range(0, GlobalModifications.Length), i))
                {
                    GlycanBox glycanBox = new GlycanBox(idCombine.ToArray());
                    glycanBox.TargetDecoy = true;
                    glycanBox.ChildBoxes = BuildChildBoxes(glycanBox.NumberOfMods, glycanBox.ModIds, glycanBox.TargetDecoy).ToArray();

                    yield return glycanBox;

                    if (buildDecoy)
                    {
                        GlycanBox glycanBox_decoy = new GlycanBox(idCombine.ToArray(), false); // decoy glycanBox
                        glycanBox_decoy.TargetDecoy = false;
                        glycanBox_decoy.ChildBoxes = BuildChildBoxes(glycanBox.NumberOfMods, glycanBox.ModIds, glycanBox.TargetDecoy).ToArray();
                        yield return glycanBox_decoy;
                    }
                }
            }
        }


        public static IEnumerable<GlycanBox> BuildChildBoxes(int maxNum, int[] ids, bool targetDecoy = true)
        {
            return ModBox.BuildChildBoxes(maxNum, ids, targetDecoy).OfType<GlycanBox>();
        }

        /// <summary>
        /// Constructor of GlycanBox.
        /// </summary>
        /// <param name="ids"> The glycanBox composition, each number represent one glycan index in the database</param>
        /// <param name="targetDecoy"></param>
        public GlycanBox(int[] ids, bool Istarget = true):base(ids)
        {
            if (Istarget)
            {
                Mass = (double)Glycan.GetMass(Kind) / 1E5;
            }
            else // For decoy glycanBox, the mass is shifted by a random value.
            {
                Random random = new Random();
                int shiftInd = random.Next(SugarShift.Length);
                Mass = (double)(Glycan.GetMass(Kind) + SugarShift[shiftInd]) / 1E5;
            }
        }
        
        public string GlycanIdString // the composition of glycanBox. Example: [1,2,3] means glycan1 + glycan2 + glycan3 are on the peptide.
        {
            get
            {
                return string.Join(",", ModIds.Select(p => p.ToString()));
            }
        }
    }
}
