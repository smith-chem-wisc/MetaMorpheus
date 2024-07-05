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
        public static Glycan[] GlobalOGlycans { get; set; } // The glycan list in the database file

        public GlycanBox[] ChildGlycanBoxes { get; set; }   // all possible glycan combinations in the glycanBox

        public static Modification[] GlobalOGlycanModifications { get; set; }

        public static GlycanBox[] OGlycanBoxes { get; set; } // all possible glycan boxes

        public byte[] Kind { get; private set; } 

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
        /// Use the glycan from database to create all possible combination glycan set into GlycanBox. 
        /// </summary>
        /// <param name="maxNum"> The maxNum is maximum glycans allowed on one peptides </param>
        /// <returns> The glycanBox collection, glycanBox[]</returns>
        public static IEnumerable<GlycanBox> BuildOGlycanBoxes(int maxNum)
        {
            return BuildOGlycanBoxes(maxNum, false);
        }
        public static IEnumerable<GlycanBox> BuildOGlycanBoxes(int maxNum, bool buildDecoy)
        {

            for (int i = 1; i <= maxNum; i++)
            {
                foreach (var idCombine in Glycan.GetKCombsWithRept(Enumerable.Range(0, GlobalOGlycans.Length), i))
                {
                    GlycanBox glycanBox = new GlycanBox(idCombine.ToArray());
                    glycanBox.TargetDecoy = true;
                    glycanBox.ChildGlycanBoxes = BuildChildOGlycanBoxes(glycanBox.NumberOfMods, glycanBox.ModIds, glycanBox.TargetDecoy).ToArray();

                    yield return glycanBox;

                    if (buildDecoy)
                    {
                        GlycanBox glycanBox_decoy = new GlycanBox(idCombine.ToArray(),false); // decoy glycanBox
                        glycanBox_decoy.TargetDecoy = false;
                        glycanBox_decoy.ChildGlycanBoxes = BuildChildOGlycanBoxes(glycanBox_decoy.NumberOfMods, glycanBox_decoy.ModIds, glycanBox_decoy.TargetDecoy).ToArray();
                        yield return glycanBox_decoy;
                    }
                }
            }
        }

        /// <summary>
        /// Convert the glycan into Modification type for MetaMorpheus to manipulate sequences. In the future we may able to combine the two type together.
        /// </summary>
        /// <param name="globalOGlycans"></param>
        /// <returns></returns>
        public static Modification[] BuildGlobalOGlycanModifications(Glycan[] globalOGlycans)
        {
            Modification[] globalOGlycanModifications = new Modification[globalOGlycans.Length];

            for (int i = 0; i < GlobalOGlycans.Length; i++)
            {
                globalOGlycanModifications[i] = Glycan.OGlycanToModification(globalOGlycans[i]);
            }
            return globalOGlycanModifications;
        }


        /// <summary>
        /// Generate all possible child/fragment box of the specific glycanBox. The childBoxes is uesd for LocalizationGraph.
        /// </summary>
        /// <param name="maxNum"></param>
        /// <param name="glycanIds"> The glycanBox, ex. [0,0,1] means glycan0 + glycan0 + glycan1 </param>
        /// <param name="targetDecoy"></param>
        /// <returns> The ChildBox collection, ChildBox[] </returns>
        public static IEnumerable<GlycanBox> BuildChildOGlycanBoxes(int maxNum, int[] glycanIds, bool targetDecoy = true)
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
                        ids.Add(glycanIds[id]);      
                    }

                    if (!seen.Contains(string.Join(",", ids.Select(p => p.ToString()))))
                    {
                        seen.Add(string.Join(",", ids.Select(p => p.ToString())));

                        GlycanBox glycanBox = new GlycanBox(ids.ToArray(), targetDecoy);

                        yield return glycanBox;
                    }

                }
            }
        }

        /// <summary>
        /// Constructor of GlycanBox.
        /// </summary>
        /// <param name="ids"> The glycanBox composition, each number represent one glycan index in the database</param>
        /// <param name="targetDecoy"></param>
        public GlycanBox(int[] ids, bool targetDecoy = true):base(ids)
        {
            byte[] kind = new byte[] { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
            foreach (var id in ModIds) //ModIds is the same as ids.
            {
                for (int i = 0; i < kind.Length; i++)   
                {
                    kind[i] += GlobalOGlycans[id].Kind[i]; //kind is the sum of all glycan Kind in the Box.
                }
            }
            Kind = kind;

            if (targetDecoy)
            {
                Mass = (double)Glycan.GetMass(Kind) / 1E5;
            }
            else
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
