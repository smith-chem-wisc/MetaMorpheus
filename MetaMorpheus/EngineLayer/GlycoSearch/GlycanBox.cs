using System.Collections.Generic;
using System.IO;
using System.Linq;
using Chemistry;
using System;
using Proteomics;
using MassSpectrometry;

namespace EngineLayer
{
    //One peptide can have several o-glycans. The combined glycans are grouped as a glycan box. Used for localization. 
    //GlycanBox -- A defined combination of glycans will be considered to modify on one peptide. 
    //The GlycanBoxMass is the total mass of all glycans on the peptide
    public class GlycanBox:ModBox
    {
        public static Glycan[] GlobalOGlycans { get; set; }

        public static Modification[] GlobalOGlycanModifications { get; set; }

        public static GlycanBox[] OGlycanBoxes { get; set; }

        //TO DO: Decoy O-glycan can be created, but the results need to be reasoned.
        //public static int[] SugarShift = new int[]{ -16205282, -20307937, -29109542, -14605791, -30709033, -15005282, -36513219, -40615874, 16205282, 20307937, 29109542, 14605791, 30709033, 15005282, 36513219, 40615874 };
        private readonly static int[] SugarShift = new int[] 
        {
            7103710, 10300920, 11502690, 12904260, 14706840, 5702150, 13705890, 12809500, 11308410, 13104050,
            11404290, 9705280, 12805860, 15610110, 8703200, 10104770, 9906840, 18607930, 16306330,
            -7103710, -10300920, -11502690, -12904260, -14706840, -5702150, -13705890, -12809500, -11308410, -13104050,
            -11404290, -9705280, -12805860, -15610110, -8703200, -10104770, -9906840, -18607930, -16306330,

        };

        //After O-glycans are read in from database, we build combinations of glycans into GlycanBox. The maxNum is maximum glycans allowed on one peptides.
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
                        GlycanBox glycanBox_decoy = new GlycanBox(idCombine.ToArray());
                        glycanBox_decoy.TargetDecoy = false;
                        glycanBox_decoy.ChildGlycanBoxes = BuildChildOGlycanBoxes(glycanBox_decoy.NumberOfMods, glycanBox_decoy.ModIds, glycanBox_decoy.TargetDecoy).ToArray();
                        yield return glycanBox_decoy;
                    }
                }
            }
        }

        //After O-glycans are read in from database, we transfer the glycans into 'Modification' class type for MetaMorpheus to manipulate sequences.
        //In the future we may able to combine the two type together. 
        public static Modification[] BuildGlobalOGlycanModifications(Glycan[] globalOGlycans)
        {
            Modification[] globalOGlycanModifications = new Modification[globalOGlycans.Length];

            for (int i = 0; i < GlobalOGlycans.Length; i++)
            {
                globalOGlycanModifications[i] = Glycan.OGlycanToModification(globalOGlycans[i]);
            }
            return globalOGlycanModifications;
        }

        //The function here is to build GlycanBoxes used for LocalizationGraph. 
        //In LocalizationGraph matrix, for each AdjNode, it represent a ChildOGlycanBox here at certain glycosite.
        public static IEnumerable<GlycanBox> BuildChildOGlycanBoxes(int maxNum, int[] glycanIds, bool targetDecoy = true)
        {
            yield return new GlycanBox(new int[0], targetDecoy);
            HashSet<string> seen = new HashSet<string>();
            for (int i = 1; i <= maxNum; i++)
            {
                foreach (var idCombine in Glycan.GetKCombs(Enumerable.Range(0, maxNum), i))
                {
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

        public GlycanBox(int[] ids, bool targetDecoy = true):base(ids)
        {
            byte[] kind = new byte[] { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
            foreach (var id in ModIds)
            {
                for (int i = 0; i < kind.Length; i++)
                {
                    kind[i] += GlobalOGlycans[id].Kind[i];
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

        public GlycanBox[] ChildGlycanBoxes { get; set; }

        public string GlycanIdString
        {
            get
            {
                return string.Join(",", ModIds.Select(p => p.ToString()));
            }
        }

        public byte[] Kind{ get; private set; }

    }
}
