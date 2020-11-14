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
        //Static members
        public static Glycan[] GlobalOGlycans { get; set; }

        public static Modification[] GlobalOGlycanModifications { get; set; }

        public static GlycanBox[] OGlycanBoxes { get; set; }

        public static Glycan[] Global_NGlycans { get; set; }

        public static Modification[] Global_NGlycanModifications { get; set; }

        public static GlycanBox[] NGlycanBoxes { get; set; }

        public static Glycan[] Global_NOGlycans { get; set; }

        public static Modification[] AllModifications;

        public static GlycanBox[] AllModBoxes;

        //
        public GlycanBox(int[] ids, Glycan[] glycans) : base(ids)
        {
            byte[] kind = new byte[Glycan.SugarLength];
            foreach (var id in ModIds)
            {
                for (int i = 0; i < kind.Length; i++)
                {
                    kind[i] += glycans[id].Kind[i];
                }
            }
            Kind = kind;
            Mass = (double)Glycan.GetMass(Kind) / 1E5;

        }

        public GlycanBox(int[] ids, Glycan[] glycans, int glycoFlag) : base(ids)
        {
            byte[] kind = new byte[Glycan.SugarLength];
            foreach (var id in ModIds)
            {
                for (int i = 0; i < kind.Length; i++)
                {
                    kind[i] += glycans[id].Kind[i];
                }
            }
            Kind = kind;

            Mass = (double)Glycan.GetMass(Kind) / 1E5;

            GlycoFlag = glycoFlag;
        }

        public GlycanBox[] ChildGlycanBoxes { get; set; }

        public string GlycanIdString
        {
            get
            {
                return string.Join(",", ModIds.Select(p => p.ToString()));
            }
        }

        public byte[] Kind { get; private set; }

        public int GlycoFlag { get; set; }

        public int OGlycanCount { get; set; }

        public int NGlycanCount { get; set; }

        //key: motif, value: all ids for this motif
        public Dictionary<string, List<int>> MotifNeeded { get; set; }

        public static void SetMotifNeeded(GlycanBox glycanBox, Modification[] allModifications)
        {
            Dictionary<string, List<int>> aa = new Dictionary<string, List<int>>();
            foreach (var id in glycanBox.ModIds)
            {
                var mod = allModifications[id];
                if (aa.ContainsKey(mod.FeatureType))
                {
                    aa[mod.FeatureType].Add(id);
                }
                else
                {
                    aa.Add(mod.FeatureType, new List<int> { id });
                }
            }
            glycanBox.MotifNeeded = aa;
        }

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
        public static IEnumerable<GlycanBox> BuildGlycanBoxes(int maxNum, Glycan[] glycans, Modification[] allModifications)
        {

            for (int i = 1; i <= maxNum; i++)
            {
                foreach (var idCombine in Glycan.GetKCombsWithRept(Enumerable.Range(0, glycans.Length), i))
                {
                    GlycanBox glycanBox = new GlycanBox(idCombine.ToArray(), glycans);
                    SetMotifNeeded(glycanBox, allModifications);
                    glycanBox.TargetDecoy = true;
                    glycanBox.ChildGlycanBoxes = BuildChildGlycanBoxes(glycanBox.NumberOfMods, glycanBox.ModIds, glycans).ToArray();

                    yield return glycanBox;
                }
            }
        }

        //After O-glycans are read in from database, we transfer the glycans into 'Modification' class type for MetaMorpheus to manipulate sequences.
        //In the future we may able to combine the two type together. 
        public static IEnumerable<Modification> BuildGlobalOGlycanModifications(Glycan[] globalOGlycans)
        {
            for (int i = 0; i < globalOGlycans.Length; i++)
            {
                var ogMod = Glycan.OGlycanToModification(globalOGlycans[i]);
                yield return ogMod;
            }

        }

        //The function here is to build GlycanBoxes used for LocalizationGraph. 
        //In LocalizationGraph matrix, for each AdjNode, it represent a ChildOGlycanBox here at certain glycosite.
        public static IEnumerable<GlycanBox> BuildChildGlycanBoxes(int maxNum, int[] glycanIds, Glycan[] glycans)
        {
            yield return new GlycanBox(new int[0], glycans);
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

                        GlycanBox glycanBox = new GlycanBox(ids.ToArray(), glycans);

                        yield return glycanBox;
                    }

                }
            }
        }

        public static IEnumerable<Modification> BuildGlobal_NGlycanModifications(Glycan[] globalNGlycans)
        {
            for (int i = 0; i < globalNGlycans.Length; i++)
            {
                var ogMod = Glycan.NGlycanToModification(globalNGlycans[i]);
                yield return ogMod;
            }

        }

        public static IEnumerable<GlycanBox> Build_NOGlycanBoxes(Glycan[] glycans, int OGlycoMaxNum, int NGlycoMaxNum, int GlobalOGlycanNumber, int GlobalNGlycoNumber)
        {
            for (int i = 1; i <= OGlycoMaxNum; i++)
            {
                foreach (var idCombine in Glycan.GetKCombsWithRept(Enumerable.Range(0, GlobalOGlycanNumber), i))
                {
                    GlycanBox o_modBox = new GlycanBox(idCombine.ToArray(), glycans, 1);

                    yield return o_modBox;

                    for (int j = 1; j <= NGlycoMaxNum; j++)
                    {
                        foreach (var jdCombine in Glycan.GetKCombsWithRept(Enumerable.Range(GlobalOGlycanNumber, GlobalNGlycoNumber), j))
                        {
                            GlycanBox n_modBox = new GlycanBox(jdCombine.ToArray(), glycans, 2);

                            yield return n_modBox;

                            var ijdCombine = idCombine.Concat(jdCombine);

                            GlycanBox no_modBox = new GlycanBox(ijdCombine.ToArray(), glycans, 3);

                            yield return no_modBox;

                        }
                    }
                }
            }

        }

    }
}
