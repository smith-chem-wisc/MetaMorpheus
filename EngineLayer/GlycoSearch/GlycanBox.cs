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

        public static Modification[] Global_NOGlycanMods;

        public static GlycanBox[] AllModBoxes;

        //
        public GlycanBox(int[] ids, string[] motifs, Glycan[] glycans, Modification[] modifications) : base(ids, motifs)
        {
            byte[] kind = new byte[Glycan.SugarLength];
            foreach (var id in ids)
            {
                for (int i = 0; i < kind.Length; i++)
                {
                    kind[i] += glycans[id].Kind[i];
                }

                if (modifications[id].FeatureType == "Nxs/t")
                {
                    NGlycanCount++;
                }
            }
            Kind = kind;
            Mass = (double)Glycan.GetMass(Kind) / 1E5;

            
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

        public int OGlycanCount { get; set; }

        public int NGlycanCount { get; set; }

        public static void SetMotifNeeded(GlycanBox glycanBox, Modification[] allModifications)
        {
            string[] motifs = new string[glycanBox.ModIds.Length];
            for (int i = 0; i < glycanBox.ModIds.Length; i++)
            {
                motifs[i] = allModifications[glycanBox.ModIds[i]].FeatureType;
            }
            glycanBox.ModMotfis = motifs;
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

        public static string[] GetGlycanBoxMotifs(int[] ids, Modification[] modifications)
        {
            string[] motifs = new string[ids.Length];
            for (int i = 0; i < ids.Length; i++)
            {
                motifs[i] = modifications[ids[i]].FeatureType;
            }
            return motifs;
        }

        //After O-glycans are read in from database, we build combinations of glycans into GlycanBox. The maxNum is maximum glycans allowed on one peptides.
        public static IEnumerable<GlycanBox> BuildGlycanBoxes(int maxNum, Glycan[] glycans, Modification[] modifications)
        {

            for (int i = 1; i <= maxNum; i++)
            {
                foreach (var idCombine in Glycan.GetKCombsWithRept(Enumerable.Range(0, glycans.Length), i))
                {
                    var motifs = GetGlycanBoxMotifs(idCombine.ToArray(), modifications);
                    GlycanBox glycanBox = new GlycanBox(idCombine.ToArray(), motifs, glycans, modifications);
                    SetMotifNeeded(glycanBox, modifications);
                    glycanBox.ChildGlycanBoxes = BuildChildGlycanBoxes(glycanBox.ModIds, glycans, modifications).ToArray();

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
        public static IEnumerable<GlycanBox> BuildChildGlycanBoxes(int[] glycanIds, Glycan[] glycans, Modification[] modifications)
        {
            yield return new GlycanBox(new int[0], new string[0], glycans, modifications);
            HashSet<string> seen = new HashSet<string>();
            for (int i = 1; i <= glycanIds.Length; i++)
            {
                foreach (var idCombine in Glycan.GetKCombs(Enumerable.Range(0, glycanIds.Length), i))
                {
                    List<int> ids = new List<int>();
                    foreach (var id in idCombine)
                    {
                        ids.Add(glycanIds[id]);
                    }

                    if (!seen.Contains(string.Join(",", ids.Select(p => p.ToString()))))
                    {
                        seen.Add(string.Join(",", ids.Select(p => p.ToString())));
                        var motifs = GetGlycanBoxMotifs(idCombine.ToArray(), modifications);
                        GlycanBox glycanBox = new GlycanBox(ids.ToArray(), motifs, glycans, modifications);
                        SetMotifNeeded(glycanBox, modifications);
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

        public static IEnumerable<GlycanBox> Build_NOGlycanBoxes(Glycan[] glycans, Modification[] modifications, int OGlycoMaxNum, int NGlycoMaxNum, int GlobalOGlycanNumber, int GlobalNGlycoNumber)
        {
            for (int i = 1; i <= OGlycoMaxNum; i++)
            {
                foreach (var idCombine in Glycan.GetKCombsWithRept(Enumerable.Range(0, GlobalOGlycanNumber), i))
                {
                    var motifs = GetGlycanBoxMotifs(idCombine.ToArray(), modifications);

                    GlycanBox o_modBox = new GlycanBox(idCombine.ToArray(), motifs, glycans, modifications);

                    yield return o_modBox;

                    for (int j = 1; j <= NGlycoMaxNum; j++)
                    {
                        foreach (var jdCombine in Glycan.GetKCombsWithRept(Enumerable.Range(GlobalOGlycanNumber, GlobalNGlycoNumber), j))
                        {
                            var n_motifs = GetGlycanBoxMotifs(jdCombine.ToArray(), modifications);

                            GlycanBox n_modBox = new GlycanBox(jdCombine.ToArray(), n_motifs, glycans, modifications);

                            yield return n_modBox;

                            var ijdCombine = idCombine.Concat(jdCombine);

                            var no_motifs = GetGlycanBoxMotifs(ijdCombine.ToArray(), modifications);

                            GlycanBox no_modBox = new GlycanBox(ijdCombine.ToArray(), no_motifs, glycans, modifications);

                            yield return no_modBox;

                        }
                    }
                }
            }

        }

    }
}
