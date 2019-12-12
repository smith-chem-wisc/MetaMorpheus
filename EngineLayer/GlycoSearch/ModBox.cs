using System.Collections.Generic;
using System.IO;
using System.Linq;
using Chemistry;
using System;
using Proteomics;
using MassSpectrometry;
using Proteomics.ProteolyticDigestion;

namespace EngineLayer.GlycoSearch
{
    public class ModBox
    {
        public static Modification[] SelectedModifications;

        public static IEnumerable<ModBox> BuildModBoxes(int maxNum)
        {
            for (int i = 1; i <= maxNum; i++)
            {
                foreach (var idCombine in Glycan.GetKCombsWithRept(Enumerable.Range(0, SelectedModifications.Length), i))
                {
                    ModBox modBox = new ModBox(idCombine.ToArray());

                    yield return modBox;
                }
            }
        }

        public static List<Tuple<int, int>[]> GetPossibleModSites(PeptideWithSetModifications peptide, ModBox modBox)
        {
            List<Tuple<int, int>[]> allPossibleMods = new List<Tuple<int, int>[]>();

            //Possible positions For each motif
            List<List<Tuple<int, int>[]>> permutations = new List<List<Tuple<int, int>[]>>();

            foreach (var mn in modBox.MotifNeeded)
            {
                List<int> possibleModSites = new List<int>();

                Modification modWithMotif = new Modification(_target: mn.Key, _locationRestriction: "Anywhere.");

                for (int r = 0; r < peptide.Length; r++)
                {
                    if (peptide.AllModsOneIsNterminus.Keys.Contains(r + 2))
                    {
                        continue;
                    }

                    if (ModificationLocalization.ModFits(modWithMotif, peptide.BaseSequence, r + 1, peptide.Length, r + 1))
                    {
                        possibleModSites.Add(r + 2);
                    }
                }

                if (possibleModSites.Count >= mn.Value.Count)
                {
                    var per = GetPermutations(possibleModSites, mn.Value.ToArray());
                    permutations.Add(per);
                }
                else
                {
                    return allPossibleMods;
                }
            }

            int[] combines = permutations.Select(p => p.Count).ToArray();

            List<int>[] combinesList = new List<int>[combines.Length];

            for (int i = 0; i < combines.Length; i++)
            {
                int num = permutations[i].Count;
                combinesList[i] = Enumerable.Range(0, num).ToList();              
            }

            var kcoms = AllCombinationsOf(combinesList);

            foreach (var k in kcoms)
            {
                List<Tuple<int, int>> tuple = new List<Tuple<int, int>>();
                for (int i = 0; i < k.Count; i++)
                {

                    tuple.AddRange(permutations.ElementAt(i).ElementAt(k[i]));
                }
                allPossibleMods.Add(tuple.ToArray());
            }


            return allPossibleMods;
        }

        public static List<Tuple<int, int>[]> GetPermutations(List<int> allModPos, int[] modBoxId)
        {
            var length = modBoxId.Length;
            var indexes = Enumerable.Range(0, length).ToArray();
            int[] orderMod = new int[length];

            List<Tuple<int, int>[]> permutateModPositions = new List<Tuple<int, int>[]>();

            var combinations = Glycan.GetKCombs(allModPos, length);

            foreach (var com in combinations)
            {
                var permutation = Glycan.GetPermutations(com, length);

                HashSet<string> keys = new HashSet<string>();

                foreach (var per in permutation)
                {
                    Array.Sort(indexes);

                    var orderedPer = per.ToArray();
                    Array.Sort(orderedPer, indexes);

                    for (int i = 0; i < length; i++)
                    {
                        orderMod[i] = modBoxId[indexes[i]];
                    }
                    var key = string.Join(",", orderMod.Select(p => p.ToString()));
                    if (!keys.Contains(key))
                    {
                        keys.Add(key);

                        Tuple<int, int>[] t = new Tuple<int, int>[length];
                        for (int i = 0; i < length; i++)
                        {
                            t[i] = new Tuple<int, int>(per.ElementAt(i), modBoxId[i]);
                        }
                        permutateModPositions.Add(t);
                    }
                }
            }

            return permutateModPositions;
        }

        public static List<List<T>> AllCombinationsOf<T>(params List<T>[] sets)
        {
            // need array bounds checking etc for production
            var combinations = new List<List<T>>();

            // prime the data
            foreach (var value in sets[0])
                combinations.Add(new List<T> { value });

            foreach (var set in sets.Skip(1))
                combinations = AddExtraSet(combinations, set);

            return combinations;
        }

        private static List<List<T>> AddExtraSet<T>
             (List<List<T>> combinations, List<T> set)
        {
            var newCombinations = from value in set
                                  from combination in combinations
                                  select new List<T>(combination) { value };

            return newCombinations.ToList();
        }

        public static PeptideWithSetModifications GetTheoreticalPeptide(Tuple<int, int>[] theModPositions, PeptideWithSetModifications peptide, ModBox modBox)
        {
            Dictionary<int, Modification> testMods = new Dictionary<int, Modification>();
            foreach (var mod in peptide.AllModsOneIsNterminus)
            {
                testMods.Add(mod.Key, mod.Value);
            }

            for (int i = 0; i < theModPositions.Count(); i++)
            {
                testMods.Add(theModPositions.ElementAt(i).Item1, SelectedModifications[theModPositions.ElementAt(i).Item2]);
            }

            var testPeptide = new PeptideWithSetModifications(peptide.Protein, peptide.DigestionParams, peptide.OneBasedStartResidueInProtein,
                peptide.OneBasedEndResidueInProtein, peptide.CleavageSpecificityForFdrCategory, peptide.PeptideDescription, peptide.MissedCleavages, testMods, peptide.NumFixedMods);

            return testPeptide;
        }

        public ModBox(int[] ids)
        {
            ModIds = ids;
            NumberOfMods = ids.Length;
        }

        public int[] ModIds { get;  }

        public int NumberOfMods { get; }

        public double Mass
        {
            get
            {
                double mass = 0;
                foreach (var id in ModIds)
                {
                    mass += SelectedModifications[id].MonoisotopicMass.Value;
                }
                return mass;
            }
        }

        public Dictionary<ModificationMotif, List<int>> MotifNeeded
        {
            get
            {
                Dictionary<ModificationMotif, List<int>> aa = new Dictionary<ModificationMotif, List<int>>();
                foreach (var id in ModIds)
                {
                    var mod = SelectedModifications[id];
                    if (aa.ContainsKey(mod.Target))
                    {
                        aa[mod.Target].Add(id);
                    }
                    else
                    {
                        aa.Add(mod.Target, new List<int>(id));
                    }
                }
                return aa;
            }
        }
    }
}
