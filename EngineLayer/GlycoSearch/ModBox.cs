using System.Collections.Generic;
using System.IO;
using System.Linq;
using Chemistry;
using System;
using Proteomics;
using MassSpectrometry;
using Proteomics.ProteolyticDigestion;
using Proteomics.Fragmentation;

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

        public static int[] GetAllPossibleModSites(PeptideWithSetModifications peptide, ModBox modBox)
        {
            List<int> possibleModSites = new List<int>();

            foreach (var mn in modBox.MotifNeeded)
            {
               
                ModificationMotif.TryGetMotif(mn.Key, out ModificationMotif motif);
                Modification modWithMotif = new Modification(_target: motif, _locationRestriction: "Anywhere.");

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

                if (possibleModSites.Count < mn.Value.Count)
                {
                    return null;
                }
            }
            return possibleModSites.OrderBy(p=>p).ToArray();
        }

        //Tuple<int modpos, int modId>[] is one possible way the peptide is modified
        public static List<Tuple<int, int>[]> GetPossibleModSites(PeptideWithSetModifications peptide, ModBox modBox)
        {
            List<Tuple<int, int>[]> allPossibleMods = new List<Tuple<int, int>[]>();

            //Possible positions For each motif
            List<List<Tuple<int, int>[]>> permutations = new List<List<Tuple<int, int>[]>>();

            foreach (var mn in modBox.MotifNeeded)
            {
                List<int> possibleModSites = new List<int>();

                ModificationMotif.TryGetMotif(mn.Key, out ModificationMotif motif);
                Modification modWithMotif = new Modification(_target: motif, _locationRestriction: "Anywhere.");

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

        //Tuple<int modpos, int modId>
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

        //Combination of List<List<int>>. See: https://stackoverflow.com/questions/545703/combination-of-listlistint
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

        public static PeptideWithSetModifications GetTheoreticalPeptide(int[] theModPositions, PeptideWithSetModifications peptide, ModBox modBox)
        {
            Modification[] modifications = new Modification[modBox.NumberOfMods];
            for (int i = 0; i < modBox.NumberOfMods; i++)
            {
                modifications[i] = ModBox.SelectedModifications[modBox.ModIds.ElementAt(i)];
            }

            Dictionary<int, Modification> testMods = new Dictionary<int, Modification>();
            foreach (var mod in peptide.AllModsOneIsNterminus)
            {
                testMods.Add(mod.Key, mod.Value);
            }

            for (int i = 0; i < theModPositions.Count(); i++)
            {
                testMods.Add(theModPositions.ElementAt(i), modifications[i]);
            }

            var testPeptide = new PeptideWithSetModifications(peptide.Protein, peptide.DigestionParams, peptide.OneBasedStartResidueInProtein,
                peptide.OneBasedEndResidueInProtein, peptide.CleavageSpecificityForFdrCategory, peptide.PeptideDescription, peptide.MissedCleavages, testMods, peptide.NumFixedMods);

            return testPeptide;
        }

        public static int[] GetFragmentHash(List<Product> products, Tuple<int, int>[] keyValuePair, int FragmentBinsPerDalton)
        {
            double[] newFragments = products.Select(p => p.NeutralMass).ToArray();
            var len = products.Count / 2;

            for (int i = 0; i < keyValuePair.Length; i++)
            {
                var j = keyValuePair[i].Item1;
                var modj = keyValuePair[i].Item2;
                while (j <= len + 1)
                {
                    newFragments[j - 2] += (double)ModBox.SelectedModifications[modj].MonoisotopicMass;
                    j++;
                }
                j = keyValuePair[i].Item1;
                while (j >= 3)
                {
                    newFragments[len * 2 - j + 2] += (double)ModBox.SelectedModifications[modj].MonoisotopicMass;
                    j--;
                }
            }

            int[] fragmentHash = new int[products.Count];
            for (int i = 0; i < products.Count; i++)
            {
                fragmentHash[i] = (int)Math.Round(newFragments[i] * FragmentBinsPerDalton);
            }
            return fragmentHash;
        }

        public static int[] GetLocalFragmentHash(List<Product> products, int peptideLength, int[] modPoses, int modInd, ModBox TotalBox, ModBox localBox, int FragmentBinsPerDalton)
        {
            List<double> newFragments = new List<double>();
            var local_c_fragments = products.Where(p => p.ProductType == ProductType.b && p.TerminusFragment.AminoAcidPosition >= modPoses[modInd] - 1 && p.TerminusFragment.AminoAcidPosition < modPoses[modInd + 1] - 1).ToList();

            foreach (var c in local_c_fragments)
            {
                var newMass = c.NeutralMass + localBox.Mass;
                newFragments.Add(newMass);
            }

            var local_z_fragments = products.Where(p => p.ProductType == ProductType.y && p.TerminusFragment.AminoAcidPosition >= modPoses[modInd] && p.TerminusFragment.AminoAcidPosition < modPoses[modInd + 1]).ToList();

            foreach (var z in local_z_fragments)
            {
                var newMass = z.NeutralMass + (TotalBox.Mass - localBox.Mass);
                newFragments.Add(newMass);
            }


            int[] fragmentHash = new int[newFragments.Count];
            for (int i = 0; i < newFragments.Count; i++)
            {
                fragmentHash[i] = (int)Math.Round(newFragments[i] * FragmentBinsPerDalton);
            }
            return fragmentHash;
        }

        public static IEnumerable<ModBox> BuildChildModBoxes(int maxNum, int[] modIds)
        {
            yield return new ModBox(new int[0]);
            HashSet<string> seen = new HashSet<string>();
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

        //key: motif, value: all ids for this motif
        public Dictionary<string, List<int>> MotifNeeded
        {
            get
            {
                Dictionary<string, List<int>> aa = new Dictionary<string, List<int>>();
                foreach (var id in ModIds)
                {
                    var mod = SelectedModifications[id];
                    if (aa.ContainsKey(mod.Target.ToString()))
                    {
                        aa[mod.Target.ToString()].Add(id);
                    }
                    else
                    {
                        aa.Add(mod.Target.ToString(), new List<int> { id });
                    }
                }
                return aa;
            }
        }
    }
}
