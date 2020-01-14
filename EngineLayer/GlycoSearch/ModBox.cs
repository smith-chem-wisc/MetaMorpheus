using System.Collections.Generic;
using System.IO;
using System.Linq;
using Chemistry;
using System;
using Proteomics;
using MassSpectrometry;
using Proteomics.ProteolyticDigestion;
using Proteomics.Fragmentation;

namespace EngineLayer
{
    public class ModBox
    {
        public static Modification[] SelectedModifications;

        public static ModBox[] ModBoxes;

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
            List<int> allPossibleModSites = new List<int>();

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

                if (possibleModSites.Count < mn.Value.Count)
                {
                    return null;
                }

                allPossibleModSites.AddRange(possibleModSites);
            }
            return allPossibleModSites.OrderBy(p=>p).ToArray();
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
