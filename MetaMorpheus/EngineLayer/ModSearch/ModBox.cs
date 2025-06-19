using Easy.Common.Extensions;
using Omics.Modifications;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Net.Http.Headers;
using System.Runtime.InteropServices;
using System.Text;
using System.Threading.Tasks;

namespace EngineLayer.ModSearch
{
    /// <summary>
    /// ModBox is a container for modifications, including glycans and regular modifications.
    /// </summary>
    public class ModBox
    {
        public List<Modification> Mods { get; } // Total modifications in this ModBox, including glycans and regular modifications
        public double? Mass { get; set; } // The total mass of the modifications in this ModBox, including glycans and regular modifications
        public bool TargetDecoy { get; set; } // Whether this ModBox is a target or decoy box. Will be use for FDR calculation.
        public int NumberOfMods { get;   } // The number of mods
        public int MaxModNum { get; set; } // Max number of modifications allowed on one peptide. This is used to limit the number of modifications in the childBox building.
        public int[] ModIds { get; } //The id is the index of the modification in the GlobalMods. 
        public bool IsChild { get; set; } 
        public List<ModBox> ChildBoxes { get; set; } // ParentBox means the original ModBox from database. The ChildBox is the subset combination of the ParentBox.
        public int[] SugarKind { get; set; } // The sugar kind of the glycan modifications in this ModBox, used for O-glycan search. The length is 11, each index represents a specific sugar kind.


        public ModBox(int[] ids) // The constructor just for glycan constructor in O-Pair
        {
            ModIds = ids;
            NumberOfMods = ids.Length;
            TargetDecoy = true;
        }

        public ModBox(List<Modification> mods, int[] ids = null, int maxModNum = 3, bool isChild = false) 
        {
            IsChild = isChild;
            Mods = mods;
            Mass = mods.Sum(b => b.MonoisotopicMass);
            MaxModNum = maxModNum;
            if (ids != null) 
            {
                ModIds = ids;
            }

            NumberOfMods = mods.Count;
            if (!isChild) //Any parent box should generate its child boxes
            {
                BuildChildBox();
            }
            BuildSugarKind(mods);
        }

        /// <summary>
        /// Builds child boxes for this ModBox based on the modifications it contains.
        /// </summary>
        public void BuildChildBox() 
        {
            ChildBoxes = new List<ModBox> ();

            // Add an  ModBox as the first child box
            var nullBox = new ModBox(new List<Modification>(), new int[0], MaxModNum, true)
            {
                TargetDecoy = true
            };
            ChildBoxes.Add(nullBox);

            var modCombinations = Enumerable.Range(0, ModIds.Length)
                                    .SelectMany(i => GetCombination(ModIds, i + 1));

            foreach (var combination in modCombinations)
            {
                var childMods = combination.Select(id => Mods[Array.IndexOf(ModIds, id)]).ToList();
                var childBox = new ModBox(childMods, combination.ToArray(), MaxModNum, true)
                {
                    TargetDecoy = true
                };
                ChildBoxes.Add(childBox);
            }

            ChildBoxes = ChildBoxes
            .GroupBy(p => string.Join(",", p.ModIds.OrderBy(id => id)))
            .Select(g => g.First())
            .ToList();
        }

        /// <summary>
        /// Build the sugar kind by their mods information
        /// </summary>
        /// <param name="modForModBox"></param>
        public void BuildSugarKind(List<Modification> modInBox)
        {
            SugarKind = new int[11];
            foreach (var mod in modInBox)
            {
                if (mod.ModificationType == "O-Glycosylation") // For now we only applied for O-glycan, will be served for N-glycan in the future
                {
                    byte[] kinds = Glycan.StringToKind(mod.OriginalId);
                    for (int i = 0; i < kinds.Length; i++)
                    {
                        SugarKind[i] += kinds[i];
                    }
                }
                
            }
        }

        private static IEnumerable<IEnumerable<int>> GetCombination(IEnumerable<int> elements, int k)
        {
            return k == 0 ? new[] { new int[0] } :
              elements.SelectMany((e, i) =>
                GetCombination(elements.Skip(i + 1), k - 1).Select(c => (new[] { e }).Concat(c)));
        }

        public static IEnumerable<IEnumerable<T>> GetCombinationWithRept<T>(IEnumerable<T> list, int length) where T : IComparable
        {
            if (length == 1) return list.Select(t => new T[] { t });
            return GetCombinationWithRept(list, length - 1).SelectMany(t => list.Where(o => o.CompareTo(t.Last()) >= 0), (t1, t2) => t1.Concat(new T[] { t2 }));
        }


    }
}
