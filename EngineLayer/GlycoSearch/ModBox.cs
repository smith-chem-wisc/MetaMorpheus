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

        public virtual double Mass { get; set; }

    }
}
