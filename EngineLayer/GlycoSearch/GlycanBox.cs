using System.Collections.Generic;
using System.IO;
using System.Linq;
using Chemistry;
using System;
using Proteomics;
using MassSpectrometry;

namespace EngineLayer
{
    public class GlycanBox:ModBox
    {
        public static Glycan[] GlobalOGlycans;

        public static Modification[] GlobalOGlycanModifications;

        public static GlycanBox[] OGlycanBoxes;

        public static IEnumerable<GlycanBox> BuildOGlycanBoxes(int maxNum, bool buildDecoy = false)
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

        public static Modification[] BuildGlobalOGlycanModifications(Glycan[] globalOGlycans)
        {
            Modification[] globalOGlycanModifications = new Modification[globalOGlycans.Length];

            for (int i = 0; i < GlobalOGlycans.Length; i++)
            {
                globalOGlycanModifications[i] = Glycan.OGlycanToModification(globalOGlycans[i]);
            }
            return globalOGlycanModifications;
        }

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
            byte[] kind = new byte[] { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
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
                var childShift = random.Next(-3000000, 3000000); //Based on pGlyco [1, 30] and GlycoPAT [-50, 50].
                Mass = (double)(Glycan.GetMass(Kind) + childShift) / 1E5;
            }

            Random random_decoyMass = new Random();
            var decoyMassShift = random_decoyMass.Next(-3000000, 3000000); //Based on pGlyco [1, 30] and GlycoPAT [-50, 50].
            DecoyMass = (double)(Glycan.GetMass(Kind) + decoyMassShift) / 1E5;
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
