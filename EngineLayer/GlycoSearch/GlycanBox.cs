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

        public static IEnumerable<GlycanBox> BuildOGlycanBoxes(int maxNum)
        {

            for (int i = 1; i <= maxNum; i++)
            {
                foreach (var idCombine in Glycan.GetKCombsWithRept(Enumerable.Range(0, GlobalOGlycans.Length), i))
                {
                    GlycanBox glycanBox = new GlycanBox(idCombine.ToArray());

                    yield return glycanBox;
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

        public GlycanBox(int[] ids):base(ids)
        {

        }

        public string GlycanIdString
        {
            get
            {
                return string.Join(",", ModIds.Select(p => p.ToString()));
            }
        }

        public byte[] Kind
        {
            get
            {
                {
                    byte[] kind = new byte[8] { 0, 0, 0, 0, 0, 0, 0, 0 };
                    foreach (var id in ModIds)
                    {
                        for (int i = 0; i < 8; i++)
                        {
                            kind[i] += GlobalOGlycans[id].Kind[i];
                        }
                    }
                    return kind;
                }
            }
        }

        public string Structure
        {
            get
            {
                return Glycan.GetKindString(Kind);
            }

        }

        public override double Mass
        {
            get
            {
                return (double)Glycan.GetMass(Kind)/1E5;
            }
        }

    }
}
