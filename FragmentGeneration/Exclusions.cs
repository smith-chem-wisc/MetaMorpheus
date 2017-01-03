using Proteomics;
using System;
using System.Collections.Generic;
using System.Linq;

namespace FragmentGeneration
{
    public class Exclusions
    {
        public static double[] PopulateExcludeList()
        {
            // Do not exclude Lysine + Anything
            // Do not exclude Lysine, Arginine, Glycine, Asparagine, Alanine, Methionine

            // Exclude -76.134779 and -48.128629 - these are real reversed phosphorylations

            HashSet<double> exclude = new HashSet<double>();

            exclude.Add(-76.134779);
            exclude.Add(-48.128629);

            HashSet<AminoAcid> doNotExcludeEvenCombos = new HashSet<AminoAcid>() { AminoAcid.GetResidue('K') };

            HashSet<AminoAcid> doNotExclude = new HashSet<AminoAcid>() {
                AminoAcid.GetResidue('K'),
                AminoAcid.GetResidue('R'),
                AminoAcid.GetResidue('G'),
                AminoAcid.GetResidue('N'),
                AminoAcid.GetResidue('A'),
                AminoAcid.GetResidue('M'),
            };

            for (char c = 'A'; c <= 'Z'; c++)
            {
                AminoAcid residue;
                if (AminoAcid.TryGetResidue(c, out residue))
                {
                    if (!doNotExclude.Contains(residue))
                        exclude.Add(residue.MonoisotopicMass);
                    for (char cc = 'A'; cc <= 'Z'; cc++)
                    {
                        AminoAcid residueCC;
                        if (AminoAcid.TryGetResidue(cc, out residueCC))
                        {
                            if (!doNotExcludeEvenCombos.Contains(residueCC) && !doNotExcludeEvenCombos.Contains(residue))
                                exclude.Add(residue.MonoisotopicMass + residueCC.MonoisotopicMass);
                        }
                    }
                }
            }
            return exclude.ToArray();
        }

        public static bool DoNotExclude(double a, double tolExclude, double[] exclude)
        {
            foreach (var heh in exclude)
                if (Math.Abs(heh - a) < tolExclude)
                    return false;
            return true;
        }
    }
}