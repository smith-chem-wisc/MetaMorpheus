using Proteomics;
using System;
using System.Collections.Generic;
using System.Linq;

namespace InternalLogicEngineLayer
{
    public static class Exclusions
    {
        #region Public Methods

        public static double[] PopulateExcludeList()
        {
            // Do not exclude Lysine + Anything
            // Do not exclude Lysine, Arginine, Glycine, Asparagine, Alanine, Methionine

            // Exclude -76.134779 and -48.128629 - these are real reversed phosphorylations

            var exclude = new HashSet<double>();

            exclude.Add(-76.134779);
            exclude.Add(-48.128629);

            var doNotExcludeEvenCombos = new HashSet<Residue> { Residue.GetResidue('K') };

            var doNotExclude = new HashSet<Residue> {
                Residue.GetResidue('K'),
                Residue.GetResidue('R'),
                Residue.GetResidue('G'),
                Residue.GetResidue('N'),
                Residue.GetResidue('A'),
                Residue.GetResidue('M')
            };

            for (char c = 'A'; c <= 'Z'; c++)
            {
                Residue residue;
                if (Residue.TryGetResidue(c, out residue))
                {
                    if (!doNotExclude.Contains(residue))
                        exclude.Add(residue.MonoisotopicMass);
                    for (char cc = 'A'; cc <= 'Z'; cc++)
                    {
                        Residue residueCC;
                        if (Residue.TryGetResidue(cc, out residueCC))
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

        #endregion Public Methods
    }
}