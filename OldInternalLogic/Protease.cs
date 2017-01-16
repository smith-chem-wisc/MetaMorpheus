using System;
using System.Collections.Generic;

namespace OldInternalLogic
{
    public class Protease
    {
        public string Name { get; private set; }

        public Terminus CleavageTerminus { get; private set; }
        public IEnumerable<string> SequencesInducingCleavage { get; private set; }
        public IEnumerable<string> SequencesPreventingCleavage { get; private set; }
        public CleavageSpecificity CleavageSpecificity { get; private set; }
        public string PsiMsAccessionNumber { get; private set; }
        public string PsiMsName { get; private set; }
        public string SiteRegexp { get; private set; }

        public Protease(string name, IEnumerable<string> sequencesInducingCleavage, IEnumerable<string> sequencesPreventingCleavage, Terminus cleavageTerminus, CleavageSpecificity cleavageSpecificity, string psiMsAccessionNumber, string psiMsName, string siteRegexp)
        {
            Name = name;
            SequencesInducingCleavage = sequencesInducingCleavage;
            SequencesPreventingCleavage = sequencesPreventingCleavage;
            CleavageTerminus = cleavageTerminus;
            CleavageSpecificity = cleavageSpecificity;
            PsiMsAccessionNumber = psiMsAccessionNumber;
            PsiMsName = psiMsName;
            SiteRegexp = siteRegexp;
        }

        public override string ToString()
        {
            return Name;
        }

        public List<int> GetDigestionSiteIndices(string sequence)
        {
            var indices = new List<int>();

            for (int i = 0; i < sequence.Length - 1; i++)
            {
                foreach (string c in SequencesInducingCleavage)
                {
                    if ((CleavageTerminus != Terminus.N && i - c.Length + 1 >= 0 && sequence.Substring(i - c.Length + 1, c.Length).Equals(c, StringComparison.InvariantCultureIgnoreCase))
                        || (CleavageTerminus == Terminus.N && i + 1 + c.Length <= sequence.Length && sequence.Substring(i + 1, c.Length).Equals(c, StringComparison.InvariantCultureIgnoreCase)))
                    {
                        bool cleave = true;
                        foreach (string nc in SequencesPreventingCleavage)
                        {
                            if ((CleavageTerminus != Terminus.N && i + 1 + nc.Length <= sequence.Length && sequence.Substring(i + 1, nc.Length).Equals(nc, StringComparison.InvariantCultureIgnoreCase))
                                || (CleavageTerminus == Terminus.N && i - nc.Length + 1 >= 0 && sequence.Substring(i - nc.Length + 1, nc.Length).Equals(nc, StringComparison.InvariantCultureIgnoreCase)))
                            {
                                cleave = false;
                                break;
                            }
                        }
                        if (cleave)
                        {
                            indices.Add(i + 1);
                        }
                    }
                }
            }

            return indices;
        }
    }
}