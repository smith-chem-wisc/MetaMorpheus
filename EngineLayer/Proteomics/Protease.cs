using System;
using System.Collections.Generic;

namespace EngineLayer
{
    public class Protease
    {
        #region Public Constructors

        public Protease(string name, IEnumerable<string> sequencesInducingCleavage, IEnumerable<string> sequencesPreventingCleavage, TerminusType cleavageTerminus, CleavageSpecificity cleavageSpecificity, string psiMSAccessionNumber, string psiMSName, string siteRegexp)
        {
            Name = name;
            SequencesInducingCleavage = sequencesInducingCleavage;
            SequencesPreventingCleavage = sequencesPreventingCleavage;
            CleavageTerminus = cleavageTerminus;
            CleavageSpecificity = cleavageSpecificity;
            PsiMsAccessionNumber = psiMSAccessionNumber;
            PsiMsName = psiMSName;
            SiteRegexp = siteRegexp;
        }

        #endregion Public Constructors

        #region Public Properties

        public string Name { get; }

        public TerminusType CleavageTerminus { get; }
        public IEnumerable<string> SequencesInducingCleavage { get; }
        public IEnumerable<string> SequencesPreventingCleavage { get; }
        public CleavageSpecificity CleavageSpecificity { get; }
        public string PsiMsAccessionNumber { get; }
        public string PsiMsName { get; }
        public string SiteRegexp { get; }

        #endregion Public Properties

        #region Public Methods

        public override string ToString()
        {
            return Name;
        }

        public override bool Equals(object obj)
        {
            var a = obj as Protease;
            return a != null
                && a.Name.Equals(Name);
        }

        public override int GetHashCode()
        {
            return Name.GetHashCode();
        }

        #endregion Public Methods

        #region Internal Methods

        internal List<int> GetDigestionSiteIndices(string sequence)
        {
            var indices = new List<int>();

            for (int i = 0; i < sequence.Length - 1; i++)
            {
                foreach (string c in SequencesInducingCleavage)
                {
                    if ((CleavageTerminus != TerminusType.N && i - c.Length + 1 >= 0 && sequence.Substring(i - c.Length + 1, c.Length).Equals(c, StringComparison.OrdinalIgnoreCase))
                        || (CleavageTerminus == TerminusType.N && i + 1 + c.Length <= sequence.Length && sequence.Substring(i + 1, c.Length).Equals(c, StringComparison.OrdinalIgnoreCase)))
                    {
                        bool cleave = true;
                        foreach (string nc in SequencesPreventingCleavage)
                        {
                            if ((CleavageTerminus != TerminusType.N && i + 1 + nc.Length <= sequence.Length && sequence.Substring(i + 1, nc.Length).Equals(nc, StringComparison.OrdinalIgnoreCase))
                                || (CleavageTerminus == TerminusType.N && i - nc.Length + 1 >= 0 && sequence.Substring(i - nc.Length + 1, nc.Length).Equals(nc, StringComparison.OrdinalIgnoreCase)))
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

        #endregion Internal Methods
    }
}