using System;
using System.Collections.Generic;

namespace OldInternalLogic
{
    public class Protease
    {

        #region Public Constructors

        public Protease(string name, IEnumerable<string> sequencesInducingCleavage, IEnumerable<string> sequencesPreventingCleavage, OldLogicTerminus cleavageTerminus, CleavageSpecificity cleavageSpecificity, string psiMSAccessionNumber, string psiMSName, string siteRegexp)
        {
            Name = name;
            SequencesInducingCleavage = sequencesInducingCleavage;
            SequencesPreventingCleavage = sequencesPreventingCleavage;
            CleavageTerminus = cleavageTerminus;
            CleavageSpecificity = cleavageSpecificity;
            PsiMSAccessionNumber = psiMSAccessionNumber;
            PsiMSName = psiMSName;
            SiteRegexp = siteRegexp;
        }

        #endregion Public Constructors

        #region Public Properties

        public string Name { get; private set; }

        public OldLogicTerminus CleavageTerminus { get; private set; }
        public IEnumerable<string> SequencesInducingCleavage { get; private set; }
        public IEnumerable<string> SequencesPreventingCleavage { get; private set; }
        public CleavageSpecificity CleavageSpecificity { get; private set; }
        public string PsiMSAccessionNumber { get; private set; }
        public string PsiMSName { get; private set; }
        public string SiteRegexp { get; private set; }

        #endregion Public Properties

        #region Public Methods

        public override string ToString()
        {
            return Name;
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
                    if ((CleavageTerminus != OldLogicTerminus.N && i - c.Length + 1 >= 0 && sequence.Substring(i - c.Length + 1, c.Length).Equals(c, StringComparison.OrdinalIgnoreCase))
                        || (CleavageTerminus == OldLogicTerminus.N && i + 1 + c.Length <= sequence.Length && sequence.Substring(i + 1, c.Length).Equals(c, StringComparison.OrdinalIgnoreCase)))
                    {
                        bool cleave = true;
                        foreach (string nc in SequencesPreventingCleavage)
                        {
                            if ((CleavageTerminus != OldLogicTerminus.N && i + 1 + nc.Length <= sequence.Length && sequence.Substring(i + 1, nc.Length).Equals(nc, StringComparison.OrdinalIgnoreCase))
                                || (CleavageTerminus == OldLogicTerminus.N && i - nc.Length + 1 >= 0 && sequence.Substring(i - nc.Length + 1, nc.Length).Equals(nc, StringComparison.OrdinalIgnoreCase)))
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