using System;
using System.Collections.Generic;

namespace EngineLayer
{
    public class Protease
    {
        #region Private Fields

        private static readonly Dictionary<CleavageSpecificity, TerminusType> semiSpecificProteaseAlterationDicitionary = new Dictionary<CleavageSpecificity, TerminusType>
        {
            {CleavageSpecificity.Full, TerminusType.None },
            {CleavageSpecificity.FullMaxN, TerminusType.N },
            {CleavageSpecificity.FullMaxC, TerminusType.C }
        };

        #endregion Private Fields

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

        public Protease(Protease protease, TerminusType terminusType)
        {
            Name = protease.Name;
            SequencesInducingCleavage = protease.SequencesInducingCleavage;
            SequencesPreventingCleavage = protease.SequencesPreventingCleavage;
            CleavageTerminus = protease.CleavageTerminus;
            PsiMsAccessionNumber = protease.PsiMsAccessionNumber;
            PsiMsName = protease.PsiMsName;
            SiteRegexp = protease.SiteRegexp;

            if (terminusType == TerminusType.N)
                CleavageSpecificity = CleavageSpecificity.FullMaxN;
            else if (terminusType == TerminusType.C)
                CleavageSpecificity = CleavageSpecificity.FullMaxC;
            else if (terminusType == TerminusType.None)
                CleavageSpecificity = CleavageSpecificity.Full;
            else
                throw new MetaMorpheusException("Terminus obtained for NonSpecific search has not been implemented.");
        }

        #endregion Public Constructors

        #region Public Properties

        public string Name { get; private set; }

        public TerminusType CleavageTerminus { get; private set; }
        public IEnumerable<string> SequencesInducingCleavage { get; private set; }
        public IEnumerable<string> SequencesPreventingCleavage { get; private set; }
        public CleavageSpecificity CleavageSpecificity { get; private set; }
        public string PsiMsAccessionNumber { get; private set; }
        public string PsiMsName { get; private set; }
        public string SiteRegexp { get; private set; }

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

        public bool ProteaseMustBeUpdated(TerminusType terminusType)
        {
            if (semiSpecificProteaseAlterationDicitionary.TryGetValue(CleavageSpecificity, out TerminusType value))
                return value == terminusType ? false : true;
            else
                return false;
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