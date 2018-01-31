using System;
using System.Collections.Generic;

namespace EngineLayer.Neo
{
    public class FusionCandidate
    {
        #region Public Constructors

        public FusionCandidate(String seq)
        {
            this.seq = seq;
            this.junctionIndexes = new List<int>();
            //this.foundIons = new bool[seq.Length];
            fusionType = FusionType.TS; //default
            this.parentInfo = new List<ParentInfo>();
            this.translatedParents = new List<TranslatedParent>();
            this.cisParents = new List<CisParent>();
            this.transParents = new List<TransParent>();
        }

        #endregion Public Constructors

        #region Public Enums

        //ordered by priority //translated, normalCis, reverseCis, trans
        public enum FusionType
        {
            TL, NC, RC, TS
        }

        #endregion Public Enums

        #region Public Properties

        public string seq { get; set; }
        public List<int> junctionIndexes { get; private set; }
        public bool[] foundIons { get; private set; }
        public List<ParentInfo> parentInfo { get; set; } //Will have a ParentInfo object for every possible fragment* (*every fragment with length of 6 or more OR just the first length that was found (<6))
        public List<TranslatedParent> translatedParents { get; set; }
        public List<CisParent> cisParents { get; set; }
        public List<TransParent> transParents { get; set; }
        public FusionType fusionType { get; set; }

        #endregion Public Properties

        #region Public Methods

        //private List<FusionCandidate> fragSources;
        public void addJunctionIndex(int index)
        {
            this.junctionIndexes.Add(index);
        }

        public void setFoundIons(bool[] foundIons)
        {
            this.foundIons = foundIons;
        }

        public void makeFoundIons()
        {
            this.foundIons = new bool[this.seq.Length]; //|A|B|C|D|E|F|K where the whole peptide peak is always placed arbitrarily at the n term
            for (int i = 0; i < this.foundIons.Length; i++)
            {
                this.foundIons[i] = false;
            }
        }

        public void deepCopyFoundIons(FusionCandidate original)
        {
            this.foundIons = new bool[original.foundIons.Length];
            for (int index = 0; index < this.foundIons.Length; index++)
            {
                this.foundIons[index] = original.foundIons[index];
            }
        }

        #endregion Public Methods
    }
}