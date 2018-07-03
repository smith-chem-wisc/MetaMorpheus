using System;
using System.Collections.Generic;

namespace EngineLayer.Neo
{
    public class FusionCandidate
    {
        public FusionCandidate(String seq)
        {
            Seq = seq;
            JunctionIndexes = new List<int>();
            //FoundIons = new bool[seq.Length];
            FusionType = FusionType.TS; //default
            ParentInfo = new List<ParentInfo>();
            TranslatedParents = new List<TranslatedParent>();
            CisParents = new List<CisParent>();
            TransParents = new List<TransParent>();
        }

        public string Seq { get; set; }
        public List<int> JunctionIndexes { get; private set; }
        public bool[] FoundIons { get; private set; }
        public List<ParentInfo> ParentInfo { get; set; } //Will have a ParentInfo object for every possible fragment* (*every fragment with length of 6 or more OR just the first length that was found (<6))
        public List<TranslatedParent> TranslatedParents { get; set; }
        public List<CisParent> CisParents { get; set; }
        public List<TransParent> TransParents { get; set; }
        public FusionType FusionType { get; set; }

        //private List<FusionCandidate> fragSources;
        public void addJunctionIndex(int index)
        {
            JunctionIndexes.Add(index);
        }

        public void setFoundIons(bool[] foundIons)
        {
            FoundIons = foundIons;
        }

        public void makeFoundIons()
        {
            FoundIons = new bool[Seq.Length]; //|A|B|C|D|E|F|K where the whole peptide peak is always placed arbitrarily at the n term
            for (int i = 0; i < FoundIons.Length; i++)
            {
                FoundIons[i] = false;
            }
        }

        public void deepCopyFoundIons(FusionCandidate original)
        {
            FoundIons = new bool[original.FoundIons.Length];
            for (int index = 0; index < FoundIons.Length; index++)
            {
                FoundIons[index] = original.FoundIons[index];
            }
        }
    }
}