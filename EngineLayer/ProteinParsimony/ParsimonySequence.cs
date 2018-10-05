using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Text;

namespace EngineLayer.ProteinParsimony
{
    class ParsimonySequence
    {
        public readonly string Sequence;
        public  readonly Protease Protease;
        public ParsimonySequence(PeptideWithSetModifications pwsm, bool TreatModPeptidesAsDifferentPeptides)
        {
            if (TreatModPeptidesAsDifferentPeptides)
            {
                Sequence = pwsm.FullSequence;
            }
            else
            {
                Sequence = pwsm.BaseSequence;
            }            
            Protease = pwsm.DigestionParams.Protease;
        }
        public override bool Equals(object obj)
        {
            ParsimonySequence other = (ParsimonySequence)obj;
            if (this.Sequence == other.Sequence && this.Protease == other.Protease)
            {
                return true;
            }
            else
            {
                return false;
            }
        }
        public override int GetHashCode()
        {
            return (Sequence.GetHashCode() + Protease.GetHashCode());
        }
    }

}
