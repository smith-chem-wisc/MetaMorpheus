using Proteomics.ProteolyticDigestion;

namespace EngineLayer.ProteinParsimony
{
    internal class ParsimonySequence
    {
        public ParsimonySequence(PeptideWithSetModifications pwsm, bool TreatModPeptidesAsDifferentPeptides)
        {
            Sequence = TreatModPeptidesAsDifferentPeptides ? pwsm.FullSequence : pwsm.BaseSequence;
            Protease = pwsm.DigestionParams.Protease;
        }

        public string Sequence { get; }
        public Protease Protease { get; }

        public override bool Equals(object obj)
        {
            ParsimonySequence other = (ParsimonySequence)obj;
            return other != null
                && (Sequence == null && other.Sequence == null || Sequence.Equals(other.Sequence))
                && (Protease == null && other.Protease == null || Protease.Equals(other.Protease));
        }

        public override int GetHashCode()
        {
            return Sequence.GetHashCode() ^ Protease.GetHashCode();
        }
    }
}