using Omics;
using Omics.Digestion;

namespace EngineLayer.ProteinParsimony
{
    internal class ParsimonySequence
    {
        public ParsimonySequence(IBioPolymerWithSetMods pwsm, bool TreatModPeptidesAsDifferentPeptides)
        {
            Sequence = TreatModPeptidesAsDifferentPeptides ? pwsm.FullSequence : pwsm.BaseSequence;
            DigestionAgent = pwsm.DigestionParams.DigestionAgent;
        }

        public string Sequence { get; }
        public DigestionAgent DigestionAgent { get; }

        public override bool Equals(object obj)
        {
            ParsimonySequence other = (ParsimonySequence)obj;
            return other != null
                && (Sequence == null && other.Sequence == null || Sequence.Equals(other.Sequence))
                && (DigestionAgent == null && other.DigestionAgent == null || DigestionAgent.Equals(other.DigestionAgent));
        }

        public override int GetHashCode()
        {
            return Sequence.GetHashCode() ^ DigestionAgent.GetHashCode();
        }
    }
}