using System.Collections.Generic;
using System.Linq;
using Chromatography;
using Chromatography.RetentionTimePrediction;
using Chromatography.RetentionTimePrediction.SSRCalc;
using Omics;

namespace Test
{
    /// <summary>
    /// SSRCalc3 for every peptidoform except the listed full sequences, for which it reports
    /// <see cref="RetentionTimeFailureReason.InvalidAminoAcid"/> and returns no value.
    ///
    /// Where <see cref="NeverPredictsRetentionTime"/> fails for everything, this fails for a chosen few, so a test
    /// can check that a failed prediction is left out of a distribution the successful ones still build.
    /// </summary>
    internal sealed class DeclinesSelectedRetentionTimes : IRetentionTimePredictor
    {
        private readonly SSRCalc3RetentionTimePredictor _inner = new();
        private readonly HashSet<string> _declinedFullSequences;

        public DeclinesSelectedRetentionTimes(IEnumerable<string> declinedFullSequences)
        {
            _declinedFullSequences = declinedFullSequences.ToHashSet();
        }

        public string PredictorName => "DeclinesSelected";
        public SeparationType SeparationType => _inner.SeparationType;

        private bool Declines(IRetentionPredictable peptide)
            => peptide is IBioPolymerWithSetMods bpwsm && _declinedFullSequences.Contains(bpwsm.FullSequence);

        public double? PredictRetentionTimeEquivalent(IRetentionPredictable peptide, out RetentionTimeFailureReason? failureReason)
        {
            if (Declines(peptide))
            {
                failureReason = RetentionTimeFailureReason.InvalidAminoAcid;
                return null;
            }
            return _inner.PredictRetentionTimeEquivalent(peptide, out failureReason);
        }

        public IReadOnlyList<(double? PredictedValue, IRetentionPredictable Peptide, RetentionTimeFailureReason? FailureReason)>
            PredictRetentionTimeEquivalents(IEnumerable<IRetentionPredictable> peptides, int maxThreads = 1)
            => peptides
                .Select(p => (PredictRetentionTimeEquivalent(p, out var reason), p, reason))
                .ToList();

        public string GetFormattedSequence(IRetentionPredictable peptide, out RetentionTimeFailureReason? failureReason)
        {
            if (Declines(peptide))
            {
                failureReason = RetentionTimeFailureReason.InvalidAminoAcid;
                return null;
            }
            return _inner.GetFormattedSequence(peptide, out failureReason);
        }

        public void Dispose() => _inner.Dispose();
    }
}
