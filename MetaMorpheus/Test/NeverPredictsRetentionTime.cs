using System;
using System.Collections.Generic;
using System.Linq;
using Chromatography;
using Chromatography.RetentionTimePrediction;
using Omics;

namespace Test
{
    /// <summary>
    /// A retention time predictor that can never produce a value, reporting
    /// <see cref="RetentionTimeFailureReason.InvalidAminoAcid"/> for every peptide.
    ///
    /// Stands in for the real cases: Chronologer declines a sequence longer than 50 residues, shorter than 7,
    /// or carrying a non-canonical amino acid such as selenocysteine, and any predictor can fail outright.
    /// What matters is that such a peptidoform is reported as UNPREDICTABLE rather than being silently scored
    /// as if its hydrophobicity were zero.
    /// </summary>
    internal sealed class NeverPredictsRetentionTime : IRetentionTimePredictor
    {
        public string PredictorName => "NeverPredicts";
        public SeparationType SeparationType => SeparationType.HPLC;

        public double? PredictRetentionTimeEquivalent(IRetentionPredictable peptide, out RetentionTimeFailureReason? failureReason)
        {
            failureReason = RetentionTimeFailureReason.InvalidAminoAcid;
            return null;
        }

        public IReadOnlyList<(double? PredictedValue, IRetentionPredictable Peptide, RetentionTimeFailureReason? FailureReason)>
            PredictRetentionTimeEquivalents(IEnumerable<IRetentionPredictable> peptides, int maxThreads = 1)
            => peptides
                .Select(p => ((double?)null, p, (RetentionTimeFailureReason?)RetentionTimeFailureReason.InvalidAminoAcid))
                .ToList();

        public string GetFormattedSequence(IRetentionPredictable peptide, out RetentionTimeFailureReason? failureReason)
        {
            failureReason = RetentionTimeFailureReason.InvalidAminoAcid;
            return null;
        }

        public void Dispose() { }
    }
}
