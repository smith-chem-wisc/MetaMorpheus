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
    /// <see cref="RetentionTimeFailureReason.IncompatibleModifications"/> for every peptide.
    ///
    /// Stands in for the real case this models: Chronologer supports roughly twenty modifications, while a
    /// G-PTM-D search can enable many times that, so a peptidoform carrying an unsupported modification cannot
    /// be represented. What matters is that such a peptidoform is reported as UNPREDICTABLE rather than being
    /// silently scored as if its hydrophobicity were zero.
    /// </summary>
    internal sealed class NeverPredictsRetentionTime : IRetentionTimePredictor
    {
        public string PredictorName => "NeverPredicts";
        public SeparationType SeparationType => SeparationType.HPLC;

        public double? PredictRetentionTimeEquivalent(IRetentionPredictable peptide, out RetentionTimeFailureReason? failureReason)
        {
            failureReason = RetentionTimeFailureReason.IncompatibleModifications;
            return null;
        }

        public IReadOnlyList<(double? PredictedValue, IRetentionPredictable Peptide, RetentionTimeFailureReason? FailureReason)>
            PredictRetentionTimeEquivalents(IEnumerable<IRetentionPredictable> peptides, int maxThreads = 1)
            => peptides
                .Select(p => ((double?)null, p, (RetentionTimeFailureReason?)RetentionTimeFailureReason.IncompatibleModifications))
                .ToList();

        public string GetFormattedSequence(IRetentionPredictable peptide, out RetentionTimeFailureReason? failureReason)
        {
            failureReason = RetentionTimeFailureReason.IncompatibleModifications;
            return null;
        }

        public void Dispose() { }
    }
}
