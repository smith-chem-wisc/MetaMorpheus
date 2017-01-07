using Spectra;
using System.Collections.Generic;

namespace IndexSearchAndAnalyze
{
    public class IntervalSearchMode : SearchMode
    {
        private List<double> acceptableSortedMassShifts;
        private Tolerance tol;

        public IntervalSearchMode(string FileNameAddition, List<double> acceptableSortedMassShifts, Tolerance tol) : base(FileNameAddition)
        {
            this.acceptableSortedMassShifts = acceptableSortedMassShifts;
            this.tol = tol;
        }

        internal IEnumerable<DoubleRange> GetAllowedPrecursorMassIntervals(double peptideMonoisotopicMass)
        {
            foreach (double huh in acceptableSortedMassShifts)
                yield return new DoubleRange(peptideMonoisotopicMass - huh, tol);
        }
    }
}