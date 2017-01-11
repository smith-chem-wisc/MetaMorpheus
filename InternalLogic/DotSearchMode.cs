using Spectra;
using System;
using System.Collections.Generic;

namespace InternalLogic
{
    public class DotSearchMode : SearchMode
    {
        private double[] acceptableSortedMassShifts;
        private Tolerance tol;

        public DotSearchMode(string FileNameAddition, double[] acceptableSortedMassShifts, Tolerance tol) : base(FileNameAddition)
        {
            this.acceptableSortedMassShifts = acceptableSortedMassShifts;
            this.tol = tol;
        }

        public override bool Accepts(double scanPrecursorMass, double peptideMass)
        {
            throw new NotImplementedException();
        }

        internal override IEnumerable<DoubleRange> GetAllowedPrecursorMassIntervals(double peptideMonoisotopicMass)
        {
            var thisDoubleRange = new DoubleRange(peptideMonoisotopicMass - acceptableSortedMassShifts[0], tol);
            double prevRangeMaximum = thisDoubleRange.Minimum;
            double prevRangeMinimum = thisDoubleRange.Maximum;
            foreach (double huh in acceptableSortedMassShifts)
            {
                thisDoubleRange = new DoubleRange(peptideMonoisotopicMass - huh, tol);
                if (thisDoubleRange.Minimum > prevRangeMaximum)
                {
                    yield return new DoubleRange(prevRangeMinimum, prevRangeMaximum);
                    prevRangeMinimum = thisDoubleRange.Minimum;
                }
                prevRangeMaximum = thisDoubleRange.Maximum;
            }
            yield return new DoubleRange(prevRangeMinimum, prevRangeMaximum);
        }
    }
}