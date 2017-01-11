using Spectra;
using System;
using System.Collections.Generic;

namespace InternalLogicEngineLayer
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
            foreach (double huh in acceptableSortedMassShifts)
            {
                yield return new DoubleRange(peptideMonoisotopicMass - huh, tol);
            }
        }
    }
}