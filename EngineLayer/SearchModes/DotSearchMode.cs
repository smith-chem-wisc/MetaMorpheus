using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Linq;

namespace EngineLayer
{
    public class DotSearchMode : SearchMode
    {

        #region Private Fields

        private readonly double[] acceptableSortedMassShifts;
        private readonly Tolerance tol;

        #endregion Private Fields

        #region Public Constructors

        public DotSearchMode(string FileNameAddition, IEnumerable<double> acceptableMassShifts, Tolerance tol) : base(FileNameAddition)
        {
            this.acceptableSortedMassShifts = acceptableMassShifts.OrderBy(b => b).ToArray();
            this.tol = tol;
            this.NumNotches = acceptableSortedMassShifts.Length;
        }

        #endregion Public Constructors

        #region Public Methods

        public override int Accepts(double scanPrecursorMass, double peptideMass)
        {
            // index of the first element that is larger than or equal to value
            int index = Array.BinarySearch(acceptableSortedMassShifts, scanPrecursorMass - peptideMass);
            if (index < 0)
                index = ~index;

            if (index < acceptableSortedMassShifts.Length && tol.Within(scanPrecursorMass, acceptableSortedMassShifts[index] + peptideMass))
                return index;

            if (index > 0 && tol.Within(scanPrecursorMass, acceptableSortedMassShifts[index - 1] + peptideMass))
                return index - 1;
            return -1;
        }

        public override IEnumerable<Tuple<DoubleRange, int>> GetAllowedPrecursorMassIntervals(double peptideMonoisotopicMass)
        {
            for (int j = 0; j < acceptableSortedMassShifts.Length; j++)
            {
                yield return new Tuple<DoubleRange, int>(new DoubleRange(peptideMonoisotopicMass + acceptableSortedMassShifts[j], tol), j);
            }
        }

        public override string ToString()
        {
            return FileNameAddition + " dot " + tol.ToString() + " " + string.Join(",", acceptableSortedMassShifts);
        }

        #endregion Public Methods

    }
}