using OldInternalLogic;
using Spectra;
using System;
using System.Collections.Generic;
using System.Linq;

namespace InternalLogicEngineLayer
{
    public class DotSearchMode : SearchMode
    {

        #region Private Fields

        private List<double> acceptableSortedMassShifts;
        private Tolerance tol;

        #endregion Private Fields

        #region Public Constructors

        public DotSearchMode(string FileNameAddition, double[] acceptableSortedMassShifts, Tolerance tol) : base(FileNameAddition)
        {
            this.acceptableSortedMassShifts = acceptableSortedMassShifts.ToList();
            this.tol = tol;
        }

        public DotSearchMode(string v, List<MorpheusModification> gptmdModifications, IEnumerable<Tuple<double, double>> combos, Tolerance tolerance) : base(v)
        {
            List<double> ok = new List<double>() { 0 };
            ok.AddRange(gptmdModifications.Select(b => b.MonoisotopicMassShift));
            ok.AddRange(combos.Select(b => b.Item1 + b.Item2));
            IEqualityComparer<double> f = new ffff(3);
            acceptableSortedMassShifts = ok.Distinct(f).OrderBy(b => b).ToList();
            tol = tolerance;
        }

        #endregion Public Constructors

        #region Public Methods

        public override bool Accepts(double scanPrecursorMass, double peptideMass)
        {
            foreach (double huh in acceptableSortedMassShifts)
            {
                if (tol.Within(scanPrecursorMass - peptideMass, huh))
                    return true;
            }
            return false;
        }

        #endregion Public Methods

        #region Internal Methods

        internal override IEnumerable<DoubleRange> GetAllowedPrecursorMassIntervals(double peptideMonoisotopicMass)
        {
            foreach (double huh in acceptableSortedMassShifts)
            {
                yield return new DoubleRange(peptideMonoisotopicMass + huh, tol);
            }
        }

        internal override string SearchModeString()
        {
            return "Tolerance of " + tol.ToString() + " around mass diffs: " + string.Join(",", acceptableSortedMassShifts);
        }

        #endregion Internal Methods

        #region Private Classes

        private class ffff : IEqualityComparer<double>
        {

            #region Private Fields

            private int i;

            #endregion Private Fields

            #region Public Constructors

            public ffff(int i)
            {
                this.i = i;
            }

            #endregion Public Constructors

            #region Public Methods

            public bool Equals(double x, double y)
            {
                return Math.Round(x, i) == Math.Round(y, i);
            }

            public int GetHashCode(double obj)
            {
                return Math.Round(obj, i).GetHashCode();
            }

            #endregion Public Methods

        }

        #endregion Private Classes

    }
}