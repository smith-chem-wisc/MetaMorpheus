using OldInternalLogic;
using Spectra;
using System;
using System.Collections.Generic;
using System.Linq;

namespace InternalLogicEngineLayer
{
    public class DotSearchMode : SearchMode
    {
        private List<double> acceptableSortedMassShifts;
        Tolerance tol;

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

        public override bool Accepts(double scanPrecursorMass, double peptideMass)
        {
            throw new NotImplementedException();
        }

        internal override IEnumerable<DoubleRange> GetAllowedPrecursorMassIntervals(double peptideMonoisotopicMass)
        {
            foreach (double huh in acceptableSortedMassShifts)
            {
                yield return new DoubleRange(peptideMonoisotopicMass + huh, tol);
            }
        }

        private class ffff : IEqualityComparer<double>
        {
            int i;
            public ffff(int i)
            {
                this.i = i;
            }
            public bool Equals(double x, double y)
            {
                return Math.Round(x, i) == Math.Round(y, i);
            }

            public int GetHashCode(double obj)
            {
                return Math.Round(obj, i).GetHashCode();
            }
        }
    }
}