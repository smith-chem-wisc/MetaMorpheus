using Spectra;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;

namespace EngineLayer
{
    public class DotSearchMode : SearchMode
    {

        #region Private Fields

        private readonly List<double> acceptableSortedMassShifts;
        private readonly Tolerance tol;

        #endregion Private Fields

        #region Public Constructors

        public DotSearchMode(IEnumerable<double> acceptableMassShifts, Tolerance tol) : this(tol.Value.ToString("F3", CultureInfo.InvariantCulture) + "around" + string.Join("-", acceptableMassShifts.Select(b => b.ToString("F3", CultureInfo.InvariantCulture))), acceptableMassShifts, tol)
        {
        }

        public DotSearchMode(string FileNameAddition, IEnumerable<double> acceptableMassShifts, Tolerance tol) : base(FileNameAddition)
        {
            HashSet<double> hs = new HashSet<double>(acceptableMassShifts.Select(b => Math.Round(b, 5)));
            this.acceptableSortedMassShifts = hs.OrderBy(b => b).ToList();
            this.tol = tol;
            this.NumNotches = acceptableSortedMassShifts.Count;
        }

        #endregion Public Constructors

        #region Public Methods

        public override int Accepts(double scanPrecursorMass, double peptideMass)
        {
            for (int j = 0; j < acceptableSortedMassShifts.Count; j++)
                if (tol.Within(scanPrecursorMass, acceptableSortedMassShifts[j] + peptideMass))
                    return j;
            return -1;
        }

        public override IEnumerable<Tuple<DoubleRange, int>> GetAllowedPrecursorMassIntervals(double peptideMonoisotopicMass)
        {
            for (int j = 0; j < acceptableSortedMassShifts.Count; j++)
            {
                yield return new Tuple<DoubleRange, int>(new DoubleRange(peptideMonoisotopicMass + acceptableSortedMassShifts[j], tol), j);
            }
        }

        #endregion Public Methods

    }
}