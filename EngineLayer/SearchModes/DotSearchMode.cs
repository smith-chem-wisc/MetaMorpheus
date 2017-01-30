using Spectra;
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

        public DotSearchMode(IEnumerable<double> acceptableSortedMassShifts, Tolerance tol) : base(tol.Value.ToString("F3", CultureInfo.InvariantCulture) + "around" + string.Join("-", acceptableSortedMassShifts.Select(b => b.ToString("F3", CultureInfo.InvariantCulture))))
        {
            this.acceptableSortedMassShifts = acceptableSortedMassShifts.ToList();
            this.tol = tol;
        }

        public DotSearchMode(string FileNameAddition, IEnumerable<double> acceptableSortedMassShifts, Tolerance tol) : base(FileNameAddition)
        {
            this.acceptableSortedMassShifts = acceptableSortedMassShifts.ToList();
            this.tol = tol;
        }

        #endregion Public Constructors

        #region Public Methods

        public override bool Accepts(double scanPrecursorMass, double peptideMass)
        {
            foreach (double huh in acceptableSortedMassShifts)
            {
                // TODO: verify that the theoretical and experimental ones make sense here
                if (tol.Within(scanPrecursorMass, huh + peptideMass))
                    return true;
            }
            return false;
        }

        public override IEnumerable<DoubleRange> GetAllowedPrecursorMassIntervals(double peptideMonoisotopicMass)
        {
            foreach (double huh in acceptableSortedMassShifts)
            {
                yield return new DoubleRange(peptideMonoisotopicMass + huh, tol);
            }
        }

        #endregion Public Methods

    }
}