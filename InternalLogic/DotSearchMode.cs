using Spectra;
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

        public override string SearchModeString()
        {
            return "Tolerance of " + tol.ToString() + " around mass diffs: " + string.Join(",", acceptableSortedMassShifts);
        }

        #endregion Public Methods

    }
}