using Spectra;
using System.Collections.Generic;

namespace InternalLogicEngineLayer
{
    public abstract class SearchMode
    {
        #region Public Constructors

        public SearchMode(string fileNameAddition)
        {
            FileNameAddition = fileNameAddition;
        }

        #endregion Public Constructors

        #region Public Properties

        public string FileNameAddition { get; internal set; }

        #endregion Public Properties

        #region Public Methods

        public abstract bool Accepts(double scanPrecursorMass, double peptideMass);

        public override string ToString()
        {
            return FileNameAddition + " " + SearchModeString();
        }

        #endregion Public Methods

        #region Internal Methods

        internal abstract IEnumerable<DoubleRange> GetAllowedPrecursorMassIntervals(double peptideMonoisotopicMass);

        internal abstract string SearchModeString();

        #endregion Internal Methods
    }
}