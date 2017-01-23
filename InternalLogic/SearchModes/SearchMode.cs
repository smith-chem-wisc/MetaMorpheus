using Spectra;
using System.Collections.Generic;

namespace InternalLogicEngineLayer
{
    public abstract class SearchMode
    {

        #region Protected Constructors

        protected SearchMode(string fileNameAddition)
        {
            FileNameAddition = fileNameAddition;
        }

        #endregion Protected Constructors

        #region Public Properties

        public string FileNameAddition { get; private set; }

        #endregion Public Properties

        #region Public Methods

        public abstract bool Accepts(double scanPrecursorMass, double peptideMass);

        public abstract IEnumerable<DoubleRange> GetAllowedPrecursorMassIntervals(double peptideMonoisotopicMass);

        #endregion Public Methods

    }
}