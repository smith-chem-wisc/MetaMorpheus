using System.Collections.Generic;

namespace EngineLayer
{
    public abstract class SearchMode
    {

        #region Protected Constructors

        protected SearchMode(string fileNameAddition)
        {
            FileNameAddition = fileNameAddition;
            NumNotches = 1;
        }

        #endregion Protected Constructors

        #region Public Properties

        public int NumNotches { get; protected set; }
        public string FileNameAddition { get; }

        #endregion Public Properties

        #region Public Methods

        /// <summary>
        /// If acceptable, returns 0 or greater, negative means does not accept
        /// </summary>
        /// <param name="scanPrecursorMass"></param>
        /// <param name="peptideMass"></param>
        /// <returns></returns>
        public abstract int Accepts(double scanPrecursorMass, double peptideMass);

        public abstract IEnumerable<AllowedIntervalWithNotch> GetAllowedPrecursorMassIntervals(double peptideMonoisotopicMass);

        #endregion Public Methods

    }
}