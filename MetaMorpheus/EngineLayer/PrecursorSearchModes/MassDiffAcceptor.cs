using System.Collections.Generic;

namespace EngineLayer
{
    public abstract class MassDiffAcceptor
    {
        protected MassDiffAcceptor(string fileNameAddition)
        {
            FileNameAddition = fileNameAddition;
            NumNotches = 1;
        }

        public int NumNotches { get; protected set; }
        public string FileNameAddition { get; }

        /// <summary>
        /// If acceptable, returns 0 or greater, negative means does not accept
        /// </summary>
        /// <param name="scanPrecursorMass"></param>
        /// <param name="peptideMass"></param>
        /// <returns></returns>
        public abstract int Accepts(double scanPrecursorMass, double peptideMass);

        public abstract IEnumerable<AllowedIntervalWithNotch> GetAllowedPrecursorMassIntervalsFromTheoreticalMass(double peptideMonoisotopicMass);
        public abstract IEnumerable<AllowedIntervalWithNotch> GetAllowedPrecursorMassIntervalsFromObservedMass(double peptideMonoisotopicMass);

        public abstract string ToProseString();
    }
}