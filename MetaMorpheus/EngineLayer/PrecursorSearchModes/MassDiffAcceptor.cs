using System.Collections.Generic;

namespace EngineLayer
{
    public abstract class MassDiffAcceptor
    {
        public const double NotchScalar = 1E5;
        protected MassDiffAcceptor(string fileNameAddition)
        {
            FileNameAddition = fileNameAddition;
            NumNotches = 1;
        }

        public int NumNotches { get; protected set; }
        public string FileNameAddition { get; }

        /// <summary>
        /// True when this acceptor matches candidates against the observed MOST-ABUNDANT isotopic peak
        /// rather than the observed monoisotopic peak. The acceptor is the only object that knows which
        /// observed mass a search keys off; scans and spectral matches carry both observations and let
        /// the acceptor choose (see PrecursorMassExtensions.GetPrecursorMassForSearch).
        /// </summary>
        public virtual bool MatchesMostAbundantPeak => false;

        /// <summary>
        /// Converts the observed mass this acceptor matched on into the observed MONOISOTOPIC mass implied
        /// for a given candidate, undoing any isotope-level offset the acceptor deliberately allowed.
        /// Identity for acceptors that already match monoisotopic-to-monoisotopic.
        /// <para>
        /// Most-abundant acceptors admit candidates whose deconvoluted monoisotopic peak is off by whole
        /// isotopologues (the apex notch set). The resulting whole-isotope offsets are isotope-assignment
        /// artifacts, not chemistry: any consumer that reads (observed - theoretical) as a MASS DIFFERENCE
        /// - a modification to localize, a G-PTM-D candidate mod, an instrument m/z error to calibrate -
        /// must go through this method first, or it will interpret an isotope offset as a ~1-2 Da mod.
        /// </para>
        /// </summary>
        public virtual double ToObservedMonoisotopicMass(double observedPrecursorMass, double theoreticalMonoisotopicMass) => observedPrecursorMass;

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