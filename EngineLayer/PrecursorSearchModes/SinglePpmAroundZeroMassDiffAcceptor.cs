using MzLibUtil;
using System;
using System.Collections.Generic;

namespace EngineLayer
{
    public class SinglePpmAroundZeroSearchMode : MassDiffAcceptor
    {
        private readonly double ppmTolerance;

        public SinglePpmAroundZeroSearchMode(double ppmTolerance) : base(ppmTolerance + "ppmAroundZero")
        {
            this.ppmTolerance = ppmTolerance;
        }

        public override int Accepts(double scanPrecursorMass, double peptideMass)
        {
            return Math.Abs((scanPrecursorMass - peptideMass) / (peptideMass) * 1e6) < ppmTolerance ? 0 : -1;
        }

        public override IEnumerable<AllowedIntervalWithNotch> GetAllowedPrecursorMassIntervals(double peptideMonoisotopicMass)
        {
            var diff = ppmTolerance / 1e6 * peptideMonoisotopicMass;
            yield return new AllowedIntervalWithNotch(new DoubleRange(peptideMonoisotopicMass - diff, peptideMonoisotopicMass + diff), 0);
        }

        public override string ToProseString()
        {
            return (String.Format("{0:0.0}", ppmTolerance) + " ppm around zero");
        }

        public override string ToString()
        {
            return FileNameAddition + " ppmAroundZero " + ppmTolerance;
        }
    }
}