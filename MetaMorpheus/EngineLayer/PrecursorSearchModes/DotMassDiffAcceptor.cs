using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Linq;

namespace EngineLayer
{
    public class DotMassDiffAcceptor : MassDiffAcceptor
    {
        private readonly double[] AcceptableSortedMassShifts;
        private readonly Tolerance Tolerance;

        public DotMassDiffAcceptor(string FileNameAddition, IEnumerable<double> acceptableMassShifts, Tolerance tol) : base(FileNameAddition)
        {
            AcceptableSortedMassShifts = acceptableMassShifts.OrderBy(b => b).ToArray();
            Tolerance = tol;
            NumNotches = AcceptableSortedMassShifts.Length;
        }

        public override int Accepts(double scanPrecursorMass, double peptideMass)
        {
            // index of the first element that is larger than or equal to value
            int index = Array.BinarySearch(AcceptableSortedMassShifts, scanPrecursorMass - peptideMass);
            if (index < 0)
                index = ~index;

            if (index < AcceptableSortedMassShifts.Length && Tolerance.Within(scanPrecursorMass, AcceptableSortedMassShifts[index] + peptideMass))
                return index;

            if (index > 0 && Tolerance.Within(scanPrecursorMass, AcceptableSortedMassShifts[index - 1] + peptideMass))
                return index - 1;
            return -1;
        }

        public override IEnumerable<AllowedIntervalWithNotch> GetAllowedPrecursorMassIntervalsFromTheoreticalMass(double peptideMonoisotopicMass)
        {
            for (int j = 0; j < AcceptableSortedMassShifts.Length; j++)
            {
                yield return new AllowedIntervalWithNotch(Tolerance.GetRange(peptideMonoisotopicMass + AcceptableSortedMassShifts[j]), j);
            }
        }

        public override IEnumerable<AllowedIntervalWithNotch> GetAllowedPrecursorMassIntervalsFromObservedMass(double peptideMonoisotopicMass)
        {
            for (int j = 0; j < AcceptableSortedMassShifts.Length; j++)
            {
                yield return new AllowedIntervalWithNotch(Tolerance.GetRange(peptideMonoisotopicMass - AcceptableSortedMassShifts[j]), j);
            }
        }

        public override string ToString()
        {
            return FileNameAddition + " dot " + Tolerance.ToString() + " " + string.Join(",", AcceptableSortedMassShifts);
        }

        public override string ToProseString()
        {
            return (Tolerance.ToString() + " around " + String.Join(",", AcceptableSortedMassShifts) + " Da");
        }
    }
}