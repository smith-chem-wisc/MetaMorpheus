using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Linq;

namespace EngineLayer
{
    public class DotMassDiffAcceptor : MassDiffAcceptor
    {
        private readonly double[] AcceptableSortedMassShifts;
        private readonly int[] NotchIntegers;
        private readonly Tolerance Tolerance;

        public DotMassDiffAcceptor(string fileNameAddition, IEnumerable<double> acceptableMassShifts, Tolerance tol) : base(fileNameAddition)
        {
            AcceptableSortedMassShifts = acceptableMassShifts.OrderBy(Math.Abs).ThenBy(p => p < 0).ToArray();
            NotchIntegers = AcceptableSortedMassShifts.Select(b => (int)Math.Round(b * NotchScalar)).ToArray();
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
                return NotchIntegers[index];

            if (index > 0 && Tolerance.Within(scanPrecursorMass, AcceptableSortedMassShifts[index - 1] + peptideMass))
                return NotchIntegers[index - 1];
            return -1;
        }

        public override IEnumerable<AllowedIntervalWithNotch> GetAllowedPrecursorMassIntervalsFromTheoreticalMass(double peptideMonoisotopicMass)
        {
            for (int j = 0; j < AcceptableSortedMassShifts.Length; j++)
            {
                var notch = NotchIntegers[j];
                var mass = peptideMonoisotopicMass + notch;
                yield return new AllowedIntervalWithNotch(Tolerance.GetMinimumValue(mass), Tolerance.GetMaximumValue(mass), notch);
            }
        }

        public override IEnumerable<AllowedIntervalWithNotch> GetAllowedPrecursorMassIntervalsFromObservedMass(double peptideMonoisotopicMass)
        {
            for (int j = 0; j < AcceptableSortedMassShifts.Length; j++)
            {
                var notch = NotchIntegers[j];
                var mass = peptideMonoisotopicMass - notch;
                yield return new AllowedIntervalWithNotch(Tolerance.GetMinimumValue(mass), Tolerance.GetMaximumValue(mass), notch);
            }
        }

        public override string ToString()
        {
            return FileNameAddition + " dot " + Tolerance.ToString() + " " + string.Join(",", AcceptableSortedMassShifts);
        }

        public override string ToProseString()
        {
            return (Tolerance.ToString() + " around " + String.Join(", ", AcceptableSortedMassShifts) + " Da");
        }
    }
}
