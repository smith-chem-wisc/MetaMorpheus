using System;
using System.Collections.Generic;
using System.Linq;
using MzLibUtil;
using Omics.Modifications;

namespace EngineLayer.PrecursorSearchModes;

internal class AdductMassDifferenceAcceptor : MassDiffAcceptor
{
    private readonly double[] AcceptableSortedMassShifts;
    private readonly Tolerance Tolerance;

    public AdductMassDifferenceAcceptor(string fileNameAddition, IEnumerable<(int, Modification)> adducts, Tolerance tol ) : base(fileNameAddition)
    {
        Tolerance = tol;
        NumNotches = AcceptableSortedMassShifts.Length;
    }

    public override int Accepts(double scanPrecursorMass, double peptideMass)
    {
        throw new System.NotImplementedException();
    }

    public override IEnumerable<AllowedIntervalWithNotch> GetAllowedPrecursorMassIntervalsFromTheoreticalMass(double peptideMonoisotopicMass)
    {
        throw new System.NotImplementedException();
    }

    public override IEnumerable<AllowedIntervalWithNotch> GetAllowedPrecursorMassIntervalsFromObservedMass(double peptideMonoisotopicMass)
    {
        throw new System.NotImplementedException();
    }

    public override string ToProseString()
    {
        throw new System.NotImplementedException();
    }
}