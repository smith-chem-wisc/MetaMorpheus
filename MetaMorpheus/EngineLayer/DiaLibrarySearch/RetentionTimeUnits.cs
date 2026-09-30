#nullable enable
using System;

namespace EngineLayer.DiaLibrarySearch;

/// <summary>
/// An indexed retention time: the run-independent scale spectral libraries store (Biognosys/Prosit iRT).
/// Distinct from <see cref="RtMinutes"/> with no conversion either way; an <see cref="IIrtMap"/> is the only bridge.
/// </summary>
public readonly record struct Irt(double Value)
{
    public override string ToString() => $"{Value} iRT";
}

/// <summary>A retention time in minutes in one run, as <see cref="MassSpectrometry.MsDataScan.RetentionTime"/> reports it.</summary>
public readonly record struct RtMinutes(double Value)
{
    public override string ToString() => $"{Value} min";
}

/// <summary>Maps one run's retention times onto the library's iRT scale and back.</summary>
public interface IIrtMap
{
    Irt ToIrt(RtMinutes retentionTime);

    RtMinutes ToRtMinutes(Irt indexedRetentionTime);
}

/// <summary>iRT = <see cref="Slope"/> × minutes + <see cref="Intercept"/>.</summary>
public sealed class LinearIrtMap : IIrtMap
{
    public double Slope { get; }
    public double Intercept { get; }

    /// <exception cref="ArgumentOutOfRangeException">The slope is zero or either value is not finite.</exception>
    public LinearIrtMap(double slope, double intercept)
    {
        if (!double.IsFinite(slope) || slope == 0)
            throw new ArgumentOutOfRangeException(nameof(slope), slope, "The slope must be finite and nonzero.");
        if (!double.IsFinite(intercept))
            throw new ArgumentOutOfRangeException(nameof(intercept), intercept, "The intercept must be finite.");
        Slope = slope;
        Intercept = intercept;
    }

    /// <summary>The map through two anchors, each a retention time and the iRT eluting there.</summary>
    public static LinearIrtMap ThroughPoints((RtMinutes Rt, Irt Irt) first, (RtMinutes Rt, Irt Irt) second)
    {
        double slope = (second.Irt.Value - first.Irt.Value) / (second.Rt.Value - first.Rt.Value);
        return new LinearIrtMap(slope, first.Irt.Value - slope * first.Rt.Value);
    }

    public Irt ToIrt(RtMinutes retentionTime) => new(Slope * retentionTime.Value + Intercept);

    public RtMinutes ToRtMinutes(Irt indexedRetentionTime) => new((indexedRetentionTime.Value - Intercept) / Slope);
}
