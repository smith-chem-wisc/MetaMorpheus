#nullable enable
using System;
using System.Linq;
using Chemistry;
using MonoIsotopicMassSpectrometry;
using MzLibUtil;

namespace EngineLayer.Util;

/// <summary>
/// Each precursor refers to one isotopic envelope observed in a single MS1 scan. 
/// It is meant to hold experimentally observed data for a peptide/proteoform precursor.
/// </summary>
public record Precursor
{
    /// <summary>
    /// Each precursor refers to one isotopic envelope observed in a single MS1 scan. 
    /// It is meant to hold experimentally observed data for a peptide/proteoform precursor.
    /// </summary>
    /// <param name="MonoisotopicPeakMz"></param>
    /// <param name="Charge"></param>
    /// <param name="MonoIsotopicMass"></param>
    /// <param name="Intensity">Either the most abundant isotope intensity, or the summed intensity from all isotopes. For MM decon, this is determined by CommonParameters </param>
    /// <param name="EnvelopePeakCount"> The number of peaks observed in the Isotopic Envelope as determined by the deconvolution. 
    /// The minimum intensity for a peak to be considered is specified by the user in decon parameters, or in arguments passed to external software</param>
    /// <param name="FractionalIntensity"> The fraction of the intensity in the MS1 scan that is accounted for by this precursor</param>
    public Precursor(double MonoisotopicPeakMz, int Charge, double MonoIsotopicMass, double Intensity, int EnvelopePeakCount, double? FractionalIntensity = null)
    {
        this.MonoisotopicPeakMz = MonoisotopicPeakMz;
        this.Charge = Charge;
        this.MonoIsotopicMass = MonoIsotopicMass;
        this.Intensity = Intensity;
        this.EnvelopePeakCount = EnvelopePeakCount;
        this.FractionalIntensity = FractionalIntensity;
    }

    public Precursor(double MonoisotopicPeakMz, int Charge, double Intensity, int EnvelopePeakCount, double? FractionalIntensity = null)
    {
        this.MonoisotopicPeakMz = MonoisotopicPeakMz;
        this.Charge = Charge;
        MonoIsotopicMass = MonoisotopicPeakMz.ToMonoIsotopicMass(Charge);
        this.Intensity = Intensity;
        this.EnvelopePeakCount = EnvelopePeakCount;
        this.FractionalIntensity = FractionalIntensity;
    }

    public Precursor(int Charge, double MonoIsotopicMass, double Intensity, int EnvelopePeakCount, double? FractionalIntensity = null)
    {
        MonoisotopicPeakMz = MonoIsotopicMass.ToMz(Charge);
        this.Charge = Charge;
        this.MonoIsotopicMass = MonoIsotopicMass;
        this.Intensity = Intensity;
        this.EnvelopePeakCount = EnvelopePeakCount;
        this.FractionalIntensity = FractionalIntensity;
    }

    public Precursor(IsotopicEnvelope envelope, double? intensity = null, double? fractionalIntensity = null)
    {
        this.Envelope = envelope ?? throw new ArgumentNullException(nameof(envelope));
        this.MonoisotopicPeakMz = envelope.MonoisotopicMonoIsotopicMass.ToMz(envelope.Charge);
        this.Charge = envelope.Charge;
        this.MonoIsotopicMass = envelope.MonoisotopicMonoIsotopicMass;
        this.Intensity = intensity ?? envelope.Peaks.Sum(p => p.intensity);
        this.EnvelopePeakCount = envelope.Peaks.Count;
        this.FractionalIntensity = fractionalIntensity;
    }

    /// <summary>
    /// The isotopic envelope associated with this precursor, if available.
    /// </summary>
    /// <remarks>
    /// This may be null if the precursor was created without an envelope, such as coming from the scan header or an external source.
    /// </remarks>
    public IsotopicEnvelope? Envelope { get; init; }

    /// <summary>
    /// The m/z of the monoisotopic peak in the isotopic envelope.
    /// </summary>
    public double MonoisotopicPeakMz { get; init; }

    public int Charge { get; init; }

    /// <summary>
    /// The neutral monoisotopic mass of the precursor.
    /// </summary>
    public double MonoIsotopicMass { get; init; }

    /// <summary>
    /// Either the most abundant isotope intensity, or the summed intensity from all isotopes. For MM decon, this is determined by CommonParameters 
    /// </summary>
    public double Intensity { get; init; }

    /// <summary> 
    /// The number of peaks observed in the Isotopic Envelope. 
    /// The minimum intensity for a peak to be considered is specified by the user in decon parameters, or in arguments passed to external software
    /// </summary>
    public int EnvelopePeakCount { get; init; }

    /// <summary> 
    /// The fraction of the intensity in the MS1 scan isolation window that is accounted for by this precursor
    /// </summary>
    public double? FractionalIntensity { get; init; }

    public bool Equals(Precursor? other, Tolerance mzTolerance)
    {
        if (other == null) return false;
        if (this.Charge != other.Charge) return false;
        if (!mzTolerance.Within(this.MonoisotopicPeakMz, other.MonoisotopicPeakMz)) return false;

        return true;
    }

    public override int GetHashCode() => HashCode.Combine(MonoisotopicPeakMz.GetHashCode(), Charge.GetHashCode());
}