#nullable enable
using System.Linq;
using MassSpectrometry;
using MathNet.Numerics;
using OxyPlot;

namespace GuiFunctions;

/// <summary>
/// Represents a deconvoluted species for visualization. No need to cache results as each property is only called once with current implementation. 
/// </summary>
public class DeconvolutedSpeciesViewModel(IsotopicEnvelope envelope) : BaseViewModel
{
    public IsotopicEnvelope Envelope { get; } = envelope;
    public OxyColor Color { get; set; } = OxyColors.Automatic;
    public double MonoisotopicMass => Envelope.MonoisotopicMass;
    public int Charge => Envelope.Charge;
    public int PeakCount => Envelope.Peaks.Count;
    public double Intensity => Envelope.Peaks.Sum(p => p.intensity);
    public string Annotation => $"M={MonoisotopicMass.Round(2)}\nz={Charge}";
    public double MostAbundantMz => Envelope.Peaks.MaxBy(p => p.intensity).mz;
    public string PeakMzs => string.Join(", ", Envelope.Peaks.Select(p => p.mz.Round(2)).OrderBy(p => p));
}