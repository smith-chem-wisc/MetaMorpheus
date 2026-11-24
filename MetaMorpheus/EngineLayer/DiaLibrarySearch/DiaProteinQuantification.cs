#nullable enable
using System;
using System.Collections.Generic;
using System.Linq;
using MassSpectrometry;
using Omics;
using Omics.BioPolymerGroup;
using Omics.SpectralMatch;
using Quantification;
using Quantification.Strategies;

namespace EngineLayer.DiaLibrarySearch;

/// <summary>
/// One identified precursor in one run, as mzLib's <see cref="QuantificationEngine"/> reads a spectral match: a single
/// intensity (its quantity) attributed to exactly one peptide.
/// </summary>
/// <remarks>
/// Deliberately not a <see cref="SpectralMatch"/>. Equality is the run and the library precursor, so two charge states
/// of one peptide stay two measurements; the engine sums them into the peptide.
/// </remarks>
public sealed class DiaQuantifiedPrecursor : ISpectralMatch
{
    private readonly IBioPolymerWithSetMods _peptide;

    /// <param name="peptide">The one peptide the quantity is attributed to.</param>
    /// <exception cref="ArgumentNullException">An argument is null.</exception>
    public DiaQuantifiedPrecursor(string fullFilePath, DiaPrecursorMatch precursor, IBioPolymerWithSetMods peptide)
    {
        FullFilePath = fullFilePath ?? throw new ArgumentNullException(nameof(fullFilePath));
        Precursor = precursor ?? throw new ArgumentNullException(nameof(precursor));
        _peptide = peptide ?? throw new ArgumentNullException(nameof(peptide));
        Intensities = double.IsFinite(precursor.Quantity) && precursor.Quantity > 0 ? [precursor.Quantity] : null;
    }

    public DiaPrecursorMatch Precursor { get; }
    public string FullFilePath { get; }
    public bool IsDecoy => Precursor.IsDecoy;
    public string Accession => _peptide.Parent.Accession;
    public string FullSequence => Precursor.FullSequence;
    public string BaseSequence => _peptide.BaseSequence;

    /// <summary>The one-based library index of the precursor; a DIA identification has no single scan.</summary>
    public int OneBasedScanNumber => Precursor.PrecursorIndex + 1;

    public double Score => Precursor.Score;

    /// <summary>The precursor's quantity, or null when it was not measured (never NaN or zero).</summary>
    public double[]? Intensities { get; }

    public IEnumerable<IBioPolymerWithSetMods> GetIdentifiedBioPolymersWithSetMods() => [_peptide];

    /// <summary>Higher score first, as MetaMorpheus orders spectral matches.</summary>
    public int CompareTo(ISpectralMatch? other) => other is null ? -1 : other.Score.CompareTo(Score);

    public bool Equals(ISpectralMatch? other) => other is DiaQuantifiedPrecursor p
        && string.Equals(FullFilePath, p.FullFilePath, StringComparison.Ordinal) && Precursor.PrecursorIndex == p.Precursor.PrecursorIndex;

    public override bool Equals(object? obj) => Equals(obj as ISpectralMatch);
    public override int GetHashCode() => HashCode.Combine(FullFilePath, Precursor.PrecursorIndex);
}

/// <summary>
/// Protein quantification for the DIA search through mzLib's <see cref="QuantificationEngine"/>, with the conservative
/// strategies MetaMorpheus uses for isobaric quantification: sums, nothing normalised, unique peptides only.
/// </summary>
public static class DiaProteinQuantification
{
    /// <summary>
    /// The run's reported precursors as quantification records: targets passing <paramref name="qValueThreshold"/>
    /// whose peptide <paramref name="inference"/> mapped to a protein. Each is attributed to one peptide: the one a
    /// protein group holds as unique if any does, else one a group holds as shared, else the first the database yields.
    /// </summary>
    /// <exception cref="ArgumentNullException">An argument is null.</exception>
    public static List<DiaQuantifiedPrecursor> QuantifiedPrecursors(string fullFilePath, IReadOnlyList<DiaPrecursorMatch> precursors,
        DiaProteinInferenceResult inference, double qValueThreshold)
    {
        ArgumentNullException.ThrowIfNull(fullFilePath);
        ArgumentNullException.ThrowIfNull(precursors);
        ArgumentNullException.ThrowIfNull(inference);

        var attribution = new Dictionary<string, IBioPolymerWithSetMods>();
        var targetGroups = inference.ProteinGroups.Where(g => !g.IsDecoy).ToList();
        foreach (var peptide in targetGroups.SelectMany(g => g.UniquePeptides).Concat(targetGroups.SelectMany(g => g.AllPeptides))
                     .Concat(inference.SpectralMatches.Where(m => !m.IsDecoy)
                         .SelectMany(m => m.BestMatchingBioPolymersWithSetMods.Select(b => b.SpecificBioPolymer))))
        {
            attribution.TryAdd(peptide.FullSequence, peptide);
        }

        return precursors
            .Where(p => !p.IsDecoy && p.QValue <= qValueThreshold && attribution.ContainsKey(p.FullSequence))
            .Select(p => new DiaQuantifiedPrecursor(fullFilePath, p, attribution[p.FullSequence]))
            .ToList();
    }

    /// <summary>
    /// Quantifies peptides and protein groups across the runs the precursors come from, one label-free sample per run.
    /// </summary>
    /// <param name="outputDirectory">Where mzLib writes its quantification tables; null writes nothing.</param>
    /// <returns>mzLib's results. A failure is reported in them (Success false, with a summary), never thrown.</returns>
    /// <exception cref="ArgumentNullException">An argument is null.</exception>
    public static QuantificationResults Quantify(IReadOnlyList<DiaQuantifiedPrecursor> precursors, IReadOnlyList<ProteinGroup> proteinGroups,
        string? outputDirectory)
    {
        ArgumentNullException.ThrowIfNull(precursors);
        ArgumentNullException.ThrowIfNull(proteinGroups);

        // One biological replicate per run, so no two runs can read as the same sample
        var design = SampleExperimentalDesign.LabelFree(precursors.Select(p => p.FullFilePath)
            .Distinct(StringComparer.Ordinal).OrderBy(path => path, StringComparer.Ordinal)
            .Select((path, i) => new SpectraFileInfo(path, "", i, 0, 0)));

        bool write = outputDirectory is not null;
        var parameters = new QuantificationParameters
        {
            SpectralMatchNormalizationStrategy = new NoNormalization(),
            SpectralMatchToPeptideRollUpStrategy = new SumRollUp(),
            PeptideNormalizationStrategy = new NoNormalization(),
            CollapseStrategy = new NoCollapse(),
            CollapseAggregationStrategy = new SumAggregation(),
            PeptideToProteinRollUpStrategy = new SumRollUp(),
            ProteinNormalizationStrategy = new NoNormalization(),
            OutputDirectory = outputDirectory ?? "",
            WriteRawInformation = write,
            WritePeptideInformation = write,
            WriteProteinInformation = write,
        };

        var peptides = precursors.SelectMany(p => p.GetIdentifiedBioPolymersWithSetMods()).Distinct().ToList();
        var groups = proteinGroups.Where(g => !g.IsDecoy).Cast<IBioPolymerGroup>().ToList();
        return new QuantificationEngine(parameters, design, precursors.Cast<ISpectralMatch>().ToList(), peptides, groups).Run();
    }
}
