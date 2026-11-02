#nullable enable
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text.RegularExpressions;
using UsefulProteomicsDatabases;

namespace EngineLayer.DatabaseLoading;

/// <summary>
/// Named FASTA header layouts that MetaMorpheus can parse. <see cref="FastaHeaderFormat.Custom"/>
/// takes the regexes from <see cref="FastaHeaderParsingParameters"/> instead of a preset;
/// <see cref="FastaHeaderFormat.Auto"/> passes none and lets mzLib detect the format per file.
/// </summary>
public enum FastaHeaderFormat
{
    UniProt,
    Ensembl,
    Gencode,
    Auto,
    Custom
}

/// <summary>
/// The six regexes handed to <see cref="ProteinDbLoader.LoadProteinFasta"/>. A null member means the
/// field is not extracted.
/// </summary>
public record FastaHeaderFieldRegexSet
{
    public FastaHeaderFieldRegex? Accession { get; init; }
    public FastaHeaderFieldRegex? FullName { get; init; }
    public FastaHeaderFieldRegex? Name { get; init; }
    public FastaHeaderFieldRegex? GeneName { get; init; }
    public FastaHeaderFieldRegex? Organism { get; init; }
    public FastaHeaderFieldRegex? OrganismId { get; init; }
}

/// <summary>
/// How to pull accession, name, gene and organism out of a protein FASTA defline.
/// Serialized as a table under [CommonParameters.FastaHeaderParsing]; every member is a string or an
/// enum so Nett needs no custom converter for the fields themselves.
/// </summary>
public class FastaHeaderParsingParameters
{
    /// <summary>
    /// Ceiling on a single validation match. Applies to validation only: mzLib's
    /// <see cref="FastaHeaderFieldRegex"/> compiles the pattern with no timeout, so a pattern that is
    /// fast on the headers checked here can still stall a load. Checking real deflines in
    /// <see cref="ValidateAgainstHeaders"/> is what narrows that gap.
    /// </summary>
    public static readonly TimeSpan ValidationMatchTimeout = TimeSpan.FromMilliseconds(250);

    /// <summary>
    /// Deflines every custom pattern is matched against before any file is read. Only used to smoke out
    /// catastrophic backtracking; a pattern that matches nothing here is still accepted.
    /// </summary>
    internal static readonly string[] ValidationHeaders =
    [
        ">sp|P12345|AATM_RABIT Aspartate aminotransferase OS=Oryctolagus cuniculus OX=9986 GN=GOT2 PE=1 SV=2",
        ">gi|16128008|ref|NP_414555.1| thr operon leader peptide [Escherichia coli str. K-12 substr. MG1655]",
        ">id:1|Name:tdbR00000010|SOterm:SO:0000254|Type:tRNA|Subtype:Ala|Feature:VGC|Cellular_Localization:prokaryotic cytosol|Species:Escherichia coli"
    ];

    public FastaHeaderFormat HeaderFormat { get; set; } = FastaHeaderFormat.UniProt;

    public string CustomAccessionRegex { get; set; } = "";
    public string CustomFullNameRegex { get; set; } = "";
    public string CustomNameRegex { get; set; } = "";
    public string CustomGeneNameRegex { get; set; } = "";
    public string CustomOrganismRegex { get; set; } = "";
    public string CustomOrganismIdRegex { get; set; } = "";

    public FastaHeaderParsingParameters()
    {
    }

    public FastaHeaderParsingParameters(FastaHeaderFormat format,
        string? accessionRegex = null, string? fullNameRegex = null, string? nameRegex = null,
        string? geneNameRegex = null, string? organismRegex = null, string? organismIdRegex = null)
    {
        HeaderFormat = format;
        CustomAccessionRegex = accessionRegex ?? "";
        CustomFullNameRegex = fullNameRegex ?? "";
        CustomNameRegex = nameRegex ?? "";
        CustomGeneNameRegex = geneNameRegex ?? "";
        CustomOrganismRegex = organismRegex ?? "";
        CustomOrganismIdRegex = organismIdRegex ?? "";
    }

    /// <summary>
    /// The format every earlier version used. Nothing is sampled, gated or added to the index cache key
    /// for it, so a default run behaves exactly as before.
    /// </summary>
    // A method, not a property, so it is not written to the task TOML.
    public bool IsDefault() => HeaderFormat == FastaHeaderFormat.UniProt;

    /// <summary>
    /// The checks that need no database: run before the first task, so a typo is reported before
    /// anything has run rather than after the tasks ahead of it.
    /// </summary>
    public bool ValidateConfiguration(out List<string> errors)
    {
        errors = new List<string>();

        if (HeaderFormat != FastaHeaderFormat.Custom)
            return true;

        // Not cosmetic: a null accession regex makes mzLib auto-detect, and detection assigns its own
        // regexes over the other five with '=', so every other custom pattern would be discarded.
        if (string.IsNullOrWhiteSpace(CustomAccessionRegex))
            errors.Add("Custom FASTA header format requires an accession regex. Set 'CustomAccessionRegex' " +
                       "(the first capture group is used as the accession), or pick a preset FASTA header format.");

        foreach (var (fieldName, pattern) in EnumerateCustomPatterns())
        {
            if (string.IsNullOrWhiteSpace(pattern))
                continue;
            if (!TryValidatePattern(pattern, ValidationHeaders, out string? reason))
                errors.Add($"The custom FASTA header regex for {fieldName} is not usable: {reason} Pattern: {pattern}");
        }

        return errors.Count == 0;
    }

    /// <summary>
    /// Checks this setting against the first deflines of one FASTA database. An accession regex that
    /// matches nothing makes mzLib use the whole defline as the accession, silently, so a choice the user
    /// made explicitly is refused when it matches none of <paramref name="sampledHeaders"/>.
    /// <paramref name="warnings"/> carries problems worth reporting that need not stop the run.
    /// </summary>
    /// <param name="sampledHeaders">Null or empty when nothing could be read from the file.</param>
    public bool ValidateAgainstHeaders(string databaseName, IReadOnlyList<string>? sampledHeaders,
        out List<string> errors, out List<string> warnings)
    {
        errors = new List<string>();
        warnings = new List<string>();

        if (IsDefault())
            return true;

        if (sampledHeaders is not { Count: > 0 })
        {
            if (HeaderFormat is FastaHeaderFormat.Auto or FastaHeaderFormat.Custom)
                errors.Add($"No FASTA header could be read from {databaseName}, so the {HeaderFormat} FASTA " +
                           "header format cannot be checked against it.");
            return errors.Count == 0;
        }

        FastaHeaderFieldRegex? accession;
        switch (HeaderFormat)
        {
            case FastaHeaderFormat.Auto:
                // mzLib detects once, from the first defline, and parses the whole file with that preset.
                var detected = ProteinDbLoader.DetectFastaHeaderFormat(sampledHeaders[0]);
                if (detected == FastaHeaderType.Unknown)
                {
                    errors.Add($"The FASTA header format of {databaseName} could not be detected automatically " +
                               $"from its first header: {sampledHeaders[0]}. Pick a preset FASTA header format, " +
                               "or use Custom and supply an accession regex.");
                    return false;
                }
                accession = MzLibPreset(detected).Accession;
                break;

            case FastaHeaderFormat.Custom:
                // Timed first: the loader will apply these with no timeout at all.
                foreach (var (fieldName, pattern) in EnumerateCustomPatterns())
                {
                    if (!string.IsNullOrWhiteSpace(pattern) && !TryValidatePattern(pattern, sampledHeaders, out string? reason))
                        errors.Add($"The custom FASTA header regex for {fieldName} is not usable on {databaseName}: {reason} Pattern: {pattern}");
                }
                if (errors.Count > 0)
                    return false;
                accession = GetFieldRegexes().Accession;
                break;

            default:
                var expected = ToFastaHeaderType(HeaderFormat);
                var looksLike = ProteinDbLoader.DetectFastaHeaderFormat(sampledHeaders[0]);
                if (looksLike != FastaHeaderType.Unknown && looksLike != expected)
                    warnings.Add($"{databaseName} looks like a {looksLike} database, but the FASTA header format " +
                                 $"is set to {HeaderFormat}, which applies to every FASTA in the task.");
                accession = GetFieldRegexes().Accession;
                break;
        }

        // An empty capture falls back to the whole defline in mzLib just as a miss does.
        if (!sampledHeaders.Any(h => !string.IsNullOrEmpty(accession?.ApplyRegex(h))))
        {
            string formatLabel = HeaderFormat == FastaHeaderFormat.Custom ? "custom accession regex" : $"{HeaderFormat} accession regex";
            errors.Add($"The {formatLabel} matches none of the first {sampledHeaders.Count} headers of {databaseName}, " +
                       $"so every accession would be the whole header line. First header: {sampledHeaders[0]}");
        }

        return errors.Count == 0;
    }

    /// <summary>
    /// Builds the regexes to hand to mzLib. Call the Validate methods first; this throws on a pattern that
    /// does not compile, as a backstop for callers that bypass the runner.
    /// </summary>
    public FastaHeaderFieldRegexSet GetFieldRegexes()
    {
        switch (HeaderFormat)
        {
            case FastaHeaderFormat.UniProt:
                // Not mzLib's UniProt bundle: MetaMorpheus has always passed the full-name regex as Name
                // and never parsed OrganismId. Matching mzLib moves Name on every entry.
                return MzLibPreset(FastaHeaderType.UniProt) with
                {
                    Name = ProteinDbLoader.UniprotFullNameRegex,
                    OrganismId = null,
                };

            case FastaHeaderFormat.Ensembl:
            case FastaHeaderFormat.Gencode:
                return MzLibPreset(ToFastaHeaderType(HeaderFormat));

            // All null on purpose: LoadProteinFasta only calls DetectFastaHeaderFormat when the
            // accession regex is null, so an empty set is how detection is requested.
            case FastaHeaderFormat.Auto:
                return new FastaHeaderFieldRegexSet();

            case FastaHeaderFormat.Custom:
                return new FastaHeaderFieldRegexSet
                {
                    Accession = Build("accession", CustomAccessionRegex),
                    FullName = Build("fullName", CustomFullNameRegex),
                    Name = Build("name", CustomNameRegex),
                    GeneName = Build("geneName", CustomGeneNameRegex),
                    Organism = Build("organism", CustomOrganismRegex),
                    OrganismId = Build("organismId", CustomOrganismIdRegex),
                };

            default:
                throw new MetaMorpheusException($"Unrecognized FASTA header format: {HeaderFormat}");
        }
    }

    /// <summary>
    /// What mzLib's LoadProteinFasta assigns after detecting <paramref name="type"/>.
    /// </summary>
    // TODO: replace with mzLib's public per-format bundles, and add Ncbi, once the mzLib preset PR
    // that follows MetaMorpheus #2782 is released.
    internal static FastaHeaderFieldRegexSet MzLibPreset(FastaHeaderType type) => type switch
    {
        FastaHeaderType.UniProt => new FastaHeaderFieldRegexSet
        {
            Accession = ProteinDbLoader.UniprotAccessionRegex,
            FullName = ProteinDbLoader.UniprotFullNameRegex,
            Name = ProteinDbLoader.UniprotNameRegex,
            GeneName = ProteinDbLoader.UniprotGeneNameRegex,
            Organism = ProteinDbLoader.UniprotOrganismRegex,
            OrganismId = ProteinDbLoader.UniprotOrganismIdRegex,
        },
        FastaHeaderType.Ensembl => new FastaHeaderFieldRegexSet
        {
            Accession = ProteinDbLoader.EnsemblAccessionRegex,
            FullName = ProteinDbLoader.EnsemblFullNameRegex,
            GeneName = ProteinDbLoader.EnsemblGeneNameRegex,
        },
        FastaHeaderType.Gencode => new FastaHeaderFieldRegexSet
        {
            Accession = ProteinDbLoader.GencodeAccessionRegex,
            FullName = ProteinDbLoader.GencodeFullNameRegex,
            GeneName = ProteinDbLoader.GencodeGeneNameRegex,
        },
        _ => throw new ArgumentOutOfRangeException(nameof(type), type, "mzLib has no preset for this FASTA header type."),
    };

    private static FastaHeaderType ToFastaHeaderType(FastaHeaderFormat format) => format switch
    {
        FastaHeaderFormat.UniProt => FastaHeaderType.UniProt,
        FastaHeaderFormat.Ensembl => FastaHeaderType.Ensembl,
        FastaHeaderFormat.Gencode => FastaHeaderType.Gencode,
        _ => throw new ArgumentOutOfRangeException(nameof(format), format, "Not a preset FASTA header format."),
    };

    /// <summary>
    /// This setting as a line of the index cache key. Accessions are how a cached peptide index finds its
    /// proteins again, so an index built under one format must not be offered to another.
    /// </summary>
    public string DescribeForIndexKey() => HeaderFormat != FastaHeaderFormat.Custom
        ? HeaderFormat.ToString()
        : HeaderFormat + "; " + string.Join("; ", EnumerateCustomPatterns().Select(p => p.FieldName + "=" + p.Pattern));

    private static FastaHeaderFieldRegex? Build(string fieldName, string pattern) =>
        string.IsNullOrWhiteSpace(pattern) ? null : new FastaHeaderFieldRegex(fieldName, pattern, 0, 1);

    private IEnumerable<(string FieldName, string Pattern)> EnumerateCustomPatterns()
    {
        yield return ("accession", CustomAccessionRegex);
        yield return ("full name", CustomFullNameRegex);
        yield return ("name", CustomNameRegex);
        yield return ("gene name", CustomGeneNameRegex);
        yield return ("organism", CustomOrganismRegex);
        yield return ("organism id", CustomOrganismIdRegex);
    }

    private static bool TryValidatePattern(string pattern, IEnumerable<string> headers, out string? reason)
    {
        Regex compiled;
        try
        {
            compiled = new Regex(pattern, RegexOptions.None, ValidationMatchTimeout);
        }
        catch (ArgumentException e)
        {
            reason = e.Message;
            return false;
        }

        if (compiled.GetGroupNumbers().Length < 2)
        {
            reason = "it has no capture group, so no field value can be taken from it.";
            return false;
        }

        foreach (string header in headers)
        {
            try
            {
                compiled.Match(header);
            }
            catch (RegexMatchTimeoutException)
            {
                reason = $"matching it against a header took longer than {ValidationMatchTimeout.TotalMilliseconds} ms.";
                return false;
            }
        }

        reason = null;
        return true;
    }
}
