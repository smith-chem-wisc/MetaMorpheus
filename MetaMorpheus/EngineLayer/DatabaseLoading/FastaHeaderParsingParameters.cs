#nullable enable
using System;
using System.Collections.Generic;
using System.Text.RegularExpressions;
using UsefulProteomicsDatabases;

namespace EngineLayer.DatabaseLoading;

/// <summary>
/// Named FASTA header layouts that MetaMorpheus can parse. <see cref="FastaHeaderFormat.Custom"/>
/// takes the regexes from <see cref="FastaHeaderParsingParameters"/> instead of a preset.
/// </summary>
public enum FastaHeaderFormat
{
    UniProt,
    Ensembl,
    Gencode,
    Ncbi,
    Modomics,
    Custom
}

/// <summary>
/// The six regexes handed to <see cref="ProteinDbLoader.LoadProteinFasta"/>. A null member means the
/// field is not extracted.
/// </summary>
public class FastaHeaderFieldRegexSet
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
    /// Ceiling on a single validation match. The user supplies these patterns, so a pathological one
    /// must be rejected at configuration time rather than stalling a load.
    /// </summary>
    public static readonly TimeSpan ValidationMatchTimeout = TimeSpan.FromMilliseconds(250);

    /// <summary>
    /// Deflines the validator matches each pattern against. Only used to smoke out catastrophic
    /// backtracking; a pattern that matches nothing here is still accepted.
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

    public FastaHeaderParsingParameters Clone() => new(HeaderFormat, CustomAccessionRegex,
        CustomFullNameRegex, CustomNameRegex, CustomGeneNameRegex, CustomOrganismRegex, CustomOrganismIdRegex);

    /// <summary>
    /// Checks the custom patterns without touching a database. Returns false and fills
    /// <paramref name="errors"/> with user-facing messages when the configuration cannot be used.
    /// </summary>
    public bool Validate(out List<string> errors)
    {
        errors = new List<string>();
        if (HeaderFormat != FastaHeaderFormat.Custom)
            return true;

        // mzLib falls back to header auto-detection when the accession regex is null, and that
        // detection throws on a header it does not recognize. Require one instead.
        if (string.IsNullOrWhiteSpace(CustomAccessionRegex))
            errors.Add("Custom FASTA header format requires an accession regex. Set 'CustomAccessionRegex' " +
                       "(the first capture group is used as the accession), or pick a preset FASTA header format.");

        foreach (var (fieldName, pattern) in EnumerateCustomPatterns())
        {
            if (string.IsNullOrWhiteSpace(pattern))
                continue;
            if (!TryValidatePattern(pattern, out string? reason))
                errors.Add($"The custom FASTA header regex for {fieldName} is not usable: {reason} Pattern: {pattern}");
        }

        return errors.Count == 0;
    }

    /// <summary>
    /// Builds the regexes to hand to mzLib. Call <see cref="Validate"/> first; this throws on a
    /// pattern that does not compile, as a backstop for callers that bypass the runner.
    /// </summary>
    public FastaHeaderFieldRegexSet GetFieldRegexes()
    {
        switch (HeaderFormat)
        {
            case FastaHeaderFormat.UniProt:
                // FullName is deliberately used for both FullName and Name: that is what MetaMorpheus
                // has always passed, and OrganismId has never been parsed. Changing either moves output.
                return new FastaHeaderFieldRegexSet
                {
                    Accession = ProteinDbLoader.UniprotAccessionRegex,
                    FullName = ProteinDbLoader.UniprotFullNameRegex,
                    Name = ProteinDbLoader.UniprotFullNameRegex,
                    GeneName = ProteinDbLoader.UniprotGeneNameRegex,
                    Organism = ProteinDbLoader.UniprotOrganismRegex,
                };

            case FastaHeaderFormat.Ensembl:
                return new FastaHeaderFieldRegexSet
                {
                    Accession = ProteinDbLoader.EnsemblAccessionRegex,
                    FullName = ProteinDbLoader.EnsemblFullNameRegex,
                    GeneName = ProteinDbLoader.EnsemblGeneNameRegex,
                };

            case FastaHeaderFormat.Gencode:
                return new FastaHeaderFieldRegexSet
                {
                    Accession = ProteinDbLoader.GencodeAccessionRegex,
                    FullName = ProteinDbLoader.GencodeFullNameRegex,
                    GeneName = ProteinDbLoader.GencodeGeneNameRegex,
                };

            case FastaHeaderFormat.Ncbi:
                return new FastaHeaderFieldRegexSet
                {
                    Accession = NcbiAccessionRegex,
                    FullName = NcbiFullNameRegex,
                    Organism = NcbiOrganismRegex,
                };

            case FastaHeaderFormat.Modomics:
                return new FastaHeaderFieldRegexSet
                {
                    Accession = ModomicsAccessionRegex,
                    FullName = ModomicsFullNameRegex,
                    Name = ModomicsNameRegex,
                    Organism = ModomicsOrganismRegex,
                };

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
    /// Legacy and modern NCBI deflines: >gi|16128008|ref|NP_414555.1| ... and >ref|NP_414555.1| ...
    /// The greedy prefix pushes the opening pipe as late as possible, so the capture lands between the
    /// last two pipes rather than spanning all of them.
    /// </summary>
    public static readonly FastaHeaderFieldRegex NcbiAccessionRegex = new("accession", @">.*\|(.*)\|", 0, 1);
    public static readonly FastaHeaderFieldRegex NcbiFullNameRegex = new("fullName", @">[^ ]*\s(.*?)(\s\[|$)", 0, 1);
    public static readonly FastaHeaderFieldRegex NcbiOrganismRegex = new("organism", @"\[([^\[\]]+)\]\s*$", 0, 1);

    // Patterns copied verbatim from mzLib RnaDbLoader.ModomicsFieldRegexes; only the mapping onto the
    // protein fields is ours. Accession uses Name (unique per entry) rather than the repeating SOterm.
    public static readonly FastaHeaderFieldRegex ModomicsAccessionRegex = new("accession", @"Name:(.+?)\|", 0, 1);
    public static readonly FastaHeaderFieldRegex ModomicsNameRegex = new("name", @"SOterm:(.+?)\|", 0, 1);
    public static readonly FastaHeaderFieldRegex ModomicsFullNameRegex = new("fullName", @"Type:(.+?)\|", 0, 1);
    public static readonly FastaHeaderFieldRegex ModomicsOrganismRegex = new("organism", @"Species:(.+?)$", 0, 1);

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

    private static bool TryValidatePattern(string pattern, out string? reason)
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

        foreach (string header in ValidationHeaders)
        {
            try
            {
                compiled.Match(header);
            }
            catch (RegexMatchTimeoutException)
            {
                reason = $"matching it against a sample header took longer than {ValidationMatchTimeout.TotalMilliseconds} ms.";
                return false;
            }
        }

        reason = null;
        return true;
    }
}
