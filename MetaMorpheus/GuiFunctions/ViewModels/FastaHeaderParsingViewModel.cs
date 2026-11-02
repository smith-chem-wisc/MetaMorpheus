#nullable enable
using System;
using System.Collections.Generic;
using EngineLayer.DatabaseLoading;

namespace GuiFunctions;

/// <summary>
/// The FASTA header panel shared by the task windows. No WPF types, so it is testable without a window.
/// </summary>
public class FastaHeaderParsingViewModel : BaseViewModel
{
    private FastaHeaderFormat _headerFormat;
    private string _customAccessionRegex;
    private string _customFullNameRegex;
    private string _customNameRegex;
    private string _customGeneNameRegex;
    private string _customOrganismRegex;
    private string _customOrganismIdRegex;

    public FastaHeaderParsingViewModel(FastaHeaderParsingParameters? parameters)
    {
        parameters ??= new FastaHeaderParsingParameters();
        _headerFormat = parameters.HeaderFormat;
        _customAccessionRegex = parameters.CustomAccessionRegex;
        _customFullNameRegex = parameters.CustomFullNameRegex;
        _customNameRegex = parameters.CustomNameRegex;
        _customGeneNameRegex = parameters.CustomGeneNameRegex;
        _customOrganismRegex = parameters.CustomOrganismRegex;
        _customOrganismIdRegex = parameters.CustomOrganismIdRegex;
    }

    public IReadOnlyList<FastaHeaderFormat> AvailableFormats { get; } = Enum.GetValues<FastaHeaderFormat>();

    public FastaHeaderFormat HeaderFormat
    {
        get => _headerFormat;
        set
        {
            _headerFormat = value;
            OnPropertyChanged(nameof(HeaderFormat));
            OnPropertyChanged(nameof(IsCustom));
            OnPropertyChanged(nameof(PresetSummary));
        }
    }

    public bool IsCustom => HeaderFormat == FastaHeaderFormat.Custom;

    public string PresetSummary => DescribePreset(HeaderFormat);

    public string CustomAccessionRegex
    {
        get => _customAccessionRegex;
        set { _customAccessionRegex = value; OnPropertyChanged(nameof(CustomAccessionRegex)); }
    }

    public string CustomFullNameRegex
    {
        get => _customFullNameRegex;
        set { _customFullNameRegex = value; OnPropertyChanged(nameof(CustomFullNameRegex)); }
    }

    public string CustomNameRegex
    {
        get => _customNameRegex;
        set { _customNameRegex = value; OnPropertyChanged(nameof(CustomNameRegex)); }
    }

    public string CustomGeneNameRegex
    {
        get => _customGeneNameRegex;
        set { _customGeneNameRegex = value; OnPropertyChanged(nameof(CustomGeneNameRegex)); }
    }

    public string CustomOrganismRegex
    {
        get => _customOrganismRegex;
        set { _customOrganismRegex = value; OnPropertyChanged(nameof(CustomOrganismRegex)); }
    }

    public string CustomOrganismIdRegex
    {
        get => _customOrganismIdRegex;
        set { _customOrganismIdRegex = value; OnPropertyChanged(nameof(CustomOrganismIdRegex)); }
    }

    public FastaHeaderParsingParameters ToParameters() => new(HeaderFormat,
        CustomAccessionRegex, CustomFullNameRegex, CustomNameRegex,
        CustomGeneNameRegex, CustomOrganismRegex, CustomOrganismIdRegex);

    // Throws on a format with no entry, so a new FastaHeaderFormat member cannot ship a blank label.
    public static string DescribePreset(FastaHeaderFormat format) => format switch
    {
        FastaHeaderFormat.UniProt => "sp|P12345|AATM_RABIT ... OS= GN=",
        FastaHeaderFormat.Ensembl => "ENSP00000001 pep:... gene:ENSG00000001",
        FastaHeaderFormat.Gencode => "ENSP0000001.1|ENST...|...|gene symbol|...",
        FastaHeaderFormat.Auto => "detected per file from its first header; changes Name and taxonomy id - see tooltip",
        FastaHeaderFormat.Custom => "",
        _ => throw new ArgumentOutOfRangeException(nameof(format), format, "No description for this FASTA header format."),
    };
}
