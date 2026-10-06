#nullable enable
using System;
using System.IO;
using MassSpectrometry;
using Readers;

namespace GuiFunctions;

public sealed class FromFileDeconParamsViewModel : DeconParamsViewModel
{
    private FromFileDeconvolutionParameters _parameters;

    public override DeconvolutionParameters Parameters
    {
        get => _parameters;
        protected set
        {
            _parameters = (FromFileDeconvolutionParameters)value;
            OnPropertyChanged(nameof(Parameters));
        }
    }

    public FromFileDeconParamsViewModel(FromFileDeconvolutionParameters parameters)
    {
        Parameters = parameters;
    }

    public string? FilePath
    {
        get => _parameters.FilePath;
        set
        {
            if (string.Equals(_parameters.FilePath, value, StringComparison.Ordinal))
                return;
            _parameters.FilePath = value;
            OnPropertyChanged(nameof(FilePath));
        }
    }

    public bool UseGenericScore
    {
        get => _parameters.UseGenericScore;
        set
        {
            if (_parameters.UseGenericScore == value)
                return;
            _parameters.UseGenericScore = value;
            OnPropertyChanged(nameof(UseGenericScore));
        }
    }

    public static bool TryValidateFilePath(string? path, out string validationMessage)
    {
        if (string.IsNullOrWhiteSpace(path))
        {
            validationMessage = "A feature file path is required.";
            return false;
        }

        if (!Path.IsPathFullyQualified(path))
        {
            validationMessage = $"A fully qualified absolute path is required; '{path}' is not one.";
            return false;
        }

        if (!File.Exists(path))
        {
            validationMessage = $"File not found: {path}";
            return false;
        }

        try
        {
            var resultType = path.GetResultFileType();
            if (resultType.GetInterfaces().Contains(typeof(IMs1FeatureFile)))
            {
                validationMessage = string.Empty;
                return true;
            }

            validationMessage =
                $"'{Path.GetFileName(path)}' is not a recognised MS1 feature file " +
                $"(detected type: {resultType.Name}). " +
                "Expected a FlashDeconv/TopFD _ms1.feature or Dinosaur .feature.tsv file.";
            return false;
        }
        catch (Exception ex)
        {
            validationMessage = $"Cannot read file: {ex.Message}";
            return false;
        }
    }

    public override string ToString() => "FromFile";
}
