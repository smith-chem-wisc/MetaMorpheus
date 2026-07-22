using Chemistry;
using MassSpectrometry;
using Readers;
using System;
using EngineLayer.Deconvolution.FeatureFileMapping;

namespace EngineLayer.Deconvolution;

/// <summary>
/// Deconvolution parameters that are mapped from a feature file. Not for use in mzLib for deconvolution, but for use in the TaskLayer to pass parameters from a feature file for deconvolution. 
/// </summary>
public class FeatureMappedFromFileDeconvolutionParameters : DeconvolutionParameters
{
    public SearchFeatureFileMap FeatureFileMap { get; set; }
    public override DeconvolutionType DeconvolutionType { get; protected set; } = DeconvolutionType.FromFile;

    public FeatureMappedFromFileDeconvolutionParameters()
        : this(new SearchFeatureFileMap())
    {
    }

    public FeatureMappedFromFileDeconvolutionParameters(SearchFeatureFileMap featureFileMap, int minCharge = 0, int maxCharge = 0, Polarity polarity = Polarity.Positive, AverageResidue averageResidueModel = null, double expectedIsotopeSpacing = Constants.C13MinusC12) : base(minCharge, maxCharge, polarity, averageResidueModel, expectedIsotopeSpacing)
    {
        FeatureFileMap = featureFileMap;
    }

    public DeconvolutionParameters ToDeconvolutionParameters(string massSpecFilePath)
    {
        // Fail loudly before any lookup: an empty embedded map means the task is misconfigured.
        FeatureFileMap.ValidateNotEmpty();

        if (FeatureFileMap.TryGetFeaturePathForMassSpecFile(massSpecFilePath, out var featureFilePath))
        {
            return new FromFileDeconvolutionParameters(featureFilePath, MinAssumedChargeState, MaxAssumedChargeState, Polarity)
            {
                AverageResidueModel = AverageResidueModel,
                ExpectedIsotopeSpacing = ExpectedIsotopeSpacing,
                UseGenericScore = UseGenericScore,
            };
        }
        else 
            throw new FeatureMappingException($"No feature file mapping found for mass spec file: {massSpecFilePath}");
    }

    public override DeconvolutionParameters Clone() => new FeatureMappedFromFileDeconvolutionParameters(FeatureFileMap.Clone(), MinAssumedChargeState, MaxAssumedChargeState, Polarity, AverageResidueModel, ExpectedIsotopeSpacing);

    public override DeconvolutionParameters ToDecoyParameters() => null;

    protected override void AddHashCodes(HashCode hash) => hash.Add(FeatureFileMap);

    protected override bool EqualProperties(DeconvolutionParameters other) => other is FeatureMappedFromFileDeconvolutionParameters otherParams && Equals(FeatureFileMap, otherParams.FeatureFileMap);
}
