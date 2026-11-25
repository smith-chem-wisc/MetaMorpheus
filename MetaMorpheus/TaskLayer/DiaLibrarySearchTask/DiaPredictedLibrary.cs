#nullable enable
using System;
using System.Collections.Generic;
using System.Linq;
using System.Net.Http;
using EngineLayer;
using MassSpectrometry;
using MzLibUtil;
using Omics.Modifications;
using Omics.SpectralMatch.MslSpectralLibrary;
using PredictionClients.Koina.AbstractClasses;
using PredictionClients.Koina.Client;
using PredictionClients.Koina.SupportedModels.FragmentIntensityModels;
using PredictionClients.Koina.SupportedModels.RetentionTimeModels;
using PredictionClients.SpectralLibraryGeneration;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using UsefulProteomicsDatabases.EntrapmentGeneration;

namespace TaskLayer;

/// <summary>
/// The library a DIA search predicts when it is given proteins and no library: Prosit 2025 spectra and iRT through Koina
/// (mzLib's <see cref="MslLibraryBuilder"/>), with decoys, written as <see cref="FileName"/> in the task's output folder and
/// then searched. The proteins are digested as the search's settings say, with its fixed and variable modifications.
/// </summary>
internal static class DiaPredictedLibrary
{
    /// <summary>The predicted library's file name in the task's output folder.</summary>
    public const string FileName = "DiaPredictedLibrary.msl";

    /// <summary>The longest peptide Prosit 2025 predicts; the library's digestion stops there.</summary>
    public const int MaxPeptideLength = 30;

    /// <summary>The precursor charges Prosit 2025 predicts.</summary>
    private const int MinCharge = 1, MaxCharge = 6;

    /// <summary>
    /// The Prosit 2025 models on Koina, at most <see cref="PostSearchAnalysisTask.MaxConcurrentPredictionRequests"/> requests
    /// at once: Koina is a public academic service, and one failed batch fails the whole build.
    /// </summary>
    public static (FragmentIntensityModel Intensity, RetentionTimeModel Irt) KoinaModels() =>
        (new Prosit2025Intensity40PTM(maxNumberOfBatchesPerRequest: PostSearchAnalysisTask.MaxConcurrentPredictionRequests),
         new Prosit2025iRT40PTM(maxNumberOfBatchesPerRequest: PostSearchAnalysisTask.MaxConcurrentPredictionRequests));

    /// <summary>
    /// Why a library cannot be predicted for this search, or null when it can: Prosit 2025 predicts peptides only, from a
    /// protease's digestion, as HCD or CID spectra, at charges 1 to 6 and a collision energy above 0.
    /// </summary>
    /// <param name="fileDissociationTypes">The fragmentation any spectra file sets for itself; one library serves them all.</param>
    public static string? Unsupported(CommonParameters commonParameters, IEnumerable<DissociationType> fileDissociationTypes,
        DiaLibrarySearchTaskParameters parameters)
    {
        if (GlobalVariables.AnalyteType != AnalyteType.Peptide)
            return $"Prosit predicts peptides, and this is a {GlobalVariables.AnalyteType} search";
        if (commonParameters.DigestionParams is not DigestionParams)
            return $"the library is digested with a protease, and this search digests with {commonParameters.DigestionParams.DigestionAgent.Name}";
        var unsupported = fileDissociationTypes.Prepend(commonParameters.DissociationType).Distinct()
            .Where(d => d != DissociationType.HCD && d != DissociationType.CID).ToList();
        if (unsupported.Count > 0)
            return $"Prosit 2025 predicts HCD and CID spectra, and this search uses {string.Join(", ", unsupported)}";
        if (parameters.LibraryPrecursorCharges.Count == 0 || parameters.LibraryPrecursorCharges.Any(z => z < MinCharge || z > MaxCharge))
            return $"Prosit 2025 predicts precursor charges {MinCharge} to {MaxCharge}, and LibraryPrecursorCharges is [{string.Join(", ", parameters.LibraryPrecursorCharges)}]";
        if (parameters.LibraryNce <= 0)
            return $"LibraryNce, the collision energy (NCE) the spectra are predicted at, must be above 0, and it is {parameters.LibraryNce}";
        return null;
    }

    /// <summary>
    /// The library's entries, targets and decoys, predicted from the proteins. Peptides longer than
    /// <see cref="MaxPeptideLength"/> are not digested; what Koina cannot predict is counted in <paramref name="result"/>.
    /// </summary>
    /// <exception cref="MetaMorpheusException">
    /// The permuter refuses the search's digestion (with the setting to use instead), or Koina did not answer (with what a
    /// prediction needs, and the way round it: an .msl library).
    /// </exception>
    public static List<MslLibraryEntry> Build(IReadOnlyList<Protein> proteins, CommonParameters commonParameters, DiaLibrarySearchTaskParameters parameters,
        List<Modification> fixedModifications, List<Modification> variableModifications, FragmentIntensityModel intensityModel, RetentionTimeModel irtModel,
        out MslLibraryBuildResult result)
    {
        var digestion = (DigestionParams)((DigestionParams)commonParameters.DigestionParams).Clone();
        digestion.MaxLength = Math.Min(digestion.MaxLength, MaxPeptideLength);
        var builder = new MslLibraryBuilder(intensityModel, irtModel)
        {
            DecoyProteins = parameters.LibraryDecoys == DiaLibraryDecoys.Permuted ? targets => PermutedDecoys(targets, digestion) : null,
        };
        var buildParameters = new MslLibraryBuildParameters(digestion, fixedModifications, variableModifications, parameters.LibraryPrecursorCharges,
            parameters.LibraryNce, commonParameters.DissociationType);
        try
        {
            return builder.Build(proteins, buildParameters, out result);
        }
        catch (MzLibException e) when (parameters.LibraryDecoys == DiaLibraryDecoys.Permuted)
        {
            throw new MetaMorpheusException($"Permuted decoys cannot be made for this search's digestion: {e.Message} " +
                "Set LibraryDecoys = \"Reversed\" to predict the library with reversed decoys.", e);
        }
        catch (Exception e) when (e is HttpRequestException or OperationCanceledException)
        {
            // Nothing in MetaMorpheus cancels the build, so a cancellation is Koina's own deadline running out
            throw new MetaMorpheusException($"The library could not be predicted, because Koina did not answer ({e.GetType().Name}: {e.Message}). " +
                $"Predicting a library needs network access to Koina ({HTTP.ModelsURL}); without it, give the search an .msl library.", e);
        }
    }

    /// <summary>
    /// Each target's permuted decoy, the benchmark libraries' recipe: mzLib's entrapment generator, fold 1 of 3, seed 1, named
    /// DECOY_&lt;accession&gt;_f1, with no peptide equal to a target's. The generator makes entrapment, so each is marked a
    /// decoy instead, and it refuses a contaminant, so a contaminant goes in as a plain copy.
    /// </summary>
    private static IReadOnlyList<Protein> PermutedDecoys(IReadOnlyList<Protein> targets, DigestionParams digestion)
    {
        var forbidden = EntrapmentProteinGenerator.TargetPeptides(targets, digestion);
        return targets.Select(target => new Protein(EntrapmentProteinGenerator.Create(target.IsContaminant ? new Protein(target, isContaminant: false) : target,
            digestion, forbidden, fold: 1, foldCount: 3, seed: 1, entrapmentIdentifier: "DECOY"), isDecoy: true, isEntrapment: false)).ToList();
    }

    /// <summary>The results.txt line for a predicted library: what it holds, how it was made, and what could not be predicted.</summary>
    public static string Describe(MslLibraryBuildResult result, int proteinCount, DiaLibrarySearchTaskParameters parameters, FragmentIntensityModel intensityModel,
        RetentionTimeModel irtModel) =>
        $"Predicted library {FileName}: {result.TargetPrecursors} target and {result.DecoyPrecursors} decoy precursors from {proteinCount} proteins, " +
        $"{parameters.LibraryDecoys.ToString().ToLowerInvariant()} decoys; {intensityModel.ModelName} and {irtModel.ModelName} at {HTTP.ModelsURL}, " +
        $"NCE {parameters.LibraryNce}, charges {string.Join(", ", parameters.LibraryPrecursorCharges)}; {result.NotPredicted} precursors not predicted, " +
        $"{result.NoFragments} with no fragments, {result.DecoysDroppedAsTargetSequences} decoys dropped as target sequences";
}
