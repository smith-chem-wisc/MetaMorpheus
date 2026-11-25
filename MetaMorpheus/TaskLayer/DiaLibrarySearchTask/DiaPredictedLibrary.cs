#nullable enable
using System;
using System.Collections.Generic;
using EngineLayer;
using MassSpectrometry;
using Omics.Modifications;
using Omics.SpectralMatch.MslSpectralLibrary;
using PredictionClients.Koina.AbstractClasses;
using PredictionClients.SpectralLibraryGeneration;
using Proteomics;

namespace TaskLayer;

/// <summary>
/// The library a DIA search predicts when it is given proteins and no library: Prosit 2025 spectra and iRT through Koina,
/// with decoys, written as <see cref="FileName"/> in the task's output folder and then searched.
/// </summary>
internal static class DiaPredictedLibrary
{
    /// <summary>The predicted library's file name in the task's output folder.</summary>
    public const string FileName = "DiaPredictedLibrary.msl";

    /// <summary>The Prosit 2025 models on Koina, at most <see cref="PostSearchAnalysisTask.MaxConcurrentPredictionRequests"/> requests at once.</summary>
    public static (FragmentIntensityModel Intensity, RetentionTimeModel Irt) KoinaModels() => throw new NotImplementedException();

    /// <summary>Why a library cannot be predicted for this search, or null when it can.</summary>
    public static string? Unsupported(CommonParameters commonParameters, IEnumerable<DissociationType> fileDissociationTypes,
        DiaLibrarySearchTaskParameters parameters) => throw new NotImplementedException();

    /// <summary>The library's entries, targets and decoys, predicted from the proteins.</summary>
    public static List<MslLibraryEntry> Build(IReadOnlyList<Protein> proteins, CommonParameters commonParameters, DiaLibrarySearchTaskParameters parameters,
        List<Modification> fixedModifications, List<Modification> variableModifications, FragmentIntensityModel intensityModel, RetentionTimeModel irtModel,
        out MslLibraryBuildResult result) => throw new NotImplementedException();
}
