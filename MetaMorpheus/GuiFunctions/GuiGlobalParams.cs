using System;
using Nett;

namespace GuiFunctions;

/// <summary>
/// Gui Parameters that are written out to a toml file. 
/// This gets parsed in and if an exception is thrown, a new default instance is created and the file is overwritten.
/// </summary>
public class GuiGlobalParams  : IEquatable<GuiGlobalParams>
{
    public bool AskAboutUpdating { get; internal set; } = true;
    public bool AskBeforeExitingMetaMorpheus { get; internal set; } = true;

    // User can set a custom proteome directory. Be sure to use double slashes in the path, otherwise it will not be read in properly. 
    [TomlMember(Key = "UserSpecifiedProteomeDir")]
    public string ProteomeDirectory { get; internal set; }
    public string DecoyIdentifier { get; internal set; } = "DECOY";

    //Ask about protease-specific parameter recommendations
    public bool AskAboutTopDownParams { get; internal set; } = true;
    public bool AskAboutChymotrypsinParams { get; internal set; } = true;
    public bool AskAboutElastaseParams { get; internal set; } = true;
    public bool AskAboutNonSpecificParams { get; internal set; } = true;
    public bool AskAboutSemiTrypsinParams { get; internal set; } = true;
    public bool AskAboutArgCParams { get; internal set; } = true;
    public bool AskAboutSpectralRecoveryParams { get; internal set; } = true;

    //Use protease-specific parameter recommendations
    public bool UseTopDownParams { get; internal set; } = true;
    public bool UseChymotrypsinParams { get; internal set; } = true;
    public bool UseElastaseParams { get; internal set; } = true;
    public bool UseNonSpecificParams { get; internal set; } = true;
    public bool UseSemiTrypsinParams { get; internal set; } = true;
    public bool UseArgCParams { get; internal set; } = true;
    public bool UseSpectralRecoveryParams { get; internal set; } = true;

    // Rna Toggles
    public bool IsRnaMode { get; internal set; }

    // Deep equality check (can be improved for more complex types)
    public bool Equals(GuiGlobalParams obj)
    {
        if (obj is not GuiGlobalParams other)
            return false;

        return
            AskAboutUpdating == other.AskAboutUpdating &&
            AskBeforeExitingMetaMorpheus == other.AskBeforeExitingMetaMorpheus &&
            ProteomeDirectory == other.ProteomeDirectory &&
            AskAboutTopDownParams == other.AskAboutTopDownParams &&
            AskAboutChymotrypsinParams == other.AskAboutChymotrypsinParams &&
            AskAboutElastaseParams == other.AskAboutElastaseParams &&
            AskAboutNonSpecificParams == other.AskAboutNonSpecificParams &&
            AskAboutSemiTrypsinParams == other.AskAboutSemiTrypsinParams &&
            AskAboutArgCParams == other.AskAboutArgCParams &&
            AskAboutSpectralRecoveryParams == other.AskAboutSpectralRecoveryParams &&
            UseTopDownParams == other.UseTopDownParams &&
            UseChymotrypsinParams == other.UseChymotrypsinParams &&
            UseElastaseParams == other.UseElastaseParams &&
            UseNonSpecificParams == other.UseNonSpecificParams &&
            UseSemiTrypsinParams == other.UseSemiTrypsinParams &&
            UseArgCParams == other.UseArgCParams &&
            UseSpectralRecoveryParams == other.UseSpectralRecoveryParams
            && IsRnaMode == other.IsRnaMode;
    }

    // Helper for deep copy
    internal GuiGlobalParams Clone()
    {
        return new GuiGlobalParams
        {
            AskAboutUpdating = AskAboutUpdating,
            AskBeforeExitingMetaMorpheus = AskBeforeExitingMetaMorpheus,
            ProteomeDirectory = ProteomeDirectory,
            AskAboutTopDownParams = AskAboutTopDownParams,
            AskAboutChymotrypsinParams = AskAboutChymotrypsinParams,
            AskAboutElastaseParams = AskAboutElastaseParams,
            AskAboutNonSpecificParams = AskAboutNonSpecificParams,
            AskAboutSemiTrypsinParams = AskAboutSemiTrypsinParams,
            AskAboutArgCParams = AskAboutArgCParams,
            AskAboutSpectralRecoveryParams = AskAboutSpectralRecoveryParams,
            UseTopDownParams = UseTopDownParams,
            UseChymotrypsinParams = UseChymotrypsinParams,
            UseElastaseParams = UseElastaseParams,
            UseNonSpecificParams = UseNonSpecificParams,
            UseSemiTrypsinParams = UseSemiTrypsinParams,
            UseArgCParams = UseArgCParams,
            UseSpectralRecoveryParams = UseSpectralRecoveryParams,
            IsRnaMode = IsRnaMode
        };
    }
}