﻿using Nett;

namespace GuiFunctions;

/// <summary>
/// Gui Parameters that are written out to a toml file. 
/// This gets parsed in and if an exception is thrown, a new default instance is created and the file is overwritten.
/// </summary>
public class GuiGlobalParams
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
}