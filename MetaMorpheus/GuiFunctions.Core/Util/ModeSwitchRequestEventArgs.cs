using System;

namespace GuiFunctions.Util;

/// <summary>
/// Event arguments for requesting mode switch confirmation from the UI
/// </summary>
public class ModeSwitchRequestEventArgs : EventArgs
{
    public ModeSwitchResult Result { get; set; } = ModeSwitchResult.Cancel;
    public bool RememberMyDecision { get; set; } = false;
}

/// <summary>
/// Represents the user's choice when switching modes
/// </summary>
public enum ModeSwitchResult
{
    /// <summary>
    /// User chose to cancel the mode switch
    /// </summary>
    Cancel,

    /// <summary>
    /// User chose to switch modes and keep all loaded files
    /// </summary>
    SwitchKeepFiles,

    /// <summary>
    /// User chose to switch modes and remove all loaded files
    /// </summary>
    SwitchRemoveFiles
}