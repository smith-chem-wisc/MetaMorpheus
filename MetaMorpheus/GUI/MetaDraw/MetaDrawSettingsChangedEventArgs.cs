using System;

namespace MetaMorpheusGUI;

public class MetaDrawSettingsChangedEventArgs : EventArgs
{
    public bool FilterChanged { get; set; } = false;
    public bool DataVisualizationChanged { get; set; } = false;
}