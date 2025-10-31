﻿using System.Windows.Controls;
using System.Windows.Input;
using GuiFunctions;

namespace MetaMorpheusGUI;

/// <summary>
/// Interaction logic for ProteinRnaImageToggle.xaml
/// </summary>
public partial class ProteinRnaImageToggle : UserControl
{
    public ProteinRnaImageToggle() => InitializeComponent();

    private void Protein_Click(object sender, MouseButtonEventArgs e) => GuiGlobalParamsViewModel.Instance.IsRnaMode = false;

    private void Rna_Click(object sender, MouseButtonEventArgs e) => GuiGlobalParamsViewModel.Instance.IsRnaMode = true;
}
