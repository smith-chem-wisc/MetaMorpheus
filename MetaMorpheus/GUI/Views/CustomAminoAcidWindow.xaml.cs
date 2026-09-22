using GuiFunctions;
using System;
using System.ComponentModel;
using System.Windows;

namespace MetaMorpheusGUI
{
    /// <summary>
    /// Interaction logic for CustomAminoAcidWindow.xaml
    /// </summary>
    public partial class CustomAminoAcidWindow : Window
    {
        private readonly CustomResidueViewModel _viewModel;

        public CustomAminoAcidWindow()
        {
            InitializeComponent();
            _viewModel = new CustomResidueViewModel();
            _viewModel.RequestClose += ViewModel_RequestClose;
            DataContext = _viewModel;
            Closing += CustomAminoAcidWindow_Closing;
        }

        private void ViewModel_RequestClose(object sender, CustomResidueDialogResultEventArgs e)
        {
            DialogResult = e.Succeeded;
        }

        private void CustomAminoAcidWindow_Closing(object sender, CancelEventArgs e)
        {
            _viewModel.RequestClose -= ViewModel_RequestClose;
            _viewModel.Dispose();
        }
    }
}
