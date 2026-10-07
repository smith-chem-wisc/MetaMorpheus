using System;
using System.Windows;
using System.Windows.Controls;
using GuiFunctions;
using Microsoft.Win32;

namespace MetaMorpheusGUI
{
    public partial class FromFileDeconParamsControl : UserControl
    {
        public FromFileDeconParamsControl()
        {
            InitializeComponent();
        }

        public Func<string?>? FilePickerOverride { get; set; }

        public string? ValidationMessage { get; private set; }

        private void BrowseButton_Click(object sender, RoutedEventArgs e) => BrowseForFile();

        public bool BrowseForFile()
        {
            string? selected = FilePickerOverride != null
                ? FilePickerOverride()
                : ShowNativeFileDialog();

            if (selected == null)
                return false;

            if (FromFileDeconParamsViewModel.TryValidateFilePath(selected, out var message))
            {
                if (DataContext is FromFileDeconParamsViewModel vm)
                    vm.FilePath = selected;

                ValidationTextBlock.Text = string.Empty;
                ValidationTextBlock.Visibility = Visibility.Collapsed;
                ValidationMessage = null;
                return true;
            }

            ValidationTextBlock.Text = message;
            ValidationTextBlock.Visibility = Visibility.Visible;
            ValidationMessage = message;
            return false;
        }

        private static string? ShowNativeFileDialog()
        {
            var dialog = new OpenFileDialog
            {
                Title = "Select MS1 Feature File",
                Filter = "MS1 Feature Files (*_ms1.feature;*.feature.tsv)|*_ms1.feature;*.feature.tsv|All Files (*.*)|*.*",
                CheckFileExists = true
            };
            return dialog.ShowDialog() == true ? dialog.FileName : null;
        }
    }
}
