using GuiFunctions;
using OxyPlot;
using System.Globalization;
using System.IO;
using System.Windows;
using System.Windows.Controls;
using Readers;

namespace MetaMorpheusGUI
{
    /// <summary>
    /// Interaction logic for MetaDrawSettingsWindow.xaml
    /// </summary>
    public partial class MetaDrawSettingsWindow : Window
    {
        public MetaDrawSettingsWindow()
        {
            InitializeComponent();
            DataContext = MetaDrawSettingsViewModel.Instance;
            PopulateChoices();
        }

        private void PopulateChoices()
        {
            StationarySequenceCheckBox.IsChecked = MetaDrawSettings.DrawStationarySequence;
            SequencenNumbersCheckBox.IsChecked = MetaDrawSettings.DrawNumbersUnderStationary;
            ShowLegendCheckBox.IsChecked = MetaDrawSettings.ShowLegend;
            SpectrumDescriptionFontSizeBox.Text = MetaDrawSettings.SpectrumDescriptionFontSize.ToString();
        }

        private void Cancel_Click(object sender, RoutedEventArgs e)
        {
            DialogResult = false;
        }

        private void Save_Click(object sender, RoutedEventArgs e)
        {
            MetaDrawSettings.DrawStationarySequence = StationarySequenceCheckBox.IsChecked.Value;
            MetaDrawSettings.DrawNumbersUnderStationary = SequencenNumbersCheckBox.IsChecked.Value;
            MetaDrawSettings.ShowLegend = ShowLegendCheckBox.IsChecked.Value;
            MetaDrawSettings.SpectrumDescriptionFontSize = double.TryParse(SpectrumDescriptionFontSizeBox.Text, out double spectrumDescriptionFontSize) ? spectrumDescriptionFontSize : 10;
            if (!ShowInternalIonsCheckBox.IsChecked.Value)
                MetaDrawSettings.InternalIonColor = OxyColors.Transparent;

            switch (MetaDrawSettings.AnnotatedFontSize)
            {
                case <= 0:
                    MessageBox.Show("Font size must be a positive integer");
                    return;
                case > 36:
                    MessageBox.Show("Font size must be <= 36");
                    return;
            }

            switch (MetaDrawSettings.AxisTitleTextSize)
            {
                case <= 0:
                    MessageBox.Show("Font size must be a positive integer");
                    return;
                case > 36:
                    MessageBox.Show("Font size must be <= 36");
                    return;
            }

            switch (MetaDrawSettings.AxisLabelTextSize)
            {
                case <= 0:
                    MessageBox.Show("Font size must be a positive integer");
                    return;
                case > 36:
                    MessageBox.Show("Font size must be <= 36");
                    return;
            }

            if (MetaDrawSettings.StrokeThicknessAnnotated <= 0)
            {
                MessageBox.Show("Stroke thickness must be a positive number");
                return;
            }

            if (MetaDrawSettings.StrokeThicknessUnannotated <= 0)
            {
                MessageBox.Show("Stroke thickness must be a positive number");
                return;
            }

            MetaDrawSettingsViewModel.Instance.Save();
            DialogResult = true;
        }

        private void setDefaultbutton_Click(object sender, RoutedEventArgs e)
        {
            Save_Click(sender, e);
            MetaDrawSettingsViewModel.Instance.SaveAsDefault();
        }

        /// <summary>
        /// Event handler for the ion type color selection. Runs the SelectionChanged command in its respective view model
        /// </summary>
        /// <param name="sender">ComboBox which was changed</param>
        /// <param name="e"></param>
        private void ComboBox_SelectionChanged(object sender, System.Windows.Controls.SelectionChangedEventArgs e)
        {
            ((IonForTreeViewModel)((ComboBox)sender).DataContext).SelectionChanged((string)((ComboBox)sender).SelectedItem);
        }

        /// <summary>
        /// Event handler for the ptm type color selection. Runs the SelectionChanged command in its respective view model
        /// </summary>
        /// <param name="sender">ComboBox which was changed</param>
        /// <param name="e"></param>
        private void ComboBox_SelectionChanged_1(object sender, SelectionChangedEventArgs e)
        {
            ((ModForTreeViewModel)((ComboBox)sender).DataContext).SelectionChanged((string)((ComboBox)sender).SelectedItem);
        }

        /// <summary>
        /// Event handler for the sequence coverage type color selection. Runs the SelectionChanged command in its respective view model
        /// </summary>
        /// <param name="sender">ComboBox which was changed</param>
        /// <param name="e"></param>
        private void ComboBox_SelectionChanged_2(object sender, SelectionChangedEventArgs e)
        {
            ((ColorForTreeViewModel)((ComboBox)sender).DataContext).SelectionChanged((string)((ComboBox)sender).SelectedItem);
        }

        /// <summary>
        /// Event handler for when the button to restore default settings is clicked
        /// </summary>
        /// <param name="sender"></param>
        /// <param name="e"></param>
        private void RestoreDefaultButton_Click(object sender, RoutedEventArgs e)
        {
            if (MessageBox.Show("Reset to default values?", "", MessageBoxButton.YesNo, MessageBoxImage.Warning) == MessageBoxResult.Yes)
            {
                if (File.Exists(MetaDrawSettingsViewModel.SettingsPath))
                    File.Delete(MetaDrawSettingsViewModel.SettingsPath);
                MetaDrawSettings.ResetSettings();
                MetaDrawSettingsViewModel settingsViewModel = MetaDrawSettingsViewModel.Instance;
                settingsViewModel.LoadSettings();
                DataContext = settingsViewModel;
                PopulateChoices();
                DialogResult = true;
            }
        }
    }
}
