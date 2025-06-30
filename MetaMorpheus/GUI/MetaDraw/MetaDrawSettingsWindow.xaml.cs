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
            foreach (string level in System.Enum.GetNames(typeof(LocalizationLevel)))
            {
                CmbGlycanLocalizationLevelStart.Items.Add(level);
                CmbGlycanLocalizationLevelEnd.Items.Add(level);
            }

            DecoysCheckBox.IsChecked = MetaDrawSettings.ShowDecoys;
            ContaminantsCheckBox.IsChecked = MetaDrawSettings.ShowContaminants;
            PrecursorChargeCheckBox.IsChecked = MetaDrawSettings.SpectrumDescription["Precursor Charge: "];
            PrecursorMassCheckBox.IsChecked = MetaDrawSettings.SpectrumDescription["Precursor Mass: "];
            TheoreticalMassCheckBox.IsChecked = MetaDrawSettings.SpectrumDescription["Theoretical Mass: "];
            ScoreCheckBox.IsChecked = MetaDrawSettings.SpectrumDescription["Score: "];
            ProteinAccessionCheckBox.IsChecked = MetaDrawSettings.SpectrumDescription["Protein Accession: "];
            ProteinCheckBox.IsChecked = MetaDrawSettings.SpectrumDescription["Protein: "];
            DecoyContaminantTargetCheckBox.IsChecked = MetaDrawSettings.SpectrumDescription["Decoy/Contaminant/Target: "];
            QValueCheckBox.IsChecked = MetaDrawSettings.SpectrumDescription["Q-Value: "];
            SequenceLengthCheckBox.IsChecked = MetaDrawSettings.SpectrumDescription["Sequence Length: "];
            AmbiguityLevelCheckBox.IsChecked = MetaDrawSettings.SpectrumDescription["Ambiguity Level: "];
            SpectralAngleCheckBox.IsChecked = MetaDrawSettings.SpectrumDescription["Spectral Angle: "];
            PEPCheckBox.IsChecked = MetaDrawSettings.SpectrumDescription["PEP: "];
            PEPQValueCheckBox.IsChecked = MetaDrawSettings.SpectrumDescription["PEP Q-Value: "];
            StationarySequenceCheckBox.IsChecked = MetaDrawSettings.DrawStationarySequence;
            SequencenNumbersCheckBox.IsChecked = MetaDrawSettings.DrawNumbersUnderStationary;
            ShowLegendCheckBox.IsChecked = MetaDrawSettings.ShowLegend;
            qValueBox.Text = MetaDrawSettings.QValueFilter.ToString();
            AmbiguityFilteringComboBox.DataContext = MetaDrawSettings.AmbiguityTypes;
            AmbiguityFilteringComboBox.SelectedItem = MetaDrawSettings.AmbiguityFilter;
            CmbGlycanLocalizationLevelStart.SelectedItem = MetaDrawSettings.LocalizationLevelStart.ToString();
            CmbGlycanLocalizationLevelEnd.SelectedItem = MetaDrawSettings.LocalizationLevelEnd.ToString();
            SpectrumDescriptionFontSizeBox.Text = MetaDrawSettings.SpectrumDescriptionFontSize.ToString();

            ExportFileFormatComboBox.ItemsSource = MetaDrawSettings.ExportTypes;
            ExportFileFormatComboBox.SelectedItem = MetaDrawSettings.ExportType;
            IonColorExpander.ItemsSource = MetaDrawSettingsViewModel.Instance.IonGroups;
            PTMColorExpander.ItemsSource = MetaDrawSettingsViewModel.Instance.Modifications;
            SequenceCoverageColorExpander.ItemsSource = MetaDrawSettingsViewModel.Instance.CoverageColors;
        }

        private void Cancel_Click(object sender, RoutedEventArgs e)
        {
            DialogResult = false;
        }

        private void Save_Click(object sender, RoutedEventArgs e)
        {
            MetaDrawSettings.ShowDecoys = DecoysCheckBox.IsChecked.Value;
            MetaDrawSettings.ShowContaminants = ContaminantsCheckBox.IsChecked.Value;
            MetaDrawSettings.SpectrumDescription["Precursor Charge: "] = PrecursorChargeCheckBox.IsChecked.Value;
            MetaDrawSettings.SpectrumDescription["Precursor Mass: "] = PrecursorMassCheckBox.IsChecked.Value;
            MetaDrawSettings.SpectrumDescription["Theoretical Mass: "] = TheoreticalMassCheckBox.IsChecked.Value;
            MetaDrawSettings.SpectrumDescription["Score: "] = ScoreCheckBox.IsChecked.Value;
            MetaDrawSettings.SpectrumDescription["Protein Accession: "] = ProteinAccessionCheckBox.IsChecked.Value;
            MetaDrawSettings.SpectrumDescription["Protein: "] = ProteinCheckBox.IsChecked.Value;
            MetaDrawSettings.SpectrumDescription["Decoy/Contaminant/Target: "] = DecoyContaminantTargetCheckBox.IsChecked.Value;
            MetaDrawSettings.SpectrumDescription["Q-Value: "] = QValueCheckBox.IsChecked.Value;
            MetaDrawSettings.SpectrumDescription["Sequence Length: "] = SequenceLengthCheckBox.IsChecked.Value;
            MetaDrawSettings.SpectrumDescription["Ambiguity Level: "] = AmbiguityLevelCheckBox.IsChecked.Value;
            MetaDrawSettings.SpectrumDescription["Spectral Angle: "] = SpectralAngleCheckBox.IsChecked.Value;
            MetaDrawSettings.SpectrumDescription["PEP: "] = PEPCheckBox.IsChecked.Value;
            MetaDrawSettings.SpectrumDescription["PEP Q-Value: "] = PEPQValueCheckBox.IsChecked.Value;
            MetaDrawSettings.DrawStationarySequence = StationarySequenceCheckBox.IsChecked.Value;
            MetaDrawSettings.DrawNumbersUnderStationary = SequencenNumbersCheckBox.IsChecked.Value;
            MetaDrawSettings.ShowLegend = ShowLegendCheckBox.IsChecked.Value;
            MetaDrawSettings.LocalizationLevelStart = (LocalizationLevel)System.Enum.Parse(typeof(LocalizationLevel), CmbGlycanLocalizationLevelStart.SelectedItem.ToString());
            MetaDrawSettings.LocalizationLevelEnd = (LocalizationLevel)System.Enum.Parse(typeof(LocalizationLevel), CmbGlycanLocalizationLevelEnd.SelectedItem.ToString());
            MetaDrawSettings.ExportType = ExportFileFormatComboBox.SelectedItem.ToString();
            MetaDrawSettings.AmbiguityFilter = AmbiguityFilteringComboBox.SelectedItem.ToString();
            MetaDrawSettings.SpectrumDescriptionFontSize = double.TryParse(SpectrumDescriptionFontSizeBox.Text, out double spectrumDescriptionFontSize) ? spectrumDescriptionFontSize : 10;
            if (!ShowInternalIonsCheckBox.IsChecked.Value)
                MetaDrawSettings.InternalIonColor = OxyColors.Transparent;
            MetaDrawSettingsViewModel.Instance.Save();

            if (!string.IsNullOrWhiteSpace(qValueBox.Text))
            {
                if (double.TryParse(qValueBox.Text, NumberStyles.Any, CultureInfo.InvariantCulture, out double qValueFilter) && qValueFilter >= 0 && qValueFilter <= 1)
                {
                    MetaDrawSettings.QValueFilter = qValueFilter;
                }
                else
                {
                    MessageBox.Show("Could not parse q-value filter; must be number between 0 and 1 inclusive");
                    return;
                }
            }
            else
            {
                MetaDrawSettings.QValueFilter = 1;
            }

            var fontSize = MetaDrawSettings.AnnotatedFontSize;
            if (fontSize <= 0)
            {
                MessageBox.Show("Font size must be a positive integer");
                return;
            }
            if (fontSize > 36)
            {
                MessageBox.Show("Font size must be <= 36");
                return;
            }

            fontSize = MetaDrawSettings.AnnotatedFontSize;
            if (fontSize <= 0)
            {
                MessageBox.Show("Font size must be a positive integer");
                return;
            }
            if (fontSize > 36)
            {
                MessageBox.Show("Font size must be <= 36");
                return;
            }

            fontSize = MetaDrawSettings.AxisLabelTextSize;
            if (fontSize <= 0)
            {
                MessageBox.Show("Font size must be a positive integer");
                return;
            }
            if (fontSize > 36)
            {
                MessageBox.Show("Font size must be <= 36");
                return;
            }

            var strokeThickness = MetaDrawSettings.StrokeThicknessAnnotated;
            if (strokeThickness <= 0)
            {
                MessageBox.Show("Stroke thickness must be a positive number");
                return;
            }

            strokeThickness = MetaDrawSettings.StrokeThicknessUnannotated;
            if (strokeThickness <= 0)
            {
                MessageBox.Show("Stroke thickness must be a positive number");
                return;
            }

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
            ((CoverageTypeForTreeViewModel)((ComboBox)sender).DataContext).SelectionChanged((string)((ComboBox)sender).SelectedItem);
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
