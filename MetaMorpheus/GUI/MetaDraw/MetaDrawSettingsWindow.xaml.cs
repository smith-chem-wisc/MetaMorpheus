﻿using EngineLayer;
using EngineLayer.GlycoSearch;
using GuiFunctions;
using Nett;
using OxyPlot;
using Omics.Fragmentation;
using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.ComponentModel;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Windows;
using System.Windows.Controls;

namespace MetaMorpheusGUI
{
    /// <summary>
    /// Interaction logic for MetaDrawSettingsWindow.xaml
    /// </summary>
    public partial class MetaDrawSettingsWindow : Window
    {
        private readonly ObservableCollection<ModTypeForTreeViewModel> Modifications = new ObservableCollection<ModTypeForTreeViewModel>();
        private readonly ObservableCollection<IonTypeForTreeViewModel> IonGroups = new ObservableCollection<IonTypeForTreeViewModel>();
        private readonly ObservableCollection<CoverageTypeForTreeViewModel> CoverageColors = new ObservableCollection<CoverageTypeForTreeViewModel>();

        private MetaDrawSettingsViewModel SettingsView;

        public MetaDrawSettingsWindow(MetaDrawSettingsViewModel view)
        {
            InitializeComponent();
            SettingsView = view;
            DataContext = SettingsView;
            PopulateChoices();
        }

        private void PopulateChoices()
        {
            foreach (string level in System.Enum.GetNames(typeof(LocalizationLevel)))
            {
                CmbGlycanLocalizationLevelStart.Items.Add(level);
                CmbGlycanLocalizationLevelEnd.Items.Add(level);
            }

            DisplayAnnotationsCheckBox.IsChecked = MetaDrawSettings.DisplayIonAnnotations;
            MZCheckBox.IsChecked = MetaDrawSettings.AnnotateMzValues;
            ChargesCheckBox.IsChecked = MetaDrawSettings.AnnotateCharges;
            BoldTextCheckBox.IsChecked = MetaDrawSettings.AnnotationBold;
            DecoysCheckBox.IsChecked = MetaDrawSettings.ShowDecoys;
            ContaminantsCheckBox.IsChecked = MetaDrawSettings.ShowContaminants;
            ShowInternalIonAnnotationsCheckBox.IsChecked = MetaDrawSettings.DisplayInternalIonAnnotations;
            PrecursorChargeCheckBox.IsChecked = MetaDrawSettings.SpectrumDescription["Precursor Charge: "];
            ShowInternalIonsCheckBox.IsChecked = MetaDrawSettings.DisplayInternalIons;
            SubSuperScriptedIonAnnotations.IsChecked = MetaDrawSettings.SubAndSuperScriptIons;
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
            TextSizeBox.Text = MetaDrawSettings.AnnotatedFontSize.ToString();
            CmbGlycanLocalizationLevelStart.SelectedItem = MetaDrawSettings.LocalizationLevelStart.ToString();
            CmbGlycanLocalizationLevelEnd.SelectedItem = MetaDrawSettings.LocalizationLevelEnd.ToString();

            ExportFileFormatComboBox.ItemsSource = MetaDrawSettings.ExportTypes;
            ExportFileFormatComboBox.SelectedItem = MetaDrawSettings.ExportType;
            IonColorExpander.ItemsSource = SettingsView.IonGroups;
            PTMColorExpander.ItemsSource = SettingsView.Modifications;
            SequenceCoverageColorExpander.ItemsSource = SettingsView.CoverageColors;
        }

        private void Cancel_Click(object sender, RoutedEventArgs e)
        {
            DialogResult = false;
        }

        private void Save_Click(object sender, RoutedEventArgs e)
        {
            MetaDrawSettings.DisplayIonAnnotations = DisplayAnnotationsCheckBox.IsChecked.Value;
            MetaDrawSettings.AnnotateMzValues = MZCheckBox.IsChecked.Value;
            MetaDrawSettings.AnnotateCharges = ChargesCheckBox.IsChecked.Value;
            MetaDrawSettings.AnnotationBold = BoldTextCheckBox.IsChecked.Value;
            MetaDrawSettings.ShowDecoys = DecoysCheckBox.IsChecked.Value;
            MetaDrawSettings.ShowContaminants = ContaminantsCheckBox.IsChecked.Value;
            MetaDrawSettings.DisplayInternalIons = ShowInternalIonsCheckBox.IsChecked.Value;
            MetaDrawSettings.SubAndSuperScriptIons = SubSuperScriptedIonAnnotations.IsChecked.Value;
            MetaDrawSettings.DisplayInternalIonAnnotations = ShowInternalIonAnnotationsCheckBox.IsChecked.Value;
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
            if (!ShowInternalIonsCheckBox.IsChecked.Value)
                MetaDrawSettings.InternalIonColor = OxyColors.Transparent;
            SettingsView.Save();

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

            if (!string.IsNullOrWhiteSpace(TextSizeBox.Text))
            {
                if (int.TryParse(TextSizeBox.Text, out int fontSize))
                {
                    if (fontSize > 18)
                    {
                        MessageBox.Show("Font size must be <= 18");
                        return;
                    }

                    MetaDrawSettings.AnnotatedFontSize = fontSize;
                }
                else
                {
                    MessageBox.Show("Could not parse font size; must be a positive integer");
                    return;
                }
            }
            else
            {
                MetaDrawSettings.AnnotatedFontSize = 12;
            }

            DialogResult = true;
        }

        private void setDefaultbutton_Click(object sender, RoutedEventArgs e)
        {
            Save_Click(sender, e);
            SettingsView.SaveAsDefault();
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
                MetaDrawSettingsViewModel settingsViewModel = new MetaDrawSettingsViewModel();
                SettingsView = settingsViewModel;
                DataContext = SettingsView;
                PopulateChoices();
                DialogResult = true;
            }
        }
    }
}
