using EngineLayer;
using EngineLayer.GlycoSearch;
using GuiFunctions;
using Nett;
using OxyPlot;
using Proteomics.Fragmentation;
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
        private readonly ObservableCollection<ModTypeForTreeView> Modifications = new ObservableCollection<ModTypeForTreeView>();
        private readonly ObservableCollection<IonTypeForTreeViewModel> IonGroups = new ObservableCollection<IonTypeForTreeViewModel>();

        public MetaDrawSettingsWindow(ObservableCollection<ModTypeForTreeView> mods)
        {
            InitializeComponent();
            Modifications = mods;
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
            PrecursorChargeCheckBox.IsChecked = MetaDrawSettings.SpectrumDescription["Precursor Charge: "];
            PrecursorMassCheckBox.IsChecked = MetaDrawSettings.SpectrumDescription["Precursor Mass: "];
            TheoreticalMassCheckBox.IsChecked = MetaDrawSettings.SpectrumDescription["Theoretical Mass: "];
            ScoreCheckBox.IsChecked = MetaDrawSettings.SpectrumDescription["Score: "];
            ProteinAccessionCheckBox.IsChecked = MetaDrawSettings.SpectrumDescription["Protein Accession: "];
            ProteinCheckBox.IsChecked = MetaDrawSettings.SpectrumDescription["Protein: "];
            DecoyContaminantTargetCheckBox.IsChecked = MetaDrawSettings.SpectrumDescription["Decoy/Contaminant/Target: "];
            QValueCheckBox.IsChecked = MetaDrawSettings.SpectrumDescription["Q-Value: "];
            SequenceLengthCheckBox.IsChecked = MetaDrawSettings.SpectrumDescription["Sequence Length: "];
            ProFormaLevelCheckBox.IsChecked = MetaDrawSettings.SpectrumDescription["ProForma Level: "];
            PEPCheckBox.IsChecked = MetaDrawSettings.SpectrumDescription["PEP: "];
            PEPQValueCheckBox.IsChecked = MetaDrawSettings.SpectrumDescription["PEP Q-Value: "];
            qValueBox.Text = MetaDrawSettings.QValueFilter.ToString();
            TextSizeBox.Text = MetaDrawSettings.AnnotatedFontSize.ToString();
            CmbGlycanLocalizationLevelStart.SelectedItem = MetaDrawSettings.LocalizationLevelStart.ToString();
            CmbGlycanLocalizationLevelEnd.SelectedItem = MetaDrawSettings.LocalizationLevelEnd.ToString();


            var ions = ((ProductType[])Enum.GetValues(typeof(ProductType)));
            var common = ions.Where(p => p.ToString().Equals("a") || p.ToString().Equals("b") || p.ToString().Equals("c")
                                          || p.ToString().Equals("x") || p.ToString().Equals("y") || p.ToString().Equals("zDot"));
            var lessCommon = ions.Where(p => !p.ToString().Equals("a") || !p.ToString().Equals("b") || !p.ToString().Equals("c")
                                          || !p.ToString().Equals("x") || !p.ToString().Equals("y") || !p.ToString().Equals("zDot"));
            ObservableCollection<string> colors = new ObservableCollection<string>(MetaDrawSettings.PossibleColors.Values.ToList());
            IonGroups.Add(new IonTypeForTreeViewModel("Common Ions", common, false, colors));
            IonGroups.Add(new IonTypeForTreeViewModel("Less Common Ions", lessCommon, false, colors));
            IonGroups.Add(new IonTypeForTreeViewModel("Cross Linked Beta Peptide", ions, true, colors));

            IonColorExpander.ItemsSource = IonGroups;
            PTMColorExpander.ItemsSource = Modifications;

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
            MetaDrawSettings.ShowDecoys = BoldTextCheckBox.IsChecked.Value;
            MetaDrawSettings.ShowContaminants = BoldTextCheckBox.IsChecked.Value;
            MetaDrawSettings.SpectrumDescription["Precursor Charge: "] = PrecursorChargeCheckBox.IsChecked.Value;
            MetaDrawSettings.SpectrumDescription["Precursor Mass: "] = PrecursorMassCheckBox.IsChecked.Value;
            MetaDrawSettings.SpectrumDescription["Theoretical Mass: "] = TheoreticalMassCheckBox.IsChecked.Value;
            MetaDrawSettings.SpectrumDescription["Score: "] = ScoreCheckBox.IsChecked.Value;
            MetaDrawSettings.SpectrumDescription["Protein Accession: "] = ProteinAccessionCheckBox.IsChecked.Value;
            MetaDrawSettings.SpectrumDescription["Protein: "] = ProteinCheckBox.IsChecked.Value;
            MetaDrawSettings.SpectrumDescription["Decoy/Contaminant/Target: "] = DecoyContaminantTargetCheckBox.IsChecked.Value;
            MetaDrawSettings.SpectrumDescription["Q-Value: "] = QValueCheckBox.IsChecked.Value;
            MetaDrawSettings.SpectrumDescription["Sequence Length: "] = SequenceLengthCheckBox.IsChecked.Value;
            MetaDrawSettings.SpectrumDescription["ProForma Level: "] = ProFormaLevelCheckBox.IsChecked.Value;
            MetaDrawSettings.SpectrumDescription["PEP: "] = PEPCheckBox.IsChecked.Value;
            MetaDrawSettings.SpectrumDescription["PEP Q-Value: "] = PEPQValueCheckBox.IsChecked.Value;
            MetaDrawSettings.LocalizationLevelStart = (LocalizationLevel)System.Enum.Parse(typeof(LocalizationLevel), CmbGlycanLocalizationLevelStart.SelectedItem.ToString());
            MetaDrawSettings.LocalizationLevelEnd = (LocalizationLevel)System.Enum.Parse(typeof(LocalizationLevel), CmbGlycanLocalizationLevelEnd.SelectedItem.ToString());

            foreach (var group in IonGroups)
            {
                foreach (var ion in group)
                {
                    if (ion.HasChanged)
                    {
                        if (ion.IsBeta)
                            MetaDrawSettings.BetaProductTypeToColor[ion.IonType] = MetaDrawSettings.PossibleColors.Keys.Where(p => p.GetColorName().Equals(ion.SelectedColor.Replace(" ", ""))).First();
                        else
                            MetaDrawSettings.ProductTypeToColor[ion.IonType] = MetaDrawSettings.PossibleColors.Keys.Where(p => p.GetColorName().Equals(ion.SelectedColor.Replace(" ", ""))).First();
                    }
                }
            }

            foreach (var group in Modifications)
            {
                foreach (var mod in group)
                {
                    if (mod.HasChanged)
                        MetaDrawSettings.ModificationTypeToColor[mod.ModName] = MetaDrawSettings.PossibleColors.Keys.Where(p => p.GetColorName().Equals(mod.SelectedColor.Replace(" ", ""))).First();
                }
            }

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
                    if (fontSize > 15)
                    {
                        MessageBox.Show("Font size must be <= 15");
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
            MetaDrawSettingsSnapshot settings = MetaDrawSettings.MakeSnapShot();
            XmlReaderWriter.WriteToXmlFile<MetaDrawSettingsSnapshot>(Path.Combine(GlobalVariables.DataDir, "DefaultParameters", @"MetaDrawSettingsDefault.xml"), settings);
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
            ((ModForTreeView)((ComboBox)sender).DataContext).SelectionChanged((string)((ComboBox)sender).SelectedItem);
        }
    }
}
