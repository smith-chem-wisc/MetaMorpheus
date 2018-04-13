using EngineLayer;
using MzLibUtil;
using Nett;
using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Input;
using TaskLayer;

namespace MetaMorpheusGUI
{
    /// <summary>
    /// Interaction logic for ChangeParametersWindow.xaml
    /// </summary>
    public partial class FileSpecificParametersWindow : Window
    {
        #region Public Constructors

        //Window that is opened if user wishes to change file specific settings (TOML) for 
        //individual or multiple spectra files. Creates a toml file where settings can be
        //viewed, loaded, and changed from it.
        public FileSpecificParametersWindow(ObservableCollection<RawDataForDataGrid> selectedSpectraFiles)
        {
            SelectedRaw = selectedSpectraFiles;
            InitializeComponent();
            PopulateChoices();
        }

        #endregion Public Constructors

        #region Internal Properties

        internal ObservableCollection<RawDataForDataGrid> SelectedRaw { get; private set; }

        #endregion Internal Properties

        #region Private Methods

        // write the toml settings file on clicking "save"
        private void Save_Click(object sender, RoutedEventArgs e)
        {
            // just using this to get variable names in case they've been renamed...
            var temp = new CommonParameters();
            var parametersToWrite = new FileSpecificParameters();

            // parse the file-specific parameters to text
            int paramsToSaveCount = 0;
            if (fileSpecificPrecursorMassTolEnabled.IsChecked.Value)
            {
                paramsToSaveCount++;
                double value = 0;
                if(double.TryParse(precursorMassToleranceTextBox.Text, out value))
                {
                    if (precursorMassToleranceComboBox.SelectedIndex == 0)
                    {
                        parametersToWrite.PrecursorMassTolerance = new AbsoluteTolerance(value);
                    }
                    else
                    {
                        parametersToWrite.PrecursorMassTolerance = new PpmTolerance(value);
                    }
                }
                else
                {
                    MessageBox.Show("Precursor tolerance must be a number", "Error", MessageBoxButton.OK, MessageBoxImage.Hand);
                    return;
                }
            }
            if (fileSpecificProductMassTolEnabled.IsChecked.Value)
            {
                paramsToSaveCount++;
                double value = 0;
                if (double.TryParse(productMassToleranceTextBox.Text, out value))
                {
                    if (productMassToleranceComboBox.SelectedIndex == 0)
                    {
                        parametersToWrite.ProductMassTolerance = new AbsoluteTolerance(value);
                    }
                    else
                    {
                        parametersToWrite.ProductMassTolerance = new PpmTolerance(value);
                    }
                }
                else
                {
                    MessageBox.Show("Product tolerance must be a number", "Error", MessageBoxButton.OK, MessageBoxImage.Hand);
                    return;
                }
            }
            if (fileSpecificProteaseEnabled.IsChecked.Value)
            {
                paramsToSaveCount++;
                parametersToWrite.Protease = (Protease)fileSpecificProtease.SelectedItem;
            }
            if (fileSpecificMinPeptideLengthEnabled.IsChecked.Value)
            {
                paramsToSaveCount++;
                if (int.TryParse(txtMinPeptideLength.Text, out int i) && i > 0)
                {
                    parametersToWrite.MinPeptideLength = i;
                }
                else
                {
                    MessageBox.Show("Min peptide length must be an integer larger than 0", "Error", MessageBoxButton.OK, MessageBoxImage.Hand);
                    return;
                }
            }
            if (fileSpecificMaxPeptideLengthEnabled.IsChecked.Value)
            {
                paramsToSaveCount++;
                if (string.IsNullOrEmpty(txtMaxPeptideLength.Text))
                {
                    parametersToWrite.MaxPeptideLength = int.MaxValue;
                }
                else if (int.TryParse(txtMaxPeptideLength.Text, out int i) && i > 0)
                {
                    parametersToWrite.MaxPeptideLength = i;
                }
                else
                {
                    MessageBox.Show("Max peptide length must be an integer larger than 0, or left blank", "Error", MessageBoxButton.OK, MessageBoxImage.Hand);
                    return;
                }
            }
            if (fileSpecificMissedCleavagesEnabled.IsChecked.Value)
            {
                paramsToSaveCount++;
                if (int.TryParse(missedCleavagesTextBox.Text, out int i) && i >= 0)
                {
                    parametersToWrite.MaxMissedCleavages = i;
                }
                else
                {
                    MessageBox.Show("Missed cleavages must be an integer greater than or equal to 0", "Error", MessageBoxButton.OK, MessageBoxImage.Hand);
                    return;
                }
            }
            if (fileSpecificMaxModNumEnabled.IsChecked.Value)
            {
                paramsToSaveCount++;
                if (int.TryParse(txtMaxModNum.Text, out int i) && i >= 0)
                {
                    parametersToWrite.MaxModsForPeptide = i;
                }
                else
                {
                    MessageBox.Show("Mods per peptide must be an integer greater than or equal to 0", "Error", MessageBoxButton.OK, MessageBoxImage.Hand);
                    return;
                }
            }
            if (fileSpecificIonTypesEnabled.IsChecked.Value)
            {
                paramsToSaveCount++;

                // don't think there's any way to mess up checkboxes... no error needed
                parametersToWrite.BIons = bCheckBox.IsChecked;
                parametersToWrite.YIons = yCheckBox.IsChecked;
                parametersToWrite.CIons = cCheckBox.IsChecked;
                parametersToWrite.ZdotIons = zdotCheckBox.IsChecked;
            }

            // write parameters to toml files for the selected spectra files
            var tomlPathsForSelectedFiles = SelectedRaw.Select(p => Path.Combine(Directory.GetParent(p.FilePath).ToString(), Path.GetFileNameWithoutExtension(p.FileName)) + ".toml");
            foreach (var tomlToWrite in tomlPathsForSelectedFiles)
            {
                if (paramsToSaveCount > 0)
                {
                    Toml.WriteFile(parametersToWrite, tomlToWrite, MetaMorpheusTask.tomlConfig);

                    // make sure the settings are able to be parsed...
                    var tempTomlTable = Toml.ReadFile(tomlToWrite, MetaMorpheusTask.tomlConfig);
                    FileSpecificParameters tempParams = new FileSpecificParameters(tempTomlTable);
                }
                else
                {
                    // user has specified that no file-specific settings should be used; delete the file-specific toml if it exists
                    File.Delete(tomlToWrite);
                }
            }

            // done
            DialogResult = true;
        }

        // exits dialog; nothing is written
        private void Cancel_Click(object sender, RoutedEventArgs e)
        {
            DialogResult = false;
        }

        private void PopulateChoices()
        {
            // use default settings to populate
            var temp = new CommonParameters();
            foreach (Protease protease in GlobalVariables.ProteaseDictionary.Values)
            {
                fileSpecificProtease.Items.Add(protease);
            }
            fileSpecificProtease.SelectedIndex = 11;

            productMassToleranceComboBox.Items.Add("Da");
            productMassToleranceComboBox.Items.Add("ppm");
            productMassToleranceComboBox.SelectedIndex = 1;

            precursorMassToleranceComboBox.Items.Add("Da");
            precursorMassToleranceComboBox.Items.Add("ppm");
            precursorMassToleranceComboBox.SelectedIndex = 1;

            precursorMassToleranceTextBox.Text = temp.PrecursorMassTolerance.Value.ToString();
            productMassToleranceTextBox.Text = temp.ProductMassTolerance.Value.ToString();
            txtMinPeptideLength.Text = temp.DigestionParams.MinPeptideLength.ToString();

            if (int.MaxValue != temp.DigestionParams.MaxPeptideLength)
            {
                txtMaxPeptideLength.Text = temp.DigestionParams.MaxPeptideLength.ToString();
            }

            txtMaxModNum.Text = temp.DigestionParams.MaxModsForPeptide.ToString();
            missedCleavagesTextBox.Text = temp.DigestionParams.MaxMissedCleavages.ToString();

            yCheckBox.IsChecked = temp.YIons;
            bCheckBox.IsChecked = temp.BIons;
            cCheckBox.IsChecked = temp.CIons;
            zdotCheckBox.IsChecked = temp.ZdotIons;
        }

        #endregion Private Methods
    }
}