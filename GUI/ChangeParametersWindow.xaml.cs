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
    public partial class ChangeFileSpecificParametersWindow : Window
    {
        #region Public Constructors

        //Window that is opened if user wishes to change file specific settings (TOML) for 
        //individual or multiple spectra files. Creates a toml file where settings can be
        //viewed, loaded, and changed from it.
        public ChangeFileSpecificParametersWindow(ObservableCollection<RawDataForDataGrid> selectedSpectraFiles)
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

        private void Save_Click(object sender, RoutedEventArgs e)
        {
            // just using this to get variable names in case they've been renamed...
            var temp = new CommonParameters();

            // parse the file-specific parameters to text
            List<string> fileSpecificOutput = new List<string>();
            if (fileSpecificPrecursorMassTolEnabled.IsChecked.Value)
            {
                Tolerance tol = null;
                if (precursorMassToleranceComboBox.SelectedIndex == 0)
                {
                    tol = new AbsoluteTolerance(double.Parse(precursorMassToleranceTextBox.Text, CultureInfo.InvariantCulture));
                }
                else
                {
                    tol = new PpmTolerance(double.Parse(precursorMassToleranceTextBox.Text, CultureInfo.InvariantCulture));
                }

                fileSpecificOutput.Add(nameof(temp.PrecursorMassTolerance) + " = \"" + tol.ToString() + "\"");
            }
            if (fileSpecificProductMassTolEnabled.IsChecked.Value)
            {
                Tolerance tol = null;
                if (productMassToleranceComboBox.SelectedIndex == 0)
                {
                    tol = new AbsoluteTolerance(double.Parse(productMassToleranceTextBox.Text, CultureInfo.InvariantCulture));
                }
                else
                {
                    tol = new PpmTolerance(double.Parse(productMassToleranceTextBox.Text, CultureInfo.InvariantCulture));
                }

                fileSpecificOutput.Add(nameof(temp.ProductMassTolerance) + " = \"" + tol.ToString() + "\"");
            }
            if (fileSpecificProteaseEnabled.IsChecked.Value)
            {
                fileSpecificOutput.Add(nameof(temp.DigestionParams.Protease) + " = \"" + fileSpecificProtease.Text + "\"");
            }
            if (fileSpecificMinPeptideLengthEnabled.IsChecked.Value)
            {
                fileSpecificOutput.Add(nameof(temp.DigestionParams.MinPeptideLength) + " = " + txtMinPeptideLength.Text);
            }
            if (fileSpecificMaxPeptideLengthEnabled.IsChecked.Value)
            {
                fileSpecificOutput.Add(nameof(temp.DigestionParams.MaxPeptideLength) + " = " + txtMaxPeptideLength.Text);
            }
            if (fileSpecificMissedCleavagesEnabled.IsChecked.Value)
            {
                fileSpecificOutput.Add(nameof(temp.DigestionParams.MaxMissedCleavages) + " = " + missedCleavagesTextBox.Text);
            }
            if (fileSpecificMaxModNumEnabled.IsChecked.Value)
            {
                fileSpecificOutput.Add(nameof(temp.DigestionParams.MaxModsForPeptide) + " = " + txtMaxModNum.Text);
            }
            if (fileSpecificIonTypesEnabled.IsChecked.Value)
            {
                fileSpecificOutput.Add(nameof(temp.BIons) + " = " + bCheckBox.IsChecked);
                fileSpecificOutput.Add(nameof(temp.YIons) + " = " + yCheckBox.IsChecked);
                fileSpecificOutput.Add(nameof(temp.CIons) + " = " + cCheckBox.IsChecked);
                fileSpecificOutput.Add(nameof(temp.ZdotIons) + " = " + zdotCheckBox.IsChecked);
            }

            // write parameters to toml files for the selected spectra files
            var tomlPathsForSelectedFiles = SelectedRaw.Select(p => Path.Combine(Directory.GetParent(p.FilePath).ToString(), Path.GetFileNameWithoutExtension(p.FileName)) + ".toml");
            foreach (var tomlToWrite in tomlPathsForSelectedFiles)
            {
                if (fileSpecificOutput.Any())
                {
                    File.WriteAllLines(tomlToWrite, fileSpecificOutput);
                }
                else
                {
                    File.Delete(tomlToWrite);
                }
            }

            // done
            DialogResult = true;
        }

        private void Cancel_Click(object sender, RoutedEventArgs e)
        {
            // exits dialog; nothing is written
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

            productMassToleranceComboBox.Items.Add("Absolute");
            productMassToleranceComboBox.Items.Add("ppm");
            productMassToleranceComboBox.SelectedIndex = 1;

            precursorMassToleranceComboBox.Items.Add("Absolute");
            precursorMassToleranceComboBox.Items.Add("ppm");
            precursorMassToleranceComboBox.SelectedIndex = 1;

            precursorMassToleranceTextBox.Text = temp.PrecursorMassTolerance.Value.ToString();
            productMassToleranceTextBox.Text = temp.ProductMassTolerance.Value.ToString();
            txtMinPeptideLength.Text = temp.DigestionParams.MinPeptideLength.ToString();

            if (!(int.MaxValue == temp.DigestionParams.MaxPeptideLength))
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