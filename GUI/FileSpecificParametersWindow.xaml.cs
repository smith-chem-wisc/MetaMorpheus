using EngineLayer;
using MzLibUtil;
using Nett;
using System;
using System.Collections.ObjectModel;
using System.IO;
using System.Linq;
using System.Windows;
using System.Windows.Input;
using TaskLayer;
using Proteomics.ProteolyticDigestion;

namespace MetaMorpheusGUI
{
    /// <summary>
    /// Interaction logic for ChangeParametersWindow.xaml
    /// </summary>
    public partial class FileSpecificParametersWindow : Window
    {

        //Window that is opened if user wishes to change file specific settings (TOML) for 
        //individual or multiple spectra files. Creates a toml file where settings can be
        //viewed, loaded, and changed from it.
        public FileSpecificParametersWindow(ObservableCollection<RawDataForDataGrid> selectedSpectraFiles)
        {
            SelectedSpectra = selectedSpectraFiles;
            InitializeComponent();
            PopulateChoices();
        }

        internal ObservableCollection<RawDataForDataGrid> SelectedSpectra { get; private set; }

        // write the toml settings file on clicking "save"
        private void Save_Click(object sender, RoutedEventArgs e)
        {
            var parametersToWrite = new FileSpecificParameters();

            // parse the file-specific parameters to text
            int paramsToSaveCount = 0;
            if (fileSpecificPrecursorMassTolEnabled.IsChecked.Value)
            {
                paramsToSaveCount++;
                if (GlobalGuiSettings.CheckPrecursorMassTolerance(precursorMassToleranceTextBox.Text))
                {
                    double value = double.Parse(precursorMassToleranceTextBox.Text);
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
                    return;
                }
            }
            if (fileSpecificProductMassTolEnabled.IsChecked.Value)
            {
                paramsToSaveCount++;
                if (GlobalGuiSettings.CheckProductMassTolerance(productMassToleranceTextBox.Text))
                {
                    double value = double.Parse(productMassToleranceTextBox.Text);
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
                if (int.TryParse(MinPeptideLengthTextBox.Text, out int i) && i > 0)
                {
                    parametersToWrite.MinPeptideLength = i;
                }
                else
                {
                    MessageBox.Show("The minimum peptide length must be a positive integer");
                    return;
                }
            }
            if (fileSpecificMaxPeptideLengthEnabled.IsChecked.Value)
            {
                paramsToSaveCount++;
                string lengthMaxPeptide = GlobalGuiSettings.MaxValueConversion(MaxPeptideLengthTextBox.Text);
                if (GlobalGuiSettings.CheckPeptideLength(MinPeptideLengthTextBox.Text, lengthMaxPeptide))
                {
                    parametersToWrite.MaxPeptideLength = int.Parse(lengthMaxPeptide);
                }
                else
                {
                    return;
                }
            }
            if (fileSpecificMissedCleavagesEnabled.IsChecked.Value)
            {
                paramsToSaveCount++;
                string lengthCleavage = GlobalGuiSettings.MaxValueConversion(missedCleavagesTextBox.Text);
                if (GlobalGuiSettings.CheckMaxMissedCleavages(lengthCleavage))
                {
                    parametersToWrite.MaxMissedCleavages = int.Parse(lengthCleavage);
                }
                else
                {
                    return;
                }
            }
            if (fileSpecificMaxModNumEnabled.IsChecked.Value)
            {
                paramsToSaveCount++;
                if (GlobalGuiSettings.CheckMaxModsPerPeptide(MaxModNumTextBox.Text))
                {
                    parametersToWrite.MaxModsForPeptide = int.Parse(MaxModNumTextBox.Text);
                }
                else
                {
                    return;
                }
            }
            //if (fileSpecificIonTypesEnabled.IsChecked.Value)
            //{
            //    paramsToSaveCount++;

            //    // don't think there's any way to mess up checkboxes... no error message needed
            //    parametersToWrite.BIons = bCheckBox.IsChecked;
            //    parametersToWrite.YIons = yCheckBox.IsChecked;
            //    parametersToWrite.CIons = cCheckBox.IsChecked;
            //    parametersToWrite.ZdotIons = zdotCheckBox.IsChecked;
            //}

            // write parameters to toml files for the selected spectra files


            var tomlPathsForSelectedFiles = SelectedSpectra.Select(p => Path.Combine(Directory.GetParent(p.FilePath).ToString(), Path.GetFileNameWithoutExtension(p.FileName)) + ".toml");
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
            var tempCommonParams = new CommonParameters();
            Protease tempProtease = tempCommonParams.DigestionParams.Protease;
            int tempMinPeptideLength = tempCommonParams.DigestionParams.MinPeptideLength;
            int tempMaxPeptideLength = tempCommonParams.DigestionParams.MaxPeptideLength;
            int tempMaxMissedCleavages = tempCommonParams.DigestionParams.MaxMissedCleavages;
            int tempMaxModsForPeptide = tempCommonParams.DigestionParams.MaxModsForPeptide;

            // do any of the selected files already have file-specific parameters specified?
            var spectraFiles = SelectedSpectra.Select(p => p.FilePath);
            foreach (string file in spectraFiles)
            {
                string tomlPath = Path.Combine(Directory.GetParent(file).ToString(), Path.GetFileNameWithoutExtension(file)) + ".toml";

                if (File.Exists(tomlPath))
                {
                    TomlTable tomlTable = Toml.ReadFile(tomlPath, MetaMorpheusTask.tomlConfig);
                    FileSpecificParameters fileSpecificParams = new FileSpecificParameters(tomlTable);

                    if (fileSpecificParams.PrecursorMassTolerance != null)
                    {
                        tempCommonParams.SetPrecursorMassTolerance(fileSpecificParams.PrecursorMassTolerance);
                        fileSpecificPrecursorMassTolEnabled.IsChecked = true;
                    }
                    if (fileSpecificParams.ProductMassTolerance != null)
                    {
                        tempCommonParams.SetProductMassTolerance(fileSpecificParams.ProductMassTolerance);
                        fileSpecificProductMassTolEnabled.IsChecked = true;
                    }
                    if (fileSpecificParams.Protease != null)
                    {
                        tempProtease = (fileSpecificParams.Protease);
                        fileSpecificProteaseEnabled.IsChecked = true;
                    }
                    if (fileSpecificParams.MinPeptideLength != null)
                    {
                        tempMinPeptideLength = (fileSpecificParams.MinPeptideLength.Value);
                        fileSpecificMinPeptideLengthEnabled.IsChecked = true;
                    }
                    if (fileSpecificParams.MaxPeptideLength != null)
                    {
                        tempMaxPeptideLength = (fileSpecificParams.MaxPeptideLength.Value);
                        fileSpecificMaxPeptideLengthEnabled.IsChecked = true;
                    }
                    if (fileSpecificParams.MaxMissedCleavages != null)
                    {
                        tempMaxMissedCleavages = (fileSpecificParams.MaxMissedCleavages.Value);
                        fileSpecificMissedCleavagesEnabled.IsChecked = true;
                    }
                    if (fileSpecificParams.MaxModsForPeptide != null)
                    {
                        tempMaxModsForPeptide = (fileSpecificParams.MaxMissedCleavages.Value);
                        fileSpecificMaxModNumEnabled.IsChecked = true;
                    }

                    //if (tempFileSpecificParams.BIons != null || tempFileSpecificParams.CIons != null
                    //    || tempFileSpecificParams.YIons != null || tempFileSpecificParams.ZdotIons != null)
                    //{
                    //    tempCommonParams.BIons = tempFileSpecificParams.BIons.Value;
                    //    tempCommonParams.YIons = tempFileSpecificParams.YIons.Value;
                    //    tempCommonParams.CIons = tempFileSpecificParams.CIons.Value;
                    //    tempCommonParams.ZdotIons = tempFileSpecificParams.ZdotIons.Value;

                    //    fileSpecificIonTypesEnabled.IsChecked = true;
                    //}
                }
            }

            DigestionParams digestParams = new DigestionParams(
                protease: tempProtease.Name,
                maxMissedCleavages: tempMaxMissedCleavages,
                minPeptideLength: tempMinPeptideLength,
                maxPeptideLength: tempMaxPeptideLength,
                maxModsForPeptides: tempMaxModsForPeptide);

            tempCommonParams.SetDigestionParams(digestParams);

            // populate the GUI
            foreach (Protease protease in ProteaseDictionary.Dictionary.Values)
            {
                fileSpecificProtease.Items.Add(protease);
            }

            fileSpecificProtease.SelectedItem = tempCommonParams.DigestionParams.Protease;

            productMassToleranceComboBox.Items.Add("Da");
            productMassToleranceComboBox.Items.Add("ppm");
            productMassToleranceComboBox.SelectedIndex = 1;

            precursorMassToleranceComboBox.Items.Add("Da");
            precursorMassToleranceComboBox.Items.Add("ppm");
            precursorMassToleranceComboBox.SelectedIndex = 1;

            precursorMassToleranceTextBox.Text = tempCommonParams.PrecursorMassTolerance.Value.ToString();
            productMassToleranceTextBox.Text = tempCommonParams.ProductMassTolerance.Value.ToString();
            MinPeptideLengthTextBox.Text = tempCommonParams.DigestionParams.MinPeptideLength.ToString();

            if (int.MaxValue != tempCommonParams.DigestionParams.MaxPeptideLength)
            {
                MaxPeptideLengthTextBox.Text = tempCommonParams.DigestionParams.MaxPeptideLength.ToString();
            }

            MaxModNumTextBox.Text = tempCommonParams.DigestionParams.MaxModsForPeptide.ToString();
            if (int.MaxValue != tempCommonParams.DigestionParams.MaxMissedCleavages)
            {
                missedCleavagesTextBox.Text = tempCommonParams.DigestionParams.MaxMissedCleavages.ToString();
            }
            //yCheckBox.IsChecked = tempCommonParams.YIons;
            //bCheckBox.IsChecked = tempCommonParams.BIons;
            //cCheckBox.IsChecked = tempCommonParams.CIons;
            //zdotCheckBox.IsChecked = tempCommonParams.ZdotIons;
        }

        private void CheckIfNumber(object sender, TextCompositionEventArgs e)
        {
            e.Handled = !GlobalGuiSettings.CheckIsNumber(e.Text);
        }
    }
}