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
using System.Globalization;

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
            bool fileSpecificParameterExists = false;
            // parse the file-specific parameters to text
            if (fileSpecificPrecursorMassTolEnabled.IsChecked.Value)
            {
                fileSpecificParameterExists = true;
                if (GlobalGuiSettings.CheckPrecursorMassTolerance(precursorMassToleranceTextBox.Text))
                {
                    double value = double.Parse(precursorMassToleranceTextBox.Text, CultureInfo.InvariantCulture);
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
                fileSpecificParameterExists = true;
                if (GlobalGuiSettings.CheckProductMassTolerance(productMassToleranceTextBox.Text))
                {
                    double value = double.Parse(productMassToleranceTextBox.Text, CultureInfo.InvariantCulture);
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
                fileSpecificParameterExists = true;
                parametersToWrite.Protease = (Protease)fileSpecificProtease.SelectedItem;
            }
            if (fileSpecificSeparationTypesEnabled.IsChecked.Value)
            {
                fileSpecificParameterExists = true;
                parametersToWrite.SeparationType = (string)fileSpecificSeparationType.SelectedItem;
            }
            if (fileSpecificMinPeptideLengthEnabled.IsChecked.Value)
            {
                fileSpecificParameterExists = true;
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
                fileSpecificParameterExists = true;
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
                fileSpecificParameterExists = true;
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
                fileSpecificParameterExists = true;
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

            foreach (var spectra in SelectedSpectra)
            {
                string directoryForThisMsFile = Directory.GetParent(spectra.FilePath).ToString();
                string filename = Path.GetFileNameWithoutExtension(spectra.FileName) + ".toml";
                string tomlToWrite = Path.Combine(directoryForThisMsFile, filename);

                //check if a toml file already exists
                if (File.Exists(tomlToWrite))
                {
                    //store the previous file-specific tomls (if any) in a folder for reproducibility
                    AccomodateNewFileSpecificToml(directoryForThisMsFile, filename);
                    //check that an old toml doesn't already exist. If it does, move it to another nested folder.
                }

                if (fileSpecificParameterExists)
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

        //recursive function that ensures no tomls are deleted when new ones are generated
        private void AccomodateNewFileSpecificToml(string directoryForThisMsFile, string filename)
        {
            //make a directory to store the old files. If the directory already exists, this method does not create a new directory
            string oldFileSpecificTomlDirectory = Path.Combine(directoryForThisMsFile, "OldFileSpecificTomls");
            System.IO.Directory.CreateDirectory(oldFileSpecificTomlDirectory);
            string fullFilePath = Path.Combine(oldFileSpecificTomlDirectory, filename);
            //if an old version already exists, we need to move it
            if(File.Exists(fullFilePath))
            {
                AccomodateNewFileSpecificToml(oldFileSpecificTomlDirectory, filename);
            }
            //move the old file into the sub-directory
            System.IO.File.Copy(Path.Combine(directoryForThisMsFile, filename), Path.Combine(oldFileSpecificTomlDirectory, filename),true);
        }

        // exits dialog; nothing is written
        private void Cancel_Click(object sender, RoutedEventArgs e)
        {
            DialogResult = false;
        }

        private void PopulateChoices()
        {
            // use default settings to populate
            var defaultParams = new CommonParameters();
            Protease tempProtease = defaultParams.DigestionParams.Protease;
            int tempMinPeptideLength = defaultParams.DigestionParams.MinPeptideLength;
            int tempMaxPeptideLength = defaultParams.DigestionParams.MaxPeptideLength;
            int tempMaxMissedCleavages = defaultParams.DigestionParams.MaxMissedCleavages;
            int tempMaxModsForPeptide = defaultParams.DigestionParams.MaxModsForPeptide;
            var tempPrecursorMassTolerance = defaultParams.PrecursorMassTolerance;
            var tempProductMassTolerance = defaultParams.ProductMassTolerance;
            string tempSeparationType = defaultParams.SeparationType;

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
                        tempPrecursorMassTolerance = fileSpecificParams.PrecursorMassTolerance;
                        fileSpecificPrecursorMassTolEnabled.IsChecked = true;
                    }
                    if (fileSpecificParams.ProductMassTolerance != null)
                    {
                        tempProductMassTolerance = fileSpecificParams.ProductMassTolerance;
                        fileSpecificProductMassTolEnabled.IsChecked = true;
                    }
                    if (fileSpecificParams.Protease != null)
                    {
                        tempProtease = (fileSpecificParams.Protease);
                        fileSpecificProteaseEnabled.IsChecked = true;
                    }
                    if (fileSpecificParams.SeparationType != null)
                    {
                        tempSeparationType = (fileSpecificParams.SeparationType);
                        fileSpecificSeparationTypesEnabled.IsChecked = true;
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
                }
            }

            DigestionParams digestParams = new DigestionParams(
                protease: tempProtease.Name,
                maxMissedCleavages: tempMaxMissedCleavages,
                minPeptideLength: tempMinPeptideLength,
                maxPeptideLength: tempMaxPeptideLength,
                maxModsForPeptides: tempMaxModsForPeptide);

            // populate the GUI
            foreach (Protease protease in ProteaseDictionary.Dictionary.Values)
            {
                fileSpecificProtease.Items.Add(protease);
            }

            fileSpecificProtease.SelectedItem = digestParams.Protease;

            fileSpecificSeparationType.Items.Add("HPLC");
            fileSpecificSeparationType.Items.Add("CZE");
            fileSpecificSeparationType.SelectedIndex = fileSpecificSeparationType.Items.IndexOf(tempSeparationType);

            productMassToleranceComboBox.Items.Add("Da");
            productMassToleranceComboBox.Items.Add("ppm");
            productMassToleranceComboBox.SelectedIndex = 1;

            precursorMassToleranceComboBox.Items.Add("Da");
            precursorMassToleranceComboBox.Items.Add("ppm");
            precursorMassToleranceComboBox.SelectedIndex = 1;

            precursorMassToleranceTextBox.Text = tempPrecursorMassTolerance.Value.ToString();
            productMassToleranceTextBox.Text = tempProductMassTolerance.Value.ToString();
            MinPeptideLengthTextBox.Text = digestParams.MinPeptideLength.ToString();

            if (int.MaxValue != digestParams.MaxPeptideLength)
            {
                MaxPeptideLengthTextBox.Text = digestParams.MaxPeptideLength.ToString();
            }

            MaxModNumTextBox.Text = digestParams.MaxModsForPeptide.ToString();
            if (int.MaxValue != digestParams.MaxMissedCleavages)
            {
                missedCleavagesTextBox.Text = digestParams.MaxMissedCleavages.ToString();
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

        private void KeyPressed(object sender, KeyEventArgs e)
        {
            if (e.Key == Key.Return)
            {
                Save_Click(sender, e);
            }
            else if (e.Key == Key.Escape)
            {
                Cancel_Click(sender, e);
            }
        }
    }
}