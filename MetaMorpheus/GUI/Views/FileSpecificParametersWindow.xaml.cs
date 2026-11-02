using EngineLayer;
using GuiFunctions;
using MassSpectrometry;
using MzLibUtil;
using Nett;
using Omics.Digestion;
using Proteomics.ProteolyticDigestion;
using Readers;
using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Windows;
using System.Windows.Input;
using TaskLayer;
using Transcriptomics.Digestion;

namespace MetaMorpheusGUI
{
    public partial class FileSpecificParametersWindow : Window
    {
        private DeconHostViewModel _deconVm;
        private bool _clearDeconRequested;
        private readonly bool _isMultiFileSelection;

        public FileSpecificParametersWindow(ObservableCollection<RawDataForDataGrid> selectedSpectraFiles)
        {
            SelectedSpectra = selectedSpectraFiles;
            _isMultiFileSelection = selectedSpectraFiles.Count > 1;

            _deconVm = new DeconHostViewModel();
            _deconVm.DoPrecursorDeconvolution = true;

            InitializeComponent();

            DeconControl.DataContext = _deconVm;

            PopulateChoices();
        }

        internal ObservableCollection<RawDataForDataGrid> SelectedSpectra { get; private set; }

        public DeconHostViewModel DeconViewModel => _deconVm;

        public string? LastValidationMessage { get; private set; }

        private Func<string?>? _featureFilePickerOverride;
        public Func<string?>? FeatureFilePickerOverride
        {
            get => _featureFilePickerOverride;
            set
            {
                _featureFilePickerOverride = value;
                if (DeconControl != null)
                    DeconControl.FilePickerOverride = value;
            }
        }

        public Func<FileSpecificParameters, string>? CandidateSerializerOverride { get; set; }

        private void Save_Click(object sender, RoutedEventArgs e)
        {
            if (ExecuteSave())
                DialogResult = true;
            else if (LastValidationMessage != null)
                MessageBox.Show(LastValidationMessage, "Invalid Feature File",
                    MessageBoxButton.OK, MessageBoxImage.Warning);
        }

        public bool BrowseForFeatureFile()
        {
            var child = DeconControl?.FromFilePrecursorControl;
            if (child == null)
                return false;

            child.FilePickerOverride = FeatureFilePickerOverride;
            bool succeeded = child.BrowseForFile();
            if (!succeeded)
                LastValidationMessage = child.ValidationMessage;
            return succeeded;
        }

        public bool ExecuteSave()
        {
            LastValidationMessage = null;
            var sharedParams = new FileSpecificParameters();
            bool sharedParamExists = false;

            if (fileSpecificPrecursorMassTolEnabled.IsChecked.Value)
            {
                sharedParamExists = true;
                if (TaskValidator.CheckPrecursorMassTolerance(precursorMassToleranceTextBox.Text))
                {
                    double value = double.Parse(precursorMassToleranceTextBox.Text, CultureInfo.InvariantCulture);
                    sharedParams.PrecursorMassTolerance = precursorMassToleranceComboBox.SelectedIndex == 0
                        ? (Tolerance)new AbsoluteTolerance(value)
                        : new PpmTolerance(value);
                }
                else
                    return false;
            }
            if (fileSpecificProductMassTolEnabled.IsChecked.Value)
            {
                sharedParamExists = true;
                if (TaskValidator.CheckProductMassTolerance(productMassToleranceTextBox.Text))
                {
                    double value = double.Parse(productMassToleranceTextBox.Text, CultureInfo.InvariantCulture);
                    sharedParams.ProductMassTolerance = productMassToleranceComboBox.SelectedIndex == 0
                        ? (Tolerance)new AbsoluteTolerance(value)
                        : new PpmTolerance(value);
                }
                else
                    return false;
            }
            if (fileSpecificProteaseEnabled.IsChecked.Value)
            {
                sharedParamExists = true;
                sharedParams.DigestionAgent = (DigestionAgent)fileSpecificProtease.SelectedItem;
            }
            if (fileSpecificDissociationTypesEnabled.IsChecked.Value)
            {
                sharedParamExists = true;
                sharedParams.DissociationType = Enum.Parse<DissociationType>(fileSpecificDissociationType.SelectedItem.ToString());
            }
            if (fileSpecificSeparationTypesEnabled.IsChecked.Value)
            {
                sharedParamExists = true;
                sharedParams.SeparationType = (string)fileSpecificSeparationType.SelectedItem;
            }
            if (fileSpecificMinPeptideLengthEnabled.IsChecked.Value)
            {
                sharedParamExists = true;
                if (int.TryParse(MinPeptideLengthTextBox.Text, out int i) && i > 0)
                    sharedParams.MinPeptideLength = i;
                else
                {
                    MessageBox.Show("The minimum peptide length must be a positive integer");
                    return false;
                }
            }
            if (fileSpecificMaxPeptideLengthEnabled.IsChecked.Value)
            {
                sharedParamExists = true;
                string lengthMaxPeptide = TaskValidator.MaxValueConversion(MaxPeptideLengthTextBox.Text);
                if (TaskValidator.CheckPeptideLength(MinPeptideLengthTextBox.Text, lengthMaxPeptide))
                    sharedParams.MaxPeptideLength = int.Parse(lengthMaxPeptide);
                else
                    return false;
            }
            if (fileSpecificMissedCleavagesEnabled.IsChecked.Value)
            {
                sharedParamExists = true;
                string lengthCleavage = TaskValidator.MaxValueConversion(missedCleavagesTextBox.Text);
                if (TaskValidator.CheckMaxMissedCleavages(lengthCleavage))
                    sharedParams.MaxMissedCleavages = int.Parse(lengthCleavage);
                else
                    return false;
            }
            if (fileSpecificMaxModNumEnabled.IsChecked.Value)
            {
                sharedParamExists = true;
                if (TaskValidator.CheckMaxModsPerPeptide(MaxModNumTextBox.Text))
                    sharedParams.MaxModsForPeptide = int.Parse(MaxModNumTextBox.Text);
                else
                    return false;
            }

            if (!_isMultiFileSelection && !_clearDeconRequested && UseFileSpecificDeconCheckBox.IsChecked == true)
            {
                var selectedVm = _deconVm.PrecursorDeconvolutionParameters;
                if (selectedVm.DeconvolutionType == DeconvolutionType.FromFile)
                {
                    var fromFileVm = (FromFileDeconParamsViewModel)selectedVm;
                    if (!FromFileDeconParamsViewModel.TryValidateFilePath(fromFileVm.FilePath, out string msg))
                    {
                        LastValidationMessage = $"MS1 feature file validation failed:\n{msg}";
                        return false;
                    }
                }
            }

            // Phase 1 — build and serialize every candidate in memory; surface any error
            var candidates = new List<(string dir, string filename, string tomlPath,
                                       string? serializedToml, bool hasParams)>();

            foreach (var spectra in SelectedSpectra)
            {
                string dir = Directory.GetParent(spectra.FilePath)!.ToString();
                string filename = Path.GetFileNameWithoutExtension(spectra.FileName) + ".toml";
                string tomlPath = Path.Combine(dir, filename);

                DeconvolutionParameters? existingDecon = null;
                if (File.Exists(tomlPath))
                {
                    try
                    {
                        var table = Toml.ReadFile(tomlPath, MetaMorpheusTask.tomlConfig);
                        existingDecon = new FileSpecificParameters(table).PrecursorDeconvolutionParameters;
                    }
                    catch (Exception ex)
                    {
                        LastValidationMessage =
                            $"Cannot read existing settings for '{spectra.FileName}': {ex.Message}\n" +
                            "Fix or remove the companion TOML before saving.";
                        return false;
                    }
                }

                DeconvolutionParameters? deconToWrite;
                if (_isMultiFileSelection)
                    deconToWrite = existingDecon;
                else if (_clearDeconRequested)
                    deconToWrite = null;
                else if (UseFileSpecificDeconCheckBox.IsChecked == true)
                    deconToWrite = _deconVm.PrecursorDeconvolutionParameters.Parameters;
                else
                    deconToWrite = existingDecon;

                var merged = sharedParams.Clone();
                merged.PrecursorDeconvolutionParameters = deconToWrite;
                bool hasParams = sharedParamExists || deconToWrite != null;

                string? serializedToml = null;
                if (hasParams)
                {
                    try
                    {
                        serializedToml = CandidateSerializerOverride != null
                            ? CandidateSerializerOverride(merged)
                            : Toml.WriteString(merged, MetaMorpheusTask.tomlConfig);
                        var roundTripTable = Toml.ReadString<TomlTable>(serializedToml,
                                                MetaMorpheusTask.tomlConfig);
                        _ = new FileSpecificParameters(roundTripTable);
                    }
                    catch (Exception ex)
                    {
                        LastValidationMessage =
                            $"Cannot serialize settings for '{spectra.FileName}': {ex.Message}";
                        return false;
                    }
                }

                candidates.Add((dir, filename, tomlPath, serializedToml, hasParams));
            }

            // Phase 2 — all candidates valid; archive existing files then commit writes
            foreach (var (dir, filename, tomlPath, serializedToml, hasParams) in candidates)
            {
                if (File.Exists(tomlPath))
                    AccomodateNewFileSpecificToml(dir, filename);

                if (hasParams)
                    File.WriteAllText(tomlPath, serializedToml!);
                else
                    File.Delete(tomlPath);
            }

            return true;
        }

        private void ClearDeconButton_Click(object sender, RoutedEventArgs e)
        {
            _clearDeconRequested = true;
            UseFileSpecificDeconCheckBox.IsChecked = false;
            ClearDeconButton.Content = "Deconvolution override will be cleared on save";
            ClearDeconButton.IsEnabled = false;
        }

        private void UseFileSpecificDeconCheckBox_Checked(object sender, RoutedEventArgs e)
        {
            _clearDeconRequested = false;
            if (!_isMultiFileSelection)
                ClearDeconButton.IsEnabled = true;
            ClearDeconButton.Content = "Clear file-specific deconvolution override";
        }

        private void AccomodateNewFileSpecificToml(string directoryForThisMsFile, string filename)
        {
            string oldFileSpecificTomlDirectory = Path.Combine(directoryForThisMsFile, "OldFileSpecificTomls");
            Directory.CreateDirectory(oldFileSpecificTomlDirectory);
            string fullFilePath = Path.Combine(oldFileSpecificTomlDirectory, filename);
            if (File.Exists(fullFilePath))
            {
                AccomodateNewFileSpecificToml(oldFileSpecificTomlDirectory, filename);
            }
            System.IO.File.Copy(Path.Combine(directoryForThisMsFile, filename),
                Path.Combine(oldFileSpecificTomlDirectory, filename), true);
        }

        private void Cancel_Click(object sender, RoutedEventArgs e)
        {
            DialogResult = false;
        }

        private void PopulateChoices()
        {
            var defaultParams = new CommonParameters();
            IDigestionParams digestionParams = GuiGlobalParamsViewModel.Instance.IsRnaMode
                ? new RnaDigestionParams("RNase T1")
                : new DigestionParams("trypsin");

            DigestionAgent tempProtease = digestionParams.DigestionAgent;
            int tempMinPeptideLength = digestionParams.MinLength;
            int tempMaxPeptideLength = digestionParams.MaxLength;
            int tempMaxMissedCleavages = digestionParams.MaxMissedCleavages;
            int tempMaxModsForPeptide = digestionParams.MaxMods;
            var tempPrecursorMassTolerance = defaultParams.PrecursorMassTolerance;
            var tempProductMassTolerance = defaultParams.ProductMassTolerance;
            DissociationType tempDissociationType = defaultParams.DissociationType;
            string tempSeparationType = defaultParams.SeparationType;

            var spectraFiles = SelectedSpectra.Select(p => p.FilePath);
            DeconvolutionParameters? loadedSingleFileDecon = null;

            foreach (string file in spectraFiles)
            {
                string tomlPath = Path.Combine(Directory.GetParent(file)!.ToString(),
                    Path.GetFileNameWithoutExtension(file)) + ".toml";

                if (File.Exists(tomlPath))
                {
                    TomlTable tomlTable = Toml.ReadFile(tomlPath, MetaMorpheusTask.tomlConfig);
                    FileSpecificParameters fileSpecificParams = new(tomlTable);

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
                    if (fileSpecificParams.DigestionAgent != null)
                    {
                        tempProtease = fileSpecificParams.DigestionAgent;
                        fileSpecificProteaseEnabled.IsChecked = true;
                    }
                    if (fileSpecificParams.DissociationType != null)
                    {
                        tempDissociationType = fileSpecificParams.DissociationType.Value;
                        fileSpecificDissociationTypesEnabled.IsChecked = true;
                    }
                    if (fileSpecificParams.SeparationType != null)
                    {
                        tempSeparationType = fileSpecificParams.SeparationType;
                        fileSpecificSeparationTypesEnabled.IsChecked = true;
                    }
                    if (fileSpecificParams.MinPeptideLength != null)
                    {
                        tempMinPeptideLength = fileSpecificParams.MinPeptideLength.Value;
                        fileSpecificMinPeptideLengthEnabled.IsChecked = true;
                    }
                    if (fileSpecificParams.MaxPeptideLength != null)
                    {
                        tempMaxPeptideLength = fileSpecificParams.MaxPeptideLength.Value;
                        fileSpecificMaxPeptideLengthEnabled.IsChecked = true;
                    }
                    if (fileSpecificParams.MaxMissedCleavages != null)
                    {
                        tempMaxMissedCleavages = fileSpecificParams.MaxMissedCleavages.Value;
                        fileSpecificMissedCleavagesEnabled.IsChecked = true;
                    }
                    if (fileSpecificParams.MaxModsForPeptide != null)
                    {
                        tempMaxModsForPeptide = fileSpecificParams.MaxMissedCleavages.Value;
                        fileSpecificMaxModNumEnabled.IsChecked = true;
                    }

                    if (!_isMultiFileSelection && fileSpecificParams.PrecursorDeconvolutionParameters != null)
                        loadedSingleFileDecon = fileSpecificParams.PrecursorDeconvolutionParameters;
                }
            }

            if (GuiGlobalParamsViewModel.Instance.IsRnaMode)
            {
                foreach (Rnase rnase in RnaseDictionary.Dictionary.Values)
                    fileSpecificProtease.Items.Add(rnase);
            }
            else
            {
                foreach (Protease protease in ProteaseDictionary.Dictionary.Values)
                    fileSpecificProtease.Items.Add(protease);
            }

            fileSpecificProtease.SelectedItem = tempProtease;

            foreach (DissociationType dissociationType in Enum.GetValues(typeof(DissociationType)))
                fileSpecificDissociationType.Items.Add(dissociationType);

            fileSpecificDissociationType.SelectedItem = DissociationType.HCD;

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
            MinPeptideLengthTextBox.Text = tempMinPeptideLength.ToString();

            if (int.MaxValue != tempMaxPeptideLength)
                MaxPeptideLengthTextBox.Text = tempMaxPeptideLength.ToString();

            MaxModNumTextBox.Text = tempMaxModsForPeptide.ToString();
            if (int.MaxValue != tempMaxMissedCleavages)
                missedCleavagesTextBox.Text = tempMaxMissedCleavages.ToString();

            if (_isMultiFileSelection)
            {
                _deconVm.EnableFromFilePrecursorDeconvolution();
                UseFileSpecificDeconCheckBox.IsEnabled = false;
                ClearDeconButton.IsEnabled = false;
                MultiFileDeconNote.Visibility = Visibility.Visible;
            }
            else if (loadedSingleFileDecon != null)
            {
                ApplyExistingDeconToVm(loadedSingleFileDecon);
                UseFileSpecificDeconCheckBox.IsChecked = true;
                ClearDeconButton.IsEnabled = true;
            }
            else
            {
                _deconVm.EnableFromFilePrecursorDeconvolution();
            }
        }

        private void ApplyExistingDeconToVm(DeconvolutionParameters existingDecon)
        {
            if (existingDecon is FromFileDeconvolutionParameters fromFileParams)
            {
                _deconVm.EnableFromFilePrecursorDeconvolution(fromFileParams);
            }
            else
            {
                _deconVm = new DeconHostViewModel(initialPrecursorParameters: existingDecon);
                _deconVm.DoPrecursorDeconvolution = true;
                _deconVm.EnableFromFilePrecursorDeconvolution();
                DeconControl.DataContext = _deconVm;
            }
        }

        private void KeyPressed(object sender, KeyEventArgs e)
        {
            if (e.Key == Key.Return)
                Save_Click(sender, e);
            else if (e.Key == Key.Escape)
                Cancel_Click(sender, e);
        }
    }
}
