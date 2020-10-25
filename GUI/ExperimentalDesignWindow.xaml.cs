using EngineLayer;
using FlashLFQ;
using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.IO;
using System.Linq;
using System.Windows;
using System.Windows.Input;

namespace MetaMorpheusGUI
{
    /// <summary>
    /// Interaction logic for ExperimentalDesignWindow.xaml
    /// </summary>
    public partial class ExperimentalDesignWindow : Window
    {
        private readonly ObservableCollection<ExperimentalDesignForDataGrid> SpectraFileExperimentalDesign;
        private readonly ObservableCollection<RawDataForDataGrid> SpectraFiles;
        private string outputPath;

        public ExperimentalDesignWindow(ObservableCollection<RawDataForDataGrid> spectraFilesObservableCollection)
        {
            InitializeComponent();
            SpectraFileExperimentalDesign = new ObservableCollection<ExperimentalDesignForDataGrid>();
            this.SpectraFiles = spectraFilesObservableCollection;
            InitializeExperDesign();
        }

        private void InitializeExperDesign()
        {
            SpectraFileExperimentalDesign.Clear();

            foreach (var item in SpectraFiles.Where(p => p.Use).Select(p => p.FilePath))
            {
                SpectraFileExperimentalDesign.Add(new ExperimentalDesignForDataGrid(item));
            }

            if (SpectraFileExperimentalDesign.Any())
            {
                outputPath = Directory.GetParent(SpectraFileExperimentalDesign.First().FullFilePathWithExtension).FullName;
                outputPath = Path.Combine(outputPath, GlobalVariables.ExperimentalDesignFileName);
            }

            DgQuant.DataContext = SpectraFileExperimentalDesign;

            // read/set existing experimental design
            var filePaths = SpectraFileExperimentalDesign.Select(p => p.FullFilePathWithExtension).ToList();
            var existingDesign = ExperimentalDesign.ReadExperimentalDesign(outputPath, filePaths, out var errors);
            foreach (SpectraFileInfo fileInfo in existingDesign)
            {
                ExperimentalDesignForDataGrid match = SpectraFileExperimentalDesign.FirstOrDefault(p => p.FullFilePathWithExtension == fileInfo.FullFilePathWithExtension);

                if (match != null)
                {
                    match.Condition = fileInfo.Condition;
                    match.Biorep = (fileInfo.BiologicalReplicate + 1).ToString();
                    match.Techrep = (fileInfo.TechnicalReplicate + 1).ToString();
                    match.Fraction = (fileInfo.Fraction + 1).ToString();
                }
            }
        }

        private void BtnSaveQuant_Click(object sender, RoutedEventArgs e)
        {
            if (outputPath == null)
            {
                // no spectra files
                DialogResult = true;
                return;
            }

            if (CheckForExperimentalDesignErrors() != null)
            {
                MessageBox.Show(CheckForExperimentalDesignErrors());
                return;
            }

            try
            {
                var fileInfos = SpectraFileExperimentalDesign
                    .Select(p => new SpectraFileInfo(p.FullFilePathWithExtension, p.Condition, int.Parse(p.Biorep) - 1, int.Parse(p.Techrep) - 1, int.Parse(p.Fraction) - 1))
                    .ToList();
                ExperimentalDesign.WriteExperimentalDesignToFile(fileInfos);
            }
            catch (Exception ex)
            {
                MessageBox.Show("Could not save experimental design!\n\n" + ex.Message);
                return;
            }

            DialogResult = true;
        }

        private void BtnCancelQuant_Click(object sender, RoutedEventArgs e)
        {
            DialogResult = false;
        }

        public string CheckForExperimentalDesignErrors()
        {
            // check for basic parsing
            foreach (var item in SpectraFileExperimentalDesign)
            {
                if (string.IsNullOrEmpty(item.Condition))
                {
                    return "Condition: " + item.Condition +
                        "\nBiorep: " + item.Biorep +
                        "\nFraction: " + item.Fraction +
                        "\nTechrep: " + item.Techrep +
                        "\n\nCondition cannot be blank";
                }

                if (!int.TryParse(item.Biorep, out int b) || b < 1)
                {
                    return "Condition: " + item.Condition +
                        "\nBiorep: " + item.Biorep +
                        "\nFraction: " + item.Fraction +
                        "\nTechrep: " + item.Techrep +
                        "\n\nBiorep must be an integer >= 1";
                }

                if (!int.TryParse(item.Fraction, out int f) || f < 1)
                {
                    return "Condition: " + item.Condition +
                        "\nBiorep: " + item.Biorep +
                        "\nFraction: " + item.Fraction +
                        "\nTechrep: " + item.Techrep +
                        "\n\nFraction must be an integer >= 1";
                }

                if (!int.TryParse(item.Techrep, out int t) || t < 1)
                {
                    return "Condition: " + item.Condition +
                        "\nBiorep: " + item.Biorep +
                        "\nFraction: " + item.Fraction +
                        "\nTechrep: " + item.Techrep +
                        "\n\nTechrep must be an integer >= 1";
                }
            }

            // check for correct iteration of integer values and duplicates
            var items = SpectraFileExperimentalDesign
                .Select(p => new SpectraFileInfo(p.FullFilePathWithExtension, p.Condition, int.Parse(p.Biorep) - 1, int.Parse(p.Techrep) - 1, int.Parse(p.Fraction) - 1))
                .ToList();
            var error = ExperimentalDesign.GetErrorsInExperimentalDesign(items);

            return error;
        }

        private void KeyPressed(object sender, KeyEventArgs e)
        {
            if (e.Key == Key.Escape)
            {
                BtnCancelQuant_Click(sender, e);
            }
        }

        private void DeleteExperDesignButton_Click(object sender, RoutedEventArgs e)
        {
            if (File.Exists(outputPath))
            {
                File.Delete(outputPath);
            }

            InitializeExperDesign();
        }

        private void DgQuant_KeyDown(object sender, KeyEventArgs e)
        {

        }
    }
}