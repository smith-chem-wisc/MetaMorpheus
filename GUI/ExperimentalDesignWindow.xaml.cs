using EngineLayer;
using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.IO;
using System.Linq;
using System.Windows;

namespace MetaMorpheusGUI
{
    /// <summary>
    /// Interaction logic for ExperimentalDesignWindow.xaml
    /// </summary>
    public partial class ExperimentalDesignWindow : Window
    {
        private readonly ObservableCollection<ExperimentalDesignForDataGrid> SpectraFilesQuantSets = new ObservableCollection<ExperimentalDesignForDataGrid>();
        private string OutputPath;

        public ExperimentalDesignWindow(ObservableCollection<RawDataForDataGrid> spectraFilesObservableCollection)
        {
            InitializeComponent();

            foreach (var item in spectraFilesObservableCollection.Where(p => p.Use).Select(p => p.FilePath))
            {
                SpectraFilesQuantSets.Add(new ExperimentalDesignForDataGrid(item));
            }

            if (spectraFilesObservableCollection.Any())
            {
                OutputPath = Directory.GetParent(spectraFilesObservableCollection.Where(p => p.Use).First().FilePath).FullName;
                OutputPath = Path.Combine(OutputPath, GlobalVariables.ExperimentalDesignFileName);
            }

            DgQuant.DataContext = SpectraFilesQuantSets;

            try
            {
                ReadExperDesignFromTsv(OutputPath);
            }
            catch (Exception ex)
            {
                MessageBox.Show("Could not read existing experimental design file!\n\n" + ex.Message);
            }
        }

        private void BtnSaveQuant_Click(object sender, RoutedEventArgs e)
        {
            if (OutputPath == null)
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
                WriteExperDesignToTsv(OutputPath);
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

        private void WriteExperDesignToTsv(string filePath)
        {
            using (StreamWriter output = new StreamWriter(filePath))
            {
                output.WriteLine("FileName\tCondition\tBiorep\tFraction\tTechrep");
                foreach (var spectraFile in SpectraFilesQuantSets)
                {
                    output.WriteLine(spectraFile.FileName +
                        "\t" + spectraFile.Condition +
                        "\t" + spectraFile.Biorep +
                        "\t" + spectraFile.Fraction +
                        "\t" + spectraFile.Techrep);
                }
            }
        }

        private void ReadExperDesignFromTsv(string filePath)
        {
            if (!File.Exists(filePath))
            {
                return;
            }

            var lines = File.ReadAllLines(filePath);
            Dictionary<string, int> typeToIndex = new Dictionary<string, int>();

            for (int l = 0; l < lines.Length; l++)
            {
                var split = lines[l].Split('\t');
                if (l == 0)
                {
                    foreach (var type in split)
                    {
                        typeToIndex.Add(type, Array.IndexOf(split, type));
                    }
                }
                else
                {
                    ExperimentalDesignForDataGrid file = SpectraFilesQuantSets.Where(p => p.FileName == split[typeToIndex["FileName"]]).FirstOrDefault();

                    if (file == null)
                    {
                        continue;
                    }

                    file.Condition = split[typeToIndex["Condition"]];
                    file.Biorep = split[typeToIndex["Biorep"]];
                    file.Fraction = split[typeToIndex["Fraction"]];
                    file.Techrep = split[typeToIndex["Techrep"]];
                }
            }
        }

        public string CheckForExperimentalDesignErrors()
        {
            // check for basic parsing
            foreach (var item in SpectraFilesQuantSets)
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
            var conditions = SpectraFilesQuantSets.GroupBy(p => p.Condition);

            foreach (var condition in conditions)
            {
                var temp = condition.OrderBy(p => p.Biorep).ThenBy(p => p.Fraction).ThenBy(p => p.Techrep);
                int numB = temp.Max(p => int.Parse(p.Biorep));

                // check bioreps are in order
                for (int b = 1; b <= numB; b++)
                {
                    var biorepFiles = temp.Where(p => int.Parse(p.Biorep) == b);

                    if (!biorepFiles.Any())
                    {
                        return "Condition \"" + condition.Key + "\" biorep " + b + " is missing!";
                    }

                    // check fractions are in order
                    int numF = biorepFiles.Max(p => int.Parse(p.Fraction));

                    for (int f = 1; f <= numF; f++)
                    {
                        var fractionFiles = biorepFiles.Where(p => int.Parse(p.Fraction) == f);

                        if (!fractionFiles.Any())
                        {
                            return "Condition \"" + condition.Key + "\" biorep " + b + " fraction " + f + " is missing!";
                        }

                        // check techreps are in order
                        int numT = fractionFiles.Max(p => int.Parse(p.Techrep));

                        for (int t = 1; t <= numT; t++)
                        {
                            var techrepFiles = fractionFiles.Where(p => int.Parse(p.Techrep) == t);

                            if (!techrepFiles.Any())
                            {
                                return "Condition \"" + condition.Key + "\" biorep " + b + " fraction " + f + " techrep " + t + " is missing!";
                            }

                            if (techrepFiles.Count() > 1)
                            {
                                return "Duplicates are not allowed:\n" +
                                    "Condition \"" + condition.Key + "\" biorep " + b + " fraction " + f + " techrep " + t;
                            }
                        }
                    }
                }
            }

            return null;
        }
    }
}