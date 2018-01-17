using EngineLayer;
using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.ComponentModel;
using System.Globalization;
using System.Linq;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Input;
using TaskLayer;
using UsefulProteomicsDatabases;

namespace MetaMorpheusGUI
{
    /// <summary>
    /// Interaction logic for Window1.xaml
    /// </summary>
    public partial class NeoSearchTaskWindow : Window
    {
        private readonly DataContextForSearchTaskWindow dataContextForSearchTaskWindow;
        private readonly ObservableCollection<SearchModeForDataGrid> SearchModesForThisTask = new ObservableCollection<SearchModeForDataGrid>();


        public NeoSearchTaskWindow()
        {
            InitializeComponent();

            TheTask = new NeoSearchTask();

            UpdateFieldsFromTask(TheTask);

            dataContextForSearchTaskWindow = new DataContextForSearchTaskWindow
            {
                ExpanderTitle = string.Join(", ", SearchModesForThisTask.Where(b => b.Use).Select(b => b.Name)),
                //ModExpanderTitle =
                //"fixed: "
                //+ string.Join(",", ModFileListInWindow.Where(b => b.Fixed).Select(b => b.FileName))
                //+ " variable: "
                //+ string.Join(",", ModFileListInWindow.Where(b => b.Variable).Select(b => b.FileName))
                //+ " localize: "
                //+ string.Join(",", ModFileListInWindow.Where(b => b.Localize).Select(b => b.FileName)),
                AnalysisExpanderTitle = "Some analysis properties...",
                SearchModeExpanderTitle = "Some search properties..."
            };
            this.DataContext = dataContextForSearchTaskWindow;
            this.saveButton.Content = "Add the Search Tasks";
        }

        internal NeoSearchTask TheTask { get; private set; }
        internal CalibrationTask CalibrationTask { get; private set; }
        internal GptmdTask GptmdTask { get; private set; }
        internal SearchTask TargetTask { get; private set; }
        internal SearchTask DecoyTask { get; private set; }
        internal SearchTask NTerminusTask { get; private set; }
        internal SearchTask CTerminusTask { get; private set; }


        private void UpdateFieldsFromTask(NeoSearchTask task)
        {
            calibrate.IsChecked = task.NeoParameters.Calibrate;
            gptmd.IsChecked = task.NeoParameters.GPTMD;
            searchTarget.IsChecked = task.NeoParameters.TargetSearch;
            targetPath.Text = task.NeoParameters.TargetFilePath != null ? task.NeoParameters.TargetFilePath : "";
            searchDecoy.IsChecked = task.NeoParameters.DecoySearch;
            decoyPath.Text = task.NeoParameters.TargetFilePath != null ? task.NeoParameters.DecoyFilePath : "";
            searchN.IsChecked = task.NeoParameters.SearchNTerminus;
            NPath.Text = task.NeoParameters.NFilePath != null ? task.NeoParameters.NFilePath : "";
            searchC.IsChecked = task.NeoParameters.SearchCTerminus;
            CPath.Text = task.NeoParameters.CFilePath != null ? task.NeoParameters.CFilePath : "";
            bCheckBox.IsChecked = task.NeoParameters.bIons;
            yCheckBox.IsChecked = task.NeoParameters.yIons;
            cCheckBox.IsChecked = task.NeoParameters.cIons;
            zdotCheckBox.IsChecked = task.NeoParameters.zdotIons;
            maxMissedConsecutiveTextBox.Text = task.NeoParameters.MaxMissedConsecutiveCleavages.ToString();
            maxMissedTextBox.Text = task.NeoParameters.MaxMissedCleavages.ToString();
            maxCandidatesPerSpectrumTextBox.Text = task.NeoParameters.MaxCandidatesPerSpectrum.ToString();
            searchNormalCis.IsChecked = task.NeoParameters.NormalCis;
            searchReverseCis.IsChecked = task.NeoParameters.ReverseCis;
        }

        private void ApmdExpander_Collapsed(object sender, RoutedEventArgs e)
        {
            dataContextForSearchTaskWindow.ExpanderTitle = string.Join(", ", SearchModesForThisTask.Where(b => b.Use).Select(b => b.Name));
            //dataContextForSearchTaskWindow.ModExpanderTitle =
            //    "fixed: "
            //    + string.Join(",", ModFileListInWindow.Where(b => b.Fixed).Select(b => b.FileName))
            //    + " variable: "
            //    + string.Join(",", ModFileListInWindow.Where(b => b.Variable).Select(b => b.FileName))
            //    + " localize: "
            //    + string.Join(",", ModFileListInWindow.Where(b => b.Localize).Select(b => b.FileName));
            dataContextForSearchTaskWindow.AnalysisExpanderTitle = "Some analysis properties...";
            dataContextForSearchTaskWindow.SearchModeExpanderTitle = "Some search properties...";
        }

        private void CancelButton_Click(object sender, RoutedEventArgs e)
        {
            DialogResult = false;
        }

        private void SaveButton_Click(object sender, RoutedEventArgs e)
        {
            #region Check Task Validity
       
            if (precursorTextBox.Text.Length != 0 && (!double.TryParse(precursorTextBox.Text, out double pre) || pre <= 0))
            {
                MessageBox.Show("The precursor tolerance contains unrecognized characters. \n You entered " + '"' + precursorTextBox.Text + '"' + "\n Please enter a positive number.");
                return;
            }
            if (productTextBox.Text.Length != 0 && (!double.TryParse(productTextBox.Text, out double pro) || pro <= 0))
            {
                MessageBox.Show("The product tolerance contains unrecognized characters. \n You entered " + '"' + productTextBox.Text + '"' + "\n Please enter a positive number.");
                return;
            }
            if (!int.TryParse(maxCandidatesPerSpectrumTextBox.Text, out int mcps) || mcps <= 0)
            {
                MessageBox.Show("The number of maximum candidates per spectra contains unrecognized characters. \n You entered " + '"' + maxCandidatesPerSpectrumTextBox.Text + '"' + "\n Please enter a positive number.");
                return;
            }
            if (!int.TryParse(maxMissedConsecutiveTextBox.Text, out int mmc) || mmc < 0)
            {
                MessageBox.Show("The number of maximum missed consecutive cleavages contains unrecognized characters. \n You entered " + '"' + maxMissedConsecutiveTextBox.Text + '"' + "\n Please enter a positive number.");
                return;
            }
            if (!int.TryParse(maxMissedTextBox.Text, out int mm) || mm <= 0)
            {
                MessageBox.Show("The number of maximum missed cleavages contains unrecognized characters. \n You entered " + '"' + maxMissedTextBox.Text + '"' + "\n Please enter a positive number.");
                return;
            }          

            #endregion Check Task Validity

            #region Save Parameters

            CommonParameters CommonParamsToSave = new CommonParameters();

            if (OutputFileNameTextBox.Text != "")
                CommonParamsToSave.TaskDescriptor = OutputFileNameTextBox.Text;
            else
                CommonParamsToSave.TaskDescriptor = "NeoSearchTask";

            //Code for determining SemiSpecific
            NeoParameters neoParameters = new NeoParameters();

            neoParameters.Calibrate = calibrate.IsChecked.Value;
            neoParameters.GPTMD = gptmd.IsChecked.Value;
            neoParameters.TargetSearch = searchTarget.IsChecked.Value;
            neoParameters.DecoySearch = searchDecoy.IsChecked.Value;
            neoParameters.SearchNTerminus = searchN.IsChecked.Value;
            neoParameters.SearchCTerminus = searchC.IsChecked.Value;
            neoParameters.bIons = bCheckBox.IsChecked.Value;
            neoParameters.yIons = yCheckBox.IsChecked.Value;
            neoParameters.cIons = cCheckBox.IsChecked.Value;
            neoParameters.zdotIons = zdotCheckBox.IsChecked.Value;
            neoParameters.MaxMissedConsecutiveCleavages = int.Parse(maxMissedConsecutiveTextBox.Text);
            neoParameters.MaxMissedCleavages = int.Parse(maxMissedConsecutiveTextBox.Text);
            neoParameters.MaxCandidatesPerSpectrum = int.Parse(maxMissedConsecutiveTextBox.Text);
            neoParameters.NormalCis = searchNormalCis.IsChecked.Value;
            neoParameters.ReverseCis = searchReverseCis.IsChecked.Value;

            #endregion Save Parameters

            DialogResult = true;
        }


        private void PreviewIfInt(object sender, TextCompositionEventArgs e)
        {
            e.Handled = !TextBoxIntAllowed(e.Text);
        }

        private static Boolean TextBoxIntAllowed(String Text2)
        {
            return Array.TrueForAll<Char>(Text2.ToCharArray(),
                delegate (Char c) { return Char.IsDigit(c) || Char.IsControl(c); });
        }

        private void addTargetSearch_Click(object sender, RoutedEventArgs e)
        {
            Microsoft.Win32.OpenFileDialog openPicker = new Microsoft.Win32.OpenFileDialog()
            {
                Filter = "Database Files|*.psmtsv",
                FilterIndex = 1,
                RestoreDirectory = true,
                Multiselect = false
            };
            if (openPicker.ShowDialog() == true)
                targetPath.Text = openPicker.FileName;
        }

        private void addDecoySearch_Click(object sender, RoutedEventArgs e)
        {
            Microsoft.Win32.OpenFileDialog openPicker = new Microsoft.Win32.OpenFileDialog()
            {
                Filter = "Database Files|*.psmtsv",
                FilterIndex = 1,
                RestoreDirectory = true,
                Multiselect = false
            };
            if (openPicker.ShowDialog() == true)
                decoyPath.Text = openPicker.FileName;
        }

        private void addN_Click(object sender, RoutedEventArgs e)
        {
            Microsoft.Win32.OpenFileDialog openPicker = new Microsoft.Win32.OpenFileDialog()
            {
                Filter = "Database Files|*.psmtsv",
                FilterIndex = 1,
                RestoreDirectory = true,
                Multiselect = false
            };
            if (openPicker.ShowDialog() == true)
                NPath.Text = openPicker.FileName;
        }

        private void addC_Click(object sender, RoutedEventArgs e)
        {
            Microsoft.Win32.OpenFileDialog openPicker = new Microsoft.Win32.OpenFileDialog()
            {
                Filter = "Database Files|*.psmtsv",
                FilterIndex = 1,
                RestoreDirectory = true,
                Multiselect = false
            };
            if (openPicker.ShowDialog() == true)
                CPath.Text = openPicker.FileName;
        }

        private void maxCisLengthTextBox_TextChanged(object sender, TextChangedEventArgs e)
        {

        }
    }
}
