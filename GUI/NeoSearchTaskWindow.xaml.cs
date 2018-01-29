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
        private readonly ObservableCollection<ModTypeForTreeView> fixedModTypeForTreeViewObservableCollection = new ObservableCollection<ModTypeForTreeView>();
        private readonly ObservableCollection<ModTypeForTreeView> variableModTypeForTreeViewObservableCollection = new ObservableCollection<ModTypeForTreeView>();
        private readonly ObservableCollection<ModTypeForTreeView> localizeModTypeForTreeViewObservableCollection = new ObservableCollection<ModTypeForTreeView>();

        public NeoSearchTaskWindow()
        {
            InitializeComponent();

            TheTask = new NeoSearchTask();

            PopulateChoices();

            UpdateFieldsFromTask(TheTask);

            dataContextForSearchTaskWindow = new DataContextForSearchTaskWindow
            {
                ExpanderTitle = string.Join(", ", SearchModesForThisTask.Where(b => b.Use).Select(b => b.Name))
            };
            this.DataContext = dataContextForSearchTaskWindow;
            this.saveButton.Content = "Add the Search Tasks";
        }

        public NeoSearchTaskWindow(NeoSearchTask task)
        {
            InitializeComponent();

            TheTask = task;

            UpdateFieldsFromTask(TheTask);

            dataContextForSearchTaskWindow = new DataContextForSearchTaskWindow
            {
                ExpanderTitle = string.Join(", ", SearchModesForThisTask.Where(b => b.Use).Select(b => b.Name))
            };
            this.DataContext = dataContextForSearchTaskWindow;
            this.saveButton.Content = "Add the Search Tasks";
        }

        internal NeoSearchTask TheTask { get; private set; }

        private void PopulateChoices()
        {
            foreach (Protease protease in GlobalVariables.ProteaseDictionary.Values)
                proteaseComboBox.Items.Add(protease);
            proteaseComboBox.SelectedIndex = 12;

            foreach (string initiatior_methionine_behavior in Enum.GetNames(typeof(InitiatorMethionineBehavior)))
                initiatorMethionineBehaviorComboBox.Items.Add(initiatior_methionine_behavior);

            productMassToleranceComboBox.Items.Add("Absolute");
            productMassToleranceComboBox.Items.Add("ppm");

            precursorMassToleranceComboBox.Items.Add("Absolute");
            precursorMassToleranceComboBox.Items.Add("ppm");

            foreach (var hm in GlobalVariables.AllModsKnown.GroupBy(b => b.modificationType))
            {
                var theModType = new ModTypeForTreeView(hm.Key, false);
                fixedModTypeForTreeViewObservableCollection.Add(theModType);
                foreach (var uah in hm)
                    theModType.Children.Add(new ModForTreeView(uah.ToString(), false, uah.id, false, theModType));
            }
            fixedModsTreeView.DataContext = fixedModTypeForTreeViewObservableCollection;
            foreach (var hm in GlobalVariables.AllModsKnown.GroupBy(b => b.modificationType))
            {
                var theModType = new ModTypeForTreeView(hm.Key, false);
                variableModTypeForTreeViewObservableCollection.Add(theModType);
                foreach (var uah in hm)
                    theModType.Children.Add(new ModForTreeView(uah.ToString(), false, uah.id, false, theModType));
            }
            variableModsTreeView.DataContext = variableModTypeForTreeViewObservableCollection;
            foreach (var hm in GlobalVariables.AllModsKnown.GroupBy(b => b.modificationType))
            {
                var theModType = new ModTypeForTreeView(hm.Key, false);
                localizeModTypeForTreeViewObservableCollection.Add(theModType);
                foreach (var uah in hm)
                    theModType.Children.Add(new ModForTreeView(uah.ToString(), false, uah.id, false, theModType));
            }
            localizeModsTreeView.DataContext = localizeModTypeForTreeViewObservableCollection;
        }

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
            maxMissedConsecutiveTextBox.Text = task.NeoParameters.MaxMissedConsecutiveCleavages.ToString();
            maxMissedTextBox.Text = task.NeoParameters.MaxMissedCleavages.ToString();
            maxCandidatesPerSpectrumTextBox.Text = task.NeoParameters.MaxCandidatesPerSpectrum.ToString();
            minCisLengthTextBox.Text = task.NeoParameters.MinDistanceAllowed.ToString();
            maxCisLengthTextBox.Text = task.NeoParameters.MaxDistanceAllowed.ToString();
            searchNormalCis.IsChecked = task.NeoParameters.NormalCis;
            searchReverseCis.IsChecked = task.NeoParameters.ReverseCis;
            proteaseComboBox.SelectedItem = task.CommonParameters.DigestionParams.Protease;
          
            missedCleavagesTextBox.Text = task.CommonParameters.DigestionParams.MaxMissedCleavages.ToString(CultureInfo.InvariantCulture);
            txtMinPeptideLength.Text = task.CommonParameters.DigestionParams.MinPeptideLength.HasValue ? task.CommonParameters.DigestionParams.MinPeptideLength.Value.ToString(CultureInfo.InvariantCulture) : "";
            txtMaxPeptideLength.Text = task.CommonParameters.DigestionParams.MaxPeptideLength.HasValue ? task.CommonParameters.DigestionParams.MaxPeptideLength.Value.ToString(CultureInfo.InvariantCulture) : "";
            proteaseComboBox.SelectedItem = task.CommonParameters.DigestionParams.Protease;
            maxModificationIsoformsTextBox.Text = task.CommonParameters.DigestionParams.MaxModificationIsoforms.ToString(CultureInfo.InvariantCulture);
            txtMaxModNum.Text = task.CommonParameters.DigestionParams.MaxModsForPeptide.ToString(CultureInfo.InvariantCulture);
            initiatorMethionineBehaviorComboBox.SelectedIndex = (int)task.CommonParameters.DigestionParams.InitiatorMethionineBehavior;
            if (task.CommonParameters.ProductMassTolerance != null)
                productMassToleranceTextBox.Text = task.CommonParameters.ProductMassTolerance.Value.ToString(CultureInfo.InvariantCulture);
            productMassToleranceComboBox.SelectedIndex = task.CommonParameters.ProductMassTolerance is AbsoluteTolerance ? 0 : 1;
            if (task.CommonParameters.PrecursorMassTolerance != null)
                precursorMassToleranceTextBox.Text = task.CommonParameters.PrecursorMassTolerance.Value.ToString(CultureInfo.InvariantCulture);
            precursorMassToleranceComboBox.SelectedIndex = task.CommonParameters.PrecursorMassTolerance is AbsoluteTolerance ? 0 : 1;
            bCheckBox.IsChecked = task.CommonParameters.BIons;
            yCheckBox.IsChecked = task.CommonParameters.YIons;
            cCheckBox.IsChecked = task.CommonParameters.CIons;
            zdotCheckBox.IsChecked = task.CommonParameters.ZdotIons;         
            OutputFileNameTextBox.Text = task.CommonParameters.TaskDescriptor;

            foreach (var mod in task.CommonParameters.ListOfModsFixed)
            {
                var theModType = fixedModTypeForTreeViewObservableCollection.FirstOrDefault(b => b.DisplayName.Equals(mod.Item1));
                if (theModType != null)
                {
                    var theMod = theModType.Children.FirstOrDefault(b => b.DisplayName.Equals(mod.Item2));
                    if (theMod != null)
                        theMod.Use = true;
                    else
                        theModType.Children.Add(new ModForTreeView("UNKNOWN MODIFICATION!", true, mod.Item2, true, theModType));
                }
                else
                {
                    theModType = new ModTypeForTreeView(mod.Item1, true);
                    fixedModTypeForTreeViewObservableCollection.Add(theModType);
                    theModType.Children.Add(new ModForTreeView("UNKNOWN MODIFICATION!", true, mod.Item2, true, theModType));
                }
            }
            foreach (var mod in task.CommonParameters.ListOfModsVariable)
            {
                var theModType = variableModTypeForTreeViewObservableCollection.FirstOrDefault(b => b.DisplayName.Equals(mod.Item1));
                if (theModType != null)
                {
                    var theMod = theModType.Children.FirstOrDefault(b => b.DisplayName.Equals(mod.Item2));
                    if (theMod != null)
                        theMod.Use = true;
                    else
                        theModType.Children.Add(new ModForTreeView("UNKNOWN MODIFICATION!", true, mod.Item2, true, theModType));
                }
                else
                {
                    theModType = new ModTypeForTreeView(mod.Item1, true);
                    variableModTypeForTreeViewObservableCollection.Add(theModType);
                    theModType.Children.Add(new ModForTreeView("UNKNOWN MODIFICATION!", true, mod.Item2, true, theModType));
                }
            }
            localizeAllCheckBox.IsChecked = task.CommonParameters.LocalizeAll;
            if (task.CommonParameters.LocalizeAll)
            {
                foreach (var heh in localizeModTypeForTreeViewObservableCollection)
                    heh.Use = true;
            }
            else
            {
                foreach (var mod in task.CommonParameters.ListOfModTypesLocalize)
                {
                    var theModType = localizeModTypeForTreeViewObservableCollection.FirstOrDefault(b => b.DisplayName.Equals(mod));
                    theModType.Use = true;
                }
            }
            foreach (var ye in variableModTypeForTreeViewObservableCollection)
                ye.VerifyCheckState();
            foreach (var ye in fixedModTypeForTreeViewObservableCollection)
                ye.VerifyCheckState();
            foreach (var ye in localizeModTypeForTreeViewObservableCollection)
                ye.VerifyCheckState();
        }

        private void ApmdExpander_Collapsed(object sender, RoutedEventArgs e)
        {
            dataContextForSearchTaskWindow.ExpanderTitle = string.Join(", ", SearchModesForThisTask.Where(b => b.Use).Select(b => b.Name));
        }

        private void CancelButton_Click(object sender, RoutedEventArgs e)
        {
            DialogResult = false;
        }

        private void SaveButton_Click(object sender, RoutedEventArgs e)
        {
            #region Check Task Validity
       
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
           
            if (missedCleavagesTextBox.Text.Length == 0)
            {
                MessageBox.Show("The number of missed cleavages was left empty. For no missed cleavages, please enter zero.");
                return;
            }
            if ((!double.TryParse(productMassToleranceTextBox.Text, out double pmt) || pmt <= 0)&&productMassToleranceTextBox.Text.Length!=0)
            {
                MessageBox.Show("The product mass tolerance contains unrecognized characters. \n You entered " + '"' + productMassToleranceTextBox.Text + '"' + "\n Please enter a positive number.");
                return;
            }
            if ((!double.TryParse(precursorMassToleranceTextBox.Text, out double premt) || premt <= 0)&& precursorMassToleranceTextBox.Text.Length != 0)
            {
                MessageBox.Show("The precursor mass tolerance contains unrecognized characters. \n You entered " + '"' + precursorMassToleranceTextBox.Text + '"' + "\n Please enter a positive number.");
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
            neoParameters.MaxMissedConsecutiveCleavages = int.Parse(maxMissedConsecutiveTextBox.Text);
            neoParameters.MaxMissedCleavages = int.Parse(maxMissedConsecutiveTextBox.Text);
            neoParameters.MaxCandidatesPerSpectrum = int.Parse(maxMissedConsecutiveTextBox.Text);
            neoParameters.MinDistanceAllowed = int.Parse(minCisLengthTextBox.Text);
            neoParameters.MaxDistanceAllowed = int.Parse(maxCisLengthTextBox.Text);
            neoParameters.NormalCis = searchNormalCis.IsChecked.Value;
            neoParameters.ReverseCis = searchReverseCis.IsChecked.Value;
            if (!searchTarget.IsChecked.Value)
                neoParameters.TargetFilePath = targetPath.Text;
            if (!searchDecoy.IsChecked.Value)
                neoParameters.DecoyFilePath = decoyPath.Text;
            if (!searchN.IsChecked.Value)
                neoParameters.NFilePath = NPath.Text;
            if (!searchC.IsChecked.Value)
                neoParameters.CFilePath = CPath.Text;

            DigestionParams digestionParamsToSave = new DigestionParams();
            digestionParamsToSave.MaxMissedCleavages = int.Parse(missedCleavagesTextBox.Text, CultureInfo.InvariantCulture);
            digestionParamsToSave.MinPeptideLength = int.TryParse(txtMinPeptideLength.Text, NumberStyles.Any, CultureInfo.InvariantCulture, out int temp) ? (int?)temp : null;
            digestionParamsToSave.MaxPeptideLength = int.TryParse(txtMaxPeptideLength.Text, NumberStyles.Any, CultureInfo.InvariantCulture, out temp) ? (int?)temp : null;
            digestionParamsToSave.Protease = (Protease)proteaseComboBox.SelectedItem;
            digestionParamsToSave.MaxModificationIsoforms = int.Parse(maxModificationIsoformsTextBox.Text, CultureInfo.InvariantCulture);
            digestionParamsToSave.MaxModsForPeptide = int.Parse(txtMaxModNum.Text, CultureInfo.InvariantCulture);
            digestionParamsToSave.InitiatorMethionineBehavior = (InitiatorMethionineBehavior)initiatorMethionineBehaviorComboBox.SelectedIndex;
            CommonParamsToSave.DigestionParams = digestionParamsToSave;

            CommonParamsToSave.BIons = bCheckBox.IsChecked.Value;
            CommonParamsToSave.YIons = yCheckBox.IsChecked.Value;
            CommonParamsToSave.CIons = cCheckBox.IsChecked.Value;
            CommonParamsToSave.ZdotIons = zdotCheckBox.IsChecked.Value;
            if (productMassToleranceComboBox.SelectedIndex == 0)
                CommonParamsToSave.ProductMassTolerance = new AbsoluteTolerance(double.Parse(productMassToleranceTextBox.Text, CultureInfo.InvariantCulture));
            else if (productMassToleranceTextBox.Text.Length == 0 || precursorMassToleranceTextBox.Text.Length==0)
                CommonParamsToSave.ProductMassTolerance = new PpmTolerance(25);
            else
                CommonParamsToSave.ProductMassTolerance = new PpmTolerance(double.Parse(productMassToleranceTextBox.Text, CultureInfo.InvariantCulture));
            if (precursorMassToleranceComboBox.SelectedIndex == 0)
                CommonParamsToSave.PrecursorMassTolerance = new AbsoluteTolerance(double.Parse(precursorMassToleranceTextBox.Text, CultureInfo.InvariantCulture));
            else if (productMassToleranceTextBox.Text.Length == 0 || precursorMassToleranceTextBox.Text.Length == 0)
                CommonParamsToSave.PrecursorMassTolerance = new PpmTolerance(15);
            else
                CommonParamsToSave.PrecursorMassTolerance = new PpmTolerance(double.Parse(precursorMassToleranceTextBox.Text, CultureInfo.InvariantCulture));

            var listOfModsVariable = new List<(string, string)>();
            foreach (var heh in variableModTypeForTreeViewObservableCollection)
                listOfModsVariable.AddRange(heh.Children.Where(b => b.Use).Select(b => (b.Parent.DisplayName, b.DisplayName)));
            CommonParamsToSave.ListOfModsVariable = listOfModsVariable;

            var listOfModsFixed = new List<(string, string)>();
            foreach (var heh in fixedModTypeForTreeViewObservableCollection)
                listOfModsFixed.AddRange(heh.Children.Where(b => b.Use).Select(b => (b.Parent.DisplayName, b.DisplayName)));
            CommonParamsToSave.ListOfModsFixed = listOfModsFixed;

            if (localizeAllCheckBox.IsChecked.Value)
            {
                CommonParamsToSave.ListOfModTypesLocalize = null;
                CommonParamsToSave.LocalizeAll = true;
            }
            else
            {
                CommonParamsToSave.LocalizeAll = false;
                CommonParamsToSave.ListOfModTypesLocalize = localizeModTypeForTreeViewObservableCollection.Where(b => b.Use.HasValue && b.Use.Value).Select(b => b.DisplayName).ToList();
            }

            TheTask.NeoParameters = neoParameters;
            TheTask.CommonParameters = CommonParamsToSave;
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
            {
                targetPath.Text = openPicker.FileName;
                searchTarget.IsChecked = false;
            }
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
            {
                decoyPath.Text = openPicker.FileName;
                searchDecoy.IsChecked = false;
            }
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
            {
                NPath.Text = openPicker.FileName; searchDecoy.IsChecked = false;
                searchN.IsChecked = false;
            }
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
            {
                CPath.Text = openPicker.FileName;
                searchC.IsChecked = false;
            }
        }

        private void maxCisLengthTextBox_TextChanged(object sender, TextChangedEventArgs e)
        {

        }

        public void ReadAllTomls()
        {

        }

        public void ReadAllTomls(string filePath)
        {

        }

        private void targetPath_TextChanged(object sender, TextChangedEventArgs e)
        {
            searchTarget.IsChecked = true;
        }

        private void decoyPath_TextChanged(object sender, TextChangedEventArgs e)
        {
            searchDecoy.IsChecked = true;
        }

        private void NPath_TextChanged(object sender, TextChangedEventArgs e)
        {
            searchN.IsChecked = true;
        }

        private void CPath_TextChanged(object sender, TextChangedEventArgs e)
        {
            searchC.IsChecked = true;
        }

        private void precursorMassToleranceTextBox_TextChanged(object sender, TextChangedEventArgs e)
        {
            calibrate.IsChecked = (precursorMassToleranceTextBox.Text.Length != 0 && productMassToleranceTextBox.Text.Length != 0) ? false : true;
        }

        private void productMassToleranceTextBox_TextChanged(object sender, TextChangedEventArgs e)
        {
            calibrate.IsChecked = (precursorMassToleranceTextBox.Text.Length != 0 && productMassToleranceTextBox.Text.Length != 0) ? false : true;
        }
    }
}
