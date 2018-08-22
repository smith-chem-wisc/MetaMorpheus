using EngineLayer;
using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.Globalization;
using System.Linq;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Input;
using TaskLayer;
using Proteomics.ProteolyticDigestion;

namespace MetaMorpheusGUI
{
    /// <summary>
    /// Interaction logic for Window1.xaml
    /// </summary>
    public partial class NeoSearchTaskWindow : Window
    {
        private readonly DataContextForSearchTaskWindow dataContextForSearchTaskWindow;

        private readonly ObservableCollection<SearchModeForDataGrid> SearchModesForThisTask = new ObservableCollection<SearchModeForDataGrid>();
        private readonly ObservableCollection<ModTypeForTreeView> FixedModTypeForTreeViewObservableCollection = new ObservableCollection<ModTypeForTreeView>();
        private readonly ObservableCollection<ModTypeForTreeView> VariableModTypeForTreeViewObservableCollection = new ObservableCollection<ModTypeForTreeView>();
        private readonly ObservableCollection<ModTypeForTreeView> LocalizeModTypeForTreeViewObservableCollection = new ObservableCollection<ModTypeForTreeView>();

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

        public void ReadAllTomls()
        {
        }

        public void ReadAllTomls(string filePath)
        {
        }

        private static Boolean TextBoxIntAllowed(String Text2)
        {
            return Array.TrueForAll<Char>(Text2.ToCharArray(),
                delegate (Char c) { return Char.IsDigit(c) || Char.IsControl(c); });
        }

        private void PopulateChoices()
        {
            foreach (Protease protease in ProteaseDictionary.Dictionary.Values)
                proteaseComboBox.Items.Add(protease);
            proteaseComboBox.SelectedIndex = 12;

            foreach (string initiatior_methionine_behavior in Enum.GetNames(typeof(InitiatorMethionineBehavior)))
                initiatorMethionineBehaviorComboBox.Items.Add(initiatior_methionine_behavior);

            productMassToleranceComboBox.Items.Add("Da");
            productMassToleranceComboBox.Items.Add("ppm");

            precursorMassToleranceComboBox.Items.Add("Da");
            precursorMassToleranceComboBox.Items.Add("ppm");

            foreach (var hm in GlobalVariables.AllModsKnown.GroupBy(b => b.modificationType))
            {
                var theModType = new ModTypeForTreeView(hm.Key, false);
                FixedModTypeForTreeViewObservableCollection.Add(theModType);
                foreach (var uah in hm)
                    theModType.Children.Add(new ModForTreeView(uah.ToString(), false, uah.id, false, theModType));
            }
            fixedModsTreeView.DataContext = FixedModTypeForTreeViewObservableCollection;
            foreach (var hm in GlobalVariables.AllModsKnown.GroupBy(b => b.modificationType))
            {
                var theModType = new ModTypeForTreeView(hm.Key, false);
                VariableModTypeForTreeViewObservableCollection.Add(theModType);
                foreach (var uah in hm)
                    theModType.Children.Add(new ModForTreeView(uah.ToString(), false, uah.id, false, theModType));
            }
            variableModsTreeView.DataContext = VariableModTypeForTreeViewObservableCollection;
            foreach (var hm in GlobalVariables.AllModsKnown.GroupBy(b => b.modificationType))
            {
                var theModType = new ModTypeForTreeView(hm.Key, false);
                LocalizeModTypeForTreeViewObservableCollection.Add(theModType);
                foreach (var uah in hm)
                    theModType.Children.Add(new ModForTreeView(uah.ToString(), false, uah.id, false, theModType));
            }
            localizeModsTreeView.DataContext = LocalizeModTypeForTreeViewObservableCollection;
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
            txtMinPeptideLength.Text = task.CommonParameters.DigestionParams.MinPeptideLength.ToString(CultureInfo.InvariantCulture);
            txtMaxPeptideLength.Text = task.CommonParameters.DigestionParams.MaxPeptideLength.ToString(CultureInfo.InvariantCulture);
            proteaseComboBox.SelectedItem = task.CommonParameters.DigestionParams.Protease;
            maxModificationIsoformsTextBox.Text = task.CommonParameters.DigestionParams.MaxModificationIsoforms.ToString(CultureInfo.InvariantCulture);
            txtMaxModNum.Text = task.CommonParameters.DigestionParams.MaxModsForPeptide.ToString(CultureInfo.InvariantCulture);
            initiatorMethionineBehaviorComboBox.SelectedIndex = (int)task.CommonParameters.DigestionParams.InitiatorMethionineBehavior;
            if (task.CommonParameters.ProductMassTolerance != null)
            {
                productMassToleranceTextBox.Text = task.CommonParameters.ProductMassTolerance.Value.ToString(CultureInfo.InvariantCulture);
            }
            productMassToleranceComboBox.SelectedIndex = task.CommonParameters.ProductMassTolerance is AbsoluteTolerance ? 0 : 1;
            if (task.CommonParameters.PrecursorMassTolerance != null)
            {
                precursorMassToleranceTextBox.Text = task.CommonParameters.PrecursorMassTolerance.Value.ToString(CultureInfo.InvariantCulture);
            }
            precursorMassToleranceComboBox.SelectedIndex = task.CommonParameters.PrecursorMassTolerance is AbsoluteTolerance ? 0 : 1;
            bCheckBox.IsChecked = task.CommonParameters.BIons;
            yCheckBox.IsChecked = task.CommonParameters.YIons;
            cCheckBox.IsChecked = task.CommonParameters.CIons;
            zdotCheckBox.IsChecked = task.CommonParameters.ZdotIons;
            OutputFileNameTextBox.Text = task.CommonParameters.TaskDescriptor;

            foreach (var mod in task.CommonParameters.ListOfModsFixed)
            {
                var theModType = FixedModTypeForTreeViewObservableCollection.FirstOrDefault(b => b.DisplayName.Equals(mod.Item1));
                if (theModType != null)
                {
                    var theMod = theModType.Children.FirstOrDefault(b => b.DisplayName.Equals(mod.Item2));
                    if (theMod != null)
                    {
                        theMod.Use = true;
                    }
                    else
                    {
                        theModType.Children.Add(new ModForTreeView("UNKNOWN MODIFICATION!", true, mod.Item2, true, theModType));
                    }
                }
                else
                {
                    theModType = new ModTypeForTreeView(mod.Item1, true);
                    FixedModTypeForTreeViewObservableCollection.Add(theModType);
                    theModType.Children.Add(new ModForTreeView("UNKNOWN MODIFICATION!", true, mod.Item2, true, theModType));
                }
            }
            foreach (var mod in task.CommonParameters.ListOfModsVariable)
            {
                var theModType = VariableModTypeForTreeViewObservableCollection.FirstOrDefault(b => b.DisplayName.Equals(mod.Item1));
                if (theModType != null)
                {
                    var theMod = theModType.Children.FirstOrDefault(b => b.DisplayName.Equals(mod.Item2));
                    if (theMod != null)
                    {
                        theMod.Use = true;
                    }
                    else
                    {
                        theModType.Children.Add(new ModForTreeView("UNKNOWN MODIFICATION!", true, mod.Item2, true, theModType));
                    }
                }
                else
                {
                    theModType = new ModTypeForTreeView(mod.Item1, true);
                    VariableModTypeForTreeViewObservableCollection.Add(theModType);
                    theModType.Children.Add(new ModForTreeView("UNKNOWN MODIFICATION!", true, mod.Item2, true, theModType));
                }
            }

            foreach (var heh in LocalizeModTypeForTreeViewObservableCollection)
            {
                heh.Use = true;
            }

            foreach (var ye in VariableModTypeForTreeViewObservableCollection)
            {
                ye.VerifyCheckState();
            }
            foreach (var ye in FixedModTypeForTreeViewObservableCollection)
            {
                ye.VerifyCheckState();
            }
            foreach (var ye in LocalizeModTypeForTreeViewObservableCollection)
            {
                ye.VerifyCheckState();
            }
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
            if ((!double.TryParse(productMassToleranceTextBox.Text, out double pmt) || pmt <= 0) && productMassToleranceTextBox.Text.Length != 0)
            {
                MessageBox.Show("The product mass tolerance contains unrecognized characters. \n You entered " + '"' + productMassToleranceTextBox.Text + '"' + "\n Please enter a positive number.");
                return;
            }
            if ((!double.TryParse(precursorMassToleranceTextBox.Text, out double premt) || premt <= 0) && precursorMassToleranceTextBox.Text.Length != 0)
            {
                MessageBox.Show("The precursor mass tolerance contains unrecognized characters. \n You entered " + '"' + precursorMassToleranceTextBox.Text + '"' + "\n Please enter a positive number.");
                return;
            }

            //Code for determining SemiSpecific
            NeoParameters neoParameters = new NeoParameters
            {
                Calibrate = calibrate.IsChecked.Value,
                GPTMD = gptmd.IsChecked.Value,
                TargetSearch = searchTarget.IsChecked.Value,
                DecoySearch = searchDecoy.IsChecked.Value,
                SearchNTerminus = searchN.IsChecked.Value,
                SearchCTerminus = searchC.IsChecked.Value,
                MaxMissedConsecutiveCleavages = int.Parse(maxMissedConsecutiveTextBox.Text),
                MaxMissedCleavages = int.Parse(maxMissedConsecutiveTextBox.Text),
                MaxCandidatesPerSpectrum = int.Parse(maxMissedConsecutiveTextBox.Text),
                MinDistanceAllowed = int.Parse(minCisLengthTextBox.Text),
                MaxDistanceAllowed = int.Parse(maxCisLengthTextBox.Text),
                NormalCis = searchNormalCis.IsChecked.Value,
                ReverseCis = searchReverseCis.IsChecked.Value
            };

            if (!searchTarget.IsChecked.Value)
                neoParameters.TargetFilePath = targetPath.Text;
            if (!searchDecoy.IsChecked.Value)
                neoParameters.DecoyFilePath = decoyPath.Text;
            if (!searchN.IsChecked.Value)
                neoParameters.NFilePath = NPath.Text;
            if (!searchC.IsChecked.Value)
                neoParameters.CFilePath = CPath.Text;

            Protease protease = (Protease)proteaseComboBox.SelectedItem;
            int maxMissedCleavages = (int.Parse(missedCleavagesTextBox.Text, CultureInfo.InvariantCulture));
            int minPeptideLength = (int.Parse(txtMinPeptideLength.Text, NumberStyles.Any, CultureInfo.InvariantCulture));
            int maxPeptideLength = (int.Parse(txtMaxPeptideLength.Text, NumberStyles.Any, CultureInfo.InvariantCulture));
            int maxModificationIsoforms = (int.Parse(maxModificationIsoformsTextBox.Text, CultureInfo.InvariantCulture));
            int maxModsForPeptide = (int.Parse(txtMaxModNum.Text, CultureInfo.InvariantCulture));
            InitiatorMethionineBehavior initiatorMethionineBehavior = (InitiatorMethionineBehavior)initiatorMethionineBehaviorComboBox.SelectedIndex;
            DigestionParams digestionParamsToSave = new DigestionParams(
                protease: protease.Name,
                maxMissedCleavages: maxMissedCleavages,
                minPeptideLength: minPeptideLength, 
                maxPeptideLength: maxPeptideLength, 
                maxModificationIsoforms: maxModificationIsoforms, 
                initiatorMethionineBehavior: initiatorMethionineBehavior, 
                maxModsForPeptides: maxModsForPeptide);

            Tolerance ProductMassTolerance;
            if (productMassToleranceComboBox.SelectedIndex == 0)
            {
                ProductMassTolerance = new AbsoluteTolerance(double.Parse(productMassToleranceTextBox.Text, CultureInfo.InvariantCulture));
            }
            else if (productMassToleranceTextBox.Text.Length == 0 || precursorMassToleranceTextBox.Text.Length == 0)
            {
                ProductMassTolerance = new PpmTolerance(25);
            }
            else
            {
                ProductMassTolerance = new PpmTolerance(double.Parse(productMassToleranceTextBox.Text, CultureInfo.InvariantCulture));
            }

            Tolerance PrecursorMassTolerance;
            if (precursorMassToleranceComboBox.SelectedIndex == 0)
            {
                PrecursorMassTolerance = new AbsoluteTolerance(double.Parse(precursorMassToleranceTextBox.Text, CultureInfo.InvariantCulture));
            }
            else if (productMassToleranceTextBox.Text.Length == 0 || precursorMassToleranceTextBox.Text.Length == 0)
            {
                PrecursorMassTolerance = new PpmTolerance(15);
            }
            else
            {
                PrecursorMassTolerance = new PpmTolerance(double.Parse(precursorMassToleranceTextBox.Text, CultureInfo.InvariantCulture));
            }

            var listOfModsVariable = new List<(string, string)>();
            foreach (var heh in VariableModTypeForTreeViewObservableCollection)
            {
                listOfModsVariable.AddRange(heh.Children.Where(b => b.Use).Select(b => (b.Parent.DisplayName, b.DisplayName)));
            }

            var listOfModsFixed = new List<(string, string)>();
            foreach (var heh in FixedModTypeForTreeViewObservableCollection)
            {
                listOfModsFixed.AddRange(heh.Children.Where(b => b.Use).Select(b => (b.Parent.DisplayName, b.DisplayName)));
            }

            CommonParameters commonParamsToSave = new CommonParameters(
                taskDescriptor: OutputFileNameTextBox.Text != "" ? OutputFileNameTextBox.Text : "NeoSearchTask",
                digestionParams: digestionParamsToSave,
                bIons: bCheckBox.IsChecked.Value,
                yIons: yCheckBox.IsChecked.Value,
                cIons: cCheckBox.IsChecked.Value,
                zDotIons: zdotCheckBox.IsChecked.Value,
                productMassTolerance: ProductMassTolerance,
                precursorMassTolerance: PrecursorMassTolerance,
                listOfModsFixed: listOfModsFixed,
                listOfModsVariable: listOfModsVariable);

            TheTask.NeoParameters = neoParameters;
            TheTask.CommonParameters = commonParamsToSave;

            DialogResult = true;
        }

        private void PreviewIfInt(object sender, TextCompositionEventArgs e)
        {
            e.Handled = !TextBoxIntAllowed(e.Text);
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