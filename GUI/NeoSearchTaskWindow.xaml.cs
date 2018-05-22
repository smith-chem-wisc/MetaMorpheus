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

namespace MetaMorpheusGUI
{
    /// <summary>
    /// Interaction logic for Window1.xaml
    /// </summary>
    public partial class NeoSearchTaskWindow : Window
    {
        #region Private Fields

        private readonly DataContextForSearchTaskWindow dataContextForSearchTaskWindow;

        private readonly ObservableCollection<SearchModeForDataGrid> SearchModesForThisTask = new ObservableCollection<SearchModeForDataGrid>();
        private readonly ObservableCollection<ModTypeForTreeView> fixedModTypeForTreeViewObservableCollection = new ObservableCollection<ModTypeForTreeView>();
        private readonly ObservableCollection<ModTypeForTreeView> variableModTypeForTreeViewObservableCollection = new ObservableCollection<ModTypeForTreeView>();
        private readonly ObservableCollection<ModTypeForTreeView> localizeModTypeForTreeViewObservableCollection = new ObservableCollection<ModTypeForTreeView>();

        #endregion Private Fields

        #region Public Constructors

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

        #endregion Public Constructors

        #region Internal Properties

        internal NeoSearchTask TheTask { get; private set; }

        #endregion Internal Properties

        #region Public Methods

        public void ReadAllTomls()
        {
        }

        public void ReadAllTomls(string filePath)
        {
        }

        #endregion Public Methods

        #region Private Methods

        private static Boolean TextBoxIntAllowed(String Text2)
        {
            return Array.TrueForAll<Char>(Text2.ToCharArray(),
                delegate (Char c) { return Char.IsDigit(c) || Char.IsControl(c); });
        }

        private void PopulateChoices()
        {
            foreach (Protease protease in GlobalVariables.ProteaseDictionary.Values)
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
            proteaseComboBox.SelectedItem = task.CommonParams.DigestionParams.Protease;

            missedCleavagesTextBox.Text = task.CommonParams.DigestionParams.MaxMissedCleavages.ToString(CultureInfo.InvariantCulture);
            txtMinPeptideLength.Text = task.CommonParams.DigestionParams.MinPeptideLength.ToString(CultureInfo.InvariantCulture);
            txtMaxPeptideLength.Text = task.CommonParams.DigestionParams.MaxPeptideLength.ToString(CultureInfo.InvariantCulture);
            proteaseComboBox.SelectedItem = task.CommonParams.DigestionParams.Protease;
            maxModificationIsoformsTextBox.Text = task.CommonParams.DigestionParams.MaxModificationIsoforms.ToString(CultureInfo.InvariantCulture);
            txtMaxModNum.Text = task.CommonParams.DigestionParams.MaxModsForPeptide.ToString(CultureInfo.InvariantCulture);
            initiatorMethionineBehaviorComboBox.SelectedIndex = (int)task.CommonParams.DigestionParams.InitiatorMethionineBehavior;
            if (task.CommonParams.ProductMassTolerance != null)
            {
                productMassToleranceTextBox.Text = task.CommonParams.ProductMassTolerance.Value.ToString(CultureInfo.InvariantCulture);
            }
            productMassToleranceComboBox.SelectedIndex = task.CommonParams.ProductMassTolerance is AbsoluteTolerance ? 0 : 1;
            if (task.CommonParams.PrecursorMassTolerance != null)
            {
                precursorMassToleranceTextBox.Text = task.CommonParams.PrecursorMassTolerance.Value.ToString(CultureInfo.InvariantCulture);
            }
            precursorMassToleranceComboBox.SelectedIndex = task.CommonParams.PrecursorMassTolerance is AbsoluteTolerance ? 0 : 1;
            bCheckBox.IsChecked = task.CommonParams.BIons;
            yCheckBox.IsChecked = task.CommonParams.YIons;
            cCheckBox.IsChecked = task.CommonParams.CIons;
            zdotCheckBox.IsChecked = task.CommonParams.ZdotIons;
            OutputFileNameTextBox.Text = task.CommonParams.TaskDescriptor;

            foreach (var mod in task.CommonParams.ListOfModsFixed)
            {
                var theModType = fixedModTypeForTreeViewObservableCollection.FirstOrDefault(b => b.DisplayName.Equals(mod.Item1));
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
                    fixedModTypeForTreeViewObservableCollection.Add(theModType);
                    theModType.Children.Add(new ModForTreeView("UNKNOWN MODIFICATION!", true, mod.Item2, true, theModType));
                }
            }
            foreach (var mod in task.CommonParams.ListOfModsVariable)
            {
                var theModType = variableModTypeForTreeViewObservableCollection.FirstOrDefault(b => b.DisplayName.Equals(mod.Item1));
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
                    variableModTypeForTreeViewObservableCollection.Add(theModType);
                    theModType.Children.Add(new ModForTreeView("UNKNOWN MODIFICATION!", true, mod.Item2, true, theModType));
                }
            }
            
            foreach (var heh in localizeModTypeForTreeViewObservableCollection)
            {
                heh.Use = true;
            }
            
            foreach (var ye in variableModTypeForTreeViewObservableCollection)
            {
                ye.VerifyCheckState();
            }
            foreach (var ye in fixedModTypeForTreeViewObservableCollection)
            {
                ye.VerifyCheckState();
            }
            foreach (var ye in localizeModTypeForTreeViewObservableCollection)
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

            #endregion Check Task Validity

            #region Save Parameters

            
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
            int MaxMissedCleavages = (int.Parse(missedCleavagesTextBox.Text, CultureInfo.InvariantCulture));
            int MinPeptideLength = (int.Parse(txtMinPeptideLength.Text, NumberStyles.Any, CultureInfo.InvariantCulture));
            int MaxPeptideLength = (int.Parse(txtMaxPeptideLength.Text, NumberStyles.Any, CultureInfo.InvariantCulture));
            int MaxModificationIsoforms = (int.Parse(maxModificationIsoformsTextBox.Text, CultureInfo.InvariantCulture));
            int MaxModsForPeptide = (int.Parse(txtMaxModNum.Text, CultureInfo.InvariantCulture));
            InitiatorMethionineBehavior InitiatorMethionineBehavior = (InitiatorMethionineBehavior)initiatorMethionineBehaviorComboBox.SelectedIndex;
            DigestionParams digestionParamsToSave = new DigestionParams(protease: protease.Name, MaxMissedCleavages: MaxMissedCleavages, MinPeptideLength: MinPeptideLength, MaxPeptideLength: MaxPeptideLength, MaxModificationIsoforms: MaxModificationIsoforms, InitiatorMethionineBehavior: InitiatorMethionineBehavior, MaxModsForPeptides: MaxModsForPeptide);

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
            double prodMassTolerance = ProductMassTolerance.Value;

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
            double preMassTolerance = PrecursorMassTolerance.Value;

            var listOfModsVariable = new List<(string, string)>();
            foreach (var heh in variableModTypeForTreeViewObservableCollection)
            {
                listOfModsVariable.AddRange(heh.Children.Where(b => b.Use).Select(b => (b.Parent.DisplayName, b.DisplayName)));
            }

            var listOfModsFixed = new List<(string, string)>();
            foreach (var heh in fixedModTypeForTreeViewObservableCollection)
            {
                listOfModsFixed.AddRange(heh.Children.Where(b => b.Use).Select(b => (b.Parent.DisplayName, b.DisplayName)));
            }

            CommonParameters CommonParamsToSave = new CommonParameters(DigestionParams: digestionParamsToSave, BIons: bCheckBox.IsChecked.Value, YIons: yCheckBox.IsChecked.Value,
                CIons: cCheckBox.IsChecked.Value,ZdotIons: zdotCheckBox.IsChecked.Value, prodMassTol: prodMassTolerance, preMassTol:preMassTolerance, ListOfModsFixed: listOfModsFixed, ListOfModsVariable: listOfModsVariable)
            {
                TaskDescriptor = (OutputFileNameTextBox.Text != "") ? OutputFileNameTextBox.Text : "NeoSearchTask"
            };


            

           
            
            
            TheTask.NeoParameters = neoParameters;
            TheTask.CommonParams = CommonParamsToSave;

            #endregion Save Parameters

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

        #endregion Private Methods
    }
}