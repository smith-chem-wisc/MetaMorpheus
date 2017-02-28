using EngineLayer;
using MzLibUtil;
using Spectra;
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
    /// Interaction logic for SearchTaskWindow.xaml
    /// </summary>
    public partial class SearchTaskWindow : Window
    {

        #region Private Fields

        // Always create a new one, even if updating an existing task
        private ObservableCollection<ModListForSearchTask> ModFileListInWindow = new ObservableCollection<ModListForSearchTask>();

        private ObservableCollection<SearchModeForDataGrid> SearchModesForThisTask = new ObservableCollection<SearchModeForDataGrid>();

        #endregion Private Fields

        #region Public Constructors

        public SearchTaskWindow()
        {
            InitializeComponent();
            PopulateChoices();

            TheTask = new SearchTask();
            UpdateFieldsFromTask(TheTask);

            this.saveButton.Content = "Add the Search Task";
        }

        public SearchTaskWindow(SearchTask task)
        {
            InitializeComponent();
            PopulateChoices();

            TheTask = task;
            UpdateFieldsFromTask(TheTask);
        }

        #endregion Public Constructors

        #region Internal Properties

        internal SearchTask TheTask { get; private set; }

        #endregion Internal Properties

        #region Private Methods

        private void Row_DoubleClick(object sender, MouseButtonEventArgs e)
        {
            var ye = sender as DataGridCell;
            var hm = ye.Content as TextBlock;
            if (hm != null && !string.IsNullOrEmpty(hm.Text))
            {
                System.Diagnostics.Process.Start(Path.Combine(@"Mods", hm.Text));
            }
        }

        private void PopulateChoices()
        {
            foreach (Protease protease in ProteaseDictionary.Instance.Values)
                proteaseComboBox.Items.Add(protease);
            proteaseComboBox.SelectedIndex = 12;

            foreach (string initiatior_methionine_behavior in Enum.GetNames(typeof(InitiatorMethionineBehavior)))
                initiatorMethionineBehaviorComboBox.Items.Add(initiatior_methionine_behavior);

            foreach (string toleranceUnit in Enum.GetNames(typeof(ToleranceUnit)))
                productMassToleranceComboBox.Items.Add(toleranceUnit);

            // Always create new ModFileList
            foreach (var uu in MetaMorpheusTask.AllModLists)
                ModFileListInWindow.Add(new ModListForSearchTask(uu));
            modificationsDataGrid.DataContext = ModFileListInWindow;

            // Always create new ModFileList
            foreach (var uu in MyEngine.SearchModesKnown)
                SearchModesForThisTask.Add(new SearchModeForDataGrid(uu));
            searchModesDataGrid.DataContext = SearchModesForThisTask;
        }

        private void UpdateFieldsFromTask(SearchTask task)
        {
            classicSearchRadioButton.IsChecked = task.ClassicSearch;
            modernSearchRadioButton.IsChecked = !task.ClassicSearch;
            checkBoxParsimony.IsChecked = task.DoParsimony;
            checkBoxNoOneHitWonders.IsChecked = task.NoOneHitWonders;
            checkBoxQuantification.IsChecked = task.Quantify;
            quantRtTolerance.Text = task.QuantifyRtTol.ToString(CultureInfo.InvariantCulture);
            quantPpmTolerance.Text = task.QuantifyPpmTol.ToString(CultureInfo.InvariantCulture);
            checkBoxHistogramAnalysis.IsChecked = task.DoHistogramAnalysis;
            checkBoxDecoy.IsChecked = task.SearchDecoy;
            missedCleavagesTextBox.Text = task.MaxMissedCleavages.ToString(CultureInfo.InvariantCulture);
            proteaseComboBox.SelectedItem = task.Protease;
            maxModificationIsoformsTextBox.Text = task.MaxModificationIsoforms.ToString(CultureInfo.InvariantCulture);
            initiatorMethionineBehaviorComboBox.SelectedIndex = (int)task.InitiatorMethionineBehavior;
            productMassToleranceTextBox.Text = task.ProductMassTolerance.Value.ToString(CultureInfo.InvariantCulture);
            productMassToleranceComboBox.SelectedIndex = (int)task.ProductMassTolerance.Unit;
            bCheckBox.IsChecked = task.BIons;
            yCheckBox.IsChecked = task.YIons;
            cCheckBox.IsChecked = task.CIons;
            zdotCheckBox.IsChecked = task.ZdotIons;

            foreach (var modList in task.ListOfModListsFixed)
                ModFileListInWindow.First(b => b.FileName.Equals(modList.FileName)).Fixed = true;
            foreach (var modList in task.ListOfModListsVariable)
                ModFileListInWindow.First(b => b.FileName.Equals(modList.FileName)).Variable = true;
            foreach (var modList in task.ListOfModListsLocalize)
                ModFileListInWindow.First(b => b.FileName.Equals(modList.FileName)).Localize = true;

            modificationsDataGrid.Items.Refresh();

            foreach (var cool in task.SearchModes)
                SearchModesForThisTask.First(b => b.searchMode.FileNameAddition.Equals(cool.FileNameAddition)).Use = true;

            searchModesDataGrid.Items.Refresh();
        }

        private void cancelButton_Click(object sender, RoutedEventArgs e)
        {
            DialogResult = false;
        }

        private void saveButton_Click(object sender, RoutedEventArgs e)
        {
            TheTask.ClassicSearch = classicSearchRadioButton.IsChecked.Value;
            TheTask.DoParsimony = checkBoxParsimony.IsChecked.Value;
            TheTask.NoOneHitWonders = checkBoxNoOneHitWonders.IsChecked.Value;
            TheTask.Quantify = checkBoxQuantification.IsChecked.Value;
            TheTask.QuantifyRtTol = double.Parse(quantRtTolerance.Text, CultureInfo.InvariantCulture);
            TheTask.QuantifyPpmTol = double.Parse(quantPpmTolerance.Text, CultureInfo.InvariantCulture);
            TheTask.SearchDecoy = checkBoxDecoy.IsChecked.Value;
            TheTask.MaxMissedCleavages = int.Parse(missedCleavagesTextBox.Text, CultureInfo.InvariantCulture);
            TheTask.Protease = (Protease)proteaseComboBox.SelectedItem;
            TheTask.MaxModificationIsoforms = int.Parse(maxModificationIsoformsTextBox.Text, CultureInfo.InvariantCulture);
            TheTask.InitiatorMethionineBehavior = (InitiatorMethionineBehavior)initiatorMethionineBehaviorComboBox.SelectedIndex;
            TheTask.ProductMassTolerance.Value = double.Parse(productMassToleranceTextBox.Text, CultureInfo.InvariantCulture);
            TheTask.ProductMassTolerance.Unit = (ToleranceUnit)productMassToleranceComboBox.SelectedIndex;
            TheTask.BIons = bCheckBox.IsChecked.Value;
            TheTask.YIons = yCheckBox.IsChecked.Value;
            TheTask.CIons = cCheckBox.IsChecked.Value;
            TheTask.ZdotIons = zdotCheckBox.IsChecked.Value;

            TheTask.ListOfModListsFixed = ModFileListInWindow.Where(b => b.Fixed).Select(b => b.ModList).ToList();
            TheTask.ListOfModListsVariable = ModFileListInWindow.Where(b => b.Variable).Select(b => b.ModList).ToList();
            TheTask.ListOfModListsLocalize = ModFileListInWindow.Where(b => b.Localize).Select(b => b.ModList).ToList();

            TheTask.SearchModes = SearchModesForThisTask.Where(b => b.Use).Select(b => b.searchMode).ToList();
            TheTask.DoHistogramAnalysis = checkBoxHistogramAnalysis.IsChecked.Value;

            DialogResult = true;
        }

        private void addNewAllowedPrecursorMassDiffsButton_Click(object sender, RoutedEventArgs e)
        {
            // Format: name, "interval", intervals
            // Format: name, "dot", num, "ppm" or "da", dots

            try
            {
                var split = newAllowedPrecursorMassDiffsTextBox.Text.Split(' ');

                switch (split[1])
                {
                    case "dot":
                        ToleranceUnit tu = ToleranceUnit.PPM;
                        if (split[3].Equals("ppm"))
                            tu = ToleranceUnit.PPM;
                        else if (split[3].Equals("da"))
                            tu = ToleranceUnit.Absolute;
                        else
                            break;
                        DotSearchMode dsm = new DotSearchMode(split[0], Array.ConvertAll(split[4].Split(','), Double.Parse), new Tolerance(tu, double.Parse(split[2], CultureInfo.InvariantCulture)));
                        MyEngine.SearchModesKnown.Add(dsm);
                        SearchModesForThisTask.Add(new SearchModeForDataGrid(dsm));
                        searchModesDataGrid.Items.Refresh();
                        break;

                    case "interval":
                        IEnumerable<DoubleRange> doubleRanges = Array.ConvertAll(split[2].Split(','), b => new DoubleRange(double.Parse(b.Trim(new char[] { '[', ']' }).Split(';')[0], CultureInfo.InvariantCulture), double.Parse(b.Trim(new char[] { '[', ']' }).Split(';')[1], CultureInfo.InvariantCulture)));
                        IntervalSearchMode ism = new IntervalSearchMode(split[0], doubleRanges);
                        MyEngine.SearchModesKnown.Add(ism);
                        SearchModesForThisTask.Add(new SearchModeForDataGrid(ism));
                        searchModesDataGrid.Items.Refresh();
                        break;
                }
            }
            catch (Exception)
            {
                MessageBox.Show("Examples:" + Environment.NewLine + "name dot 5 ppm 0,1.003,2.006" + Environment.NewLine + "name interval [-4;-3],[-0.5;0.5],[101;102]", "Error parsing search mode text box", MessageBoxButton.OK, MessageBoxImage.Warning);
            }
        }

        #endregion Private Methods

    }
}