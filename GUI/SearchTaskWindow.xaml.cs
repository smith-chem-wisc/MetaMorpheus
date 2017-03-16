using EngineLayer;
using MzLibUtil;
using System;
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
                System.Diagnostics.Process.Start(hm.Text);
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
            conserveMemoryCheckBox.IsChecked = task.ConserveMemory;

            foreach (var modList in task.ListOfModListsFixed)
                ModFileListInWindow.First(b => b.FileName.Equals(modList)).Fixed = true;
            foreach (var modList in task.ListOfModListsVariable)
                ModFileListInWindow.First(b => b.FileName.Equals(modList)).Variable = true;
            foreach (var modList in task.ListOfModListsLocalize)
                ModFileListInWindow.First(b => b.FileName.Equals(modList)).Localize = true;

            modificationsDataGrid.Items.Refresh();

            foreach (var cool in task.SearchModes)
                SearchModesForThisTask.First(b => b.searchMode.FileNameAddition.Equals(cool.FileNameAddition)).Use = true;

            searchModesDataGrid.Items.Refresh();
        }

        private void CancelButton_Click(object sender, RoutedEventArgs e)
        {
            DialogResult = false;
        }

        private void SaveButton_Click(object sender, RoutedEventArgs e)
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
            TheTask.ConserveMemory = conserveMemoryCheckBox.IsChecked.Value;

            TheTask.ListOfModListsFixed = ModFileListInWindow.Where(b => b.Fixed).Select(b => b.FileName).ToList();
            TheTask.ListOfModListsVariable = ModFileListInWindow.Where(b => b.Variable).Select(b => b.FileName).ToList();
            TheTask.ListOfModListsLocalize = ModFileListInWindow.Where(b => b.Localize).Select(b => b.FileName).ToList();

            TheTask.SearchModes = SearchModesForThisTask.Where(b => b.Use).Select(b => b.searchMode).ToList();
            TheTask.DoHistogramAnalysis = checkBoxHistogramAnalysis.IsChecked.Value;

            DialogResult = true;
        }

        private void AddNewAllowedPrecursorMassDiffsButton_Click(object sender, RoutedEventArgs e)
        {
            try
            {
                var ye = MetaMorpheusTask.ParseSearchMode(newAllowedPrecursorMassDiffsTextBox.Text);
                MyEngine.SearchModesKnown.Add(ye);
                SearchModesForThisTask.Add(new SearchModeForDataGrid(ye));
                searchModesDataGrid.Items.Refresh();
            }
            catch
            {
                MessageBox.Show("Examples:" + Environment.NewLine + "name dot 5 ppm 0,1.003,2.006" + Environment.NewLine + "name interval [-4;-3],[-0.5;0.5],[101;102]", "Error parsing search mode text box", MessageBoxButton.OK, MessageBoxImage.Warning);
            }
        }

        #endregion Private Methods

    }
}