using InternalLogicEngineLayer;
using InternalLogicTaskLayer;
using OldInternalLogic;
using Spectra;
using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.Globalization;
using System.Linq;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Input;

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

        private ObservableCollection<SearchModeFoSearch> SearchModes = new ObservableCollection<SearchModeFoSearch>();

        #endregion Private Fields

        #region Public Constructors

        public SearchTaskWindow(IEnumerable<ModList> modList, IEnumerable<SearchMode> searchModes)
        {
            InitializeComponent();
            PopulateChoices(modList, searchModes);

            TheTask = new SearchTask(modList, searchModes);
            UpdateFieldsFromTask(TheTask);

            this.saveButton.Content = "Add the Search Task";
        }

        public SearchTaskWindow(SearchTask task, IEnumerable<ModList> modList, IEnumerable<SearchMode> searchModes)
        {
            InitializeComponent();
            PopulateChoices(modList, searchModes);

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

        private void PopulateChoices(IEnumerable<ModList> modList, IEnumerable<SearchMode> searchModes)
        {
            foreach (Protease protease in ProteaseDictionary.Instance.Values)
                proteaseComboBox.Items.Add(protease);
            proteaseComboBox.SelectedIndex = 12;

            foreach (string initiatior_methionine_behavior in Enum.GetNames(typeof(InitiatorMethionineBehavior)))
                initiatorMethionineBehaviorComboBox.Items.Add(initiatior_methionine_behavior);

            foreach (string toleranceUnit in Enum.GetNames(typeof(ToleranceUnit)))
                productMassToleranceComboBox.Items.Add(toleranceUnit);

            // Always create new ModFileList
            foreach (var uu in modList)
                ModFileListInWindow.Add(new ModListForSearchTask(uu));
            modificationsDataGrid.DataContext = ModFileListInWindow;

            // Always create new ModFileList
            foreach (var uu in searchModes)
                SearchModes.Add(new SearchModeFoSearch(uu));
            allowedPrecursorMassDiffsDataGrid.DataContext = SearchModes;
        }

        private void UpdateFieldsFromTask(SearchTask task)
        {
            classicSearchRadioButton.IsChecked = task.ClassicSearch;
            modernSearchRadioButton.IsChecked = !task.ClassicSearch;
            checkBoxParsimony.IsChecked = task.DoParsimony;
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
            for (int i = 0; i < ModFileListInWindow.Count; i++)
            {
                if (task.ListOfModListsForSearch[i].Fixed)
                    ModFileListInWindow[i].Fixed = true;
                if (task.ListOfModListsForSearch[i].Variable)
                    ModFileListInWindow[i].Variable = true;
                if (task.ListOfModListsForSearch[i].Localize)
                    ModFileListInWindow[i].Localize = true;
            }

            modificationsDataGrid.Items.Refresh();

            for (int i = 0; i < SearchModes.Count; i++)
            {
                if (task.SearchModes[i].Use)
                    SearchModes[i].Use = true;
            }
            allowedPrecursorMassDiffsDataGrid.Items.Refresh();
        }

        private void cancelButton_Click(object sender, RoutedEventArgs e)
        {
            DialogResult = false;
        }

        private void saveButton_Click(object sender, RoutedEventArgs e)
        {
            TheTask.ClassicSearch = classicSearchRadioButton.IsChecked.Value;
            TheTask.DoParsimony = checkBoxParsimony.IsChecked.Value;
            TheTask.SearchDecoy = checkBoxDecoy.IsChecked.Value;
            TheTask.MaxMissedCleavages = int.Parse(missedCleavagesTextBox.Text, CultureInfo.InvariantCulture);
            TheTask.Protease = (Protease)proteaseComboBox.SelectedItem;
            TheTask.MaxModificationIsoforms = int.Parse(maxModificationIsoformsTextBox.Text, CultureInfo.InvariantCulture);
            TheTask.InitiatorMethionineBehavior = (InitiatorMethionineBehavior)initiatorMethionineBehaviorComboBox.SelectedIndex;
            TheTask.ProductMassTolerance.Value = double.Parse(productMassToleranceTextBox.Text, CultureInfo.InvariantCulture);
            TheTask.ProductMassTolerance.Unit = (ToleranceUnit)productMassToleranceComboBox.SelectedIndex;
            TheTask.BIons = bCheckBox.IsChecked.Value;
            TheTask.YIons = yCheckBox.IsChecked.Value;
            TheTask.ListOfModListsForSearch = ModFileListInWindow.ToList();
            TheTask.SearchModes = SearchModes.ToList();
            TheTask.DoHistogramAnalysis = checkBoxHistogramAnalysis.IsChecked.Value;

            DialogResult = true;
        }

        private void addNewAllowedPrecursorMassDiffsButton_Click(object sender, RoutedEventArgs e)
        {
            // Format: name, "interval", intervals
            // Format: name, "dot", num, "ppm" or "da", dots

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
                    SearchModes.Add(new SearchModeFoSearch(dsm));
                    allowedPrecursorMassDiffsDataGrid.Items.Refresh();
                    break;

                case "interval":
                    IEnumerable<DoubleRange> doubleRanges = Array.ConvertAll(split[2].Split(','), b => new DoubleRange(double.Parse(b.Trim(new char[] { '[', ']' }).Split('-')[0], CultureInfo.InvariantCulture), double.Parse(b.Trim(new char[] { '[', ']' }).Split('-')[1], CultureInfo.InvariantCulture)));
                    IntervalSearchMode ism = new IntervalSearchMode(split[0], doubleRanges);
                    SearchModes.Add(new SearchModeFoSearch(ism));
                    allowedPrecursorMassDiffsDataGrid.Items.Refresh();
                    break;
            }
        }

        #endregion Private Methods

    }
}