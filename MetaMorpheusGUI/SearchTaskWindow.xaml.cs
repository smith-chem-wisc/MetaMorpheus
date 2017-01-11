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

namespace MetaMorpheusGUI
{
    /// <summary>
    /// Interaction logic for SearchTaskWindow.xaml
    /// </summary>
    public partial class SearchTaskWindow : Window
    {
        // Always create a new one, even if updating an existing task
        private ObservableCollection<ModListForSearch> ModFileListInWindow = new ObservableCollection<ModListForSearch>();

        private ObservableCollection<SearchModeFoSearch> SearchModes = new ObservableCollection<SearchModeFoSearch>();

        public SearchTaskWindow(IEnumerable<ModList> modList, IEnumerable<SearchMode> searchModes)
        {
            InitializeComponent();
            PopulateChoices(modList, searchModes);

            TheTask = new SearchTask(modList, searchModes);
            UpdateFieldsFromTask(TheTask);
        }

        public SearchTaskWindow(SearchTask task, IEnumerable<ModList> modList, IEnumerable<SearchMode> searchModes)
        {
            InitializeComponent();
            PopulateChoices(modList, searchModes);

            TheTask = task;
            UpdateFieldsFromTask(TheTask);
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
                ModFileListInWindow.Add(new ModListForSearch(uu));
            modificationsDataGrid.DataContext = ModFileListInWindow;

            // Always create new ModFileList
            foreach (var uu in searchModes)
                SearchModes.Add(new SearchModeFoSearch(uu));
            allowedPrecursorMassDiffsDataGrid.DataContext = SearchModes;
        }

        private void UpdateFieldsFromTask(SearchTask task)
        {
            classicSearchRadioButton.IsChecked = task.classicSearch;
            checkBoxParsimony.IsChecked = task.doParsimony;
            checkBoxDecoy.IsChecked = task.searchDecoy;
            missedCleavagesTextBox.Text = task.maxMissedCleavages.ToString(CultureInfo.InvariantCulture);
            proteaseComboBox.SelectedItem = task.protease;
            maxModificationIsoformsTextBox.Text = task.maxModificationIsoforms.ToString(CultureInfo.InvariantCulture);
            initiatorMethionineBehaviorComboBox.SelectedIndex = (int)task.initiatorMethionineBehavior;
            productMassToleranceTextBox.Text = task.productMassTolerance.Value.ToString(CultureInfo.InvariantCulture);
            productMassToleranceComboBox.SelectedIndex = (int)task.productMassTolerance.Unit;
            bCheckBox.IsChecked = task.bIons;
            yCheckBox.IsChecked = task.yIons;
            for (int i = 0; i < ModFileListInWindow.Count; i++)
            {
                if (task.listOfModListsForSearch[i].Fixed)
                    ModFileListInWindow[i].Fixed = true;
                if (task.listOfModListsForSearch[i].Variable)
                    ModFileListInWindow[i].Variable = true;
                if (task.listOfModListsForSearch[i].Localize)
                    ModFileListInWindow[i].Localize = true;
            }

            modificationsDataGrid.Items.Refresh();

            for (int i = 0; i < SearchModes.Count; i++)
            {
                if (task.searchModes[i].Use)
                    SearchModes[i].Use = true;
            }
            allowedPrecursorMassDiffsDataGrid.Items.Refresh();
        }

        internal SearchTask TheTask { get; private set; }

        private void cancelButton_Click(object sender, RoutedEventArgs e)
        {
            DialogResult = false;
        }

        private void saveButton_Click(object sender, RoutedEventArgs e)
        {
            TheTask.classicSearch = classicSearchRadioButton.IsChecked.Value;
            TheTask.searchDecoy = checkBoxDecoy.IsChecked.Value;
            TheTask.maxMissedCleavages = int.Parse(missedCleavagesTextBox.Text);
            TheTask.protease = (Protease)proteaseComboBox.SelectedItem;
            TheTask.maxModificationIsoforms = int.Parse(maxModificationIsoformsTextBox.Text);
            TheTask.initiatorMethionineBehavior = (InitiatorMethionineBehavior)initiatorMethionineBehaviorComboBox.SelectedIndex;
            TheTask.productMassTolerance.Value = double.Parse(productMassToleranceTextBox.Text);
            TheTask.productMassTolerance.Unit = (ToleranceUnit)productMassToleranceComboBox.SelectedIndex;
            TheTask.bIons = bCheckBox.IsChecked.Value;
            TheTask.yIons = yCheckBox.IsChecked.Value;
            TheTask.listOfModListsForSearch = ModFileListInWindow.ToList();
            TheTask.searchModes = SearchModes.ToList();

            DialogResult = true;
        }

        private void addNewAllowedPrecursorMassDiffsButton_Click(object sender, RoutedEventArgs e)
        {
            // Format: name, "interval", intervals
            // Format: name, "dot", num, "ppm" or "da", dots

            var split = newAllowedPrecursorMassDiffsTextBox.Text.Split(' ');

            try
            {
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
                        DotSearchMode dsm = new DotSearchMode(split[0], Array.ConvertAll(split[4].Split(','), Double.Parse), new Tolerance(tu, double.Parse(split[2])));
                        SearchModes.Add(new SearchModeFoSearch(dsm));
                        allowedPrecursorMassDiffsDataGrid.Items.Refresh();
                        break;

                    case "interval":
                        IEnumerable<DoubleRange> doubleRanges = Array.ConvertAll(split[2].Split(','), b => new DoubleRange(double.Parse(b.Trim(new char[] { '[', ']' }).Split('-')[0]), double.Parse(b.Trim(new char[] { '[', ']' }).Split('-')[1])));
                        IntervalSearchMode ism = new IntervalSearchMode(split[0], doubleRanges);
                        SearchModes.Add(new SearchModeFoSearch(ism));
                        allowedPrecursorMassDiffsDataGrid.Items.Refresh();
                        break;
                }
            }
            catch (Exception ee)
            {
                MessageBoxResult result = MessageBox.Show(ee.ToString(), "Error", MessageBoxButton.OK, MessageBoxImage.Warning);
            }
        }
    }
}