using IndexSearchAndAnalyze;
using MetaMorpheus;
using Spectra;
using System;
using System.Collections.ObjectModel;
using System.Globalization;
using System.Linq;
using System.Windows;

namespace GoodGUI
{
    /// <summary>
    /// Interaction logic for CalibrateTaskWindow.xaml
    /// </summary>
    public partial class CalibrateTaskWindow : Window
    {
        // Always create a new one, even if updating an existing task
        private ObservableCollection<ModListForSearch> ModFileListInWindow = new ObservableCollection<ModListForSearch>();

        public CalibrateTaskWindow(ObservableCollection<ModList> modList)
        {
            InitializeComponent();
            PopulateChoices(modList);

            TheTask = new MyCalibrateTask(modList);
            UpdateFieldsFromTask(TheTask);
        }

        private void UpdateFieldsFromTask(MyCalibrateTask task)
        {
            missedCleavagesTextBox.Text = task.maxMissedCleavages.ToString(CultureInfo.InvariantCulture);
            proteaseComboBox.SelectedItem = task.protease;
            maxModificationIsoformsTextBox.Text = task.maxModificationIsoforms.ToString(CultureInfo.InvariantCulture);
            initiatorMethionineBehaviorComboBox.SelectedIndex = (int)task.initiatorMethionineBehavior;
            productMassToleranceTextBox.Text = task.productMassTolerance.Value.ToString(CultureInfo.InvariantCulture);
            productMassToleranceComboBox.SelectedIndex = (int)task.productMassTolerance.Unit;
            precursorMassToleranceTextBox.Text = task.precursorMassTolerance.Value.ToString(CultureInfo.InvariantCulture);
            precursorMassToleranceComboBox.SelectedIndex = (int)task.precursorMassTolerance.Unit;

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
        }

        private void PopulateChoices(ObservableCollection<ModList> modList)
        {
            foreach (Protease protease in ProteaseDictionary.Instance.Values)
                proteaseComboBox.Items.Add(protease);
            proteaseComboBox.SelectedIndex = 12;

            foreach (string initiatior_methionine_behavior in Enum.GetNames(typeof(InitiatorMethionineBehavior)))
                initiatorMethionineBehaviorComboBox.Items.Add(initiatior_methionine_behavior);

            foreach (string toleranceUnit in Enum.GetNames(typeof(ToleranceUnit)))
                productMassToleranceComboBox.Items.Add(toleranceUnit);

            foreach (string toleranceUnit in Enum.GetNames(typeof(ToleranceUnit)))
                productMassToleranceComboBox.Items.Add(toleranceUnit);

            // Always create new ModFileList
            foreach (var uu in modList)
                ModFileListInWindow.Add(new ModListForSearch(uu));
            modificationsDataGrid.DataContext = ModFileListInWindow;
        }

        public CalibrateTaskWindow(MyCalibrateTask myCalibrateTask, ObservableCollection<ModList> modList)
        {
            InitializeComponent();
            PopulateChoices(modList);

            TheTask = myCalibrateTask;
            UpdateFieldsFromTask(TheTask);
        }

        internal MyCalibrateTask TheTask { get; private set; }

        private void cancelButton_Click(object sender, RoutedEventArgs e)
        {
            DialogResult = false;
        }

        private void saveButton_Click(object sender, RoutedEventArgs e)
        {
            TheTask.maxMissedCleavages = int.Parse(missedCleavagesTextBox.Text);
            TheTask.protease = (Protease)proteaseComboBox.SelectedItem;
            TheTask.maxModificationIsoforms = int.Parse(maxModificationIsoformsTextBox.Text);
            TheTask.initiatorMethionineBehavior = (InitiatorMethionineBehavior)initiatorMethionineBehaviorComboBox.SelectedIndex;
            TheTask.productMassTolerance.Value = double.Parse(productMassToleranceTextBox.Text);
            TheTask.productMassTolerance.Unit = (ToleranceUnit)productMassToleranceComboBox.SelectedIndex;
            TheTask.bIons = bCheckBox.IsChecked.Value;
            TheTask.yIons = yCheckBox.IsChecked.Value;
            TheTask.listOfModListsForSearch = ModFileListInWindow.ToList();
            TheTask.productMassTolerance.Value = double.Parse(precursorMassToleranceTextBox.Text);
            TheTask.productMassTolerance.Unit = (ToleranceUnit)precursorMassToleranceComboBox.SelectedIndex;

            DialogResult = true;
        }
    }
}