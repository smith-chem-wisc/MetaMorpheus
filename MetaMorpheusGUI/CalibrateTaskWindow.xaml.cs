using InternalLogicTaskLayer;
using OldInternalLogic;
using Spectra;
using System;
using System.Collections.ObjectModel;
using System.Globalization;
using System.Linq;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Input;

namespace MetaMorpheusGUI
{
    /// <summary>
    /// Interaction logic for CalibrateTaskWindow.xaml
    /// </summary>
    public partial class CalibrateTaskWindow : Window
    {
        #region Private Fields

        // Always create a new one, even if updating an existing task
        private ObservableCollection<ModListForCalibrationTask> ModFileListInWindow = new ObservableCollection<ModListForCalibrationTask>();

        #endregion Private Fields

        #region Public Constructors

        public CalibrateTaskWindow(ObservableCollection<ModList> modList)
        {
            InitializeComponent();
            PopulateChoices(modList);

            TheTask = new CalibrationTask(modList);
            UpdateFieldsFromTask(TheTask);
        }

        public CalibrateTaskWindow(CalibrationTask myCalibrateTask, ObservableCollection<ModList> modList)
        {
            InitializeComponent();
            PopulateChoices(modList);

            TheTask = myCalibrateTask;
            UpdateFieldsFromTask(TheTask);
        }

        #endregion Public Constructors

        #region Internal Properties

        internal CalibrationTask TheTask { get; private set; }

        #endregion Internal Properties

        #region Private Methods

        private void UpdateFieldsFromTask(CalibrationTask task)
        {
            missedCleavagesTextBox.Text = task.maxMissedCleavages.ToString(CultureInfo.InvariantCulture);
            proteaseComboBox.SelectedItem = task.protease;
            maxModificationIsoformsTextBox.Text = task.maxModificationIsoforms.ToString(CultureInfo.InvariantCulture);
            initiatorMethionineBehaviorComboBox.SelectedIndex = (int)task.initiatorMethionineBehavior;
            productMassToleranceTextBox.Text = task.productMassToleranceInDaltons.ToString(CultureInfo.InvariantCulture);
            precursorMassToleranceTextBox.Text = task.precursorMassTolerance.Value.ToString(CultureInfo.InvariantCulture);
            precursorMassToleranceComboBox.SelectedIndex = (int)task.precursorMassTolerance.Unit;

            bCheckBox.IsChecked = task.bIons;
            yCheckBox.IsChecked = task.yIons;
            for (int i = 0; i < ModFileListInWindow.Count; i++)
            {
                if (task.listOfModListsForCalibration[i].Fixed)
                    ModFileListInWindow[i].Fixed = true;
                if (task.listOfModListsForCalibration[i].Variable)
                    ModFileListInWindow[i].Variable = true;
                if (task.listOfModListsForCalibration[i].Localize)
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
                precursorMassToleranceComboBox.Items.Add(toleranceUnit);

            // Always create new ModFileList
            foreach (var uu in modList)
                ModFileListInWindow.Add(new ModListForCalibrationTask(uu));
            modificationsDataGrid.DataContext = ModFileListInWindow;
        }

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
            TheTask.productMassToleranceInDaltons = double.Parse(productMassToleranceTextBox.Text);
            TheTask.bIons = bCheckBox.IsChecked.Value;
            TheTask.yIons = yCheckBox.IsChecked.Value;
            TheTask.listOfModListsForCalibration = ModFileListInWindow.ToList();
            TheTask.precursorMassTolerance.Value = double.Parse(precursorMassToleranceTextBox.Text);
            TheTask.precursorMassTolerance.Unit = (ToleranceUnit)precursorMassToleranceComboBox.SelectedIndex;

            DialogResult = true;
        }

        private void Row_DoubleClick(object sender, MouseButtonEventArgs e)
        {
            var ye = sender as DataGridCell;
            var hm = ye.Content as TextBlock;
            if (hm != null && !hm.Text.Equals(""))
            {
                System.Diagnostics.Process.Start(hm.Text);
            }
        }

        #endregion Private Methods
    }
}