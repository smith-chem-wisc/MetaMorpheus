using InternalLogicTaskLayer;
using OldInternalLogic;
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

            this.saveButton.Content = "Add the Calibration Task";
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
            missedCleavagesTextBox.Text = task.MaxMissedCleavages.ToString(CultureInfo.InvariantCulture);
            proteaseComboBox.SelectedItem = task.Protease;
            maxModificationIsoformsTextBox.Text = task.MaxModificationIsoforms.ToString(CultureInfo.InvariantCulture);
            initiatorMethionineBehaviorComboBox.SelectedIndex = (int)task.InitiatorMethionineBehavior;
            productMassToleranceTextBox.Text = task.ProductMassToleranceInDaltons.ToString(CultureInfo.InvariantCulture);
            precursorMassToleranceTextBox.Text = task.PrecursorMassToleranceInDaltons.ToString(CultureInfo.InvariantCulture);

            bCheckBox.IsChecked = task.BIons;
            yCheckBox.IsChecked = task.YIons;
            for (int i = 0; i < ModFileListInWindow.Count; i++)
            {
                if (task.ListOfModListsForCalibration[i].Fixed)
                    ModFileListInWindow[i].Fixed = true;
                if (task.ListOfModListsForCalibration[i].Variable)
                    ModFileListInWindow[i].Variable = true;
                if (task.ListOfModListsForCalibration[i].Localize)
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
            TheTask.MaxMissedCleavages = int.Parse(missedCleavagesTextBox.Text, CultureInfo.InvariantCulture);
            TheTask.Protease = (Protease)proteaseComboBox.SelectedItem;
            TheTask.MaxModificationIsoforms = int.Parse(maxModificationIsoformsTextBox.Text, CultureInfo.InvariantCulture);
            TheTask.InitiatorMethionineBehavior = (InitiatorMethionineBehavior)initiatorMethionineBehaviorComboBox.SelectedIndex;
            TheTask.ProductMassToleranceInDaltons = double.Parse(productMassToleranceTextBox.Text, CultureInfo.InvariantCulture);
            TheTask.BIons = bCheckBox.IsChecked.Value;
            TheTask.YIons = yCheckBox.IsChecked.Value;
            TheTask.ListOfModListsForCalibration = ModFileListInWindow.ToList();
            TheTask.PrecursorMassToleranceInDaltons = double.Parse(precursorMassToleranceTextBox.Text, CultureInfo.InvariantCulture);

            DialogResult = true;
        }

        private void Row_DoubleClick(object sender, MouseButtonEventArgs e)
        {
            var ye = sender as DataGridCell;
            var hm = ye.Content as TextBlock;
            if (hm != null && !string.IsNullOrEmpty(hm.Text))
            {
                System.Diagnostics.Process.Start(hm.Text);
            }
        }

        #endregion Private Methods

    }
}