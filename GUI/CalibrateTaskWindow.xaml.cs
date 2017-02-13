using EngineLayer;
using MzLibUtil;
using Spectra;
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
    /// Interaction logic for CalibrateTaskWindow.xaml
    /// </summary>
    public partial class CalibrateTaskWindow : Window
    {

        #region Private Fields

        // Always create a new one, even if updating an existing task
        private ObservableCollection<ModListForCalibrationTask> ModFileListInWindow = new ObservableCollection<ModListForCalibrationTask>();

        #endregion Private Fields

        #region Public Constructors

        public CalibrateTaskWindow()
        {
            InitializeComponent();
            PopulateChoices();

            TheTask = new CalibrationTask();
            UpdateFieldsFromTask(TheTask);

            this.saveButton.Content = "Add the Calibration Task";
        }

        public CalibrateTaskWindow(CalibrationTask myCalibrateTask)
        {
            InitializeComponent();
            PopulateChoices();

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

            productMassToleranceTextBox.Text = task.ProductMassTolerance.Value.ToString(CultureInfo.InvariantCulture);
            productMassToleranceComboBox.SelectedIndex = (int)task.ProductMassTolerance.Unit;
            precursorMassToleranceTextBox.Text = task.PrecursorMassTolerance.Value.ToString(CultureInfo.InvariantCulture);
            precursorMassToleranceComboBox.SelectedIndex = (int)task.PrecursorMassTolerance.Unit;

            bCheckBox.IsChecked = task.BIons;
            yCheckBox.IsChecked = task.YIons;
            zdotCheckBox.IsChecked = task.ZdotIons;
            cCheckBox.IsChecked = task.CIons;

            foreach (var modList in task.ListOfModListsFixed)
                ModFileListInWindow.First(b => b.FileName.Equals(modList.FileName)).Fixed = true;
            foreach (var modList in task.ListOfModListsVariable)
                ModFileListInWindow.First(b => b.FileName.Equals(modList.FileName)).Variable = true;
            foreach (var modList in task.ListOfModListsLocalize)
                ModFileListInWindow.First(b => b.FileName.Equals(modList.FileName)).Localize = true;

            modificationsDataGrid.Items.Refresh();
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

            foreach (string toleranceUnit in Enum.GetNames(typeof(ToleranceUnit)))
                precursorMassToleranceComboBox.Items.Add(toleranceUnit);

            // Always create new ModFileList
            foreach (var uu in MetaMorpheusTask.AllModLists)
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

            TheTask.BIons = bCheckBox.IsChecked.Value;
            TheTask.YIons = yCheckBox.IsChecked.Value;
            TheTask.CIons = cCheckBox.IsChecked.Value;
            TheTask.ZdotIons = zdotCheckBox.IsChecked.Value;

            TheTask.ListOfModListsFixed = ModFileListInWindow.Where(b => b.Fixed).Select(b => b.ModList).ToList();
            TheTask.ListOfModListsVariable = ModFileListInWindow.Where(b => b.Variable).Select(b => b.ModList).ToList();
            TheTask.ListOfModListsLocalize = ModFileListInWindow.Where(b => b.Localize).Select(b => b.ModList).ToList();

            TheTask.ProductMassTolerance.Value = double.Parse(productMassToleranceTextBox.Text, CultureInfo.InvariantCulture);
            TheTask.ProductMassTolerance.Unit = (ToleranceUnit)productMassToleranceComboBox.SelectedIndex;

            TheTask.PrecursorMassTolerance.Value = double.Parse(precursorMassToleranceTextBox.Text, CultureInfo.InvariantCulture);
            TheTask.PrecursorMassTolerance.Unit = (ToleranceUnit)precursorMassToleranceComboBox.SelectedIndex;

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