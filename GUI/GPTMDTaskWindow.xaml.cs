using EngineLayer;
using MzLibUtil;
using System;
using System.Collections.ObjectModel;
using System.Globalization;
using System.Linq;
using System.Text.RegularExpressions;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Input;
using TaskLayer;

namespace MetaMorpheusGUI
{
    /// <summary>
    /// Interaction logic for GptmdTaskWindow.xaml
    /// </summary>
    public partial class GptmdTaskWindow : Window
    {

        #region Private Fields

        // Always create a new one, even if updating an existing task
        private ObservableCollection<ModListForGptmdTask> ModFileListInWindow = new ObservableCollection<ModListForGptmdTask>();

        #endregion Private Fields

        #region Public Constructors

        public GptmdTaskWindow()
        {
            InitializeComponent();
            PopulateChoices();

            TheTask = new GptmdTask();
            UpdateFieldsFromTask(TheTask);

            this.saveButton.Content = "Add the GPTMD Task";
        }

        public GptmdTaskWindow(GptmdTask myGPTMDtask)
        {
            InitializeComponent();
            PopulateChoices();

            TheTask = myGPTMDtask;
            UpdateFieldsFromTask(TheTask);
        }

        #endregion Public Constructors

        #region Internal Properties

        internal GptmdTask TheTask { get; private set; }

        #endregion Internal Properties

        #region Private Methods

        private void Row_DoubleClick(object sender, MouseButtonEventArgs e)
        {
            var ye = sender as DataGridCell;
            if (ye.Content is TextBlock hm && !string.IsNullOrEmpty(hm.Text))
            {
                System.Diagnostics.Process.Start(hm.Text);
            }
        }

        private void UpdateFieldsFromTask(GptmdTask task)
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
            cCheckBox.IsChecked = task.CIons;
            zdotCheckBox.IsChecked = task.ZdotIons;

            foreach (var modList in task.ListOfModListsFixed)
                ModFileListInWindow.First(b => b.FileName.Equals(modList)).Fixed = true;
            foreach (var modList in task.ListOfModListsVariable)
                ModFileListInWindow.First(b => b.FileName.Equals(modList)).Variable = true;
            foreach (var modList in task.ListOfModListsLocalize)
                ModFileListInWindow.First(b => b.FileName.Equals(modList)).Localize = true;
            foreach (var modList in task.ListOfModListsGptmd)
                ModFileListInWindow.First(b => b.FileName.Equals(modList)).Gptmd = true;

            modificationsDataGrid.Items.Refresh();
        }

        private void PopulateChoices()
        {
            foreach (Protease protease in GlobalTaskLevelSettings.ProteaseDictionary.Values)
                proteaseComboBox.Items.Add(protease);
            proteaseComboBox.SelectedIndex = 12;

            foreach (string initiatior_methionine_behavior in Enum.GetNames(typeof(InitiatorMethionineBehavior)))
                initiatorMethionineBehaviorComboBox.Items.Add(initiatior_methionine_behavior);

            foreach (string toleranceUnit in Enum.GetNames(typeof(ToleranceUnit)))
                productMassToleranceComboBox.Items.Add(toleranceUnit);

            foreach (string toleranceUnit in Enum.GetNames(typeof(ToleranceUnit)))
                precursorMassToleranceComboBox.Items.Add(toleranceUnit);

            // Always create new ModFileList
            //foreach (var uu in MetaMorpheusTask.AllModLists)
            //    ModFileListInWindow.Add(new ModListForGptmdTask(uu));
            modificationsDataGrid.DataContext = ModFileListInWindow;
        }

        private void CancelButton_Click(object sender, RoutedEventArgs e)
        {
            DialogResult = false;
        }

        private void SaveButton_Click(object sender, RoutedEventArgs e)
        {
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

            TheTask.ListOfModListsFixed = ModFileListInWindow.Where(b => b.Fixed).Select(b => b.FileName).ToList();
            TheTask.ListOfModListsVariable = ModFileListInWindow.Where(b => b.Variable).Select(b => b.FileName).ToList();
            TheTask.ListOfModListsLocalize = ModFileListInWindow.Where(b => b.Localize).Select(b => b.FileName).ToList();
            TheTask.ListOfModListsGptmd = ModFileListInWindow.Where(b => b.Gptmd).Select(b => b.FileName).ToList();

            TheTask.PrecursorMassTolerance.Value = double.Parse(precursorMassToleranceTextBox.Text, CultureInfo.InvariantCulture);
            TheTask.PrecursorMassTolerance.Unit = (ToleranceUnit)precursorMassToleranceComboBox.SelectedIndex;

            Regex substituteWhitespace = new Regex(@"\s+");
            TheTask.DatabaseReferencesToKeep = substituteWhitespace.Replace(databaseReferencesTextBox.Text, "").Split(',');

            DialogResult = true;
        }

        #endregion Private Methods

    }
}