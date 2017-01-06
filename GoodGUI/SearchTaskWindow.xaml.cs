using MetaMorpheus;
using Spectra;
using System.Globalization;
using System.Windows;
using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;

namespace GoodGUI
{
    /// <summary>
    /// Interaction logic for SearchTaskWindow.xaml
    /// </summary>
    public partial class SearchTaskWindow : Window
    {
        private ObservableCollection<ModListForSearch> ModFileList = new ObservableCollection<ModListForSearch>();

        public SearchTaskWindow(IEnumerable<ModList> modList)
        {
            InitializeComponent();
            PopulateChoices();

            foreach (var uu in modList)
                ModFileList.Add(new ModListForSearch(uu));
            modificationsDataGrid.DataContext = ModFileList;

            TheTask = new MySearchTask(ModFileList);
            UpdateFields(TheTask);
        }

        public SearchTaskWindow(MySearchTask task)
        {
            InitializeComponent();
            PopulateChoices();

            modificationsDataGrid.DataContext = task.modFileList;

            TheTask = task;
            UpdateFields(TheTask);
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
        }

        private void UpdateFields(MySearchTask task)
        {
            checkBoxDecoy.IsChecked = task.searchDecoy;
            missedCleavagesTextBox.Text = task.maxMissedCleavages.ToString(CultureInfo.InvariantCulture);
            proteaseComboBox.SelectedItem = task.protease;
            maxModificationIsoformsTextBox.Text = task.maxModificationIsoforms.ToString(CultureInfo.InvariantCulture);
            initiatorMethionineBehaviorComboBox.SelectedIndex = (int)task.initiatorMethionineBehavior;
            productMassToleranceTextBox.Text = task.productMassTolerance.Value.ToString(CultureInfo.InvariantCulture);
            productMassToleranceComboBox.SelectedIndex = (int)task.productMassTolerance.Unit;
            bCheckBox.IsChecked = task.bIons;
            yCheckBox.IsChecked = task.yIons;
            modificationsDataGrid.Items.Refresh();
        }

        internal MySearchTask TheTask { get; private set; }

        private void cancelButton_Click(object sender, RoutedEventArgs e)
        {
            DialogResult = false;
        }
        private void saveButton_Click(object sender, RoutedEventArgs e)
        {
            TheTask.searchDecoy = checkBoxDecoy.IsChecked.Value;
            TheTask.maxMissedCleavages = int.Parse(missedCleavagesTextBox.Text);
            TheTask.protease = (Protease)proteaseComboBox.SelectedItem;
            TheTask.maxModificationIsoforms = int.Parse(maxModificationIsoformsTextBox.Text);
            TheTask.initiatorMethionineBehavior = (InitiatorMethionineBehavior)initiatorMethionineBehaviorComboBox.SelectedIndex;
            TheTask.productMassTolerance.Value = double.Parse(productMassToleranceTextBox.Text);
            TheTask.productMassTolerance.Unit = (ToleranceUnit)productMassToleranceComboBox.SelectedIndex;
            TheTask.bIons = bCheckBox.IsChecked.Value;
            TheTask.yIons = yCheckBox.IsChecked.Value;

            DialogResult = true;
        }
    }
}
