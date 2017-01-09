using IndexSearchAndAnalyze;
using MetaMorpheus;
using System.Collections.ObjectModel;
using System.Windows;

namespace GoodGUI
{
    /// <summary>
    /// Interaction logic for GPTMDTaskWindow.xaml
    /// </summary>
    public partial class GPTMDTaskWindow : Window
    {
        // Always create a new one, even if updating an existing task
        private ObservableCollection<ModListForSearch> ModFileListInWindow = new ObservableCollection<ModListForSearch>();

        public GPTMDTaskWindow(ObservableCollection<ModList> modFileList)
        {
            InitializeComponent();
        }

        public GPTMDTaskWindow(MyGPTMDtask myGPTMDtask, ObservableCollection<ModList> modFileList)
        {
            InitializeComponent();
        }

        internal MyGPTMDtask TheTask { get; private set; }

        private void cancelButton_Click(object sender, RoutedEventArgs e)
        {
            DialogResult = false;
        }

        private void saveButton_Click(object sender, RoutedEventArgs e)
        {
            //TheTask.searchDecoy = checkBoxDecoy.IsChecked.Value;
            //TheTask.maxMissedCleavages = int.Parse(missedCleavagesTextBox.Text);
            //TheTask.protease = (Protease)proteaseComboBox.SelectedItem;
            //TheTask.maxModificationIsoforms = int.Parse(maxModificationIsoformsTextBox.Text);
            //TheTask.initiatorMethionineBehavior = (InitiatorMethionineBehavior)initiatorMethionineBehaviorComboBox.SelectedIndex;
            //TheTask.productMassTolerance.Value = double.Parse(productMassToleranceTextBox.Text);
            //TheTask.productMassTolerance.Unit = (ToleranceUnit)productMassToleranceComboBox.SelectedIndex;
            //TheTask.bIons = bCheckBox.IsChecked.Value;
            //TheTask.yIons = yCheckBox.IsChecked.Value;
            //TheTask.listOfModListsForSearch = ModFileListInWindow.ToList();

            DialogResult = true;
        }
    }
}