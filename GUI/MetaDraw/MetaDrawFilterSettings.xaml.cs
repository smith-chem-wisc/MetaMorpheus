using System.ComponentModel;
using System.Windows;
using EngineLayer.GlycoSearch;

namespace MetaMorpheusGUI
{
    /// <summary>
    /// Interaction logic for MetaDrawFilterSettings.xaml
    /// </summary>
    public partial class MetaDrawFilterSettings : Window
    {
        public bool ShowDecoys { get; private set; }
        public bool ShowContaminants { get; private set; }
        public double QValueFilter { get; private set; }

        public LocalizationLevel LocalizationLevelStart { get; set; }
        public LocalizationLevel LocalizationLevelEnd { get; set; }

        public MetaDrawFilterSettings(bool showDecoys = false, bool showContaminants = true, double qValueFilter = 0.01, LocalizationLevel localizationLevelStart = LocalizationLevel.Level1, LocalizationLevel localizationLevelEnd = LocalizationLevel.Level3)
        {
            InitializeComponent();
            base.Closing += this.OnClosing;

            foreach (string level in System.Enum.GetNames(typeof(LocalizationLevel)))
            {
                CmbGlycanLocalizationLevelStart.Items.Add(level);
                CmbGlycanLocalizationLevelEnd.Items.Add(level);
            }

            QValueFilter = qValueFilter;
            ShowContaminants = showContaminants;
            ShowDecoys = showDecoys;
            LocalizationLevelStart = localizationLevelStart;
            LocalizationLevelEnd = localizationLevelEnd;

            CmbGlycanLocalizationLevelStart.SelectedItem = localizationLevelStart.ToString();
            CmbGlycanLocalizationLevelEnd.SelectedItem = localizationLevelEnd.ToString();
        }

        private void Cancel_Click(object sender, RoutedEventArgs e)
        {
            this.Visibility = Visibility.Hidden;
        }

        private void Save_Click(object sender, RoutedEventArgs e)
        {
            ShowDecoys = DecoysCheckBox.IsChecked.Value;
            ShowContaminants = ContaminantsCheckBox.IsChecked.Value;
            
            if (!string.IsNullOrWhiteSpace(qValueBox.Text))
            {
                if (double.TryParse(qValueBox.Text, out double qValueFilter))
                {
                    QValueFilter = qValueFilter;
                }
                else
                {
                    MessageBox.Show("Could not parse q-value filter");
                    return;
                }
            }
            LocalizationLevelStart = (LocalizationLevel)System.Enum.Parse(typeof(LocalizationLevel), CmbGlycanLocalizationLevelStart.SelectedItem.ToString());
            LocalizationLevelEnd = (LocalizationLevel)System.Enum.Parse(typeof(LocalizationLevel), CmbGlycanLocalizationLevelEnd.SelectedItem.ToString());
            this.Visibility = Visibility.Hidden;
        }

        private void OnClosing(object sender, CancelEventArgs e)
        {
            if (this.Visibility == Visibility.Visible)
            {
                this.Visibility = Visibility.Hidden;
                e.Cancel = true;
            }
        }

    }
}
