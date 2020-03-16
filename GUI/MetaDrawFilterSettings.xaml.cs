using System.ComponentModel;
using System.Windows;

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

        public int? GlycanLocalizationLevelFilter { get; set; }

        public MetaDrawFilterSettings(bool showDecoys = false, bool showContaminants = true, double qValueFilter = 0.01, int? glycanLocalizationLevelFilter = null)
        {
            InitializeComponent();
            base.Closing += this.OnClosing;

            QValueFilter = qValueFilter;
            ShowContaminants = showContaminants;
            ShowDecoys = showDecoys;
            GlycanLocalizationLevelFilter = glycanLocalizationLevelFilter;
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

            if (!string.IsNullOrWhiteSpace(glycanLocalizationLelveBox.Text))
            {
                if (int.TryParse(glycanLocalizationLelveBox.Text, out int glycanLocalizationLevelFilter))
                {
                    GlycanLocalizationLevelFilter = glycanLocalizationLevelFilter;
                }
                else
                {
                    MessageBox.Show("Could not parse Glycan Localization Level filter");
                    return;
                }
            }
            else
            {
                GlycanLocalizationLevelFilter = null;
            }

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
