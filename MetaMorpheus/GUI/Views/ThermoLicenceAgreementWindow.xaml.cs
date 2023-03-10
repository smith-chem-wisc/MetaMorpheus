using System.Windows;

namespace MetaMorpheusGUI
{
    /// <summary>
    /// Interaction logic for ThermoLicenceAgreementWindow.xaml
    /// </summary>
    public partial class ThermoLicenceAgreementWindow : Window
    {
        public ThermoLicenceAgreementWindow()
        {
            InitializeComponent();
        }

        private void AgreeButton_Click(object sender, RoutedEventArgs e)
        {
            DialogResult = true;
            return;
        }

        private void DisagreeButton_Click(object sender, RoutedEventArgs e)
        {
            DialogResult = false;
            return;
        }
    }
}
