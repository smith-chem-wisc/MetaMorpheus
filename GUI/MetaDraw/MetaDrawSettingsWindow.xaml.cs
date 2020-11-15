using EngineLayer;
using System.ComponentModel;
using System.Windows;

namespace MetaMorpheusGUI
{
    /// <summary>
    /// Interaction logic for MetaDrawSettingsWindow.xaml
    /// </summary>
    public partial class MetaDrawSettingsWindow : Window
    {
        public MetaDrawSettingsWindow()
        {
            InitializeComponent();
        }

        private void Cancel_Click(object sender, RoutedEventArgs e)
        {
            DialogResult = false;
        }

        private void Save_Click(object sender, RoutedEventArgs e)
        {
            MetaDrawSettings.ShowMzValues = MZCheckBox.IsChecked.Value;
            MetaDrawSettings.ShowAnnotationCharges = ChargesCheckBox.IsChecked.Value;
            MetaDrawSettings.BoldText = BoldTextCheckBox.IsChecked.Value;

            if (!string.IsNullOrWhiteSpace(TextSizeBox.Text))
            {
                if (int.TryParse(TextSizeBox.Text, out int fontSize))
                {
                    if (fontSize > 15)
                    {
                        MessageBox.Show("Font size must be <= 15");
                        return;
                    }

                    MetaDrawSettings.AnnotatedFontSize = fontSize;
                }
                else
                {
                    MessageBox.Show("Could not parse font size");
                    return;
                }
            }

            DialogResult = true;
        }
    }
}
