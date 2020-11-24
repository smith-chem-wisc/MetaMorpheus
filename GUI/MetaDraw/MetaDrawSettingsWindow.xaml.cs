using EngineLayer;
using EngineLayer.GlycoSearch;
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
            PopulateChoices();
        }

        private void PopulateChoices()
        {
            foreach (string level in System.Enum.GetNames(typeof(LocalizationLevel)))
            {
                CmbGlycanLocalizationLevelStart.Items.Add(level);
                CmbGlycanLocalizationLevelEnd.Items.Add(level);
            }

            DisplayAnnotationsCheckBox.IsChecked = MetaDrawSettings.DisplayIonAnnotations;
            MZCheckBox.IsChecked = MetaDrawSettings.AnnotateMzValues;
            ChargesCheckBox.IsChecked = MetaDrawSettings.AnnotateCharges;
            BoldTextCheckBox.IsChecked = MetaDrawSettings.AnnotationBold;
            DecoysCheckBox.IsChecked = MetaDrawSettings.ShowDecoys;
            ContaminantsCheckBox.IsChecked = MetaDrawSettings.ShowContaminants;
            qValueBox.Text = MetaDrawSettings.QValueFilter.ToString();
            TextSizeBox.Text = MetaDrawSettings.AnnotatedFontSize.ToString();
            CmbGlycanLocalizationLevelStart.SelectedItem = MetaDrawSettings.LocalizationLevelStart.ToString();
            CmbGlycanLocalizationLevelEnd.SelectedItem = MetaDrawSettings.LocalizationLevelEnd.ToString();
        }

        private void Cancel_Click(object sender, RoutedEventArgs e)
        {
            DialogResult = false;
        }

        private void Save_Click(object sender, RoutedEventArgs e)
        {
            MetaDrawSettings.DisplayIonAnnotations = DisplayAnnotationsCheckBox.IsChecked.Value;
            MetaDrawSettings.AnnotateMzValues = MZCheckBox.IsChecked.Value;
            MetaDrawSettings.AnnotateCharges = ChargesCheckBox.IsChecked.Value;
            MetaDrawSettings.AnnotationBold = BoldTextCheckBox.IsChecked.Value;
            MetaDrawSettings.ShowDecoys = BoldTextCheckBox.IsChecked.Value;
            MetaDrawSettings.ShowContaminants = BoldTextCheckBox.IsChecked.Value;
            MetaDrawSettings.LocalizationLevelStart = (LocalizationLevel)System.Enum.Parse(typeof(LocalizationLevel), CmbGlycanLocalizationLevelStart.SelectedItem.ToString());
            MetaDrawSettings.LocalizationLevelEnd = (LocalizationLevel)System.Enum.Parse(typeof(LocalizationLevel), CmbGlycanLocalizationLevelEnd.SelectedItem.ToString());

            if (!string.IsNullOrWhiteSpace(qValueBox.Text))
            {
                if (double.TryParse(qValueBox.Text, out double qValueFilter) && qValueFilter >= 0 && qValueFilter <= 1)
                {
                    MetaDrawSettings.QValueFilter = qValueFilter;
                }
                else
                {
                    MessageBox.Show("Could not parse q-value filter; must be number between 0 and 1 inclusive");
                    return;
                }
            }
            else
            {
                MetaDrawSettings.QValueFilter = 1;
            }

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
                    MessageBox.Show("Could not parse font size; must be a positive integer");
                    return;
                }
            }
            else
            {
                MetaDrawSettings.AnnotatedFontSize = 12;
            }

            DialogResult = true;
        }
    }
}
