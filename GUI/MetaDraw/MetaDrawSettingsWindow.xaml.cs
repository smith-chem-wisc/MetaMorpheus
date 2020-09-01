using System.ComponentModel;
using System.Windows;

namespace MetaMorpheusGUI
{
    /// <summary>
    /// Interaction logic for MetaDrawSettingsWindow.xaml
    /// </summary>
    public partial class MetaDrawGraphicalSettings : Window
    {
        public bool ShowMzValues { get; set; }
        public bool ShowAnnotationCharges { get; set; }
        public int AnnotatedFontSize { get; set; }
        public bool BoldText { get; set; }

        public MetaDrawGraphicalSettings(bool annotateMzs = false, bool annotateCharges = false, int annotatedFontSize = 12, bool boldAnnotations = false)
        {
            InitializeComponent();
            base.Closing += this.OnClosing;
            this.ShowMzValues = annotateMzs;
            this.ShowAnnotationCharges = annotateCharges;
            this.AnnotatedFontSize = annotatedFontSize;
            this.BoldText = boldAnnotations;
        }

        private void Cancel_Click(object sender, RoutedEventArgs e)
        {
            this.Visibility = Visibility.Hidden;
        }

        private void Save_Click(object sender, RoutedEventArgs e)
        {
            ShowMzValues = MZCheckBox.IsChecked.Value;
            ShowAnnotationCharges = ChargesCheckBox.IsChecked.Value;
            BoldText = BoldTextCheckBox.IsChecked.Value;

            if (!string.IsNullOrWhiteSpace(TextSizeBox.Text))
            {
                if (int.TryParse(TextSizeBox.Text, out int fontSize))
                {
                    if (fontSize > 15)
                    {
                        MessageBox.Show("Font size must be <= 15");
                        return;
                    }

                    AnnotatedFontSize = fontSize;
                }
                else
                {
                    MessageBox.Show("Could not parse font size");
                    return;
                }
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
