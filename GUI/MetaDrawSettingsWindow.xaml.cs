using EngineLayer;
using System;
using System.ComponentModel;
using System.IO;
using System.Windows;
using System.Windows.Shapes;

namespace MetaMorpheusGUI
{
    /// <summary>
    /// Interaction logic for MetaDrawSettingsWindow.xaml
    /// </summary>
    public partial class MetaDrawGraphicalSettings : Window
    {
        public MetaDrawSettings settings = new MetaDrawSettings();

        public MetaDrawGraphicalSettings()
        {
            InitializeComponent();
            settings.Size = 12;
            base.Closing += this.OnClosing;
        }

        private void Cancel_Click(object sender, RoutedEventArgs e) {
            this.Visibility = Visibility.Hidden;
        }
        
        private void Save_Click(object sender, RoutedEventArgs e)
        {
            // get check box values and font size, and send to draw
            // have a settings file that stores information? 
            settings.ShowMzValues = MZCheckBox.IsChecked == true ? true : false;
            settings.ShowAnnotationCharges = ChargesCheckBox.IsChecked == true ? true : false;
            settings.BoldText = BoldTextCheckBox.IsChecked == true ? true : false;
            settings.Size = TextSizeBox.Text.Length > 0 ? Convert.ToInt32(TextSizeBox.Text) : 12;

            this.Visibility = Visibility.Hidden;
        }

        private void OnClosing(object sender, CancelEventArgs e)
        {
            // this window only closes when task is added
            if (this.Visibility == Visibility.Visible)
            {
                this.Visibility = Visibility.Hidden;
                e.Cancel = true;
            }
        }
    }
}
