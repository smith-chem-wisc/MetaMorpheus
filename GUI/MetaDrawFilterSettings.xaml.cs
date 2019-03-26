using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Shapes;

namespace MetaMorpheusGUI
{
    /// <summary>
    /// Interaction logic for MetaDrawFilterSettings.xaml
    /// </summary>
    public partial class MetaDrawFilterSettings : Window
    {
        public bool ShowDecoys { get; private set; }
        public bool ShowContaminants { get; private set; }
        public double QValue { get; private set; }
        public bool SelectContaminants { get; set; } = true;

        public MetaDrawFilterSettings()
        {
            InitializeComponent();
            ContaminantsCheckBox.IsChecked = true; 
            base.Closing += this.OnClosing;

            // default values
            QValue = 0.01;
            ShowContaminants = true;
            ShowDecoys = false;
        }
        
        private void Cancel_Click(object sender, RoutedEventArgs e)
        {
            this.Visibility = Visibility.Hidden;
        }

        private void Save_Click(object sender, RoutedEventArgs e)
        {
            ShowDecoys = DecoysCheckBox.IsChecked == true ? true : false;
            ShowContaminants = ContaminantsCheckBox.IsChecked == true ? true : false;
            if (qValueBox.Text.Length > 0)
            {
                QValue = Convert.ToDouble(qValueBox.Text);
            }

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
