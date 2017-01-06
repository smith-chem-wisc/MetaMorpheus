using System;
using System.Collections.Generic;
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

namespace GoodGUI
{
    /// <summary>
    /// Interaction logic for CalibrateTaskWindow.xaml
    /// </summary>
    public partial class CalibrateTaskWindow : Window
    {
        public CalibrateTaskWindow()
        {
            InitializeComponent();
        }

        internal MyTask TheTask { get; set; }

        private void cancelButton_Click(object sender, RoutedEventArgs e)
        {
            DialogResult = false;
        }
        private void saveButton_Click(object sender, RoutedEventArgs e)
        {
            DialogResult = true;
        }
    }
}
