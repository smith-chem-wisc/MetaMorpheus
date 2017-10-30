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
using System.IO;
using MetaMorpheusGUI;
namespace MetaMorpheusGUI
{
    /// <summary>
    /// Interaction logic for MetaUpdater.xaml
    /// </summary>
    public partial class MetaUpdater : Window
    {
        public MetaUpdater()
        {
            InitializeComponent();
            lbl.Text = "A newer version: " + MainWindow.VersionCheck + " is available, would you like to install?";
        }

        private void CheckBox_Checked(object sender, RoutedEventArgs e)
        {
            String LogPath = "update.txt";
            File.WriteAllText(LogPath, "false");
        }
        private void CheckBox_unChecked(object sender, RoutedEventArgs e)
        {
            String LogPath = "update.txt";
            File.WriteAllText(LogPath, "true");
        }
        private void YesClicked(object semder, RoutedEventArgs e)
        {
            this.Close();
            MetaUpdater_Yes newwind = new MetaUpdater_Yes();
            newwind.Show();
        }
        private void NoClicked(object semder, RoutedEventArgs e)
        {
            this.Close();

        }
    }
}
