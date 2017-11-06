using EngineLayer;
using System;
using System.IO;
using System.Windows;

namespace MetaMorpheusGUI
{
    /// <summary>
    /// Interaction logic for MetaUpdater.xaml
    /// </summary>
    public partial class MetaUpdater : Window
    {
        #region Public Constructors

        public MetaUpdater()
        {
            InitializeComponent();
            lbl.Text = "A newer version: " + GlobalEngineLevelSettings.NewestVersion + " is available, would you like to download it?";
        }

        #endregion Public Constructors

        #region Private Methods

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

        #endregion Private Methods
    }
}