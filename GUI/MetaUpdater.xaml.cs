using EngineLayer;
using Nett;
using System;
using System.Diagnostics;
using System.Net;
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
            lbl.Text = "A newer version: " + GlobalEngineLevelSettings.NewestVersion + " is available!";
        }

        #endregion Public Constructors

        #region Private Methods

        private void InstallerClicked(object semder, RoutedEventArgs e)
        {
            DialogResult = true;
            using (var client = new WebClient())
            {
                var uri = new Uri(@"https://github.com/smith-chem-wisc/MetaMorpheus/releases/download/" + GlobalEngineLevelSettings.NewestVersion + @"/MetaMorpheusInstaller.msi");

                try
                {
                    var tempDownloadLocation = System.IO.Path.Combine(System.IO.Path.GetTempPath(), "MetaMorpheusInstaller.msi");
                    client.DownloadFile(uri, tempDownloadLocation);
                    Process p = new Process();
                    p.StartInfo.FileName = tempDownloadLocation;
                    Application.Current.Shutdown();
                    p.Start();
                }
                catch (Exception ex)
                {
                    MessageBox.Show(ex.Message);
                }
            }
        }

        private void PortableClicked(object semder, RoutedEventArgs e)
        {
            DialogResult = true;
            using (var client = new WebClient())
            {
                var uri = new Uri(@"https://github.com/smith-chem-wisc/MetaMorpheus/releases/download/" + GlobalEngineLevelSettings.NewestVersion + @"/MetaMorpheusGuiDotNetFrameworkAppveyor.zip");

                try
                {
                    var tempDownloadLocation = System.IO.Path.Combine(System.IO.Path.GetTempPath(), "MetaMorpheusGuiDotNetFrameworkAppveyor.zip");
                    client.DownloadFile(uri, tempDownloadLocation);
                    Process p = new Process();
                    p.StartInfo.FileName = tempDownloadLocation;
                    Application.Current.Shutdown();
                    p.Start();
                }
                catch (Exception ex)
                {
                    MessageBox.Show(ex.Message);
                }
            }
        }

        private void NoClicked(object semder, RoutedEventArgs e)
        {
            DialogResult = false;
        }

        private void CheckBox_Checked(object sender, RoutedEventArgs e)
        {
            try
            {
                TomlTable obj = Toml.ReadFile(GlobalEngineLevelSettings.settingsTomlLocation);
                obj.Update("AskAboutUpdating", false);
                Toml.WriteFile(obj, GlobalEngineLevelSettings.settingsTomlLocation);
            }
            catch (Exception inner)
            {
                MessageBox.Show(inner.ToString());
            }
        }

        private void CheckBox_Unchecked(object sender, RoutedEventArgs e)
        {
            try
            {
                TomlTable obj = Toml.ReadFile(GlobalEngineLevelSettings.settingsTomlLocation);
                obj.Update("AskAboutUpdating", true);
                Toml.WriteFile(obj, GlobalEngineLevelSettings.settingsTomlLocation);
            }
            catch (Exception inner)
            {
                MessageBox.Show(inner.ToString());
            }
        }

        #endregion Private Methods
    }
}