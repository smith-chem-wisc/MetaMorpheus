using EngineLayer;
using System;
using System.Diagnostics;
using System.Net;
using System.Windows;

namespace MetaMorpheusGUI
{
    /// <summary>
    /// Interaction logic for MetaUpdater_Yes.xaml
    /// </summary>
    public partial class MetaUpdater_Yes : Window
    {
        #region Public Constructors

        public MetaUpdater_Yes()
        {
            InitializeComponent();
        }

        #endregion Public Constructors

        #region Private Methods

        private void InstallerClicked(object sender, RoutedEventArgs e)
        {
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

        private void FilesClicked(object sender, RoutedEventArgs e)
        {
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

        #endregion Private Methods
    }
}