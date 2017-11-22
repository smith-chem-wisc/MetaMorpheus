using EngineLayer;
using Nett;
using System;
using System.Diagnostics;
using System.Net;
using System.Net.Http;
using System.Windows;
using Newtonsoft.Json.Linq;
using Proteomics;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using Newtonsoft.Json;
using System.Text.RegularExpressions;
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
            ReleaseHandler();
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

        private void ReleaseHandler()
        {
            using (var client = new HttpClient())
            {
                client.DefaultRequestHeaders.Add("User-Agent", "Mozilla/5.0 (compatible; MSIE 10.0; Windows NT 6.2; WOW64; Trident/6.0)");

                using (var response = client.GetAsync("https://api.github.com/repos/smith-chem-wisc/MetaMorpheus/releases").Result)
                {
                    string json = response.Content.ReadAsStringAsync().Result;
                    JArray GitArray = JArray.Parse(json);
                    bool DebugV;
#if DEBUG
                    DebugV = true;
#else
                    DebugV = false;
#endif
                    int currV =0;
                    if (!DebugV)
                        currV = parser(GlobalEngineLevelSettings.NewestVersion);
                    currV = 190;
                    foreach (JObject obj in GitArray.Children<JObject>())
                    {
                        string str = obj.SelectToken("tag_name") + "";
                        if (!DebugV && parser(str) <= currV)
                            break;
                        releases.Text += (str).TrimEnd('\r', '\n')+ "\n";
                        str = obj.SelectToken("body") + "";
                        if (str.Trim('\n','\r',' ','\t').Equals(""))
                            releases.Text += "N/A\n";
                        else
                            releases.Text += (str).TrimEnd('\r', '\n') + "\n";
                    }
                }
            }
        }
        public int parser(string VersionNode)
        {
            string pattern = @"\d\.\d\.(\d+)";
            return Int32.Parse(Regex.Match(VersionNode, pattern).Groups[1].Value);
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