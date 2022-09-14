using EngineLayer;
using Newtonsoft.Json.Linq;
using System;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Net.Http;
using System.Text;
using System.Threading.Tasks;
using System.Windows;

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
            lbl.Text = "A newer version: " + MainWindow.NewestKnownMetaMorpheusVersion + " is available!";
            ReleaseHandler();
        }

        public static (int, int, int) GetVersionNumber(string VersionNode)
        {
            try
            {
                var split = VersionNode.Split('.');

                return (int.Parse(split[0]), int.Parse(split[1]), int.Parse(split[2]));
            }
            catch (FormatException)
            {
                return (0, 0, 0);
            }
        }

        public static string GetVersionNumbersFromWeb()
        {
            // Attempt to get current MetaMorpheus version
            using (var client = new HttpClient())
            {
                client.DefaultRequestHeaders.Add("User-Agent", "Mozilla/5.0 (compatible; MSIE 10.0; Windows NT 6.2; WOW64; Trident/6.0)");

                using (var response = client.GetAsync("https://api.github.com/repos/smith-chem-wisc/MetaMorpheus/releases/latest").Result)
                {
                    var json = response.Content.ReadAsStringAsync().Result;
                    JObject deserialized = JObject.Parse(json);
                    var assets = deserialized["assets"].Select(b => b["name"].ToString()).ToList();
                    if (!assets.Contains("MetaMorpheusInstaller.msi"))
                    {
                        throw new MetaMorpheusException("A new version of MetaMorpheus was detected, but the files haven't been" +
                            " uploaded yet. Try again in a few minutes.");
                    }

                    return deserialized["tag_name"].ToString();
                }
            }
        }

        private void InstallerClicked(object sender, RoutedEventArgs e)
        {
            DialogResult = true;

            HttpClient client = new();
            var uri = new Uri(@"https://github.com/smith-chem-wisc/MetaMorpheus/releases/download/" + MainWindow.NewestKnownMetaMorpheusVersion + @"/MetaMorpheusInstaller.msi");

            Exception exception = null;
            try
            {
                var tempDownloadLocation = Path.Combine(Path.GetTempPath(), "MetaMorpheusInstaller.msi");

                // download the installer
                HttpResponseMessage urlResponse = Task.Run(() => client.GetAsync(uri)).Result;
                using (FileStream stream = new(tempDownloadLocation, FileMode.CreateNew))
                {
                    Task.Run(() => urlResponse.Content.CopyToAsync(stream)).Wait();
                }

                // start the installer
                GlobalVariables.StartProcess(tempDownloadLocation);
            }
            catch (Exception ex)
            {
                exception = ex;
                MessageBox.Show(ex.Message);
            }

            if (exception == null)
            {
                // close metamorpheus if the installer was started successfully
                Application.Current.Shutdown();
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
                    var currV = GetVersionNumber(GlobalVariables.MetaMorpheusVersion);
                    StringBuilder allVersionsText = new StringBuilder();
                    foreach (JObject obj in GitArray.Children<JObject>())
                    {
                        var checkVersion = GetVersionNumber(obj.SelectToken("tag_name").ToString());
                        if (checkVersion.Item1 < currV.Item1 ||
                            (checkVersion.Item1 == currV.Item1 && checkVersion.Item2 < currV.Item2) ||
                            (checkVersion.Item1 == currV.Item1 && checkVersion.Item2 == currV.Item2 && checkVersion.Item3 <= currV.Item3))
                            break;
                        string body = new MarkdownSharp.Markdown().Transform(obj.SelectToken("body").ToString());
                        allVersionsText.AppendLine("<font face=\"Arial\" size=2>");
                        allVersionsText.AppendLine("<h3>" + obj.SelectToken("tag_name").ToString() + "</h3>");
                        allVersionsText.AppendLine(body);
                        allVersionsText.AppendLine();
                        allVersionsText.AppendLine("</font>");
                    }
                    releases.NavigateToString(allVersionsText.ToString());
                    releases.Navigating += Releases_Navigating;
                }
            }
        }

        public void Releases_Navigating(object sender, System.Windows.Navigation.NavigatingCancelEventArgs e)
        {
            //cancel the current event
            e.Cancel = true;

            //this opens the URL in the user's default browser
            GlobalVariables.StartProcess(e.Uri.ToString());
        }

        private void PortableClicked(object semder, RoutedEventArgs e)
        {
            DialogResult = true;
            HttpClient client = new();
            var uri = new Uri(@"https://github.com/smith-chem-wisc/MetaMorpheus/releases/download/" + MainWindow.NewestKnownMetaMorpheusVersion + @"/MetaMorpheusGuiDotNetFrameworkAppveyor.zip");

            try
            {
                var tempDownloadLocation = Path.Combine(Path.GetTempPath(), "MetaMorpheusGuiDotNetFrameworkAppveyor.zip");
                HttpResponseMessage urlResponse = Task.Run(() => client.GetAsync(uri)).Result;
                using (FileStream stream = new(tempDownloadLocation, FileMode.CreateNew))
                {
                    Task.Run(() => urlResponse.Content.CopyToAsync(stream)).Wait();
                }
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

        private void NoClicked(object semder, RoutedEventArgs e)
        {
            DialogResult = false;
        }
    }
}