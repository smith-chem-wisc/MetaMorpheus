using EngineLayer;
using GuiFunctions.Databases;
using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.IO;
using System.Linq;
using System.Net.Http;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;

namespace MetaMorpheusGUI
{
    /// <summary>
    /// Interaction logic for DownloadUniProtDatabaseWindow.xaml
    /// </summary>
    public partial class DownloadUniProtDatabaseWindow : Window
    {
        private List<string> availableProteomes = new();
        private ObservableCollection<string> filteredProteomes = new();
        private string selectedProteome;
        private ObservableCollection<ProteinDbForDataGrid> proteinDatabases;
        private string AbsolutePathToStorageDirectory = Path.Combine(GlobalVariables.DataDir, @"Proteomes");

        public DownloadUniProtDatabaseWindow(ObservableCollection<ProteinDbForDataGrid> ProteinDatabases)
        {
            proteinDatabases = ProteinDatabases; // passed by reference from the MainWindow.xaml. downloaded database created here is added to this collection
            InitializeComponent();
            availableProteomes.Clear();
            filteredProteomes.Clear();
            LoadAvailableProteomes();

            foreach (var item in availableProteomes)
            {
                filteredProteomes.Add(item);
            }

            availableProteomesListbox.ItemsSource = filteredProteomes;
        }

        private void LoadAvailableProteomes()
        {
            foreach (var item in GlobalVariables.AvailableUniProtProteomes)
            {
                availableProteomes.Add(item.Value);
            }
        }

        private void DownloadProteomeWithProgress(string DownloadLink, string TargetPath)
        {
            Uri uri = new(DownloadLink);
            HttpClient client = new();
            ProgressBarIndeterminate progress = new ProgressBarIndeterminate(client)
            {
                WindowStartupLocation = WindowStartupLocation.CenterScreen
            };
            progress.Show();
            Task.Factory.StartNew(() =>
            {
                client.Timeout = TimeSpan.FromMinutes(30); //This the time in case the process take response long time in this case have 30 min
                AppContext.SetSwitch("System.Net.Http.SocketsHttpHandler.Http2UnencryptedSupport", true);
                try
                {
                    Task.Run(() => ProgressInformation(uri, client, progress , TargetPath)).Wait();
                }
                catch
                {
                    //The task was canceled.
                }
            }).ContinueWith(task =>
            {
                progress.Close();
                if (ProteomeDownloaded(TargetPath))
                {
                    FileInfo fi = new(TargetPath);
                    if (fi.Length > 0)
                    {
                        proteinDatabases.Add(new ProteinDbForDataGrid(TargetPath));
                    }
                    File.SetCreationTime(TargetPath, DateTime.Now); //this is needed because overwriting an old file will not update the file creation time
                }
                this.Close();
            }, System.Threading.CancellationToken.None, TaskContinuationOptions.None,
            TaskScheduler.FromCurrentSynchronizationContext());
        }

        async Task ProgressInformation(Uri uri, HttpClient client, ProgressBarIndeterminate progressBar , string TargetPath)
        {

            var totalBytesRead = 0L;
            var buffer = new byte[8192];
            var isMoreToRead = true;
            var response = client.GetAsync(uri, HttpCompletionOption.ResponseHeadersRead).Result;
            using (var contentStream = await response.Content.ReadAsStreamAsync())
            {
                using (var fileStream = new FileStream(TargetPath, FileMode.Create, FileAccess.Write, FileShare.None, 8192, true))
                {
                    do
                    {
                        var bytesRead = await contentStream.ReadAsync(buffer, 0, buffer.Length);
                        if (bytesRead == 0)
                        {
                            isMoreToRead = false;
                            fileStream.Close();
                            continue;
                        }

                        if(!progressBar.open && isMoreToRead)
                        {
                            isMoreToRead = false;
                            fileStream.Close();
                            File.Delete(TargetPath);
                            break;
                        }

                        await fileStream.WriteAsync(buffer, 0, bytesRead);
                        totalBytesRead += bytesRead;

                        progressBar.Bytes = totalBytesRead < 10000000 ? (totalBytesRead / 1024).ToString("N0") + "KB" : (totalBytesRead / 1000000).ToString("N0") + "MB";
                    }
                    while (isMoreToRead);
                }
            }
        }

        private void DownloadProteomeButton_Click(object sender, RoutedEventArgs e)
        {

            if (!string.IsNullOrEmpty(selectedProteome) && availableProteomes.Contains(selectedProteome))
            {
                string filename = DownloadUniProtDatabaseFunctions.GetUniprotFilename(GetProteomeId(selectedProteome), reviewedCheckBox.IsChecked.Value, addIsoformsCheckBox.IsChecked.Value, xmlBox.IsChecked.Value, compressedCheckBox.IsChecked.Value);
                string htmlQueryString = DownloadUniProtDatabaseFunctions.GetUniProtHtmlQueryString(GetProteomeId(selectedProteome), reviewedCheckBox.IsChecked.Value, addIsoformsCheckBox.IsChecked.Value, xmlBox.IsChecked.Value, compressedCheckBox.IsChecked.Value);
                DownloadProteomeWithProgress(htmlQueryString, AbsolutePathToStorageDirectory + filename);
            }
        }

        private static bool ProteomeDownloaded(string downloadedProteomeFullPath)
        {
            return File.Exists(downloadedProteomeFullPath);
        }

        private static string GetProteomeId(string selectedProteome)
        {
            return (GlobalVariables.AvailableUniProtProteomes.FirstOrDefault(x => x.Value == selectedProteome).Key);
        }

        //any time text is changed in the search box, we want to filter the list displayed in the gridview
        private void ProteomesSearchBox_TextChanged(object sender, TextChangedEventArgs e)
        {
            string searchText = (sender as TextBox).Text;

            if (string.IsNullOrEmpty(searchText) == false)
            {
                filteredProteomes.Clear();
                foreach (var item in availableProteomes)
                {
                    if (item.ToUpper().StartsWith(searchText.ToUpper()) || item.ToUpper().Contains(searchText.ToUpper()))
                    {
                        filteredProteomes.Add(item);
                    }
                    else if (searchText == "")
                    {
                        filteredProteomes.Add(item);
                    }
                }
            }
        }

        private void AvailableProteomesListbox_SelectionChanged(object sender, SelectionChangedEventArgs e)
        {
            if (availableProteomesListbox.SelectedItem == null)
            {
                return;
            }

            selectedProteome = (string)availableProteomesListbox.SelectedItem;
            selectedProteomeBox.Text = selectedProteome;
        }
    }
}