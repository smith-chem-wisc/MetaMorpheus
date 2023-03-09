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
            foreach (var item in EngineLayer.GlobalVariables.AvailableUniProtProteomes)
            {
                availableProteomes.Add(item.Value);
            }
        }

        private void DownloadProteomeWithProgress(string DownloadLink, string TargetPath)
        {
            Uri uri = new(DownloadLink);
            HttpClient client = new();
            ProgressBarIndeterminate progress = new ProgressBarIndeterminate
            {
                WindowStartupLocation = WindowStartupLocation.CenterScreen
            };
            progress.Show();
            Task.Factory.StartNew(() =>
            {
                HttpResponseMessage urlResponse = Task.Run(() => client.GetAsync(uri)).Result;
                using (FileStream stream = new(TargetPath, FileMode.CreateNew))
                {
                    Task.Run(() => urlResponse.Content.CopyToAsync(stream)).Wait();
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
            }, System.Threading.CancellationToken.None, System.Threading.Tasks.TaskContinuationOptions.None,
            System.Threading.Tasks.TaskScheduler.FromCurrentSynchronizationContext());
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