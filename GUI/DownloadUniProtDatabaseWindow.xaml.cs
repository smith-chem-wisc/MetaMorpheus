using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.IO;
using System.Linq;
using System.Net;
using System.Windows;
using System.Windows.Controls;
using UsefulProteomicsDatabases;
using static UsefulProteomicsDatabases.ProteinDbRetriever;

namespace MetaMorpheusGUI
{
    /// <summary>
    /// Interaction logic for DownloadUniProtDatabaseWindow.xaml
    /// </summary>
    public partial class DownloadUniProtDatabaseWindow : Window
    {

        private List<string> availableProteomes = new List<string>();
        private ObservableCollection<string> filteredProteomes = new ObservableCollection<string>();
        private string selectedProteome;
        private string _downloadedProteomeFullPath;
        public string DownloadedFilepath { get { return _downloadedProteomeFullPath; } }
        public DownloadUniProtDatabaseWindow()
        {

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
            using (WebClient Client = new WebClient())
            {
                ProgressBarIndeterminate progress = new ProgressBarIndeterminate
                {
                    WindowStartupLocation = WindowStartupLocation.CenterScreen
                };
                progress.Show();
                System.Threading.Tasks.Task.Factory.StartNew(() =>
                {
                    Client.DownloadFile(DownloadLink, TargetPath);
                }).ContinueWith(task =>
                {
                    progress.Close();
                }, System.Threading.CancellationToken.None, System.Threading.Tasks.TaskContinuationOptions.None, System.Threading.Tasks.TaskScheduler.FromCurrentSynchronizationContext());
            }
            
            if (File.Exists(TargetPath))
            {
                _downloadedProteomeFullPath = TargetPath;
            }
        }

        public static string GetQueryString(string proteomeID, ProteomeFormat format, Reviewed reviewed, Compress compress, IncludeIsoforms include)
        {
            string htmlQueryString = "";
            if (format == ProteomeFormat.fasta)
            {

                if (reviewed == Reviewed.yes)
                {
                    htmlQueryString = "https://www.uniprot.org/uniprot/?query=proteome:" + proteomeID + " reviewed:" + reviewed + "&compress=" + compress + "&format=" + format + "&include:" + include;
                }
                else
                {
                    htmlQueryString = "https://www.uniprot.org/uniprot/?query=proteome:" + proteomeID + "&compress=" + compress + "&format=" + format + "&include:" + include;
                }
            }
            else if (format == ProteomeFormat.xml)
            {
                htmlQueryString = "https://www.uniprot.org/uniprot/?query=proteome:" + proteomeID + " reviewed:" + reviewed + "&compress=" + compress + "&format=" + format;
            }

            return htmlQueryString;
        }

        private void downloadProteomeButton_Click(object sender, RoutedEventArgs e)
        {
            if (string.IsNullOrEmpty(selectedProteome) == false && availableProteomes.Contains(selectedProteome))
            {
                string absolutePathToStorageDirectory = "E:\\junk";
                string filename = "\\" + GetProteomeId(selectedProteome);

                ProteinDbRetriever.Reviewed reviewed = new ProteinDbRetriever.Reviewed();
                if (reviewedCheckBox.IsChecked == true)
                {
                    reviewed = ProteinDbRetriever.Reviewed.yes;
                    filename += "_reviewed";
                }
                else
                {
                    reviewed = ProteinDbRetriever.Reviewed.no;
                    filename += "_unreviewed";
                }

                ProteinDbRetriever.IncludeIsoforms isoforms;
                if (addIsoformsCheckBox.IsChecked == true)
                {
                    isoforms = ProteinDbRetriever.IncludeIsoforms.yes;
                    filename += "_isoform";
                }
                else
                {
                    isoforms = ProteinDbRetriever.IncludeIsoforms.no;
                }
                ProteinDbRetriever.ProteomeFormat format = new ProteinDbRetriever.ProteomeFormat();
                if (xmlBox.IsChecked == true)
                {
                    format = ProteinDbRetriever.ProteomeFormat.xml;
                    filename += ".xml";
                }
                else
                {
                    format = ProteinDbRetriever.ProteomeFormat.fasta;
                    filename += ".fasta";
                }
                ProteinDbRetriever.Compress compressed;
                if (compressedCheckBox.IsChecked == true)
                {
                    compressed = ProteinDbRetriever.Compress.yes;
                    filename += ".gz";
                }
                else
                {
                    compressed = ProteinDbRetriever.Compress.no;
                }
                string htmlQueryString = GetQueryString(GetProteomeId(selectedProteome), format, reviewed, compressed, isoforms);

                DownloadProteomeWithProgress(htmlQueryString, absolutePathToStorageDirectory + filename);

                this.Close();
            }
        }

        private string GetProteomeId(string selectedProteome)
        {
            return (EngineLayer.GlobalVariables.AvailableUniProtProteomes.FirstOrDefault(x => x.Value == selectedProteome).Key);
        }

        //any time text is changed in the search box, we want to filter the list displayed in the gridview
        private void proteomesSearchBox_TextChanged(object sender, TextChangedEventArgs e)
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

        private void availableProteomesListbox_SelectionChanged(object sender, SelectionChangedEventArgs e)
        {
            if (availableProteomesListbox.SelectedItem == null)
            {
                return;
            }

            selectedProteome = (string)availableProteomesListbox.SelectedItem;
            selectedProteomeBox.Text = selectedProteome;
        }

        //private string Window_Closing(object sender, System.ComponentModel.CancelEventArgs e)
        //{
        //    return "";
        //}
    }
}