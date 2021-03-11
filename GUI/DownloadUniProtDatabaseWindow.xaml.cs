using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.IO;
using System.Linq;
using System.Net;
using System.Windows;
using System.Windows.Controls;
using UsefulProteomicsDatabases;

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
            Stream remoteStream = null;
            Stream localStream = null;
            WebResponse response = null;

            try
            {
                double TotalBytesToReceive = 0;
                var SizewebRequest = HttpWebRequest.Create(new Uri(DownloadLink));
                SizewebRequest.Method = "HEAD";

                using (var webResponse = SizewebRequest.GetResponse())
                {
                    var fileSize = webResponse.Headers.Get("Content-Length");
                    TotalBytesToReceive = Convert.ToDouble(fileSize);
                }

                ProgressBarIndeterminate progress = new ProgressBarIndeterminate
                {
                    WindowStartupLocation = WindowStartupLocation.CenterScreen
                };
                progress.Show();

                WebRequest request = WebRequest.Create(DownloadLink);
                if (request != null)
                {
                    response = request.GetResponse();
                    if (response != null)
                    {
                        remoteStream = response.GetResponseStream();
                        string filePath = TargetPath;
                        localStream = File.Create(filePath);
                        byte[] buffer = new byte[1024];
                        int bytesRead = 0;
                        do
                        {
                            bytesRead = remoteStream.Read(buffer, 0, buffer.Length);
                            localStream.Write(buffer, 0, bytesRead);

                        } while (bytesRead > 0);


                    }
                }
                progress.Close();
            }
            catch (Exception ex)
            {
                string myException = ex.ToString();
            }
            finally
            {
                if (response != null) response.Close();
                if (remoteStream != null) remoteStream.Close();
                if (localStream != null) localStream.Close();
            }
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
                string htmlQueryString = ProteinDbRetriever.RetrieveProteome(GetProteomeId(selectedProteome), @"E:\junk", format, reviewed, compressed, isoforms);

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
    }
}