using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.IO;
using System.Linq;
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

        private void downloadProteomeButton_Click(object sender, RoutedEventArgs e)
        {
            if (string.IsNullOrEmpty(selectedProteome) == false && availableProteomes.Contains(selectedProteome))
            {
                ProteinDbRetriever.ProteomeFormat format = new ProteinDbRetriever.ProteomeFormat();
                if(xmlBox.IsChecked == true)
                {
                    format = ProteinDbRetriever.ProteomeFormat.xml;
                }
                else
                {
                    format = ProteinDbRetriever.ProteomeFormat.fasta;
                }
                ProteinDbRetriever.Reviewed reviewed = new ProteinDbRetriever.Reviewed();
                if (reviewedCheckBox.IsChecked == true)
                {
                    reviewed = ProteinDbRetriever.Reviewed.yes;
                }
                else
                {
                    reviewed = ProteinDbRetriever.Reviewed.no;
                }
                ProteinDbRetriever.Compress compressed;
                if(compressedCheckBox.IsChecked == true)
                {
                    compressed = ProteinDbRetriever.Compress.yes;
                }
                else
                {
                    compressed = ProteinDbRetriever.Compress.no;
                }
                ProteinDbRetriever.IncludeIsoforms isoforms;
                if(addIsoformsCheckBox.IsChecked == true)
                {
                    isoforms = ProteinDbRetriever.IncludeIsoforms.yes;
                }
                else
                {
                    isoforms = ProteinDbRetriever.IncludeIsoforms.no;
                }

                ProteinDbRetriever.RetrieveProteome(GetProteomeId(selectedProteome), @"E:\junk", format, reviewed, compressed, isoforms);
                this.Close();
            }
        }


        public string GetProteomeId(string selectedProteome)
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