using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.ComponentModel;
using System.Data;
using System.IO;
using System.Text;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Shapes;
using UsefulProteomicsDatabases;

namespace MetaMorpheusGUI
{
    /// <summary>
    /// Interaction logic for DownloadUniProtDatabaseWindow.xaml
    /// </summary>
    public partial class DownloadUniProtDatabaseWindow : Window
    {
        List<string> filteredStrings = new List<string>();
        public object ThreadLocker;
        public ObservableCollection<string> FilteredListOfPsms { get; private set; }
        public ICollectionView PeptideSpectralMatchesView;
        private readonly DataTable propertyView;
        public DownloadUniProtDatabaseWindow()
        {
            InitializeComponent();
            BindingOperations.EnableCollectionSynchronization(filteredStrings, ThreadLocker);

            propertyView = new DataTable();
            propertyView.Columns.Add("Name", typeof(string));
            propertyView.Columns.Add("Value", typeof(string));
            //dataGridProperties.DataContext = propertyView.DefaultView;

            dataGridScanNums.DataContext = PeptideSpectralMatchesView;
        }

        private void downloadProteomeButton_Click(object sender, RoutedEventArgs e)
        {
            
        }

        private void lookUpProteome_Click(object sender, RoutedEventArgs e)
        {

        }

        //any time text is changed in the search box, we want to filter the list displayed in the gridview
        private void proteomesSearchBox_TextChanged(object sender, TextChangedEventArgs e)
        {
            string searchText = (sender as TextBox).Text;
            FilterProteomesByString(searchText);
        }

        public void FilterProteomesByString(string searchText)
        {
            if (searchText == "")
            {
                //PeptideSpectralMatchesView.Filter = null;
            }
            else
            {
                //PeptideSpectralMatchesView.Filter = obj =>
                //{
                //    PsmFromTsv psm = obj as PsmFromTsv;
                //    return ((psm.Ms2ScanNumber.ToString()).StartsWith(searchString) || psm.FullSequence.ToUpper().Contains(searchString.ToUpper()));
                //};
            }
        }

        private void dataGridScanNums_SelectedCellsChanged(object sender, SelectedCellsChangedEventArgs e)
        {

        }
    }
}
