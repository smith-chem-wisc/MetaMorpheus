using MzLibUtil;
using System;
using System.Collections.Generic;
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
        public DownloadUniProtDatabaseWindow()
        {
            InitializeComponent();
            //string filepath = System.IO.Path.Combine(DataDir, @"Proteomes", @"availableUniProtProteomes.txt");
            //    //string downloadedFilePath = ProteinDbRetriever.DownloadAvailableUniProtProteomes(filepath);

            //    //string downloadedFilePath = @"E:\junk\availableUniProtProteomes.txt.gz";
            ////Dictionary<string, string> uniprotProteoms = ProteinDbRetriever.UniprotProteomesList(downloadedFilePath);
            //Dictionary<string, string> uniprotProteoms = ProteinDbRetriever.UniprotProteomesList(filepath);
        }
    }
}
