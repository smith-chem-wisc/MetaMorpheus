using System;
using System.Linq;
using System.Windows;
using System.Collections.ObjectModel;
using System.IO;

namespace MetaMorpheusGUI
{
    /// <summary>
    /// Interaction logic for QuantSettingWindow.xaml
    /// </summary>
    public partial class QuantSettingWindow : Window
    {
        private readonly ObservableCollection<QuantForDataGrid> spectraFilesQuantSets = new ObservableCollection<QuantForDataGrid>();

        public ObservableCollection<QuantForDataGrid> SpectraFilesQuantSets { get; set; }
        public bool DoQuant { get; set; }

        public QuantSettingWindow(ObservableCollection<RawDataForDataGrid> spectraFilesObservableCollection)
        {
            InitializeComponent();
            foreach (var item in spectraFilesObservableCollection.Where(p => p.Use == true).Select(p => p.FilePath))
            {
                spectraFilesQuantSets.Add(new QuantForDataGrid(item));
            }
            DgQuant.DataContext = spectraFilesQuantSets;

        }

        private void BtnSaveQuant_Click(object sender, RoutedEventArgs e)
        {
            foreach (var item in spectraFilesQuantSets)
            {
                //Add more conditions
                if (item.FileName == null || item.FilePath ==null || 
                    item.QbioRep == null || item.Qcondition == null || item.Qfraction == null ||item.QtechRep ==null)
                {
                    MessageBox.Show("Please check quantification setting");
                    return;
                }
            }

            //Do something in GUI.Main
            SpectraFilesQuantSets = spectraFilesQuantSets;
            DoQuant = CbQuant.IsChecked.Value;

            //Write Quan Setting
            try
            {
                WriteQuantSetToTsv(System.IO.Path.Combine(txtOuantSetOut.Text, "QuantSet.txt"), SpectraFilesQuantSets, CbQuant.IsChecked.Value);
            }
            catch (Exception)
            {

                throw;
            }
           
            
            DialogResult = true;
        }

        private void BtnCancelQuant_Click(object sender, RoutedEventArgs e)
        {
            DialogResult = false;
        }

        public static void WriteQuantSetToTsv(string filePath, ObservableCollection<QuantForDataGrid> SpectraFilesQuantSets, bool doQuant)
        {
            using (StreamWriter output = new StreamWriter(filePath))
            {
                output.WriteLine("Quantification" + doQuant.ToString());
                foreach (var heh in SpectraFilesQuantSets)
                {
                    output.WriteLine(heh.FileName + 
                        "\t" + heh.QbioRep + 
                        "\t" + heh.Qcondition + 
                        "\t" + heh.Qfraction + 
                        "\t" + heh.QtechRep + 
                        "\t" + heh.FilePath);
                }
            }
        }
    }
}
