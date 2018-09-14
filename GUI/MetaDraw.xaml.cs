using System.Linq;
using System.Windows;
using System.Windows.Controls;
using System.Collections.ObjectModel;
using System.IO;
using ViewModels;
using EngineLayer;
using System.Collections.Generic;
using MassSpectrometry;
using TaskLayer;
using System.ComponentModel;
using System.Threading.Tasks;
using System.Windows.Threading;
using System;
using System.Data;
using System.Windows.Data;
using System.Windows.Media;
using System.Windows.Media.Imaging;

namespace MetaMorpheusGUI
{
    /// <summary>
    /// Interaction logic for MetaDraw.xaml
    /// </summary>
    public partial class MetaDraw : Window
    {
        private PsmAnnotationViewModel mainViewModel;
        private MyFileManager spectraFileManager;
        private MsDataFile MsDataFile;
        private readonly ObservableCollection<MetaDrawPsm> peptideSpectralMatches;
        ICollectionView peptideSpectralMatchesView;
        private readonly DataTable propertyView;
        private string spectraFilePath;
        private string tsvResultsFilePath;

        public MetaDraw()
        {
            InitializeComponent();

            mainViewModel = new PsmAnnotationViewModel();
            plotView.DataContext = mainViewModel;
            peptideSpectralMatches = new ObservableCollection<MetaDrawPsm>();
            propertyView = new DataTable();
            propertyView.Columns.Add("Name", typeof(string));
            propertyView.Columns.Add("Value", typeof(string));
            peptideSpectralMatchesView = CollectionViewSource.GetDefaultView(peptideSpectralMatches);
            dataGridScanNums.DataContext = peptideSpectralMatchesView;
            dataGridProperties.DataContext = propertyView.DefaultView;
            Title = "MetaDraw: version " + GlobalVariables.MetaMorpheusVersion;
            spectraFileManager = new MyFileManager(true);

        }

        private void Window_Drop(object sender, DragEventArgs e)
        {
            string[] files = ((string[])e.Data.GetData(DataFormats.FileDrop)).OrderBy(p => p).ToArray();

            if (files != null)
            {
                foreach (var draggedFilePath in files)
                {
                    if (File.Exists(draggedFilePath))
                    {
                        LoadFile(draggedFilePath);
                    }
                }
            }
        }

        private void LoadFile(string filePath)
        {
            var theExtension = Path.GetExtension(filePath).ToLower();

            switch (theExtension)
            {
                case ".raw":
                case ".mzml":
                    spectraFilePath = filePath;
                    spectraFileNameLabel.Text = filePath;
                    break;
                case ".psmtsv":
                    tsvResultsFilePath = filePath;
                    psmFileNameLabel.Text = filePath;
                    break;
                default:
                    MessageBox.Show("Cannot read file type: " + theExtension);
                    break;
            }
        }

        private void LoadPsms(string filename)
        {
            string fileNameWithExtension = Path.GetFileName(spectraFilePath);
            string fileNameWithoutExtension = Path.GetFileNameWithoutExtension(spectraFilePath);

            foreach (var psm in TsvResultReader.ReadTsv(filename))
            {
                if (psm.FileName == fileNameWithExtension || psm.FileName == fileNameWithoutExtension)
                {
                    Dispatcher.BeginInvoke(new Action(() =>
                    {
                        peptideSpectralMatches.Add(psm);
                    }));
                }
            }
        }

        private void DrawPsm(int oneBasedScanNumber, string fullSequence = null)
        {
            MsDataScan msDataScanToDraw = MsDataFile.GetOneBasedScan(oneBasedScanNumber);
            IEnumerable<MetaDrawPsm> scanPsms = peptideSpectralMatches.Where(p => p.ScanNum == oneBasedScanNumber);

            if (fullSequence != null)
            {
                scanPsms = scanPsms.Where(p => p.FullSequence == fullSequence);
            }

            MetaDrawPsm psmToDraw = scanPsms.FirstOrDefault();

            // draw annotated spectrum
            mainViewModel.DrawPeptideSpectralMatch(msDataScanToDraw, psmToDraw);

            // draw annotated base sequence
            DrawAnnotatedBaseSequence(psmToDraw);
        }

        /// <summary>
        /// Event triggers when a different cell is selected in the PSM data grid
        /// </summary>
        private void dataGridScanNums_SelectedCellsChanged(object sender, SelectedCellsChangedEventArgs e)
        {
            if (dataGridScanNums.SelectedItem == null)
            {
                return;
            }

            // draw the selected PSM
            propertyView.Clear();
            MetaDrawPsm row = (MetaDrawPsm)dataGridScanNums.SelectedItem;
            System.Reflection.PropertyInfo[] temp = row.GetType().GetProperties();

            for (int i = 4; i < temp.Length; i++)
            {
                propertyView.Rows.Add(temp[i].Name, temp[i].GetValue(row, null));
            }
            dataGridProperties.Items.Refresh();
            DrawPsm(row.ScanNum, row.FullSequence);
        }

        private void selectSpectraFileButton_Click(object sender, RoutedEventArgs e)
        {
            Microsoft.Win32.OpenFileDialog openFileDialog1 = new Microsoft.Win32.OpenFileDialog
            {
                Filter = "Spectra Files(*.raw;*.mzML)|*.raw;*.mzML",
                FilterIndex = 1,
                RestoreDirectory = true,
                Multiselect = false
            };
            if (openFileDialog1.ShowDialog() == true)
            {
                foreach (var filePath in openFileDialog1.FileNames.OrderBy(p => p))
                {
                    LoadFile(filePath);
                }
            }
        }

        private void selectPsmFileButton_Click(object sender, RoutedEventArgs e)
        {
            Microsoft.Win32.OpenFileDialog openFileDialog1 = new Microsoft.Win32.OpenFileDialog
            {
                Filter = "Result Files(*.psmtsv)|*.psmtsv",
                FilterIndex = 1,
                RestoreDirectory = true,
                Multiselect = false
            };
            if (openFileDialog1.ShowDialog() == true)
            {
                foreach (var filePath in openFileDialog1.FileNames.OrderBy(p => p))
                {
                    LoadFile(filePath);
                }
            }
        }

        private async void loadFilesButton_Click(object sender, RoutedEventArgs e)
        {
            // check for validity
            propertyView.Clear();
            if (spectraFilePath == null)
            {
                MessageBox.Show("Please add a spectra file.");
                return;
            }

            if (tsvResultsFilePath == null)
            {
                MessageBox.Show("Please add a search result file.");
                return;
            }

            // load the spectra file
            (sender as Button).IsEnabled = false;
            selectSpectraFileButton.IsEnabled = false;
            selectPsmFileButton.IsEnabled = false;
            prgsFeed.IsOpen = true;
            prgsText.Content = "Loading spectra file...";

            var slowProcess = Task<MsDataFile>.Factory.StartNew(() => spectraFileManager.LoadFile(spectraFilePath, null, null, false, false, new CommonParameters()));
            await slowProcess;
            MsDataFile = slowProcess.Result;

            // load the PSMs
            this.prgsText.Content = "Loading PSMs...";
            await Task.Run(() => LoadPsms(tsvResultsFilePath));

            // done loading - restore controls
            this.prgsFeed.IsOpen = false;
            (sender as Button).IsEnabled = true;
            selectSpectraFileButton.IsEnabled = true;
            selectPsmFileButton.IsEnabled = true;
        }

        private void TextBox_TextChanged(object sender, TextChangedEventArgs e)
        {
            string txt = (sender as TextBox).Text;
            if (txt == "")
            {
                peptideSpectralMatchesView.Filter = null;
            }
            else
            {
                peptideSpectralMatchesView.Filter = obj =>
                {
                    MetaDrawPsm psm = obj as MetaDrawPsm;
                    return ((psm.ScanNum.ToString()).StartsWith(txt) || psm.FullSequence.ToUpper().Contains(txt.ToUpper()));
                };
            }
        }

        private void DrawAnnotatedBaseSequence(MetaDrawPsm psm)
        {
            double spacing = 22;
            BaseDraw.clearCanvas(canvas);

            // draw base sequence
            for (int r = 0; r < psm.BaseSequence.Length; r++)
            {
                BaseDraw.txtDrawing(canvas, new Point(r * spacing + 10, 10), psm.BaseSequence[r].ToString(), Brushes.Black);
            }

            // draw b ions
            foreach (var bIon in psm.FragmentIons.Where(p => p.ProductType == ProductType.B))
            {
                int residue = bIon.IonNumber;
                BaseDraw.botSplittingDrawing(canvas, new Point(residue * spacing + 8, 50), Colors.Blue, bIon.ProductType.ToString().ToLower() + bIon.IonNumber);
            }

            // draw c ions
            foreach (var cIon in psm.FragmentIons.Where(p => p.ProductType == ProductType.C))
            {
                int residue = psm.BaseSequence.Length - cIon.IonNumber;
                BaseDraw.botSplittingDrawing(canvas, new Point(residue * spacing + 8, 50), Colors.Gold, cIon.ProductType.ToString().ToLower() + cIon.IonNumber);
            }

            // draw y ions
            foreach (var yIon in psm.FragmentIons.Where(p => p.ProductType == ProductType.Y))
            {
                int residue = psm.BaseSequence.Length - yIon.IonNumber;
                BaseDraw.topSplittingDrawing(canvas, new Point(residue * spacing + 8, 0), Colors.Purple, yIon.ProductType.ToString().ToLower() + yIon.IonNumber);
            }

            // draw zdot ions
            foreach (var zDotIon in psm.FragmentIons.Where(p => p.ProductType == ProductType.Zdot))
            {
                int residue = zDotIon.IonNumber;
                BaseDraw.topSplittingDrawing(canvas, new Point(residue * spacing + 8, 0), Colors.Orange, zDotIon.ProductType.ToString().ToLower() + zDotIon.IonNumber);
            }

            // draw modifications
            int aa = 0;
            bool currentlyReadingMod = false;
            for (int c = 0; c < psm.FullSequence.Length; c++)
            {
                switch (psm.FullSequence[c])
                {
                    case '[':
                        currentlyReadingMod = true;
                        BaseDraw.circledTxtDraw(canvas, new Point(aa * spacing - 17, 12), Brushes.Yellow);
                        break;
                    case ']':
                        currentlyReadingMod = false;
                        break;
                    default:
                        if (!currentlyReadingMod)
                        {
                            aa++;
                        }
                        break;
                }
            }
        }

        private void dataGridProperties_SelectedCellsChanged(object sender, SelectedCellsChangedEventArgs e)
        {
            (sender as DataGrid).UnselectAll();
        }

        private void PDFButton_Click(object sender, RoutedEventArgs e)
        {
            if (dataGridScanNums.SelectedCells.Count == 0)
            {
                MessageBox.Show("Please select at least one scan to export");
            }

            int num = dataGridScanNums.SelectedItems.Count;
            string writeDirectory = Path.Combine(Directory.GetParent(tsvResultsFilePath).FullName, "PDF");

            if (!Directory.Exists(writeDirectory))
            {
                Directory.CreateDirectory(writeDirectory);
            }

            foreach (object selectedItem in dataGridScanNums.SelectedItems)
            {
                MetaDrawPsm psm = (MetaDrawPsm)selectedItem;
                ExportToPdf(psm, Path.Combine(writeDirectory, psm.ScanNum + "_" + psm.FullSequence + ".pdf"));
            }

            dataGridScanNums.SelectedItem = dataGridScanNums.SelectedItem;
            MessageBox.Show(string.Format("{0} PDFs exported to " + writeDirectory, num));
        }

        private void ExportToPdf(MetaDrawPsm psm, string path)
        {
            System.Reflection.PropertyInfo[] temp = psm.GetType().GetProperties();

            for (int i = 4; i < temp.Length; i++)
            {
                propertyView.Rows.Add(temp[i].Name, temp[i].GetValue(psm, null));
            }
            dataGridProperties.Items.Refresh();
            DrawPsm(psm.ScanNum, psm.FullSequence);
            
            double wid = 0;
            dataGridProperties.HorizontalScrollBarVisibility = ScrollBarVisibility.Disabled;
            dataGridProperties.VerticalScrollBarVisibility = ScrollBarVisibility.Disabled;
            foreach (DataGridColumn col in dataGridProperties.Columns)
            {
                wid += col.ActualWidth;
            }
            PDFOutPut.Background = Brushes.White;
            PDFOutPut.ColumnDefinitions[0].Width = new GridLength(wid + 10);
            PDFOutPut.Measure(new Size(wid + gbPSM.ActualWidth + 10, 600));
            PDFOutPut.Arrange(new Rect(new Size(wid + gbPSM.ActualWidth + 10, 600)));
            dataGridProperties.Measure(new Size(wid + 22, 600));
            dataGridProperties.Arrange(new Rect(new Size(wid + 5, 600)));

            dataGridProperties.Arrange(new Rect(new Size(wid + 5, 600)));
            var rtb = new RenderTargetBitmap((int)(wid + gbPSM.ActualWidth) + 11, 600, 96, 96, PixelFormats.Pbgra32);

            rtb.Render(PDFOutPut);
            BitmapFrame bf = BitmapFrame.Create(rtb);

            var encoder = new BmpBitmapEncoder();
            encoder.Frames.Add(bf);
            using (var stream = new MemoryStream())
            {
                encoder.Save(stream);
                var img = System.Drawing.Bitmap.FromStream(stream);
                PdfWriter.WriteToPdf(img, (int)(wid + gbPSM.ActualWidth) + 11, 600, path);
            }

            dataGridProperties.HorizontalScrollBarVisibility = ScrollBarVisibility.Visible;
            dataGridProperties.VerticalScrollBarVisibility = ScrollBarVisibility.Visible;
            PDFOutPut.ColumnDefinitions[0].Width = new GridLength(1, GridUnitType.Star);
        }
    }
}
