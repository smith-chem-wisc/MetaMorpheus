using EngineLayer;
using FlashLFQ;
using System;
using System.Collections.ObjectModel;
using System.ComponentModel;
using System.IO;
using System.Linq;
using System.Windows;
using System.Windows.Input;
using MetaMorpheusGUI.ForDisplayingInDataGrids;
using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Drawing;
using System.Runtime.CompilerServices;
using System.Windows.Forms;
using System.Windows.Controls;

namespace MetaMorpheusGUI.Views
{
    /// <summary>
    /// Interaction logic for SelectedPrecursorWindow.xaml
    /// </summary>
    public partial class SelectedPrecursorWindow : Window
    {
        private string outputFolder = null; // path to the experimental design file
        private string dataFileName = null;
        private string precursorFilePath = null; // path to the precursor file
        private const string precursorFileName = "_selected_precursors.tsv";
        private BindingList<SelectedPrecursorForDataGrid> PrecursorsFromDataGrid;


        public SelectedPrecursorWindow(string dataFilePath)
        {
            InitializeComponent();
            outputFolder = Path.GetDirectoryName(dataFilePath);
            dataFileName = Path.GetFileNameWithoutExtension(dataFilePath);
            precursorFilePath = Path.Combine(outputFolder, dataFileName + precursorFileName);
            if (File.Exists(Path.Combine(precursorFilePath)))
            {
                // read the existing precursor file and populate the data grid
                List<string> errors;
                var precursors = SelectedPrecursors.ReadSelectedPrecursors(precursorFilePath, out errors);
                if (errors.Count > 0)
                {
                    MessageBox.Show("Errors reading the selected precursors file:\n" + string.Join("\n", errors));
                }
                PrecursorsFromDataGrid = new BindingList<SelectedPrecursorForDataGrid>(precursors.Select(p => new SelectedPrecursorForDataGrid(p.Mz, p.Charge, p.RtStartInMinutes, p.RtEndInMinutes)).ToList());
            }
            else
            {
                PrecursorsFromDataGrid = new();// create a new empty list
            }
           
            PrecursorDataGrid.ItemsSource = PrecursorsFromDataGrid;
        }

        private void AddRow_Click(object sender, RoutedEventArgs e)
        {
            // add a new row to the data grid
            PrecursorsFromDataGrid.Add(new SelectedPrecursorForDataGrid());
        }

        private void Paste(object sender, ExecutedRoutedEventArgs e)
        {
            var selectedCell = PrecursorDataGrid.SelectedCells.FirstOrDefault(); // never null

             // this is the row index of the selected cell
            int columnIndex = PrecursorDataGrid.Columns.IndexOf(selectedCell.Column);

            // get data from clipboard in text format
            // clipboardRawData will be null if the data on the clipboard is not text
            object clipboardRawData = System.Windows.Clipboard.GetDataObject().GetData(DataFormats.Text);

            if (clipboardRawData != null)
            {
                string pastedText = clipboardRawData as string;

                // each line is delimited by a newline and/or return character
                var pastedLines = pastedText.Split(new char[] { '\r', '\n' }).Where(p => !string.IsNullOrWhiteSpace(p)).ToList();

                for (int i = 0; i < pastedLines.Count; i++)
                {
                    var pastedLine = pastedLines[i];

                    // each element in the line is delimited by tabs or commas
                    var pastedCells = pastedLine.Split(new char[] { '\t', ',' });

                    try
                    {
                        double mz = double.Parse(pastedCells[0]);
                        int charge = int.Parse(pastedCells[1]);
                        double rtStart = double.Parse(pastedCells[2]);
                        double rtEnd = double.Parse(pastedCells[3]);
                        PrecursorsFromDataGrid.Add(new SelectedPrecursorForDataGrid(mz, charge, rtStart, rtEnd));
                    }
                    catch (Exception)
                    {
                        // don't really need to print a warning
                    }
                }
            }
        }

        private void SavePrecursors_Click(object sender, RoutedEventArgs e)
        {
            if (outputFolder == null)
            {
                // no spectra files
                DialogResult = true;
                return;
            }

            try
            {
                var precursors = PrecursorsFromDataGrid.Select(p => new SelectedPrecursors.PrecursorInfo(p.Mz, p.Charge, p.RtStartInMinutes, p.RtEndInMinutes)).ToList();
                SelectedPrecursors.WriteSelectedPrecursorsToFile(precursors, outputFolder, dataFileName + precursorFileName);
            }
            catch (Exception ex)
            {
                MessageBox.Show("Could not save the precursor list!\n\n" + ex.Message);
                return;
            }

            DialogResult = true;
        }

        private void CancelPrecursors_Click(object sender, RoutedEventArgs e)
        {
            DialogResult = false;
        }

        private void KeyPressed(object sender, KeyEventArgs e)
        {
            if (e.Key == Key.Escape)
            {
                CancelPrecursors_Click(sender, e);
            }
        }


    }
}
