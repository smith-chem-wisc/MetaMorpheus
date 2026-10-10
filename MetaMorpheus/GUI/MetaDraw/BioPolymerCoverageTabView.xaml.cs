using System;
using System.Linq;
using EngineLayer;
using System.Windows;
using System.Windows.Controls;
using GuiFunctions.MetaDraw;
using Readers;

namespace MetaMorpheusGUI;

/// <summary>
/// Interaction logic for BioPolymerCoverageTabView.xaml
/// </summary>
public partial class BioPolymerCoverageTabView : UserControl
{
    public BioPolymerCoverageTabView()
    {
        InitializeComponent();
    }

    private void SelectDatabaseButton_OnClick(object sender, RoutedEventArgs e)
    {
        var dataContext = DataContext as BioPolymerTabViewModel;
        if (sender == null || dataContext == null)
        {
            return;
        }

        string filterString = string.Join(";", GlobalVariables.AcceptedDatabaseFormats.Select(p => "*" + p));

        Microsoft.Win32.OpenFileDialog openFileDialog1 = new Microsoft.Win32.OpenFileDialog
        {
            Filter = "Database Files(" + filterString + ")|" + filterString,
            FilterIndex = 1,
            RestoreDirectory = true,
            Multiselect = true
        };

        if (openFileDialog1.ShowDialog() == true)
        {
            foreach (var filePath in openFileDialog1.FileNames)
            {
                if (GlobalVariables.AcceptedDatabaseFormats.Contains(GlobalVariables.GetFileExtension(filePath).ToLowerInvariant()))
                {
                    if (!dataContext.DatabasePaths.Contains(filePath))
                    {
                        dataContext.DatabasePaths.Add(filePath);
                    }
                }
                else
                {
                    MessageBox.Show("Cannot read file type: " + GlobalVariables.GetFileExtension(filePath).ToLowerInvariant());
                }
            }
            
            dataContext.OnPropertyChanged(nameof(dataContext.DatabaseName));
            dataContext.OnPropertyChanged(nameof(dataContext.DatabasePathsTooltip));
        }
    }

    private void SelectProteinGroupsButton_OnClick(object sender, RoutedEventArgs e)
    {

        // ==================================================
        // TODO: Change this behavior once Mzlib supports all of our protein/transcript group files (Working on that PR)
        // ==================================================



        var dataContext = DataContext as BioPolymerTabViewModel;
        if (sender == null || dataContext == null)
            return;

        var openFileDialog = new Microsoft.Win32.OpenFileDialog
        {
            Filter = "Protein Group Files|*.tsv|All Files|*.*",
            FilterIndex = 1,
            RestoreDirectory = true,
            Multiselect = false
        };

        if (openFileDialog.ShowDialog() != true)
            return;

        string filePath = openFileDialog.FileName;
        try
        {
            if (!BioPolymerTabViewModel.SupportedProteinGroupFileTypes.Contains(filePath.ParseFileType()))
            {
                MessageBox.Show("Cannot read file type: " + GlobalVariables.GetFileExtension(filePath).ToLowerInvariant());
                return;
            }

            dataContext.ProteinGroupsFilePath = filePath;
        }
        catch (Exception)
        {
            MessageBox.Show("Cannot read file type: " + GlobalVariables.GetFileExtension(filePath).ToLowerInvariant());
        }
    }

    private void BioPolymerDataGrid_OnSelectionChanged(object sender, SelectionChangedEventArgs e)
    {
        //throw new System.NotImplementedException();
    }
}
