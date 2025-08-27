using System.Windows;
using System.Windows.Controls;
using GuiFunctions;

namespace MetaMorpheusGUI;

/// <summary>
/// Interaction logic for BioPolymerCoverageMapView.xaml
/// </summary>
public partial class BioPolymerCoverageMapView : UserControl
{
    public BioPolymerCoverageMapView()
    {
        InitializeComponent();
        this.SizeChanged += OnSizeChanged;
    }

    private BioPolymerCoverageMapViewModel ViewModel => DataContext as BioPolymerCoverageMapViewModel;

    private void OnSizeChanged(object sender, SizeChangedEventArgs e)
    {
        ViewModel?.UpdateLettersPerRow(CoverageImage.ActualWidth);
    }
}
