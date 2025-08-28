using System.Windows;
using System.Windows.Controls;
using GuiFunctions;

namespace MetaMorpheusGUI;

/// <summary>
/// Interaction logic for BioPolymerCoverageMapView.xaml
/// </summary>
public partial class BioPolymerCoverageMapView : UserControl
{
    private double _lastWidth = -1;
    public BioPolymerCoverageMapView()
    {
        InitializeComponent();
        this.SizeChanged += OnSizeChanged;
    }

    private BioPolymerCoverageMapViewModel ViewModel => DataContext as BioPolymerCoverageMapViewModel;

    private void OnSizeChanged(object sender, SizeChangedEventArgs e)
    {
        // Only update if the UserControl's size actually changed
        if (e.WidthChanged && this.ActualWidth > 0 && this.ActualWidth != _lastWidth)
        {
            _lastWidth = this.ActualWidth;
            ViewModel?.UpdateLettersPerRow(this.ActualWidth);
        }
    }
}
