using GuiFunctions;
using GuiFunctions.MetaDraw;
using System;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Input;
using System.Windows.Media;

namespace MetaMorpheusGUI
{
    /// <summary>
    /// Interaction logic for ChimeraAnalysisTabView.xaml
    /// </summary>
    public partial class ChimeraAnalysisTabView : UserControl
    {

        public ChimeraAnalysisTabView()
        {
            InitializeComponent();

            // Update legend automatically if legend plotting settings change
            MetaDrawSettingsViewModel.Instance.PropertyChanged += (s, e) =>
            {
                if (e.PropertyName is nameof(MetaDrawSettingsViewModel.ChimeraLegendMainTextType) or nameof(MetaDrawSettingsViewModel.ChimeraLegendSubTextType) or nameof(MetaDrawSettingsViewModel.DisplayChimeraLegend) or nameof(MetaDrawSettingsViewModel.ChimeraLegendTakeFirstIfAmbiguous)or nameof(MetaDrawSettingsViewModel.ChimeraLegendMaxWidth))
                {
                    if (DataContext is ChimeraAnalysisTabViewModel { SelectedChimeraGroup: not null } context)
                        context.LegendCanvas = new(context.SelectedChimeraGroup);

                    AttachLegendCanvasEvents();
                }
            };

        }

        private void ChermicDataGrid_OnSelectedCellsChanged(object sender, SelectedCellsChangedEventArgs e)
        {
            var dataContext = DataContext as ChimeraAnalysisTabViewModel;
            if (chermicDataGrid.SelectedItem == null || sender == null || dataContext == null)
            {
                return;
            }

            ChimeraGroupViewModel chimeraGroup = dataContext.SelectedChimeraGroup;
            if (chimeraGroup == null)
            {
                MessageBox.Show("No chimera group found for this PSM");
                return;
            }
            dataContext.Ms1ChimeraPlot = new Ms1ChimeraPlot(ms1ChimeraOverlaPlot, chimeraGroup);
            dataContext.ChimeraSpectrumMatchPlot = new ChimeraSpectrumMatchPlot(ms2ChimeraPlot, chimeraGroup);
            dataContext.ChimeraDrawnSequence = dataContext.ChimeraDrawnSequence is null
                    ? new ChimeraDrawnSequence(chimeraSequenceCanvas, chimeraGroup, dataContext)
                    : dataContext.ChimeraDrawnSequence.UpdateData(chimeraGroup)
                ;
            AttachLegendCanvasEvents();

            // Reposition the initial location of the legend to be the upper right of the ms2 plot. 
            if (!legendHasBeenMoved)
                Dispatcher.BeginInvoke(new Action(() =>
                {
                    var ms2Plot = ms2ChimeraPlot;
                    var legend = ChimeraLegend;

                    // If you use a child legend canvas, get it:
                    var legendCanvas = legend.Children.Count > 0 ? legend.Children[0] as Canvas : null;
                    if (legendCanvas == null) return;

                    // Get MS2 plot position relative to the parent grid
                    var transform = ms2Plot.TransformToAncestor(MainGrid as Visual);
                    var ms2PlotTopLeft = transform.Transform(new Point(0, 0));

                    // Calculate top-right position for the legend
                    double legendLeft = ms2PlotTopLeft.X + ms2Plot.ActualWidth - SelectionDataGrid.ActualWidth - legendCanvas.Width - 16; // 10px margin
                    double legendTop = ms2PlotTopLeft.Y + 10; // 10px margin from top

                    Canvas.SetLeft(legendCanvas, legendLeft);
                    Canvas.SetTop(legendCanvas, legendTop);
                }));
        }

        #region Movable Legend

        private bool legendHasBeenMoved = false;
        private bool isDragging = false;
        private Point clickPosition; 
        private ChimeraLegendCanvas _legendCanvas;

        private void AttachLegendCanvasEvents()
        {
            var dataContext = DataContext as ChimeraAnalysisTabViewModel;
            var legendCanvas = dataContext?.LegendCanvas;

            // Preserve last position if it exists
            double lastLeft = 0, lastTop = 0;
            if (_legendCanvas != null && _legendCanvas.Parent is Canvas)
            {
                lastLeft = Canvas.GetLeft(_legendCanvas);
                lastTop = Canvas.GetTop(_legendCanvas);
            }
            ChimeraLegend.Children.Clear();

            if (!MetaDrawSettings.DisplayChimeraLegend)
                return;

            if (legendCanvas != null)
            {
                ChimeraLegend.Children.Add(dataContext?.LegendCanvas);
                legendCanvas.UpdateLayout();
                double parentWidth = ChimeraLegend.ActualWidth;
                double parentHeight = ChimeraLegend.ActualHeight;
                double legendWidth = legendCanvas.ActualWidth;
                double legendHeight = legendCanvas.ActualHeight;

                // Check if previous position is valid and within bounds
                bool validLeft = !double.IsNaN(lastLeft) && lastLeft >= 0 && lastLeft + legendWidth <= parentWidth;
                bool validTop = !double.IsNaN(lastTop) && lastTop >= 0 && lastTop + legendHeight <= parentHeight;

                if (validLeft && validTop && (lastLeft != 0 || lastTop != 0))
                {
                    Canvas.SetLeft(legendCanvas, lastLeft);
                    Canvas.SetTop(legendCanvas, lastTop);
                }
                else
                {
                    double initialLeft = Math.Max(0, parentWidth - legendWidth - 10);
                    Canvas.SetLeft(legendCanvas, initialLeft);
                    Canvas.SetTop(legendCanvas, 10);
                }

                legendCanvas.MouseLeftButtonDown += Legend_MouseLeftButtonDown;
                legendCanvas.MouseLeftButtonUp += Legend_MouseLeftButtonUp;
                legendCanvas.MouseMove += Legend_MouseMove;

                _legendCanvas = legendCanvas;
            }
        }

        private void Legend_MouseLeftButtonDown(object sender, MouseButtonEventArgs e)
        {
            isDragging = true;
            clickPosition = e.GetPosition(_legendCanvas);
            _legendCanvas.CaptureMouse();
        }

        private void Legend_MouseLeftButtonUp(object sender, MouseButtonEventArgs e)
        {
            isDragging = false;
            _legendCanvas.ReleaseMouseCapture();
        }

        private void Legend_MouseMove(object sender, MouseEventArgs e)
        {
            if (isDragging && _legendCanvas.Parent is Canvas parent)
            {
                Point mousePos = e.GetPosition(parent);
                double left = mousePos.X - clickPosition.X;
                double top = mousePos.Y - clickPosition.Y;

                // Clamp to parent bounds
                double maxLeft = parent.ActualWidth - _legendCanvas.ActualWidth;
                double maxTop = parent.ActualHeight - _legendCanvas.ActualHeight;

                if (left < 0) left = 0;
                if (top < 0) top = 0;
                if (left > maxLeft) left = maxLeft;
                if (top > maxTop) top = maxTop;

                Canvas.SetLeft(_legendCanvas, left);
                Canvas.SetTop(_legendCanvas, top);

                // Mark as moved
                legendHasBeenMoved = true;
            }
        }

        #endregion
    }
}
