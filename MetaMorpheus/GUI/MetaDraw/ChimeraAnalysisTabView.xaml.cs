using GuiFunctions;
using GuiFunctions.MetaDraw.Chimeras;
using System;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Input;

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

            MetaDrawSettingsViewModel.Instance.PropertyChanged += (s, e) =>
            {
                if (e.PropertyName is nameof(MetaDrawSettingsViewModel.ChimeraLegendMainTextType) or nameof(MetaDrawSettingsViewModel.ChimeraLegendSubTextType) or nameof(MetaDrawSettingsViewModel.DisplayChimeraLegend))
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
            dataContext.ChimeraDrawnSequence =
                  dataContext.ChimeraDrawnSequence is null ?
                 new ChimeraDrawnSequence(chimeraSequenceCanvas, chimeraGroup, dataContext)
                 : dataContext.ChimeraDrawnSequence.UpdateData(chimeraGroup)
                 ;
            AttachLegendCanvasEvents();
        }

        #region Movable Legend

        private bool isDragging = false;
        private Point clickPosition; 
        private Canvas _legendCanvas;

        private void AttachLegendCanvasEvents()
        {
            var dataContext = DataContext as ChimeraAnalysisTabViewModel;
            var legendCanvas = dataContext?.LegendCanvas;

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
                ChimeraLegend.Children.Add(legendCanvas);

                double parentWidth = ChimeraLegend.ActualWidth;
                double parentHeight = ChimeraLegend.ActualHeight;
                double legendWidth = legendCanvas.Width;
                double legendHeight = legendCanvas.Height;

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
            }
        }

        #endregion
    }
}
