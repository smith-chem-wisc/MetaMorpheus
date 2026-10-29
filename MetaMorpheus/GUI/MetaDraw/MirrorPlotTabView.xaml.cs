using GuiFunctions;
using GuiFunctions.MetaDraw;

namespace MetaMorpheusGUI
{
    public partial class MirrorPlotTabView : System.Windows.Controls.UserControl
    {
        public MirrorPlotTabView()
        {
            InitializeComponent();

            MirrorPlotView.DataContextChanged += (s, e) =>
            {
                if (DataContext is MirrorPlotTabViewModel vm)
                {
                    vm.MirrorPlotView = MirrorPlotView;
                    vm.MirrorPlotExportElement = MirrorPlotExportPanel;
                    vm.TopCanvasExport = TopSequenceCanvas;
                    vm.BottomCanvasExport = BottomSequenceCanvas;
                }
            };
        }

        private void LeftDataGrid_SelectedCellsChanged(object sender, System.Windows.Controls.SelectedCellsChangedEventArgs e)
        {
            if (DataContext is MirrorPlotTabViewModel vm && vm.SelectedLeftPsm != null && vm.SelectedRightPsm != null)
            {
                vm.RefreshSequences();
            }
        }

        private void RightDataGrid_SelectedCellsChanged(object sender, System.Windows.Controls.SelectedCellsChangedEventArgs e)
        {
            if (DataContext is MirrorPlotTabViewModel vm && vm.SelectedLeftPsm != null && vm.SelectedRightPsm != null)
            {
                vm.RefreshSequences();
            }
        }
    }
}
