using GuiFunctions;
using GuiFunctions.MetaDraw;
using System.Windows;
using System.Windows.Controls;

namespace MetaMorpheusGUI
{
    public partial class MirrorPlotTabView : UserControl
    {
        public MirrorPlotTabView()
        {
            InitializeComponent();

            MirrorPlotView.DataContextChanged += (s, e) =>
            {
                if (DataContext is MirrorPlotTabViewModel vm)
                {
                    vm.MirrorPlotView = MirrorPlotView;
                    vm.TopCanvasExport = TopSequenceCanvas;
                    vm.BottomCanvasExport = BottomSequenceCanvas;
                }
            };
        }

        private void LeftDataGrid_SelectedCellsChanged(object sender, SelectedCellsChangedEventArgs e)
        {
            if (DataContext is MirrorPlotTabViewModel vm && vm.SelectedLeftPsm != null && vm.SelectedRightPsm != null)
            {
                RebuildSequences(vm);
            }
        }

        private void RightDataGrid_SelectedCellsChanged(object sender, SelectedCellsChangedEventArgs e)
        {
            if (DataContext is MirrorPlotTabViewModel vm && vm.SelectedLeftPsm != null && vm.SelectedRightPsm != null)
            {
                RebuildSequences(vm);
            }
        }

        private void RebuildSequences(MirrorPlotTabViewModel vm)
        {
            TopSequenceCanvas.Children.Clear();
            BottomSequenceCanvas.Children.Clear();

            var topSequence = new DrawnSequence(TopSequenceCanvas, vm.SelectedLeftPsm, false);
            topSequence.AnnotateBaseSequence(
                vm.SelectedLeftPsm.BaseSeq,
                vm.SelectedLeftPsm.FullSequence,
                10,
                vm.SelectedLeftPsm.MatchedIons,
                vm.SelectedLeftPsm);

            var bottomSequence = new DrawnSequence(BottomSequenceCanvas, vm.SelectedRightPsm, false);
            bottomSequence.AnnotateBaseSequence(
                vm.SelectedRightPsm.BaseSeq,
                vm.SelectedRightPsm.FullSequence,
                10,
                vm.SelectedRightPsm.MatchedIons,
                vm.SelectedRightPsm);
        }
    }
}
