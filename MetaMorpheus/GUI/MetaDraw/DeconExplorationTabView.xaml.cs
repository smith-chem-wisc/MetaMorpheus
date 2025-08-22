using GuiFunctions;
using System.Windows.Controls;

namespace MetaMorpheusGUI
{
    /// <summary>
    /// Interaction logic for DeconExplorationTabView.xaml
    /// </summary>
    public partial class DeconExplorationTabView : UserControl
    {
        public DeconExplorationTabView()
        {
            InitializeComponent();
        }

        private void DataGrid_OnSelectedCellsChanged(object sender, SelectedCellsChangedEventArgs e)
        {
            var dc = DataContext as DeconExplorationViewModel;
            dc!.RunDeconvolutionCommand.Execute(DeconPlot);
        }
    }
}
