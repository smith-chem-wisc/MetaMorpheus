using System.Linq;
using System.Windows;
using System.Windows.Controls;
using GuiFunctions;

namespace MetaMorpheusGUI
{
    /// <summary>
    /// Interaction logic for FragmentationParamsControl.xaml
    /// </summary>
    public partial class FragmentationParamsControl : UserControl
    {
        public FragmentationParamsControl()
        {
            InitializeComponent();
        }

        private void AddCustomMIonLoss_Click(object sender, RoutedEventArgs e)
        {
            var window = new CustomMIonLossWindow();
            if (window.ShowDialog() == true)
            {
                // Refresh the data context to reload loss
                if (DataContext is FragmentationParamsViewModel viewModel)
                {
                    // Trigger a reload of the M-Ion loss
                    viewModel.ReloadMIonLosses();

                    foreach (var loss in window.CreatedLosses)
                    {
                        var vm = viewModel.AvailableMIonLosses.FirstOrDefault(p => p.Name == loss.Name);
                        if (vm is not null)
                            vm.IsSelected = true;
                    }
                }
            }
        }
    }
}
