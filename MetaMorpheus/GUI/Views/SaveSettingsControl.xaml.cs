using GuiFunctions;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Navigation;
using System.Windows.Shapes;

namespace MetaMorpheusGUI
{
    /// <summary>
    /// Interaction logic for SaveSettingsControl.xaml
    /// </summary>
    public partial class SaveSettingsControl : UserControl
    {
        public SaveSettingsControl()
        {
            InitializeComponent();
        }

        /// <summary>
        /// Open the SaveSettings Window
        /// </summary>
        /// <param name="sender"></param>
        /// <param name="e"></param>
        private void SaveAs_Click(object sender, RoutedEventArgs e)
        {
            var taskSettingsViewModel = DataContext as TaskSettingViewModel 
                ?? throw new ArgumentException($"DataContext is of type {DataContext.GetType()}, TaskSettingsViewModel is required.");


            SaveSettingsWindow newWindow = new SaveSettingsWindow(taskSettingsViewModel);

            // if user hit save button
            if (newWindow.ShowDialog() == true)
            {
                taskSettingsViewModel.SaveSettingsFromWindow();
            }
            // if user hit cancel
            else
            {

            }
        }
    }
}
