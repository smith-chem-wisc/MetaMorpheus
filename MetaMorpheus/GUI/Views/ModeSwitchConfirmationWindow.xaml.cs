using System.ComponentModel;
using System.Configuration;
using System.Windows;
using GuiFunctions.Util;

namespace MetaMorpheusGUI
{
    /// <summary>
    /// Interaction logic for ModeSwitchConfirmationWindow.xaml
    /// </summary>
    public partial class ModeSwitchConfirmationWindow : Window
    {
        public ModeSwitchRequestEventArgs Result { get; init; }

        public ModeSwitchConfirmationWindow(ModeSwitchRequestEventArgs args)
        {
            Result = args;
            InitializeComponent();
        }

        private void YesKeepFiles_Click(object sender, RoutedEventArgs e)
        {
            Result.Result = ModeSwitchResult.SwitchKeepFiles;
            Result.RememberMyDecision = RememberCheckBox.IsChecked.Value;
            Close();
        }

        private void YesRemoveFiles_Click(object sender, RoutedEventArgs e)
        {
            Result.Result = ModeSwitchResult.SwitchRemoveFiles;
            Result.RememberMyDecision = RememberCheckBox.IsChecked.Value;
            Close();
        }

        private void No_Click(object sender, RoutedEventArgs e)
        {
            Result.Result = ModeSwitchResult.Cancel;
            Result.RememberMyDecision = RememberCheckBox.IsChecked.Value;
            Close();
        }
    }
}
