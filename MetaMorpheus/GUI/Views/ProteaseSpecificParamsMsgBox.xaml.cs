using System.Windows;

namespace MetaMorpheusGUI
{
    /// <summary>
    /// Interaction logic for CustomMsgBox.xaml
    /// </summary>
    public partial class ProteaseSpecificMsgBox : Window
    {
        public ProteaseSpecificMsgBox(string title, string caption)
        {
            InitializeComponent();
            Title = title;
            Label.Content = caption;
        }

        //the defaults are used when the window is closed via the top right corner
        bool UseSettings = false;
        bool AskAgain = true;

        public static (bool UseSettings, bool AskAgain) Show(string title, string caption)
        {
            var msgBox = new ProteaseSpecificMsgBox(title, caption);
            msgBox.ShowDialog();

            return (msgBox.UseSettings, msgBox.AskAgain);
        }

        private void YesButton_Click(object sender, RoutedEventArgs e)
        {
            UseSettings = true;
            AskAgain = !DoNotAskAgainCheckBox.IsChecked.Value;
            Close();
        }

        private void NoButton_Click(object sender, RoutedEventArgs e)
        {
            UseSettings = false;
            AskAgain = !DoNotAskAgainCheckBox.IsChecked.Value;
            Close();
        }
    }
}