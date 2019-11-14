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

            YesButton.Click += new RoutedEventHandler(YesButton_Click);
            NoButton.Click += new RoutedEventHandler(NoButton_Click);
        }

        static ProteaseSpecificMsgBox MsgBox;
        //the defaults are used when the window is closed via the top right corner
        static bool UseSettings = false;
        static bool AskAgain = true;

        public static (bool UseSettings, bool AskAgain) Show(string title, string caption)
        {
            MsgBox = new ProteaseSpecificMsgBox(title, caption);
            MsgBox.ShowDialog();

            return (UseSettings, AskAgain);
        }

        private void YesButton_Click(object sender, RoutedEventArgs e)
        {
            UseSettings = true;
            AskAgain = !DoNotAskAgainCheckBox.IsChecked.Value;
            MsgBox.Close();
        }

        private void NoButton_Click(object sender, RoutedEventArgs e)
        {
            UseSettings = false;
            AskAgain = !DoNotAskAgainCheckBox.IsChecked.Value;
            MsgBox.Close();
        }
    }
}