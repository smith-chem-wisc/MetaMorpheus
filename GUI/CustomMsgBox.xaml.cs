using System.Windows;

namespace MetaMorpheusGUI
{
    /// <summary>
    /// Interaction logic for CustomMsgBox.xaml
    /// </summary>
    public partial class CustomMsgBox : Window
    {
        public CustomMsgBox(string title, string caption, string btn1, string btn2, string btn3)
        {
            InitializeComponent();
            this.Title = title;
            this.Label.Content = caption;
            this.Button1.Content = btn1;
            this.Button2.Content = btn2;
            this.Button3.Content = btn3;

            Button1.Click += new RoutedEventHandler(Button1_Click);
            Button2.Click += new RoutedEventHandler(Button2_Click);
            Button3.Click += new RoutedEventHandler(Button3_Click);
        }

        static CustomMsgBox MsgBox;
        static MessageBoxResult result = MessageBoxResult.No;

        public static MessageBoxResult Show(string title, string caption, string btn1, string btn2, string btn3)
        {
            MsgBox = new CustomMsgBox(title, caption, btn1, btn2, btn3);
            MsgBox.ShowDialog();
            return result;
        }

        private void Button1_Click(object sender, RoutedEventArgs e)
        {
            result = MessageBoxResult.Yes;
            MsgBox.Close();
        }

        private void Button2_Click(object sender, RoutedEventArgs e)
        {
            result = MessageBoxResult.No;
            MsgBox.Close();
        }

        private void Button3_Click(object sender, RoutedEventArgs e)
        {
            result = MessageBoxResult.OK;
            MsgBox.Close();
        }
    }
}