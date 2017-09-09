using System.Windows;
using System.Windows.Input;

namespace MetaMorpheusGUI
{
    /// <summary>
    /// Interaction logic for DialogWindow.xaml
    /// </summary>
    public partial class DialogWindow : Window
    {
        #region Internal Fields

        internal string stringSuffix = "";

        #endregion Internal Fields

        #region Public Constructors

        public DialogWindow()
        {
            InitializeComponent();
        }

        #endregion Public Constructors

        #region Private Methods

        private void OnKeyDownHandler(object sender, KeyEventArgs e)
        {
            if (e.Key == Key.Return)
            {
                if (!string.IsNullOrWhiteSpace(metaMorpheusTaskSuffixTextbox.Text))
                {
                    stringSuffix = "-" + metaMorpheusTaskSuffixTextbox.Text;
                }
                this.Close();
            }
        }

        private void DialogWindowLoaded(object sender, RoutedEventArgs e)
        {
            metaMorpheusTaskSuffixTextbox.Focus();
        }

        #endregion Private Methods
    }
}