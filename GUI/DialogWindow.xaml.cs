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
using System.Windows.Shapes;

using MetaMorpheusGUI;



namespace MetaMorpheusGUI
{
    /// <summary>
    /// Interaction logic for DialogWindow.xaml
    /// </summary>
    public partial class DialogWindow : Window
    {
        internal string stringSuffix = "";

        public DialogWindow()
        {
            InitializeComponent();
        }

        private void OnKeyDownHandler(object sender, KeyEventArgs e)
        {
            if (e.Key == Key.Return)
            {
                if(!string.IsNullOrWhiteSpace(metaMorpheusTaskSuffixTextbox.Text))
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
    }
}
