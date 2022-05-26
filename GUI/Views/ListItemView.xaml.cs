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
    /// Interaction logic for ListItemView.xaml
    /// </summary>
    public partial class ListItemView : UserControl
    {
        public ListItemView()
        {
            InitializeComponent();
        }

        private void ComboBox_SelectionChanged(object sender, SelectionChangedEventArgs e)
        {
            ((IonForTreeViewModel)((ComboBox)sender).DataContext).SelectionChanged((string)((ComboBox)sender).SelectedItem);
        }
    }
}
