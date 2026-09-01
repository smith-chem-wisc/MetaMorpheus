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
    /// Interaction logic for FragmentReanalysisControl.xaml
    /// </summary>
    public partial class FragmentReanalysisControl : UserControl
    {
        private MetaDraw parent;
        public FragmentReanalysisControl()
        {
            InitializeComponent();
        }

        internal void LinkMetaDraw(MetaDraw metaDraw)
        {
            parent = metaDraw;
        }

        private void SearchWithNewIons_OnClick(object sender, RoutedEventArgs e)
        {
            parent.SearchWithNewIons_OnClick(sender, e);
        }
    }
}
