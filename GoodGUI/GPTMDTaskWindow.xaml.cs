using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
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
using IndexSearchAndAnalyze;
using MetaMorpheus;

namespace GoodGUI
{
    /// <summary>
    /// Interaction logic for GPTMDTaskWindow.xaml
    /// </summary>
    public partial class GPTMDTaskWindow : Window
    {
        private ObservableCollection<ModList> modFileList;
        private MyGPTMDtask myGPTMDtask;

        public GPTMDTaskWindow()
        {
            InitializeComponent();
        }

        public GPTMDTaskWindow(ObservableCollection<ModList> modFileList)
        {
            this.modFileList = modFileList;
        }

        public GPTMDTaskWindow(MyGPTMDtask myGPTMDtask, ObservableCollection<ModList> modFileList)
        {
            this.myGPTMDtask = myGPTMDtask;
            this.modFileList = modFileList;
        }

        public MyTask TheTask { get; private set; }
    }
}
