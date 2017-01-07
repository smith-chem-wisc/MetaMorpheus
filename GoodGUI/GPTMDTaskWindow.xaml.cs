using MetaMorpheus;
using System.Collections.ObjectModel;
using System.Windows;

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

        internal MyTask TheTask { get; set; }

        private void cancelButton_Click(object sender, RoutedEventArgs e)
        {
            DialogResult = false;
        }

        private void saveButton_Click(object sender, RoutedEventArgs e)
        {
            DialogResult = true;
        }
    }
}