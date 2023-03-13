using System;
using System.ComponentModel;
using System.Net.Http;
using System.Windows;


namespace MetaMorpheusGUI
{
    /// <summary>
    /// Interaction logic for ProgressBarIndeterminate.xaml
    /// </summary>
    public partial class ProgressBarIndeterminate : Window , INotifyPropertyChanged
    {
        public event PropertyChangedEventHandler PropertyChanged;
        private void NotifyPropertyChanged(String info)
        {
            if (PropertyChanged != null)
            {
                PropertyChanged(this, new PropertyChangedEventArgs(info));
            }
        }
        string bytes;
        public string Bytes
        {
            get
            {
                return bytes;
            }
            set
            {
                bytes = value;
                NotifyPropertyChanged(nameof(Bytes));
            }
        }

        HttpClient client;
        public bool open;

        public ProgressBarIndeterminate()
        {
            InitializeComponent();
        }

        public ProgressBarIndeterminate(HttpClient Client)
        {
            InitializeComponent();
            client = Client;
            open = true;
            Bytes = "0KB";
        }

        private void Cancel_Click(object sender, RoutedEventArgs e)
        {
            try
            {
                client.CancelPendingRequests();
            }
            catch
            {

            }
            client.Dispose();
            open = false;
            Close();
        }
        protected override void OnClosing(CancelEventArgs e)
        {
            try
            {
                client.CancelPendingRequests();
            }
            catch
            {
                
            }
            open = false;
            client.Dispose();
            e.Cancel = false;
        }

    }
}
