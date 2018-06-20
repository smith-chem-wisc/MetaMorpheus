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

namespace RealTimeGUI
{
    /// <summary>
    /// Interaction logic for MainWindow.xaml
    /// </summary>
    public partial class MainWindow : Window
    {
        public MainWindow()
        {
            InitializeComponent();
            DataReceiver.DataReceiverNotificationEventHandler += UpdateTbNotification;
        }

        private void UpdateTbNotification(object sender, NotificationEventArgs e)
        {
            TbNotifications.AppendText(e.Notification);
        }

        private void BtnConnection_Click(object sender, RoutedEventArgs e)
        {
            var x = Connection.DoJob();
            TbNotifications.AppendText(x);
        }

        private void BtnRealTimeData_Click(object sender, RoutedEventArgs e)
        {
            DataReceiver dataReceiver = new DataReceiver();
            dataReceiver.DoJob();
        }
    }
}
