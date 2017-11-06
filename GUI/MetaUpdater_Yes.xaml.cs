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
using System.Net;
using System.Diagnostics;
namespace MetaMorpheusGUI
{
    /// <summary>
    /// Interaction logic for MetaUpdater_Yes.xaml
    /// </summary>
    public partial class MetaUpdater_Yes : Window
    {
        public MetaUpdater_Yes()
        {
            InitializeComponent();
        }

        private void InstallerClicked(object sender, RoutedEventArgs e)
        {
            this.Close();
            using (var client = new WebClient())
            {
                var uri = new Uri("https://github.com/smith-chem-wisc/MetaMorpheus/releases/download/" + MainWindow.VersionCheck + "/MetaMorpheusInstaller.msi");
                client.DownloadFileAsync(uri, "MetaMorpheusInstaller.msi");//downloaded succeed!
                client.DownloadFileCompleted += (sender1, e1) =>
                {

                    Process p = new Process();
                    //p.StartInfo.FileName = @"E:\XRSheeranQ\0.0.201\MetaMorpheus-master\installer\Debug\installer.msi"; //installer version
                    p.StartInfo.FileName = "MetaMorpheusInstaller.msi"; //zip version

                    //Get user feedback on installation
                    MessageBoxResult result1 = MessageBox.Show("Download completed! Program will be terminated, are you sure to exit and install?", "Metamorpheus Auto Update", MessageBoxButton.YesNoCancel);
                    if (result1.Equals(MessageBoxResult.Yes))
                    {

                        System.Windows.Application.Current.Shutdown();
                        //auto remove or choose to remove for zip<-
                        //Or installer will take care of the rest
                        p.Start();
                    }
                    else
                    {
                       

                    }


                };
            }
        }

        private void FilesClicked(object sender, RoutedEventArgs e)
        {
            this.Close();
            using (var client = new WebClient())
            {
                var uri = new Uri("https://github.com/smith-chem-wisc/MetaMorpheus/releases/download/" + MainWindow.VersionCheck + "/MetaMorpheusGuiDotNetFrameworkAppveyor.zip");
                client.DownloadFileAsync(uri, "MetaMorpheusGuiDotNetFrameworkAppveyor.zip");//downloaded succeed!
                client.DownloadFileCompleted += (sender1, e1) =>
                {

                    Process p = new Process();
                    p.StartInfo.FileName = "MetaMorpheusGuiDotNetFrameworkAppveyor.zip"; //zip version

                    //Get user feedback on installation
                    MessageBoxResult result1 = MessageBox.Show("Download completed! Program will be terminated, are you sure to exit and install?", "Metamorpheus Auto Update", MessageBoxButton.YesNoCancel);
                    if (result1.Equals(MessageBoxResult.Yes))
                    {

                        System.Windows.Application.Current.Shutdown();
                        //auto remove or choose to remove for zip<-
                        //Or installer will take care of the rest
                        p.Start();
                    }
                    else
                    {


                    }


                };
            }
        }

        private void Button_Click_2(object sender, RoutedEventArgs e)
        {
            this.Close();
            MessageBox.Show("please check github for more information");
        }
    }
}
