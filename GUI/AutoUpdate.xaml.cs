using EngineLayer;
using Nett;
using Proteomics;
using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.IO;
using System.Linq;
using System.Threading;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Input;
using TaskLayer;
using System.Net;
using System.Text.RegularExpressions;
using System.Diagnostics;
namespace MetaMorpheusGUI
{
    /// <summary>
    /// Interaction logic for Window1.xaml
    /// </summary>
    public partial class Window1 : Window
    {
        public Window1()
        {
            InitializeComponent();
            //show version in label
            lbl.Content = "NewVersion " + MainWindow.VersionCheck + " is available, would you like to install?";
            
        }
        
        //checkbox checked event
        private void ckb_Checked(object sender, RoutedEventArgs e)
        {
            String LogPath = @"E:\XRSheeranQ\update.txt";
            File.WriteAllText(LogPath, "false");
        }
        //checkbox unchecked event
        private void ckb_Unchecked(object sender, RoutedEventArgs e)
        {
            String LogPath = @"E:\XRSheeranQ\update.txt";
            File.WriteAllText(LogPath, "true");
        }
        //Yes button
        private void Yes_Click(object sender, RoutedEventArgs e)
        {
            this.Close();
            using (var client = new WebClient())
            {
                var uri = new Uri("https://github.com/smith-chem-wisc/MetaMorpheus/releases/download/" + MainWindow.VersionCheck + "/MetaMorpheusGuiDotNetFrameworkAppveyor.zip");
                client.DownloadFileAsync(uri, @"E:\XRSheeranQ\MetaMorpheusGuiDotNetFrameworkAppveyor.zip");//downloaded succeed!
                client.DownloadFileCompleted += (sender1, e1) =>
                {
                    
                    Process p = new Process();
                    //p.StartInfo.FileName = @"E:\XRSheeranQ\0.0.201\MetaMorpheus-master\installer\Debug\installer.msi"; //installer version
                    p.StartInfo.FileName = @"E:\XRSheeranQ\MetaMorpheusGuiDotNetFrameworkAppveyor.zip"; //zip version

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
                        MessageBox.Show("Package downloaded in " + @"E:\XRSheeranQ\MetaMorpheusGuiDotNetFrameworkAppveyor.zip" + " Please update manually");
         
                    }


                };
            }

        }

        private void No_Click(object sender, RoutedEventArgs e)
        {
            this.Close();
            Window2 newwd = new Window2();
            newwd.Show();
        }
    }
}
