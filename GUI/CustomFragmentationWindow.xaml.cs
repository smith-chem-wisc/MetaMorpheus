using EngineLayer;
using Nett;
using Proteomics.Fragmentation;
using System.Collections.Generic;
using System.Windows;
using TaskLayer;
using System.IO;
using System.ComponentModel;

namespace MetaMorpheusGUI
{
    /// <summary>
    /// Interaction logic for CustomFragmentationIons.xaml
    /// </summary>
    public partial class CustomFragmentationWindow : Window
    {
        private CustomFragmentation CustomFragmentation = new CustomFragmentation();

        public CustomFragmentationWindow() : this(null)
        {
        }

        public CustomFragmentationWindow(List<ProductType> list)
        {
            InitializeComponent();
            // update checkboxes
            if (list != null)
            {
                bCheckBox.IsChecked = list.Contains(ProductType.b);
                cCheckBox.IsChecked = list.Contains(ProductType.c);
                yCheckBox.IsChecked = list.Contains(ProductType.y);
                zdotCheckBox.IsChecked = list.Contains(ProductType.zDot);
            }

            base.Closing += this.OnClosing;
        }

        private void Save_Click(object sender, RoutedEventArgs e)
        {
            CustomFragmentation.BIons = bCheckBox.IsChecked.Value;
            CustomFragmentation.CIons = cCheckBox.IsChecked.Value;
            CustomFragmentation.YIons = yCheckBox.IsChecked.Value;
            CustomFragmentation.ZDotIons = zdotCheckBox.IsChecked.Value;

            // creates temporary toml file to store custom type settings
            Toml.WriteFile(CustomFragmentation, Path.Combine(GlobalVariables.DataDir, @"customDissociationType.toml"), MetaMorpheusTask.tomlConfig);
            DialogResult = true;
        }

        private void Cancel_Click(object sender, RoutedEventArgs e)
        {
           DialogResult = true;
        }

        private void OnClosing(object sender, CancelEventArgs e)
        {
            this.Visibility = Visibility.Hidden;
            e.Cancel = true;
        }
    }
}
