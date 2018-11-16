using EngineLayer;
using Nett;
using Proteomics.Fragmentation;
using System.Collections.Generic;
using System.Windows;
using TaskLayer;
using System.IO;

namespace MetaMorpheusGUI
{
    /// <summary>
    /// Interaction logic for CustomFragmentationIons.xaml
    /// </summary>
    public partial class CustomDissociationTypeWindow : Window
    {
        private CustomDissociationType CustomDissociationType = new CustomDissociationType();

        public CustomDissociationTypeWindow() : this(null)
        {
        }

        public CustomDissociationTypeWindow(List<ProductType> list)
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
        }

        private void Save_Click(object sender, RoutedEventArgs e)
        {
            CustomDissociationType.BIons = bCheckBox.IsChecked.Value;
            CustomDissociationType.CIons = cCheckBox.IsChecked.Value;
            CustomDissociationType.YIons = yCheckBox.IsChecked.Value;
            CustomDissociationType.ZDotIons = zdotCheckBox.IsChecked.Value;

            // creates temporary toml file to store custom type settings
            Toml.WriteFile(CustomDissociationType, Path.Combine(GlobalVariables.DataDir, @"customDissociationType.toml"), MetaMorpheusTask.tomlConfig);
            DialogResult = true;
        }

        private void Cancel_Click(object sender, RoutedEventArgs e)
        {
            DialogResult = false;
        }
    }
}
