using Proteomics.Fragmentation;
using System.Collections.Generic;
using System.Windows;
using System.ComponentModel;
using MassSpectrometry;
using System.Linq;
using System.Collections.ObjectModel;
using System;

namespace MetaMorpheusGUI
{
    /// <summary>
    /// Interaction logic for CustomFragmentationIons.xaml
    /// </summary>
    public partial class CustomFragmentationWindow : Window
    {
        private ObservableCollection<BoolStringClass> TheList { get; set; }

        public CustomFragmentationWindow() : this(null)
        {
        }

        public CustomFragmentationWindow(List<ProductType> list)
        {
            InitializeComponent();
            PopulateChoices();

            // update checkboxes
            if (list != null)
            {
                var res = TheList.Where(n => list.Contains(n.Type));
                foreach(var r in res)
                {
                    r.IsSelected = true;
                }               
            }

            base.Closing += this.OnClosing;
        }

        private void PopulateChoices()
        {
            TheList = new ObservableCollection<BoolStringClass>();

            var knownProductTypes = Enum.GetValues(typeof(ProductType)).Cast<ProductType>().ToList();
            knownProductTypes.Remove(ProductType.D);
            knownProductTypes.Remove(ProductType.M);
            knownProductTypes.Remove(ProductType.zPlusOne);

            foreach (ProductType productType in knownProductTypes)
            {
                TheList.Add(new BoolStringClass {
                    IsSelected = false,
                    Type = productType,
                    ToolTip = DissociationTypeCollection.GetMassShiftFromProductType(productType).ToString("F4") + " Da; " 
                        + TerminusSpecificProductTypes.ProductTypeToFragmentationTerminus[productType] + " terminus" });
            }

            ProductTypeList.ItemsSource = TheList;
        }

        private void Save_Click(object sender, RoutedEventArgs e)
        {
            var selectedIons = TheList.Where(p => p.IsSelected).Select(p => p.Type);
            DissociationTypeCollection.ProductsFromDissociationType[DissociationType.Custom] = selectedIons.ToList();
            this.Visibility = Visibility.Hidden;
        }

        private void Cancel_Click(object sender, RoutedEventArgs e)
        {
            this.Visibility = Visibility.Hidden;
        }

        private void OnClosing(object sender, CancelEventArgs e)
        {
            // this window only closes when task is added
            if (this.Visibility == Visibility.Visible)
            { 
                this.Visibility = Visibility.Hidden;
                e.Cancel = true;
            }
        }
    }

    public class BoolStringClass 
    {
        public ProductType Type { get; set; }
        public bool IsSelected { get; set; }
        public string ToolTip { get; set; }
    }
}
