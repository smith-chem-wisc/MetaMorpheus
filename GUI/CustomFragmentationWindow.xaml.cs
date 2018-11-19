using Proteomics.Fragmentation;
using System.Collections.Generic;
using System.Windows;
using System.ComponentModel;
using MassSpectrometry;
using System.Linq;
using System.Collections.ObjectModel;

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
            var productTypes = DissociationTypeCollection.ProductsFromDissociationType.Values.SelectMany(p => p).Distinct();

            foreach (ProductType p in productTypes)
            {
                TheList.Add(new BoolStringClass { IsSelected = false, Type = p });
            }

            ProductTypeList.ItemsSource = TheList;
        }

        private void Save_Click(object sender, RoutedEventArgs e)
        {
            var selectedIons = TheList.Where(p => p.IsSelected).Select(p => p.Type);
            DissociationTypeCollection.ProductsFromDissociationType[DissociationType.Custom] = selectedIons.ToList();
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

    public class BoolStringClass 
    {
        public ProductType Type { get; set; }
        public bool IsSelected { get; set; }
    }
}
