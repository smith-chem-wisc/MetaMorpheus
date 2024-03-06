using Omics.Fragmentation;
using System.Collections.Generic;
using System.Windows;
using System.ComponentModel;
using MassSpectrometry;
using System.Linq;
using System.Collections.ObjectModel;
using System;
using Omics.Fragmentation.Peptide;

namespace MetaMorpheusGUI
{
    /// <summary>
    /// Interaction logic for CustomFragmentationIons.xaml
    /// </summary>
    public partial class CustomFragmentationWindow : Window
    {
        private bool isRna;
        private ObservableCollection<BoolStringClass> TheList { get; set; }

        public CustomFragmentationWindow() : this(null)
        {
        }

        public CustomFragmentationWindow(List<ProductType> list, bool isRna = false)
        {
            this.isRna = isRna;
            InitializeComponent();
            PopulateChoices();

            // update checkboxes
            if (list != null)
            {
                var res = TheList.Where(n => list.Contains(n.Type));
                foreach (var r in res)
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
            knownProductTypes.Remove(ProductType.Y);
            knownProductTypes.Remove(ProductType.Ycore);
            knownProductTypes.Remove(ProductType.Ycore);
            knownProductTypes.Remove(ProductType.Y);

            if (isRna)
            {
                knownProductTypes.Remove(ProductType.aStar);
                knownProductTypes.Remove(ProductType.aDegree);
                knownProductTypes.Remove(ProductType.bAmmoniaLoss);
                knownProductTypes.Remove(ProductType.yAmmoniaLoss);
                knownProductTypes.Remove(ProductType.zPlusOne);
                knownProductTypes.Remove(ProductType.zDot);
            }
            else
            {
                knownProductTypes.Remove(ProductType.aBaseLoss);
                knownProductTypes.Remove(ProductType.aWaterLoss);
                knownProductTypes.Remove(ProductType.bBaseLoss);
                knownProductTypes.Remove(ProductType.cBaseLoss);
                knownProductTypes.Remove(ProductType.cWaterLoss);
                knownProductTypes.Remove(ProductType.d);
                knownProductTypes.Remove(ProductType.dBaseLoss);
                knownProductTypes.Remove(ProductType.dWaterLoss);

                knownProductTypes.Remove(ProductType.w);
                knownProductTypes.Remove(ProductType.wBaseLoss);
                knownProductTypes.Remove(ProductType.wWaterLoss);
                knownProductTypes.Remove(ProductType.xBaseLoss);
                knownProductTypes.Remove(ProductType.xWaterLoss);
                knownProductTypes.Remove(ProductType.yBaseLoss);
                knownProductTypes.Remove(ProductType.z);
                knownProductTypes.Remove(ProductType.zBaseLoss);
                knownProductTypes.Remove(ProductType.zWaterLoss);
            }

            foreach (ProductType productType in knownProductTypes)
            {
                var tooltip = isRna
                    ? Omics.Fragmentation.Oligo.DissociationTypeCollection
                        .GetRnaMassShiftFromProductType(productType).ToString("F4") + " Da; "
                    : DissociationTypeCollection.GetMassShiftFromProductType(productType).ToString("F4") +
                      " Da; "
                      + TerminusSpecificProductTypes.ProductTypeToFragmentationTerminus[productType] +
                      " terminus";


                TheList.Add(new BoolStringClass
                {
                    IsSelected = false,
                    Type = productType,
                    ToolTip = tooltip
                });
            }

            ProductTypeList.ItemsSource = TheList;
        }

        private void Save_Click(object sender, RoutedEventArgs e)
        {
            var selectedIons = TheList.Where(p => p.IsSelected).Select(p => p.Type);
            if (isRna)
                throw new NotImplementedException("No RNA just yet");
                //Omics.Fragmentation.Oligo.DissociationTypeCollection.ProductsFromDissociationType[DissociationType.Custom] = selectedIons.ToList();
            else
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