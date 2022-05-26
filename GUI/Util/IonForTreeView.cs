using GuiFunctions;
using OxyPlot;
using Proteomics.Fragmentation;
using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.ComponentModel;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MetaMorpheusGUI
{
    public class IonForTreeView : INotifyPropertyChanged
    {
        #region Public Properties

        public ProductType ProductType { get; set; }    
        public string IonName { get { return ProductType.ToString(); } }
        public ObservableCollection<OxyColor> PossibleColors { get; set; }
        public OxyColor SelectedColor { get; set; }
        public string BoxName { get; set; }

        #endregion

        #region Constructor

        public IonForTreeView(ProductType type, OxyColor selected = new OxyColor())
        {
            ProductType = type;
            BoxName = IonName + "ComboBox";
            SelectedColor = selected;
            PossibleColors = new ObservableCollection<OxyColor>();
            foreach (var color in MetaDrawSettings.PossibleColors)
            {
                PossibleColors.Add(color.Key);
            }
        }

        #endregion

        #region Commands

        #endregion

        #region PropertyChanged

        public event PropertyChangedEventHandler PropertyChanged;

        #endregion

    }
}
