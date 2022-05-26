using GuiFunctions;
using OxyPlot;
using Proteomics.Fragmentation;
using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MetaMorpheusGUI
{
    public class IonForTreeDisplayModel : IonForTreeViewModel
    {
        /// <summary>
        /// static instance with example data for design time editing
        /// </summary>
        public static IonForTreeViewModel Instance => new IonForTreeDisplayModel();

        public IonForTreeDisplayModel() : base (((ProductType[])Enum.GetValues(typeof(ProductType))).First(), false)
        {
            IonType = ((ProductType[])Enum.GetValues(typeof(ProductType))).First();
            IonName = IonType.ToString();
            PossibleColors = new ObservableCollection<string>(MetaDrawSettings.PossibleColors.Values.ToList());
            AddSpaces(PossibleColors);
            IsBeta = false;
            HasChanged = false;
            if (IsBeta)
                SelectedColor = AddSpaces(MetaDrawSettings.BetaProductTypeToColor[IonType].GetColorName());
            else
                SelectedColor = AddSpaces(MetaDrawSettings.ProductTypeToColor[IonType].GetColorName());
        }
    }
}
