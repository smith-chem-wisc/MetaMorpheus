using GuiFunctions;
using OxyPlot;
using Proteomics.Fragmentation;
using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.ComponentModel;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Controls;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Shapes;

namespace MetaMorpheusGUI
{
    public class IonForTreeViewModel : BaseViewModel
    {
        #region Private Properties

        protected string _selectedColor;

        #endregion

        #region Public Properties

        public ProductType IonType { get; set; }    
        public string IonName { get; set; }
        public ObservableCollection<string> PossibleColors { get; set; }
        public string SelectedColor
        {
            get { return _selectedColor; }
            set 
            { 
                _selectedColor = value;
                OnPropertyChanged(nameof(SelectedColor));
            }
        }
        SolidColorBrush ColorBrush { get; set; }
            
        public bool IsBeta { get; set; }
        public bool HasChanged { get; set; } = false;

        #endregion

        #region Constructor

        public IonForTreeViewModel(ProductType type, bool beta)
        {
            IonType = type;
            IonName = IonType.ToString();
            PossibleColors = new ObservableCollection<string>(MetaDrawSettings.PossibleColors.Values.ToList());
            AddSpaces(PossibleColors);
            IsBeta = beta;
            OxyColor color = MetaDrawSettings.BetaProductTypeToColor[IonType];

            if (IsBeta)
                SelectedColor = AddSpaces(color.GetColorName());
            else
                SelectedColor = AddSpaces(color.GetColorName());

            var colorVal = color.ToByteString().Split(',');
            ColorBrush = new SolidColorBrush(System.Windows.Media.Color.FromArgb(Byte.Parse(colorVal[0]), Byte.Parse(colorVal[1]), Byte.Parse(colorVal[2]), Byte.Parse(colorVal[3])));

        }


        ///// <summary>
        ///// static instance with example data for design time editing
        ///// </summary>
        //public static IonForTreeViewModel Instance => new IonForTreeViewModel();

        //public IonForTreeViewModel()
        //{
        //    IonType = ((ProductType[])Enum.GetValues(typeof(ProductType))).First();
        //    PossibleColors = new ObservableCollection<string>(MetaDrawSettings.PossibleColors.Values.ToList());
        //    AddSpaces(PossibleColors);
        //    IsBeta = false;

        //    if (IsBeta)
        //        SelectedColor = AddSpaces(MetaDrawSettings.BetaProductTypeToColor[IonType].GetColorName());
        //    else
        //        SelectedColor = AddSpaces(MetaDrawSettings.ProductTypeToColor[IonType].GetColorName());
        //}

        #endregion

        public void SelectionChanged(string newColor)
        {
            SelectedColor = newColor;
            HasChanged = true;
        }

        protected string AddSpaces(string text)
        {
            if (string.IsNullOrWhiteSpace(text))
                return "";
            StringBuilder newText = new StringBuilder(text.Length * 2);
            newText.Append(text[0]);
            for (int i = 1; i < text.Length; i++)
            {
                if (char.IsUpper(text[i]) && text[i - 1] != ' ')
                    newText.Append(' ');
                newText.Append(text[i]);
            }
            return newText.ToString();
        }

        protected void AddSpaces(ObservableCollection<string> strings)
        {
            for (int i = 0; i < strings.Count; i++)
            {
                strings[i] = AddSpaces(strings[i]);
            }
        }

    }
}
