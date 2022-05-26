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
using System.Windows.Input;

namespace MetaMorpheusGUI
{
    public class IonForTreeViewModel : BaseViewModel
    {
        #region Private Properties

        private string _selectedColor;

        #endregion

        #region Public Properties

        public ProductType IonType { get; set; }    
        public string IonName { get { return IonType.ToString(); } }
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
            
        public bool IsBeta { get; set; }
        public bool HasChanged { get; set; } = false;

        #endregion

        #region Constructor

        public IonForTreeViewModel(ProductType type, bool beta)
        {
            IonType = type;
            PossibleColors = new ObservableCollection<string>(MetaDrawSettings.PossibleColors.Values.ToList());
            AddSpaces(PossibleColors);
            IsBeta = beta;

            if (IsBeta)
                SelectedColor = AddSpaces(MetaDrawSettings.BetaProductTypeToColor[type].GetColorName());
            else
                SelectedColor = AddSpaces(MetaDrawSettings.ProductTypeToColor[type].GetColorName());

        }

        #endregion

        public void SelectionChanged(string newColor)
        {
            SelectedColor = newColor;
            HasChanged = true;
        }

        private string AddSpaces(string text)
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

        private void AddSpaces(ObservableCollection<string> strings)
        {
            for (int i = 0; i < strings.Count; i++)
            {
                strings[i] = AddSpaces(strings[i]);
            }
        }

    }
}
