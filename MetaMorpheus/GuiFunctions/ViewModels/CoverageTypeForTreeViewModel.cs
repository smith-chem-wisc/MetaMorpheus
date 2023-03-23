using GuiFunctions;
using OxyPlot;
using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Media;

namespace GuiFunctions
{
    public class CoverageTypeForTreeViewModel : BaseViewModel
    {
        #region Private Properties

        protected string _selectedColor;
        protected SolidColorBrush _colorBrush;

        #endregion

        #region Public Properties

        public string Name { get; set; }
        public bool HasChanged { get; set; } = false;
        public string SelectedColor
        {
            get { return _selectedColor; }
            set
            {
                _selectedColor = value;
                ColorBrush = DrawnSequence.ParseColorBrushFromName(_selectedColor);
                OnPropertyChanged(nameof(SelectedColor));
            }
        }
        public SolidColorBrush ColorBrush
        {
            get { return _colorBrush; }
            set
            {
                _colorBrush = value;
                OnPropertyChanged(nameof(ColorBrush));
            }
        }

        #endregion

        #region Constructor

        public CoverageTypeForTreeViewModel(string name)
        {
            Name = name;
            OxyColor color = MetaDrawSettings.CoverageTypeToColor[name];
            SelectedColor = AddSpaces(color.GetColorName());
            ColorBrush = DrawnSequence.ParseColorBrushFromOxyColor(color);
        }

        #endregion

        public void SelectionChanged(string newColor)
        {
            SelectedColor = newColor;
            HasChanged = true;
        }
    }
}
