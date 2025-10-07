using OxyPlot;
using OxyPlot.Wpf;
using System.Windows.Media;

namespace GuiFunctions
{
    public class ColorForTreeViewModel : BaseViewModel
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

        public ColorForTreeViewModel(string name, OxyColor defaultColor)
        {
            Name = name;
            SelectedColor = AddSpaces(defaultColor.GetColorName());
            ColorBrush = DrawnSequence.ParseColorBrushFromOxyColor(defaultColor);
        }

        public ColorForTreeViewModel(string name, SolidColorBrush brush)
        {
            Name = name;
            ColorBrush = brush;
            SelectedColor = AddSpaces(brush.ToOxyColor().GetColorName());
        }

        #endregion

        public void SelectionChanged(string newColor)
        {
            SelectedColor = newColor;
            HasChanged = true;
        }
    }

    public class CoverageTypeForTreeViewModel(string name)
        : ColorForTreeViewModel(name, MetaDrawSettings.CoverageTypeToColor[name]);
}
