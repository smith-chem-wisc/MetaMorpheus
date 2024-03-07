using OxyPlot;
using System.Windows.Media;
using Omics.Fragmentation;

namespace GuiFunctions
{
    public class IonForTreeViewModel : BaseViewModel
    {
        #region Private Properties

        protected string _selectedColor;
        protected SolidColorBrush _colorBrush;

        #endregion

        #region Public Properties

        public ProductType IonType { get; set; }    
        public string IonName { get; set; }
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

        public bool IsBeta { get; set; }
        public bool HasChanged { get; set; } = false;

        #endregion

        #region Constructors

        public IonForTreeViewModel(ProductType type, bool beta)
        {
            IonType = type;
            IonName = IonType.ToString() + " - Ion";
            IsBeta = beta;
            OxyColor color;
            if (IsBeta)
                color = MetaDrawSettings.BetaProductTypeToColor[IonType];
            else
                color = MetaDrawSettings.ProductTypeToColor[IonType];
            SelectedColor = AddSpaces(color.GetColorName());
            ColorBrush = DrawnSequence.ParseColorBrushFromOxyColor(color);
        }

        // should only be used for unannotated peak or internal ion color as it is not an actual product type
        public IonForTreeViewModel(string type, bool beta)
        {
            if (type.Equals("Unannotated Peak"))
            {
                IonName = type;
                IsBeta = beta;
                OxyColor color = MetaDrawSettings.UnannotatedPeakColor;
                SelectedColor = AddSpaces(color.GetColorName());
                ColorBrush = DrawnSequence.ParseColorBrushFromOxyColor(color);
            }
            else if (type.Equals("Internal Ion"))
            {
                IonName = type;
                IsBeta = beta;
                OxyColor color = MetaDrawSettings.InternalIonColor;
                SelectedColor = AddSpaces(color.GetColorName());
                ColorBrush = DrawnSequence.ParseColorBrushFromOxyColor(color);
            }
        }

        #endregion

        public void SelectionChanged(string newColor)
        {
            SelectedColor = newColor;
            HasChanged = true;
        }

    }
}
