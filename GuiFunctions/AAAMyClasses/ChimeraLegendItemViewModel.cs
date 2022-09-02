using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Media;
using OxyPlot;

namespace GuiFunctions
{
    public class ChimeraLegendItemViewModel : BaseViewModel
    {
        #region Private Properties

        private SolidColorBrush _bColorBrush;
        private SolidColorBrush _yColorBrush;
        private SolidColorBrush _internalColorBrush;

        #endregion

        #region Public Properties

        public string Name { get; }
        public SolidColorBrush BColorBrush
        {
            get { return _bColorBrush; }
            set
            {
                _bColorBrush = value;
                OnPropertyChanged(nameof(BColorBrush));
            }
        }
        public SolidColorBrush YColorBrush
        {
            get { return _yColorBrush; }
            set
            {
                _yColorBrush = value;
                OnPropertyChanged(nameof(YColorBrush));
            }
        }
        public SolidColorBrush InternalColorBrush
        {
            get { return _internalColorBrush; }
            set
            {
                _internalColorBrush = value;
                OnPropertyChanged(nameof(InternalColorBrush));
            }
        }

        #endregion

        #region Constructor

        public ChimeraLegendItemViewModel(string mods, (OxyColor, OxyColor, OxyColor) colors)
        {
            if (mods == null || mods == "")
                Name = "No Modifications";
            else
                Name = mods;
            BColorBrush = DrawnSequence.ParseColorBrushFromOxyColor(colors.Item1);
            YColorBrush = DrawnSequence.ParseColorBrushFromOxyColor(colors.Item2);
            InternalColorBrush = DrawnSequence.ParseColorBrushFromOxyColor(colors.Item3);
        }


        #endregion

    }
}
