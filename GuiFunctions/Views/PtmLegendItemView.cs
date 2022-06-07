using OxyPlot;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Media;

namespace GuiFunctions
{
    /// <summary>
    /// View Model class for each mod in the ptm legend
    /// </summary>
    public class PtmLegendItemView : BaseView
    {
        #region Private Properties

        private SolidColorBrush _colorBrush;

        #endregion

        #region Public Properties

        public string ModName { get; }
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

        public PtmLegendItemView(string modName)
        {
            ModName = modName;
            OxyColor color = MetaDrawSettings.ModificationTypeToColor[modName];
            ColorBrush = DrawnSequence.ParseColorBrushFromOxyColor(color);
        }

        #endregion
    }
}
