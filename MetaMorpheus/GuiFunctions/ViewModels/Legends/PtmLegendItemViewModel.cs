using OxyPlot;

namespace GuiFunctions
{
    /// <summary>
    /// View Model class for each mod in the ptm legend
    /// </summary>
    public class PtmLegendItemViewModel : LegendItemViewModel
    {

        #region Constructor

        public PtmLegendItemViewModel(string modName)
        {
            Name = modName;
            OxyColor color = MetaDrawSettings.ModificationTypeToColor[modName];
            ColorBrush = DrawnSequence.ParseColorBrushFromOxyColor(color);
        }

        #endregion
    }
}
