using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using EngineLayer;
using GuiFunctions.ViewModels.Legends;

namespace GuiFunctions
{
    public class ChimeraLegendViewModel : LegendViewModel
    {
        #region Private Properties

        private Visibility sharedIonStackPanelVisibility;

        #endregion

        #region Public Properties

        public ChimeraLegendItemViewModel SharedIons { get; set; }

        public Visibility SharedIonStackPanelVisibility
        {
            get => sharedIonStackPanelVisibility;
            set { sharedIonStackPanelVisibility = value; OnPropertyChanged(nameof(SharedIonStackPanelVisibility)); }
        }

        #endregion

        public ChimeraLegendViewModel(List<PsmFromTsv> chimericIDs, string baseSequence, double offset = 0) : base()
        {
            Header = baseSequence;
            if (legendItems.Count > 2)
            {
                SharedIons = new ChimeraLegendItemViewModel("Shared Ions", ChimeraSpectrumMatchPlot.MultipleProteinSharedColor);
                SharedIonStackPanelVisibility = Visibility.Visible;
            }
            else
            {
                SharedIonStackPanelVisibility = Visibility.Collapsed;
            }
            legendItems.ForEach(p => LegendItemViewModels.Add(p));
            TopOffset = offset;
            Visibility = Visibility.Visible;
        }

        public void ParseLegendItemsFromPsms(List<PsmFromTsv> psms)
        {

        }
    }
}
