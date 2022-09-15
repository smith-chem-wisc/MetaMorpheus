using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using EngineLayer;
using GuiFunctions.ViewModels.Legends;
using Proteomics.ProteolyticDigestion;

namespace GuiFunctions
{
    public class ChimeraLegendViewModel : LegendViewModel
    {
        #region Private Properties

        private Visibility sharedIonStackPanelVisibility;
        private Dictionary<string, List<ChimeraLegendItemViewModel>> chimeraLegendItems;

        #endregion

        #region Public Properties

        public Dictionary<string, List<ChimeraLegendItemViewModel>> ChimeraLegendItems
        {
            get => chimeraLegendItems;
            set { chimeraLegendItems = value; OnPropertyChanged(nameof(ChimeraLegendItems)); }
        }

        public Visibility SharedIonStackPanelVisibility
        {
            get => sharedIonStackPanelVisibility;
            set { sharedIonStackPanelVisibility = value; OnPropertyChanged(nameof(SharedIonStackPanelVisibility)); }
        }

        #endregion

        public ChimeraLegendViewModel(List<PsmFromTsv> chimericIDs, double offset = 0) : base()
        {
            TopOffset = offset;
            Visibility = Visibility.Visible;
            ChimeraLegendItems = new();
            ParseLegendItemsFromPsms(chimericIDs);
        }

        public void ParseLegendItemsFromPsms(List<PsmFromTsv> chimericIDs)
        {
            var groupedByProtein = chimericIDs.GroupBy(p => p.BaseSeq);
            int proteinIndex = 0;
            foreach (var protein in groupedByProtein)
            {
                ChimeraLegendItems.Add(protein.Key, new List<ChimeraLegendItemViewModel>());
                if (protein.Count() > 1)
                {
                    ChimeraLegendItems[protein.Key].Add(new("Shared Ions", ChimeraSpectrumMatchPlot.ColorByProteinDictionary[proteinIndex][0]));
                }
                for (int i = 0; i < protein.Count(); i++)
                {
                    PeptideWithSetModifications peptideWithSetMods =
                        new(protein.ToList()[i].FullSequence, GlobalVariables.AllModsKnownDictionary);
                    var modsString = String.Join(", ",
                        peptideWithSetMods.AllModsOneIsNterminus.Select(p => p.Key + " - " + p.Value.IdWithMotif));
                    ChimeraLegendItems[protein.Key].Add(new(modsString, ChimeraSpectrumMatchPlot.ColorByProteinDictionary[proteinIndex][i + 1]));
                }

            }
        }
    }

    public class ChimeraLegendModel : ChimeraLegendViewModel
    {
        public static List<PsmFromTsv> chimericPsms { get; set; }
        public ChimeraLegendModel Instance => new ChimeraLegendModel();
        //public ChimeraLegendModel()
        //{

        //}

        static ChimeraLegendModel()
        {

        }
    }

}
