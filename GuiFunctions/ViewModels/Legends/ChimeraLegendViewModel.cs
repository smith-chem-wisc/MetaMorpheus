﻿using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using EngineLayer;
using GuiFunctions.ViewModels.Legends;
using OxyPlot;
using Proteomics;
using Proteomics.ProteolyticDigestion;

namespace GuiFunctions
{
    public class ChimeraLegendViewModel : LegendViewModel
    {
        #region Private Properties

        private Dictionary<string, List<ChimeraLegendItemViewModel>> chimeraLegendItems;

        #endregion

        #region Public Properties

        public Dictionary<string, List<ChimeraLegendItemViewModel>> ChimeraLegendItems
        {
            get => chimeraLegendItems;
            set { chimeraLegendItems = value; OnPropertyChanged(nameof(ChimeraLegendItems)); }
        }

        public bool DisplaySharedIonLabel
        {
            get
            {
                if (ChimeraLegendItems != null && chimeraLegendItems.Count > 1)
                    return true;
                else
                    return false;
            }
        }

        #endregion

        public ChimeraLegendViewModel(List<PsmFromTsv> chimericIDs, double offset = 0) : base()
        {
            TopOffset = offset;
            ChimeraLegendItems = new();
            ParseLegendItemsFromPsms(chimericIDs);
        }

        /// <summary>
        /// Populates legend items dictionary from a list of PsmTsv
        /// </summary>
        /// <param name="chimericIDs"></param>
        public void ParseLegendItemsFromPsms(List<PsmFromTsv> chimericIDs)
        {
            var groupedByProtein = chimericIDs.GroupBy(p => p.BaseSeq).OrderByDescending(p => p.Count());

            Queue<OxyColor> overflowColors = ChimeraSpectrumMatchPlot.OverflowColors;
            int proteinIndex = 0;
            foreach (var protein in groupedByProtein)
            {
                // more proteins than protein programmed colors
                if (proteinIndex >= ChimeraSpectrumMatchPlot.ColorByProteinDictionary.Keys.Count)
                {
                    proteinIndex = 0;
                }

                ChimeraLegendItems.Add(protein.Key, new List<ChimeraLegendItemViewModel>());
                if (protein.Count() > 1)
                {
                    ChimeraLegendItems[protein.Key].Add(new("Shared Ions", ChimeraSpectrumMatchPlot.ColorByProteinDictionary[proteinIndex][0]));
                }
                for (int i = 0; i < protein.Count(); i++)
                {
                    OxyColor color;
                    // more proteoforms than programmed colors
                    if (i + 1 >= ChimeraSpectrumMatchPlot.ColorByProteinDictionary[proteinIndex].Count)
                    {
                        color = overflowColors.Dequeue();
                    }
                    else
                    {
                        color = ChimeraSpectrumMatchPlot.ColorByProteinDictionary[proteinIndex][i + 1];
                    }

                    PeptideWithSetModifications peptideWithSetMods =
                        new(protein.ToList()[i].FullSequence.Split("|")[0], GlobalVariables.AllModsKnownDictionary);
                    var modsString = String.Join(", ",
                        peptideWithSetMods.AllModsOneIsNterminus.Select(p => p.Key + " - " + p.Value.IdWithMotif));
                    ChimeraLegendItems[protein.Key].Add(new(modsString, color));
                }
                proteinIndex++;
            }
        }
    } 
}
