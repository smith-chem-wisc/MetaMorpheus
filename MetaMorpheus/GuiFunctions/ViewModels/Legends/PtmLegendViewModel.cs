using GuiFunctions.ViewModels.Legends;
using Proteomics;
using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using EngineLayer;
using Omics.Modifications;
using Proteomics.ProteolyticDigestion;
using Readers;
using Omics;

namespace GuiFunctions
{
    /// <summary>
    /// View Model class for the ptm color legend
    /// </summary>
    public class PtmLegendViewModel : LegendViewModel
    {
        public double SpectrumDescriptionFontSize => GuiFunctions.MetaDrawSettings.SpectrumDescriptionFontSize;

        public double LegendVerticalOffset
        {
            get
            {
                // These are magic numbers that help make the legend always below the description
                var offset =  (MetaDrawSettings.SpectrumDescription.Count(p => p.Value) + 2)
                    * MetaDrawSettings.SpectrumDescriptionFontSize
                    * 1.4;
                return offset;
            }
        }

        /// <summary>
        /// Segments per row in the sequence annotation 
        /// </summary>
        public int SegmentsPerRow
        {
            get { return MetaDrawSettings.SequenceAnnotationSegmentPerRow; }
            set
            {
                if (value <= 0) throw new IndexOutOfRangeException("SegmentsPerRow cannot be less than one");
                MetaDrawSettings.SequenceAnnotationSegmentPerRow = value;
                OnPropertyChanged(nameof(SegmentsPerRow));
            }
        }

        /// <summary>
        /// Residues per segment in the sequence annotation
        /// </summary>
        public int ResiduesPerSegment
        {
            get { return MetaDrawSettings.SequenceAnnotaitonResiduesPerSegment; }
            set
            {
                if (value <= 0) throw new IndexOutOfRangeException("ResiduesPerSegment cannot be less than one");
                MetaDrawSettings.SequenceAnnotaitonResiduesPerSegment = value;
                OnPropertyChanged(nameof(ResiduesPerSegment));
            }
        }

        #region Constructor

        public PtmLegendViewModel(SpectrumMatchFromTsv psm, double offset = 0) : base()
        {
            ParseModsFromPsmTsv(psm);
            TopOffset = offset;
            Visibility = LegendItemViewModels.Count > 0 ? true : false;
        }

        #endregion

        #region Commands

        /// <summary>
        /// Command to increase the residues per segent of the sequence annotation view by one 
        /// </summary>
        public void IncreaseResiduesPerSegment()
        {
            ResiduesPerSegment += 1;
            double maxDisplayedPerRow = MetaDrawSettings.NumberOfAAOnScreen + 7;
            int segmentsPerRow = (int)Math.Floor(maxDisplayedPerRow / (MetaDrawSettings.SequenceAnnotaitonResiduesPerSegment + 1));
            SegmentsPerRow = segmentsPerRow > 0 ? segmentsPerRow : 1;
        }

        /// <summary>
        /// Command to decrease the residues per segment of the sequence annotaiton by one 
        /// </summary>
        public void DecreaseResiduesPerSegment()
        {
            ResiduesPerSegment -= 1;
            double maxDisplayedPerRow = MetaDrawSettings.NumberOfAAOnScreen + 6;
            int segmentsPerRow = (int)Math.Floor(maxDisplayedPerRow / (MetaDrawSettings.SequenceAnnotaitonResiduesPerSegment + 1));
            SegmentsPerRow = segmentsPerRow > 0 ? segmentsPerRow : 1;
        }

        /// <summary>
        /// Command to increase the segmetns per row of the sequence annotaiton by one 
        /// </summary>
        public void IncreaseSegmentsPerRow()
        {
            SegmentsPerRow += 1;
        }

        /// <summary>
        /// Command to decrease the segments per row of the sequence annotiation by one 
        /// </summary>
        public void DecreaseSegmentsPerRow()
        {
            SegmentsPerRow -= 1;
        }

        #endregion

        private void ParseModsFromPsmTsv(SpectrumMatchFromTsv psm)
        {
            IBioPolymerWithSetMods bioPolymerWithSetMods = psm.ToBioPolymerWithSetMods();
            List<Modification> mods = bioPolymerWithSetMods.AllModsOneIsNterminus.Values.ToList();
            foreach (var mod in mods.Distinct())
            {
                var modItem = new PtmLegendItemViewModel(mod.IdWithMotif);
                LegendItemViewModels.Add(modItem);
            }
        }
    }
}
