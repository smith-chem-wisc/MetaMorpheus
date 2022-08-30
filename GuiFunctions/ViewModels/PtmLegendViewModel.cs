using Proteomics;
using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;

namespace GuiFunctions
{
    /// <summary>
    /// View Model class for the ptm color legend
    /// </summary>
    public class PtmLegendViewModel : BaseViewModel
    {
        private Visibility visibility;
        private double topOffset;
        public string Header { get; set; } = "Legend";
        public int HeaderSize { get; set; } = 12;
        public ObservableCollection<PtmLegendItemViewModel> LegendItems { get; set; }

        public Visibility Visibility
        {
            get { return visibility; }
            set
            {
                visibility = value;
                OnPropertyChanged(nameof(Visibility));
            }
        }

        public double TopOffset
        {
            get => topOffset;
            set { topOffset = value; OnPropertyChanged(nameof(TopOffset)); }
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

        public PtmLegendViewModel(List<Modification> mods, double offset = 0)
        {
            LegendItems = new ObservableCollection<PtmLegendItemViewModel>();
            foreach (var mod in mods.Distinct())
            {
                var modItem = new PtmLegendItemViewModel(mod.IdWithMotif);
                LegendItems.Add(modItem);
            }

            topOffset = offset;
            Visibility = mods.Count > 0 ? Visibility.Visible : Visibility.Hidden;
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
            int segmentsPerRow = (int)Math.Floor(maxDisplayedPerRow / (double)(MetaDrawSettings.SequenceAnnotaitonResiduesPerSegment + 1));
            SegmentsPerRow = segmentsPerRow > 0 ? segmentsPerRow : 1;
        }

        /// <summary>
        /// Command to decrease the residues per segment of the sequence annotaiton by one 
        /// </summary>
        public void DecreaseResiduesPerSegment()
        {
            ResiduesPerSegment -= 1;
            double maxDisplayedPerRow = MetaDrawSettings.NumberOfAAOnScreen + 6;
            int segmentsPerRow = (int)Math.Floor(maxDisplayedPerRow / (double)(MetaDrawSettings.SequenceAnnotaitonResiduesPerSegment + 1));
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
    }
}
