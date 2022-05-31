using EngineLayer.GlycoSearch;
using OxyPlot;
using Proteomics.Fragmentation;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace GuiFunctions
{
    /// <summary>
    /// Class for exporting and saving an instance of the static class MetaDrawSettings 
    /// </summary>
    public class MetaDrawSettingsSnapshot 
    {
        public Dictionary<string, bool> SpectrumDescription { get; set; }
        public bool DisplayIonAnnotations { get; set; } = true;
        public bool AnnotateMzValues { get; set; } = false;
        public bool AnnotateCharges { get; set; } = false;
        public bool AnnotationBold { get; set; } = false;
        public static Dictionary<ProductType, OxyColor> ProductTypeToColor { get; set; }
        public static Dictionary<ProductType, OxyColor> BetaProductTypeToColor { get; set; }

        // filter settings
        public bool ShowDecoys { get; set; } = false;
        public bool ShowContaminants { get; set; } = true;
        public double QValueFilter { get; set; } = 0.01;
        public LocalizationLevel LocalizationLevelStart { get; set; } = LocalizationLevel.Level1;
        public LocalizationLevel LocalizationLevelEnd { get; set; } = LocalizationLevel.Level3;

        /// <summary>
        /// Creates an instance of the snapshot based upon default settings for design time display
        /// </summary>
        public static MetaDrawSettingsSnapshot Instance => new MetaDrawSettingsSnapshot();

        public MetaDrawSettingsSnapshot()
        {
            SpectrumDescription = MetaDrawSettings.SpectrumDescription;
            DisplayIonAnnotations = MetaDrawSettings.DisplayIonAnnotations;
            AnnotateMzValues = MetaDrawSettings.AnnotateMzValues;
            AnnotateCharges = MetaDrawSettings.AnnotateCharges;
            AnnotationBold = MetaDrawSettings.AnnotationBold;
            ShowDecoys = MetaDrawSettings.ShowDecoys;
            ShowContaminants = MetaDrawSettings.ShowContaminants;
            QValueFilter = MetaDrawSettings.QValueFilter;
            LocalizationLevelStart = MetaDrawSettings.LocalizationLevelStart;
            LocalizationLevelEnd = MetaDrawSettings.LocalizationLevelEnd;
            ProductTypeToColor = MetaDrawSettings.ProductTypeToColor;
            BetaProductTypeToColor = MetaDrawSettings.BetaProductTypeToColor;
        }
    }
}
