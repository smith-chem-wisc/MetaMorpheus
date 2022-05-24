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
        public bool DisplayIonAnnotations { get; set; } = true;
        public bool AnnotateMzValues { get; set; } = false;
        public bool AnnotateCharges { get; set; } = false;
        public bool AnnotationBold { get; set; } = false;
        public List<bool> SpectrumDescriptionValues { get; set; }
        public List<string> ProductTypeToColorValues { get; set; }
        public List<string> BetaProductTypeToColorValues { get; set; }

        // filter settings
        public bool ShowDecoys { get; set; } = false;
        public bool ShowContaminants { get; set; } = true;
        public double QValueFilter { get; set; } = 0.01;
        public LocalizationLevel LocalizationLevelStart { get; set; } = LocalizationLevel.Level1;
        public LocalizationLevel LocalizationLevelEnd { get; set; } = LocalizationLevel.Level3;
    }
}
