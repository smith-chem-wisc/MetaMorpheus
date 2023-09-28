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
        // graphic settings
        public bool DisplayIonAnnotations { get; set; } = true;
        public bool AnnotateMzValues { get; set; } = false;
        public bool AnnotateCharges { get; set; } = true;
        public bool AnnotationBold { get; set; } = false;
        public bool DisplayInternalIons { get; set; } = true;
        public bool DisplayInternalIonAnnotations { get; set; } = true;
        public bool SubAndSuperScriptIons { get; set; } = true;
        public bool DrawStationarySequence { get; set; } = true;
        public bool DrawNumbersUnderStationary { get; set; } = true;
        public bool ShowLegend { get; set; } = true;
        public List<bool> SpectrumDescriptionValues { get; set; }
        public List<string> ProductTypeToColorValues { get; set; }
        public List<string> BetaProductTypeToColorValues { get; set; }
        public List<string> ModificationTypeToColorValues { get; set; }
        public List<string> CoverageTypeToColorValues { get; set; }
        public string UnannotatedPeakColor { get; set; }
        public string InternalIonColor { get; set; }

        // filter settings
        public bool ShowDecoys { get; set; } = false;
        public bool ShowContaminants { get; set; } = true;
        public double QValueFilter { get; set; } = 0.01;
        public string AmbiguityFilter { get; set; } = "No Filter";
        public LocalizationLevel LocalizationLevelStart { get; set; } = LocalizationLevel.Level1;
        public LocalizationLevel LocalizationLevelEnd { get; set; } = LocalizationLevel.Level3;
        public string ExportType { get; set; }
    }
}
