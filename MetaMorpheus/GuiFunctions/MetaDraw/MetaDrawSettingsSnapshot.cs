using EngineLayer.GlycoSearch;
using OxyPlot;
using Omics.Fragmentation;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Readers;
using GuiFunctions.MetaDraw;

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
        public List<string> SpectrumDescriptionValues { get; set; }
        public List<string> ProductTypeToColorValues { get; set; }
        public List<string> BetaProductTypeToColorValues { get; set; }
        public List<string> ModificationTypeToColorValues { get; set; }
        public List<string> CoverageTypeToColorValues { get; set; }
        public string UnannotatedPeakColor { get; set; }
        public string InternalIonColor { get; set; }
        public int AnnotatedFontSize { get; set; } = 14;
        public int AxisTitleTextSize { get; set; } = 14;
        public int AxisLabelTextSize { get; set; } = 12;
        public double StrokeThicknessUnannotated { get; set; } = 0.7;
        public double StrokeThicknessAnnotated { get; set; } = 1.0;
        public double SpectrumDescriptionFontSize { get; set; } = 10;
        public bool DisplayChimeraLegend { get; set; } = true;
        public bool SuppressMessageBoxes { get; set; } = true;
        public bool ChimeraLegendTakeFirstIfAmbiguous { get; set; }
        public LegendDisplayProperty ChimeraLegendMainTextType { get; set; } = LegendDisplayProperty.ProteinName;
        public LegendDisplayProperty ChimeraLegendSubTextType { get; set; } = LegendDisplayProperty.Modifications;

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
