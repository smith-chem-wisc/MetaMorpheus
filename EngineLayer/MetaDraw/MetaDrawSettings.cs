using EngineLayer;
using EngineLayer.GlycoSearch;
using OxyPlot;
using Proteomics.Fragmentation;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Windows.Media;

namespace EngineLayer
{
    public static class MetaDrawSettings
    {
        // graphics settings
        public static Dictionary<ProductType, double> ProductTypeToYOffset { get; set; }
        public static Dictionary<ProductType, OxyColor> ProductTypeToColor { get; set; }
        public static Dictionary<ProductType, OxyColor> BetaProductTypeToColor { get; set; }
        public static OxyColor VariantCrossColor { get; set; } = OxyColors.Green;
        public static OxyColor UnannotatedPeakColor { get; set; } = OxyColors.LightGray;
        public static SolidColorBrush ModificationAnnotationColor { get; set; } = Brushes.Orange;
        public static double CanvasPdfExportDpi { get; set; } = 300;
        public static bool DisplayIonAnnotations { get; set; } = true;
        public static bool AnnotateMzValues { get; set; } = false;
        public static bool AnnotateCharges { get; set; } = false;
        public static int AnnotatedFontSize { get; set; } = 12;
        public static bool AnnotationBold { get; set; } = false;
        public static double StrokeThicknessUnannotated { get; set; } = 0.7;
        public static double StrokeThicknessAnnotated { get; set; } = 1.0;
        public static double AnnotatedSequenceTextSpacing { get; set; } = 22;

        // filter settings
        public static bool ShowDecoys { get; set; } = false;
        public static bool ShowContaminants { get; set; } = true;
        public static double QValueFilter { get; set; } = 0.01;
        public static LocalizationLevel LocalizationLevelStart { get; set; } = LocalizationLevel.Level1;
        public static LocalizationLevel LocalizationLevelEnd { get; set; } = LocalizationLevel.Level3;

        static MetaDrawSettings()
        {
            InitializeDictionaries();
        }

        public static bool FilterAcceptsPsm(PsmFromTsv psm)
        {
            if (psm.QValue <= QValueFilter
                 && (psm.QValueNotch == null || psm.QValueNotch <= QValueFilter)
                 && (psm.DecoyContamTarget == "T" || (psm.DecoyContamTarget == "D" && ShowDecoys) || (psm.DecoyContamTarget == "C" && ShowContaminants))
                 && (psm.GlycanLocalizationLevel == null || psm.GlycanLocalizationLevel >= LocalizationLevelStart && psm.GlycanLocalizationLevel <= LocalizationLevelEnd))
            {
                return true;
            }

            return false;
        }

        private static void InitializeDictionaries()
        {
            // colors of each fragment to annotate
            ProductTypeToColor = ((ProductType[])Enum.GetValues(typeof(ProductType))).ToDictionary(p => p, p => OxyColors.Aqua);
            ProductTypeToColor[ProductType.b] = OxyColors.Blue;
            ProductTypeToColor[ProductType.y] = OxyColors.Purple;
            ProductTypeToColor[ProductType.zDot] = OxyColors.Orange;
            ProductTypeToColor[ProductType.c] = OxyColors.Gold;
            ProductTypeToColor[ProductType.D] = OxyColors.DodgerBlue;
            ProductTypeToColor[ProductType.M] = OxyColors.Firebrick;

            // colors of each fragment to annotate
            BetaProductTypeToColor = ((ProductType[])Enum.GetValues(typeof(ProductType))).ToDictionary(p => p, p => OxyColors.Aqua);
            BetaProductTypeToColor[ProductType.b] = OxyColors.LightBlue;
            BetaProductTypeToColor[ProductType.y] = OxyColors.MediumPurple;
            BetaProductTypeToColor[ProductType.zDot] = OxyColors.LightGoldenrodYellow;
            BetaProductTypeToColor[ProductType.c] = OxyColors.OrangeRed;
            BetaProductTypeToColor[ProductType.D] = OxyColors.AliceBlue;
            BetaProductTypeToColor[ProductType.M] = OxyColors.LightCoral;

            // offset for annotation on base sequence
            ProductTypeToYOffset = ((ProductType[])Enum.GetValues(typeof(ProductType))).ToDictionary(p => p, p => 0.0);
            ProductTypeToYOffset[ProductType.b] = 40;
            ProductTypeToYOffset[ProductType.y] = -10;
            ProductTypeToYOffset[ProductType.c] = 43.6;
            ProductTypeToYOffset[ProductType.zDot] = -13.6;
        }
    }
}