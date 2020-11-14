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
        public static Dictionary<ProductType, double> productTypeToYOffset { get; set; }
        public static Dictionary<ProductType, OxyColor> productTypeToColor { get; set; }
        public static Dictionary<ProductType, OxyColor> betaProductTypeToColor { get; set; }
        public static OxyColor variantCrossColor { get; set; } = OxyColors.Green;
        public static OxyColor UnannotatedPeakColor { get; set; } = OxyColors.LightGray;
        public static SolidColorBrush ModificationAnnotationColor { get; set; } = Brushes.Orange;
        public static bool ShowMzValues { get; set; } = false;
        public static bool ShowAnnotationCharges { get; set; } = false;
        public static int AnnotatedFontSize { get; set; } = 12;
        public static bool BoldText { get; set; } = false;

        public static double StrokeThicknessUnannotated = 0.7;
        public static double StrokeThicknessAnnotated = 1.0;
        public static double AnnotatedSequenceTextSpacing = 22;

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
            productTypeToColor = ((ProductType[])Enum.GetValues(typeof(ProductType))).ToDictionary(p => p, p => OxyColors.Aqua);
            productTypeToColor[ProductType.b] = OxyColors.Blue;
            productTypeToColor[ProductType.y] = OxyColors.Purple;
            productTypeToColor[ProductType.zDot] = OxyColors.Orange;
            productTypeToColor[ProductType.c] = OxyColors.Gold;
            productTypeToColor[ProductType.D] = OxyColors.DodgerBlue;
            productTypeToColor[ProductType.M] = OxyColors.Firebrick;

            // colors of each fragment to annotate
            betaProductTypeToColor = ((ProductType[])Enum.GetValues(typeof(ProductType))).ToDictionary(p => p, p => OxyColors.Aqua);
            betaProductTypeToColor[ProductType.b] = OxyColors.LightBlue;
            betaProductTypeToColor[ProductType.y] = OxyColors.MediumPurple;
            betaProductTypeToColor[ProductType.zDot] = OxyColors.LightGoldenrodYellow;
            betaProductTypeToColor[ProductType.c] = OxyColors.OrangeRed;
            betaProductTypeToColor[ProductType.D] = OxyColors.AliceBlue;
            betaProductTypeToColor[ProductType.M] = OxyColors.LightCoral;

            // offset for annotation on base sequence
            productTypeToYOffset = ((ProductType[])Enum.GetValues(typeof(ProductType))).ToDictionary(p => p, p => 0.0);
            productTypeToYOffset[ProductType.b] = 40;
            productTypeToYOffset[ProductType.y] = -10;
            productTypeToYOffset[ProductType.c] = 43.6;
            productTypeToYOffset[ProductType.zDot] = -13.6;
        }
    }
}