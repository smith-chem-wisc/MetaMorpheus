using EngineLayer;
using EngineLayer.GlycoSearch;
using OxyPlot;
using Proteomics.Fragmentation;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Windows.Media;

namespace GuiFunctions
{
    public static class MetaDrawSettings
    {

        #region Customizable
        // graphic settings
        public static Dictionary<string, bool> SpectrumDescription { get; set; }
        public static bool DisplayIonAnnotations { get; set; } = true;
        public static bool AnnotateMzValues { get; set; } = false;
        public static bool AnnotateCharges { get; set; } = false;
        public static bool AnnotationBold { get; set; } = false;
        public static Dictionary<ProductType, OxyColor> ProductTypeToColor { get; set; }
        public static Dictionary<ProductType, OxyColor> BetaProductTypeToColor { get; set; }

        // filter settings
        public static bool ShowDecoys { get; set; } = false;
        public static bool ShowContaminants { get; set; } = true;
        public static double QValueFilter { get; set; } = 0.01;
        public static LocalizationLevel LocalizationLevelStart { get; set; } = LocalizationLevel.Level1;
        public static LocalizationLevel LocalizationLevelEnd { get; set; } = LocalizationLevel.Level3;

        #endregion

        public static Dictionary<ProductType, double> ProductTypeToYOffset { get; set; }
        public static OxyColor VariantCrossColor { get; set; } = OxyColors.Green;
        public static OxyColor UnannotatedPeakColor { get; set; } = OxyColors.LightGray;
        public static SolidColorBrush ModificationAnnotationColor { get; set; } = Brushes.Orange;
        public static double CanvasPdfExportDpi { get; set; } = 300;
        public static double StrokeThicknessUnannotated { get; set; } = 0.7;
        public static double StrokeThicknessAnnotated { get; set; } = 1.0;
        public static double AnnotatedSequenceTextSpacing { get; set; } = 22;
        public static int AnnotatedFontSize { get; set; } = 12;
        public static int NumberOfAAOnScreen { get; set; }
        public static int FirstAAonScreenIndex { get; set; }
        public static bool DrawMatchedIons { get; set; } = true;

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
            ProductTypeToColor[ProductType.y] = OxyColors.Red;
            ProductTypeToColor[ProductType.zDot] = OxyColors.Orange;
            ProductTypeToColor[ProductType.c] = OxyColors.Gold;
            ProductTypeToColor[ProductType.D] = OxyColors.DodgerBlue;
            ProductTypeToColor[ProductType.M] = OxyColors.Firebrick;

            // colors of each fragment to annotate
            BetaProductTypeToColor = ((ProductType[])Enum.GetValues(typeof(ProductType))).ToDictionary(p => p, p => OxyColors.Aqua);
            BetaProductTypeToColor[ProductType.b] = OxyColors.LightBlue;
            BetaProductTypeToColor[ProductType.y] = OxyColors.OrangeRed;
            BetaProductTypeToColor[ProductType.zDot] = OxyColors.LightGoldenrodYellow;
            BetaProductTypeToColor[ProductType.c] = OxyColors.Orange;
            BetaProductTypeToColor[ProductType.D] = OxyColors.AliceBlue;
            BetaProductTypeToColor[ProductType.M] = OxyColors.LightCoral;

            // offset for annotation on base sequence
            ProductTypeToYOffset = ((ProductType[])Enum.GetValues(typeof(ProductType))).ToDictionary(p => p, p => 0.0);
            ProductTypeToYOffset[ProductType.b] = 40;
            ProductTypeToYOffset[ProductType.y] = -10;
            ProductTypeToYOffset[ProductType.c] = 43.6;
            ProductTypeToYOffset[ProductType.zDot] = -13.6;

            // lines to be written on the spectrum
            SpectrumDescription = new Dictionary<string, bool>()
            {
                {"Precursor Charge: ", true },
                {"Precursor Mass: ", true },
                {"Theoretical Mass: ", true },
                {"Score: ", true },
                {"Protein Accession: ", true },
                {"Protein: ", true },
                {"Decoy/Contaminant/Target: ", true },
                {"Q-Value: ", true },
                {"Sequence Length: ", true },
                {"ProForma Level: ", true },
                {"PEP: ", true },
                {"PEP Q-Value: ", true }
            };
        }

        /// <summary>
        /// Create an instance of the metadraw settings to be saved
        /// </summary>
        /// <returns></returns>
        public static MetaDrawSettingsSnapshot MakeSnapShot()
        {
            return new MetaDrawSettingsSnapshot()
            {
                SpectrumDescription = SpectrumDescription,
                DisplayIonAnnotations = DisplayIonAnnotations,
                AnnotateMzValues = AnnotateMzValues,
                AnnotateCharges = AnnotateCharges,
                AnnotationBold = AnnotationBold,
                ShowDecoys = ShowDecoys,    
                ShowContaminants = ShowContaminants,
                QValueFilter = QValueFilter,    
                LocalizationLevelStart = LocalizationLevelStart,
                LocalizationLevelEnd = LocalizationLevelEnd
            };
        }

        /// <summary>
        /// Loads in settings based upon SettingsSnapshot parameter
        /// </summary>
        public static void LoadSettings(MetaDrawSettingsSnapshot settings)
        {
            SpectrumDescription = settings.SpectrumDescription;
            DisplayIonAnnotations = settings.DisplayIonAnnotations;
            AnnotateMzValues = settings.AnnotateMzValues;
            AnnotateCharges = settings.AnnotateCharges;
            AnnotationBold = settings.AnnotationBold;
            ShowDecoys = settings.ShowDecoys;    
            ShowContaminants = settings.ShowContaminants;
            QValueFilter = settings.QValueFilter;
            LocalizationLevelStart = settings.LocalizationLevelStart;
            LocalizationLevelEnd = settings.LocalizationLevelEnd;
        }
    }

}