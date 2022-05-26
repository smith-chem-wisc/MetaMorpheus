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
        public static Dictionary<OxyColor, string> PossibleColors { get; set; }
        public static Dictionary<ProductType, OxyColor> ProductTypeToColor { get; set; }
        public static Dictionary<ProductType, OxyColor> BetaProductTypeToColor { get; set; }

        // filter settings
        public static bool ShowDecoys { get; set; } = false;
        public static bool ShowContaminants { get; set; } = true;
        public static double QValueFilter { get; set; } = 0.01;
        public static LocalizationLevel LocalizationLevelStart { get; set; } = LocalizationLevel.Level1;
        public static LocalizationLevel LocalizationLevelEnd { get; set; } = LocalizationLevel.Level3;

        #endregion

        public static string[] SpectrumDescriptors { get; set; } = 
        {"Precursor Charge: ", "Precursor Mass: ", "Theoretical Mass: ", "Protein Accession: ", "Protein: ",
        "Decoy/Contaminant/Target: ", "Sequence Length: ", "ProForma Level: ", "Score: ", "Q-Value: ", "PEP: ", "PEP Q-Value: "};
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
            // default color of each fragment to annotate
            ProductTypeToColor = ((ProductType[])Enum.GetValues(typeof(ProductType))).ToDictionary(p => p, p => OxyColors.Aqua);
            ProductTypeToColor[ProductType.b] = OxyColors.Blue;
            ProductTypeToColor[ProductType.y] = OxyColors.Red;
            ProductTypeToColor[ProductType.zDot] = OxyColors.Orange;
            ProductTypeToColor[ProductType.c] = OxyColors.Gold;
            ProductTypeToColor[ProductType.D] = OxyColors.DodgerBlue;
            ProductTypeToColor[ProductType.M] = OxyColors.Firebrick;

            // default color of each fragment to annotate
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
            SpectrumDescription = SpectrumDescriptors.ToDictionary(p => p, p => true);
            PossibleColors = ((ColorEnum[])Enum.GetValues(typeof(ColorEnum))).ToDictionary(p => NameToOxyColorConverter(p.ToString()), p => p.ToString());
        }


        /// <summary>
        /// Create an instance of the metadraw settings to be saved
        /// </summary>
        /// <returns></returns>
        public static MetaDrawSettingsSnapshot MakeSnapShot()
        {
            return new MetaDrawSettingsSnapshot()
            {
                DisplayIonAnnotations = DisplayIonAnnotations,
                AnnotateMzValues = AnnotateMzValues,
                AnnotateCharges = AnnotateCharges,
                AnnotationBold = AnnotationBold,
                ShowDecoys = ShowDecoys,
                ShowContaminants = ShowContaminants,
                QValueFilter = QValueFilter,
                LocalizationLevelStart = LocalizationLevelStart,
                LocalizationLevelEnd = LocalizationLevelEnd,
                ProductTypeToColorValues = ProductTypeToColor.Values.Select(p => OxyColorToNameConverter(p)).ToList(),
                BetaProductTypeToColorValues = BetaProductTypeToColor.Values.Select(p => OxyColorToNameConverter(p)).ToList(),
                SpectrumDescriptionValues = SpectrumDescription.Values.ToList()
            };
        }

        /// <summary>
        /// Loads in settings based upon SettingsSnapshot parameter
        /// </summary>
        public static void LoadSettings(MetaDrawSettingsSnapshot settings)
        {
            DisplayIonAnnotations = settings.DisplayIonAnnotations;
            AnnotateMzValues = settings.AnnotateMzValues;
            AnnotateCharges = settings.AnnotateCharges;
            AnnotationBold = settings.AnnotationBold;
            ShowDecoys = settings.ShowDecoys;
            ShowContaminants = settings.ShowContaminants;
            QValueFilter = settings.QValueFilter;
            LocalizationLevelStart = settings.LocalizationLevelStart;
            LocalizationLevelEnd = settings.LocalizationLevelEnd;

            ProductTypeToColor = ((ProductType[])Enum.GetValues(typeof(ProductType))).ToDictionary(p => p, p => NameToOxyColorConverter(settings.ProductTypeToColorValues[Array.IndexOf(((ProductType[])Enum.GetValues(typeof(ProductType))), p)]));
            BetaProductTypeToColor = ((ProductType[])Enum.GetValues(typeof(ProductType))).ToDictionary(p => p, p => NameToOxyColorConverter(settings.BetaProductTypeToColorValues[Array.IndexOf(((ProductType[])Enum.GetValues(typeof(ProductType))), p)]));
            SpectrumDescription = SpectrumDescriptors.ToDictionary(p => p, p => settings.SpectrumDescriptionValues[Array.IndexOf(SpectrumDescriptors, p)]);

        }

        /// <summary>
        /// Converts the string representation of the color to a OxyColor
        /// </summary>
        /// <param name="name"></param>
        /// <returns></returns>
        public static OxyColor NameToOxyColorConverter(string name)
        {
            return name switch
            {
                "Blue" => OxyColors.Blue,
                "Red" => OxyColors.Red,
                "Orange" => OxyColors.Orange,
                "Violet" => OxyColors.Violet,
                "Gold" => OxyColors.Gold,
                "Black" => OxyColors.Black,
                "Green" => OxyColors.Green,
                "HotPink" => OxyColors.HotPink,
                "Indigo" => OxyColors.Indigo,
                "Lime" => OxyColors.Lime,
                "Magenta" => OxyColors.Magenta,
                "MidnightBlue" => OxyColors.MidnightBlue,
                "Olive" => OxyColors.Olive,
                "Purple" => OxyColors.Purple,
                "DodgerBlue" => OxyColors.DodgerBlue,
                "Firebrick" => OxyColors.Firebrick,
                "LightBlue" => OxyColors.LightBlue,
                "OrangeRed" => OxyColors.OrangeRed,
                "LightGoldenrodYellow" => OxyColors.LightGoldenrodYellow,
                "AliceBlue" => OxyColors.AliceBlue,
                "LightCoral" => OxyColors.LightCoral,
                "Aqua" => OxyColors.Aqua,
                "Chartreuse" => OxyColors.Chartreuse,
                "BurlyWood" => OxyColors.BurlyWood,
                _ => OxyColors.Sienna,
            };
        }

        /// <summary>
        /// Converts OxyColor to a string representation
        /// </summary>
        /// <param name="color"></param>
        /// <returns></returns>
        public static string OxyColorToNameConverter(OxyColor color)
        {
            string name = color.ToString();

            if (color == OxyColors.Blue)
                return "Blue";
            else if (color == OxyColors.Red)
                return "Red";
            else if (color == OxyColors.Orange)
                return "Orange";
            else if (color == OxyColors.Violet)
                return "Violet";
            else if (color == OxyColors.Gold)
                return "Gold";
            else if (color == OxyColors.Black)
                return "Black";
            else if (color == OxyColors.Green)
                return "Green";
            else if (color == OxyColors.HotPink)
                return "HotPink";
            else if (color == OxyColors.Indigo)
                return "Indigo";
            else if (color == OxyColors.Lime)
                return "Lime";
            else if (color == OxyColors.Magenta)
                return "Magenta";
            else if (color == OxyColors.MidnightBlue)
                return "MidnightBlue";
            else if (color == OxyColors.Olive)
                return "Olive";
            else if (color == OxyColors.Purple)
                return "Purple";
            else if (color == OxyColors.DodgerBlue)
                return "DodgerBlue";
            else if (color == OxyColors.Firebrick)
                return "Firebrick";
            else if (color == OxyColors.LightBlue)
                return "LightBlue";
            else if (color == OxyColors.OrangeRed)
                return "OrangeRed";
            else if (color == OxyColors.LightGoldenrodYellow)
                return "LightGoldenrodYellow";
            else if (color == OxyColors.AliceBlue)
                return "AliceBlue";
            else if (color == OxyColors.LightCoral)
                return "LightCoral";
            else if (color == OxyColors.Aqua)
                return "Aqua";
            else if (color == OxyColors.Chartreuse)
                return "Chartreuse";
            else if (color == OxyColors.BurlyWood)
                return "BurlyWood";
            else if (color == OxyColors.Sienna)
                return "Sienna";
            else
                return "Blue";
        }

        /// <summary>
        /// Enum full of all selectable colors
        /// </summary>
        public enum ColorEnum
        {
            Blue,
            Red,
            Orange,
            Violet,
            Gold,
            Black,
            Green,
            HotPink,
            Indigo,
            Lime,
            Magenta,
            MidnightBlue,
            Olive,
            Purple,
            DodgerBlue,
            Firebrick,
            LightBlue,
            OrangeRed,
            LightGoldenrodYellow,
            AliceBlue,
            LightCoral,
            Aqua,
            Chartreuse,
            BurlyWood
        }
    }
}