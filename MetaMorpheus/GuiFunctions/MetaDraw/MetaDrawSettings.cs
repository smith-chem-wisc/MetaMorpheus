using EngineLayer;
using EngineLayer.GlycoSearch;
using OxyPlot;
using Proteomics;
using Omics.Fragmentation;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Reflection;
using System.Text;
using System.Windows.Media;

namespace GuiFunctions
{
    public static class MetaDrawSettings
    {
        #region Constants 

        public static char[] SubScriptNumbers = {
            '\u2080', '\u2081', '\u2082', '\u2083', '\u2084',
            '\u2085', '\u2086', '\u2087', '\u2088', '\u2089'
        };

        public static char[] SuperScriptNumbers = {
            '\u2070', '\u00b9', '\u00b2', '\u00b3', '\u2074',
            '\u2075', '\u2076', '\u2077', '\u2078', '\u2079'
        }; 
        
        #endregion

        #region Customizable Settings

        // graphic settings
        public static Dictionary<string, bool> SpectrumDescription { get; set; }
        public static bool DisplayIonAnnotations { get; set; } = true;
        public static bool AnnotateMzValues { get; set; } = false;
        public static bool AnnotateCharges { get; set; } = true;
        public static bool AnnotationBold { get; set; } = false;
        public static bool DisplayInternalIons { get; set; } = true;
        public static bool DisplayInternalIonAnnotations { get; set; }= true;
        public static Dictionary<OxyColor, string> PossibleColors { get; set; }
        public static Dictionary<ProductType, OxyColor> ProductTypeToColor { get; set; }
        public static Dictionary<ProductType, OxyColor> BetaProductTypeToColor { get; set; }
        public static Dictionary<string, OxyColor> ModificationTypeToColor {get; set; }
        public static Dictionary<string, OxyColor> CoverageTypeToColor { get; set; }
        public static bool DrawStationarySequence { get; set; } = true;
        public static bool DrawNumbersUnderStationary { get; set; } = true;
        public static bool ShowLegend { get; set; } = true;

        public static bool SubAndSuperScriptIons = true;

        // filter settings
        public static bool ShowDecoys { get; set; } = false;
        public static bool ShowContaminants { get; set; } = true;
        public static double QValueFilter { get; set; } = 0.01;
        public static string AmbiguityFilter { get; set; } = "No Filter";
        public static LocalizationLevel LocalizationLevelStart { get; set; } = LocalizationLevel.Level1;
        public static LocalizationLevel LocalizationLevelEnd { get; set; } = LocalizationLevel.Level3;
        public static string ExportType { get; set; } = "Pdf"; 

        #endregion

        // used for constructing data structures and mathcing them with the saved settings
        #region Data Lists
        public static List<OxyColor> AllColors { get; set; } = new List<OxyColor>()
        {   OxyColors.Undefined, OxyColors.Automatic, OxyColors.AliceBlue, OxyColors.AntiqueWhite, OxyColors.Aqua, OxyColors.Aquamarine,
            OxyColors.Azure, OxyColors.Beige, OxyColors.Bisque, OxyColors.Black, OxyColors.BlanchedAlmond, OxyColors.Blue, OxyColors.BlueViolet,
            OxyColors.Brown, OxyColors.BurlyWood, OxyColors.CadetBlue, OxyColors.Chartreuse, OxyColors.Chocolate, OxyColors.Coral, OxyColors.CornflowerBlue,
            OxyColors.Cornsilk, OxyColors.Crimson, OxyColors.DarkBlue, OxyColors.DarkCyan, OxyColors.DarkGoldenrod, OxyColors.DarkGray, OxyColors.YellowGreen,
            OxyColors.DarkGreen, OxyColors.DarkKhaki, OxyColors.DarkMagenta, OxyColors.DarkOliveGreen, OxyColors.DarkOrange, OxyColors.DarkOrchid,
            OxyColors.DarkRed, OxyColors.DarkSalmon, OxyColors.DarkSeaGreen, OxyColors.DarkSlateBlue, OxyColors.DarkSlateGray, OxyColors.DarkTurquoise,
            OxyColors.DarkViolet, OxyColors.DeepPink, OxyColors.DeepSkyBlue, OxyColors.DimGray, OxyColors.DodgerBlue, OxyColors.Firebrick, OxyColors.FloralWhite,
            OxyColors.ForestGreen, OxyColors.Gainsboro, OxyColors.GhostWhite, OxyColors.Gold, OxyColors.Goldenrod, OxyColors.Gray, OxyColors.Green,
            OxyColors.GreenYellow, OxyColors.Honeydew, OxyColors.HotPink, OxyColors.IndianRed, OxyColors.Indigo, OxyColors.Ivory, OxyColors.Khaki, OxyColors.Lavender,
            OxyColors.LavenderBlush, OxyColors.LawnGreen, OxyColors.LemonChiffon, OxyColors.LightBlue, OxyColors.LightCoral, OxyColors.LightCyan, OxyColors.LightGoldenrodYellow,
            OxyColors.LightGray, OxyColors.LightGreen, OxyColors.LightPink, OxyColors.LightSalmon, OxyColors.LightSeaGreen, OxyColors.LightSkyBlue, OxyColors.LightSlateGray,
            OxyColors.LightSteelBlue, OxyColors.LightYellow, OxyColors.Lime, OxyColors.LimeGreen, OxyColors.Linen, OxyColors.Fuchsia, OxyColors.Maroon, OxyColors.MediumAquamarine,
            OxyColors.MediumBlue, OxyColors.MediumOrchid, OxyColors.MediumPurple, OxyColors.MediumSeaGreen, OxyColors.MediumSlateBlue, OxyColors.MediumSpringGreen, OxyColors.MediumTurquoise,
            OxyColors.MediumVioletRed, OxyColors.MidnightBlue, OxyColors.MintCream, OxyColors.MistyRose, OxyColors.Moccasin, OxyColors.NavajoWhite, OxyColors.Navy, OxyColors.OldLace,
            OxyColors.Olive, OxyColors.OliveDrab, OxyColors.Orange, OxyColors.OrangeRed, OxyColors.Orchid, OxyColors.PaleGoldenrod, OxyColors.PaleGreen, OxyColors.PaleTurquoise,
            OxyColors.PaleVioletRed, OxyColors.PapayaWhip, OxyColors.PeachPuff, OxyColors.Peru, OxyColors.Pink, OxyColors.Plum, OxyColors.PowderBlue, OxyColors.Purple, OxyColors.Red,
            OxyColors.RosyBrown, OxyColors.RoyalBlue, OxyColors.SaddleBrown, OxyColors.Salmon, OxyColors.SandyBrown, OxyColors.SeaGreen, OxyColors.SeaShell, OxyColors.Sienna,
            OxyColors.Silver, OxyColors.SkyBlue, OxyColors.SlateBlue, OxyColors.SlateGray, OxyColors.Snow, OxyColors.SpringGreen, OxyColors.SteelBlue, OxyColors.Tan, OxyColors.Teal,
            OxyColors.Thistle, OxyColors.Tomato, OxyColors.Transparent, OxyColors.Turquoise, OxyColors.Violet, OxyColors.Wheat, OxyColors.White, OxyColors.WhiteSmoke, OxyColors.Yellow
        };
        public static string[] SpectrumDescriptors { get; set; } =
        {"Precursor Charge: ", "Precursor Mass: ", "Theoretical Mass: ", "Protein Accession: ", "Protein: ",
        "Decoy/Contaminant/Target: ", "Sequence Length: ", "Ambiguity Level: ", "Spectral Angle: ", "Score: ", "Q-Value: ", "PEP: ", "PEP Q-Value: "};
        public static string[] CoverageTypes { get; set; } = { "N-Terminal Color", "C-Terminal Color", "Internal Color" };
        public static string[] ExportTypes { get; set; } = { "Pdf", "Png", "Jpeg", "Tiff", "Wmf", "Bmp" };
        public static string[] AmbiguityTypes { get; set; } = { "No Filter", "1", "2A", "2B", "2C", "2D", "3", "4", "5" };

        #endregion

        public static Dictionary<ProductType, double> ProductTypeToYOffset { get; set; }
        public static Dictionary<ProductType, double> ProductTypeToXOffset { get; set; }
        public static OxyColor VariantCrossColor { get; set; } = OxyColors.Green;
        public static OxyColor UnannotatedPeakColor { get; set; } = OxyColors.LightGray;
        public static OxyColor InternalIonColor { get; set; } = OxyColors.Purple;
        public static SolidColorBrush ModificationAnnotationColor { get; set; } = Brushes.Orange;
        public static double CanvasPdfExportDpi { get; set; } = 600;
        public static double StrokeThicknessUnannotated { get; set; } = 0.7;
        public static double StrokeThicknessAnnotated { get; set; } = 1.0;
        public static double AnnotatedSequenceTextSpacing { get; set; } = 22;
        public static int AnnotatedFontSize { get; set; } = 14;
        public static int NumberOfAAOnScreen { get; set; }
        public static int FirstAAonScreenIndex { get; set; }
        public static bool DrawMatchedIons { get; set; } = true;
        public static int SequenceAnnotationSegmentPerRow { get; set; } = 3;
        public static int SequenceAnnotaitonResiduesPerSegment { get; set; } = 10;

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
                // Ambiguity filtering conditionals, should only be hit if Ambiguity Filtering is selected
                if (AmbiguityFilter == "No Filter" || psm.AmbiguityLevel == AmbiguityFilter)
                {
                    return true;
                }
                else
                {
                    return false;
                }
            }

            return false;
        }

        private static void InitializeDictionaries()
        {
            // If no default settings are saved
            string settingsPath = Path.Combine(GlobalVariables.DataDir, "DefaultParameters", @"MetaDrawSettingsDefault.xml");
            if (!File.Exists(settingsPath))
            {
                // default color of each fragment to annotate
                ProductTypeToColor = ((ProductType[])Enum.GetValues(typeof(ProductType))).ToDictionary(p => p, p => OxyColors.Aqua);
                ProductTypeToColor[ProductType.a] = OxyColors.DarkOrange;
                ProductTypeToColor[ProductType.aBaseLoss] = OxyColors.SandyBrown;
                ProductTypeToColor[ProductType.aWaterLoss] = OxyColors.Orange;
                ProductTypeToColor[ProductType.b] = OxyColors.Blue;
                ProductTypeToColor[ProductType.bBaseLoss] = OxyColors.DarkSlateBlue;
                ProductTypeToColor[ProductType.bAmmoniaLoss] = OxyColors.DarkSlateBlue;
                ProductTypeToColor[ProductType.bWaterLoss] = OxyColors.LightBlue;
                ProductTypeToColor[ProductType.c] = OxyColors.Gold;
                ProductTypeToColor[ProductType.cBaseLoss] = OxyColors.Goldenrod;
                ProductTypeToColor[ProductType.cWaterLoss] = OxyColors.Khaki;
                ProductTypeToColor[ProductType.d] = OxyColors.Purple;
                ProductTypeToColor[ProductType.dBaseLoss] = OxyColors.DarkViolet;
                ProductTypeToColor[ProductType.dWaterLoss] = OxyColors.MediumPurple;

                ProductTypeToColor[ProductType.w] = OxyColors.Green;
                ProductTypeToColor[ProductType.wBaseLoss] = OxyColors.DarkGreen;
                ProductTypeToColor[ProductType.wWaterLoss] = OxyColors.LightGreen;
                ProductTypeToColor[ProductType.x] = OxyColors.Peru;
                ProductTypeToColor[ProductType.xBaseLoss] = OxyColors.Sienna;
                ProductTypeToColor[ProductType.xWaterLoss] = OxyColors.BurlyWood;
                ProductTypeToColor[ProductType.y] = OxyColors.Red;
                ProductTypeToColor[ProductType.yBaseLoss] = OxyColors.DarkSalmon;
                ProductTypeToColor[ProductType.yAmmoniaLoss] = OxyColors.DarkSalmon;
                ProductTypeToColor[ProductType.yWaterLoss] = OxyColors.Tomato;
                ProductTypeToColor[ProductType.z] = OxyColors.Magenta;
                ProductTypeToColor[ProductType.zBaseLoss] = OxyColors.DarkMagenta;
                ProductTypeToColor[ProductType.zWaterLoss] = OxyColors.Plum;

                ProductTypeToColor[ProductType.zDot] = OxyColors.Orange;
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

                #region Setting Color Defaults

                    ModificationTypeToColor = GlobalVariables.AllModsKnownDictionary.Values.ToDictionary(p => p.IdWithMotif, p => OxyColors.Orange);
                    
                    // setting whole groups
                    foreach (var mod in GlobalVariables.AllModsKnownDictionary.Values.Where(p => p.ModificationType == "Common Biological").Select(p => p.IdWithMotif))
                    {
                        ModificationTypeToColor[mod] = OxyColors.Plum;
                    }
                    
                    foreach (var mod in GlobalVariables.AllModsKnownDictionary.Values.Where(p => p.ModificationType == "Less Common").Select(p => p.IdWithMotif))
                    {
                        ModificationTypeToColor[mod] = OxyColors.PowderBlue;
                    }
                    
                    foreach (var mod in GlobalVariables.AllModsKnownDictionary.Values.Where(p => p.ModificationType == "Common Artifact").Select(p => p.IdWithMotif))
                    {
                        ModificationTypeToColor[mod] = OxyColors.Teal;
                    }
                    
                    foreach (var mod in GlobalVariables.AllModsKnownDictionary.Values.Where(p => p.ModificationType == "Metal").Select(p => p.IdWithMotif))
                    {
                        ModificationTypeToColor[mod] = OxyColors.Maroon;
                    }
                    
                    foreach (var mod in GlobalVariables.AllModsKnownDictionary.Values.Where(p => p.ModificationType.Contains("Glycosylation")).Select(p => p.IdWithMotif))
                    {
                        ModificationTypeToColor[mod] = OxyColors.Maroon;
                    }

                    // setting individual specific
                    foreach (var mod in ModificationTypeToColor.Where(p => p.Key.Contains("Phosphorylation")))
                    {
                        ModificationTypeToColor[mod.Key] = OxyColors.Red;
                    }

                    foreach (var mod in ModificationTypeToColor.Where(p => p.Key.Contains("Acetylation")))
                    {
                        ModificationTypeToColor[mod.Key] = OxyColors.Purple; ;
                    }

                    ModificationTypeToColor["Carbamidomethyl on C"] = OxyColors.Green;
                    ModificationTypeToColor["Carbamidomethyl on U"] = OxyColors.Green;
                    ModificationTypeToColor["Oxidation on M"] = OxyColors.HotPink;

                    CoverageTypeToColor = CoverageTypes.ToDictionary(p => p, p => OxyColors.Blue);
                    CoverageTypeToColor["C-Terminal Color"] = OxyColors.Red;
                    CoverageTypeToColor["Internal Color"] = OxyColors.Purple;

                    UnannotatedPeakColor = OxyColors.LightGray;
                    InternalIonColor = OxyColors.Purple;

                #endregion

                // lines to be written on the spectrum
                SpectrumDescription = SpectrumDescriptors.ToDictionary(p => p, p => true);
                SpectrumDescription["Spectral Angle: "] = true;
            }
            
            // offset for annotation on base sequence
            ProductTypeToYOffset = ((ProductType[])Enum.GetValues(typeof(ProductType))).ToDictionary(p => p, p => 0.0);
            ProductTypeToYOffset[ProductType.aStar] = 35;
            ProductTypeToYOffset[ProductType.aDegree] = 36;
            ProductTypeToYOffset[ProductType.a] = 37;
            ProductTypeToYOffset[ProductType.aBaseLoss] = 38;
            ProductTypeToYOffset[ProductType.aWaterLoss] = 39;
            ProductTypeToYOffset[ProductType.b] = 40;
            ProductTypeToYOffset[ProductType.bBaseLoss] = 41;
            ProductTypeToYOffset[ProductType.bAmmoniaLoss] = 41;
            ProductTypeToYOffset[ProductType.bWaterLoss] = 42;
            ProductTypeToYOffset[ProductType.c] = 43;
            ProductTypeToYOffset[ProductType.cBaseLoss] = 44;
            ProductTypeToYOffset[ProductType.cWaterLoss] = 45;
            ProductTypeToYOffset[ProductType.d] = 46;
            ProductTypeToYOffset[ProductType.dBaseLoss] = 47;
            ProductTypeToYOffset[ProductType.dWaterLoss] = 48;
            ProductTypeToYOffset[ProductType.w] = -6;
            ProductTypeToYOffset[ProductType.wBaseLoss] = -7;
            ProductTypeToYOffset[ProductType.wWaterLoss] = -8;
            ProductTypeToYOffset[ProductType.x] = -9;
            ProductTypeToYOffset[ProductType.xBaseLoss] = -10;
            ProductTypeToYOffset[ProductType.xWaterLoss] = -11;
            ProductTypeToYOffset[ProductType.y] = -12;
            ProductTypeToYOffset[ProductType.yBaseLoss] = -13;
            ProductTypeToYOffset[ProductType.yAmmoniaLoss] = -13;
            ProductTypeToYOffset[ProductType.yWaterLoss] = -14;
            ProductTypeToYOffset[ProductType.z] = -15;
            ProductTypeToYOffset[ProductType.zBaseLoss] = -16;
            ProductTypeToYOffset[ProductType.zWaterLoss] = -17;
            ProductTypeToYOffset[ProductType.zDot] = -18;
            ProductTypeToYOffset[ProductType.zPlusOne] = -19;

            ProductTypeToXOffset = ((ProductType[])Enum.GetValues(typeof(ProductType))).ToDictionary(p => p, p => 0.0);
            ProductTypeToXOffset[ProductType.a] = -1.5;
            ProductTypeToXOffset[ProductType.aStar] = -1.5;
            ProductTypeToXOffset[ProductType.aDegree] = -1.5;
            ProductTypeToXOffset[ProductType.aBaseLoss] = -1.5;
            ProductTypeToXOffset[ProductType.aWaterLoss] = -1.5;
            ProductTypeToXOffset[ProductType.b] = -0.5;
            ProductTypeToXOffset[ProductType.bBaseLoss] = -0.5;
            ProductTypeToXOffset[ProductType.bWaterLoss] = -0.5;
            ProductTypeToXOffset[ProductType.bAmmoniaLoss] = -0.5;
            ProductTypeToXOffset[ProductType.c] = 0.5;
            ProductTypeToXOffset[ProductType.cBaseLoss] = 0.5;
            ProductTypeToXOffset[ProductType.cWaterLoss] = 0.5;
            ProductTypeToXOffset[ProductType.d] = 1.5;
            ProductTypeToXOffset[ProductType.dBaseLoss] = 1.5;
            ProductTypeToXOffset[ProductType.dWaterLoss] = 1.5;
            ProductTypeToXOffset[ProductType.w] = 1.5;
            ProductTypeToXOffset[ProductType.wBaseLoss] = 1.5;
            ProductTypeToXOffset[ProductType.wWaterLoss] = 1.5;
            ProductTypeToXOffset[ProductType.x] = 0.5;
            ProductTypeToXOffset[ProductType.xBaseLoss] = 0.5;
            ProductTypeToXOffset[ProductType.xWaterLoss] = 0.5;
            ProductTypeToXOffset[ProductType.y] = -0.5;
            ProductTypeToXOffset[ProductType.yBaseLoss] = -0.5;
            ProductTypeToXOffset[ProductType.yAmmoniaLoss] = -0.5;
            ProductTypeToXOffset[ProductType.yWaterLoss] = -0.5;
            ProductTypeToXOffset[ProductType.z] = -1.5;
            ProductTypeToXOffset[ProductType.zPlusOne] = -1.5;
            ProductTypeToXOffset[ProductType.zBaseLoss] = -1.5;
            ProductTypeToXOffset[ProductType.zWaterLoss] = -1.5;
            ProductTypeToXOffset[ProductType.zDot] = -1.5;

            PossibleColors = AllColors.ToDictionary(p => p, p => p.GetColorName());
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
                DisplayInternalIons = DisplayInternalIons,
                DisplayInternalIonAnnotations = DisplayInternalIonAnnotations,
                SubAndSuperScriptIons = SubAndSuperScriptIons,
                QValueFilter = QValueFilter,
                AmbiguityFilter = AmbiguityFilter,
                DrawStationarySequence = DrawStationarySequence,
                DrawNumbersUnderStationary = DrawNumbersUnderStationary,
                ShowLegend = ShowLegend,
                LocalizationLevelStart = LocalizationLevelStart,
                LocalizationLevelEnd = LocalizationLevelEnd,
                ExportType = ExportType,
                ProductTypeToColorValues = ProductTypeToColor.Values.Select(p => p.GetColorName()).ToList(),
                BetaProductTypeToColorValues = BetaProductTypeToColor.Values.Select(p => p.GetColorName()).ToList(),
                ModificationTypeToColorValues = ModificationTypeToColor.Values.Select(p => p.GetColorName()).ToList(),
                CoverageTypeToColorValues = CoverageTypeToColor.Values.Select(p => p.GetColorName()).ToList(),
                SpectrumDescriptionValues = SpectrumDescription.Values.ToList(),
                UnannotatedPeakColor = UnannotatedPeakColor.GetColorName(),
                InternalIonColor = InternalIonColor.GetColorName(),
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
            DisplayInternalIons = settings.DisplayInternalIons;
            DisplayInternalIonAnnotations = settings.DisplayInternalIonAnnotations;
            SubAndSuperScriptIons = settings.SubAndSuperScriptIons;
            QValueFilter = settings.QValueFilter;
            AmbiguityFilter = settings.AmbiguityFilter;
            DrawStationarySequence = settings.DrawStationarySequence;
            DrawNumbersUnderStationary = settings.DrawNumbersUnderStationary;
            ShowLegend = settings.ShowLegend;
            LocalizationLevelStart = settings.LocalizationLevelStart;
            LocalizationLevelEnd = settings.LocalizationLevelEnd;
            ExportType = settings.ExportType;

            ProductTypeToColor = ((ProductType[])Enum.GetValues(typeof(ProductType))).ToDictionary(p => p, p => DrawnSequence.ParseOxyColorFromName(settings.ProductTypeToColorValues[Array.IndexOf(((ProductType[])Enum.GetValues(typeof(ProductType))), p)]));
            BetaProductTypeToColor = ((ProductType[])Enum.GetValues(typeof(ProductType))).ToDictionary(p => p, p => DrawnSequence.ParseOxyColorFromName(settings.BetaProductTypeToColorValues[Array.IndexOf(((ProductType[])Enum.GetValues(typeof(ProductType))), p)]));
            ModificationTypeToColor = GlobalVariables.AllModsKnown.Select(p => p.IdWithMotif).ToDictionary(p => p, p => DrawnSequence.ParseOxyColorFromName(settings.ModificationTypeToColorValues[Array.IndexOf(GlobalVariables.AllModsKnown.Select(p => p.IdWithMotif).ToArray(), p)]));
            CoverageTypeToColor = CoverageTypes.ToDictionary(p => p, p => DrawnSequence.ParseOxyColorFromName(settings.CoverageTypeToColorValues[Array.IndexOf(CoverageTypes, p)]));
            SpectrumDescription = SpectrumDescriptors.ToDictionary(p => p, p => settings.SpectrumDescriptionValues[Array.IndexOf(SpectrumDescriptors, p)]);
            UnannotatedPeakColor = DrawnSequence.ParseOxyColorFromName(settings.UnannotatedPeakColor);
            InternalIonColor = DrawnSequence.ParseOxyColorFromName(settings.InternalIonColor);
        }

        /// <summary>
        /// Used to reset the settings to their default value, particullary needed for unit testing
        /// </summary>
        public static void ResetSettings()
        {
            InitializeDictionaries();
            DrawMatchedIons = true;
            DisplayIonAnnotations  = true;
            AnnotateMzValues = false;
            AnnotateCharges = false;
            AnnotationBold = false;
            DisplayInternalIons = true;
            DisplayInternalIonAnnotations = true;
            SubAndSuperScriptIons = true;
            DrawStationarySequence = true;
            DrawNumbersUnderStationary = true;
            ShowLegend = true;
            ShowDecoys = false;
            ShowContaminants = true;
            QValueFilter = 0.01;
            AmbiguityFilter = "No Filter";
            LocalizationLevelStart = LocalizationLevel.Level1;
            LocalizationLevelEnd = LocalizationLevel.Level3;
            DrawMatchedIons  = true;
            AnnotatedFontSize = 12;
            SequenceAnnotaitonResiduesPerSegment = 10;
            SequenceAnnotationSegmentPerRow = 3;
            ExportType = "Pdf";
            UnannotatedPeakColor = OxyColors.LightGray;
            InternalIonColor = OxyColors.Purple;
        }

        public static TValue GetValueOrDefault<TKey, TValue>(this IDictionary<TKey, TValue> dictionary, TKey key, TValue defaultValue)
        {
            return dictionary.TryGetValue(key, out var value) ? value : defaultValue;
        }
    }
}