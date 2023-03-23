using EngineLayer;
using EngineLayer.GlycoSearch;
using OxyPlot;
using Proteomics;
using Proteomics.Fragmentation;
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

        #region Customizable Settings

        // graphic settings
        public static Dictionary<string, bool> SpectrumDescription { get; set; }
        public static bool DisplayIonAnnotations { get; set; } = true;
        public static bool AnnotateMzValues { get; set; } = false;
        public static bool AnnotateCharges { get; set; } = false;
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
        public static OxyColor VariantCrossColor { get; set; } = OxyColors.Green;
        public static OxyColor UnannotatedPeakColor { get; set; } = OxyColors.LightGray;
        public static OxyColor InternalIonColor { get; set; } = OxyColors.Purple;
        public static SolidColorBrush ModificationAnnotationColor { get; set; } = Brushes.Orange;
        public static double CanvasPdfExportDpi { get; set; } = 300;
        public static double StrokeThicknessUnannotated { get; set; } = 0.7;
        public static double StrokeThicknessAnnotated { get; set; } = 1.0;
        public static double AnnotatedSequenceTextSpacing { get; set; } = 22;
        public static int AnnotatedFontSize { get; set; } = 12;
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
                SpectrumDescription["Spectral Angle: "] = false;
            }
            
            // offset for annotation on base sequence
            ProductTypeToYOffset = ((ProductType[])Enum.GetValues(typeof(ProductType))).ToDictionary(p => p, p => 0.0);
            ProductTypeToYOffset[ProductType.b] = 40;
            ProductTypeToYOffset[ProductType.y] = -10;
            ProductTypeToYOffset[ProductType.c] = 43.6;
            ProductTypeToYOffset[ProductType.zDot] = -13.6;

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