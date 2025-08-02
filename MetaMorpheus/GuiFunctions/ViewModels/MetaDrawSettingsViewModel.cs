using EngineLayer;
using Omics.Fragmentation;
using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.IO;
using System.Linq;
using System.Runtime.CompilerServices;
using System.Text;
using System.Threading.Tasks;
using Nett;
using System.Windows.Input;
using MassSpectrometry;
using GuiFunctions.MetaDraw;
using Readers;
using System.Globalization;
using System.Windows;

namespace GuiFunctions
{
    /// <summary>
    /// Model for the metadraw settings window
    /// In a perfect world, nothing would be in the MetaDrawSettings.xaml.cs file and all would be located here
    /// This would allow for more extensive testing of the GUI Elements
    /// </summary>
    public class MetaDrawSettingsViewModel : BaseViewModel, IAsyncInitialization
    {
        private static readonly Lazy<MetaDrawSettingsViewModel> _instance =
            new(() => new MetaDrawSettingsViewModel());

        // Singleton to ensure only one instance of MetaDrawSettingsViewModel exists
        public static MetaDrawSettingsViewModel Instance => _instance.Value;

        #region Private Properties

        private ObservableCollection<ModTypeForTreeViewModel> _Modifications = new ObservableCollection<ModTypeForTreeViewModel>();
        private ObservableCollection<IonTypeForTreeViewModel> _IonGroups = new ObservableCollection<IonTypeForTreeViewModel>();
        private ObservableCollection<CoverageTypeForTreeViewModel> _CoverageColors = new ObservableCollection<CoverageTypeForTreeViewModel>();
        private bool _LoadedIons { get { return (_IonGroups.Count > 0); } }
        private bool _LoadedPTMs { get { return (_Modifications.Count > 0); } }
        private bool _LoadedSequenceCoverage { get { return (_CoverageColors.Count > 0); } }

        #endregion

        #region Public Properties

        public ObservableCollection<ModTypeForTreeViewModel> Modifications
        {
            get { return _Modifications; }
            set
            {
                _Modifications = value;
                OnPropertyChanged(nameof(Modifications));
            } 
        }

        public ObservableCollection<IonTypeForTreeViewModel> IonGroups
        {
            get { return _IonGroups; }
            set
            {
                _IonGroups = value;
                OnPropertyChanged(nameof(IonGroups));
            }
        }

        public ObservableCollection<CoverageTypeForTreeViewModel> CoverageColors
        {
            get { return _CoverageColors; }
            set
            {
                _CoverageColors = value;
                OnPropertyChanged(nameof(CoverageColors));
            }
        }

        public ObservableCollection<SpectrumDescriptorViewModel> SpectrumDescriptors { get; }
        public ObservableCollection<string> ExportTypes { get; } = [.. MetaDrawSettings.ExportTypes];

        public DeconHostViewModel DeconHostViewModel { get; set; }

        public ObservableCollection<string> PossibleColors { get; set; }
        public ObservableCollection<LegendDisplayProperty> ChimericLegendDisplayProperties { get; } = [..Enum.GetValues<LegendDisplayProperty>()];
        public ObservableCollection<string> AmbiguityFilters { get; } =  [..MetaDrawSettings.AmbiguityTypes];
        public ObservableCollection<LocalizationLevel> GlycanLocalizationLevels { get; } = [.. Enum.GetValues<LocalizationLevel>()];

        public bool HasDefaultSaved { get { return File.Exists(SettingsPath); } }
        public bool CanOpen { get { return (_LoadedIons && _LoadedPTMs && _LoadedSequenceCoverage); } }
        public Task Initialization { get; private set; }
        public static string SettingsPath = Path.Combine(GlobalVariables.DataDir, "DefaultParameters", @"MetaDrawSettingsDefault.xml");

        public bool ShowDecoys
        {
            get => MetaDrawSettings.ShowDecoys;
            set { MetaDrawSettings.ShowDecoys = value; OnPropertyChanged(nameof(ShowDecoys)); }
        }

        public bool ShowContaminants
        {
            get => MetaDrawSettings.ShowContaminants;
            set { MetaDrawSettings.ShowContaminants = value; OnPropertyChanged(nameof(ShowContaminants)); }
        }

        public double QValueFilter
        {
            get => MetaDrawSettings.QValueFilter;
            set 
            {
                MetaDrawSettings.QValueFilter = value;
                OnPropertyChanged(nameof(QValueFilter));
            }
        }

        public string AmbiguityFilter
        {
            get => MetaDrawSettings.AmbiguityFilter;
            set { MetaDrawSettings.AmbiguityFilter = value; OnPropertyChanged(nameof(AmbiguityFilter)); }
        }

        public LocalizationLevel LocalizationLevelStart
        {
            get => MetaDrawSettings.LocalizationLevelStart;
            set { MetaDrawSettings.LocalizationLevelStart = value; OnPropertyChanged(nameof(LocalizationLevelStart)); }
        }

        public LocalizationLevel LocalizationLevelEnd
        {
            get => MetaDrawSettings.LocalizationLevelEnd;
            set { MetaDrawSettings.LocalizationLevelEnd = value; OnPropertyChanged(nameof(LocalizationLevelEnd)); }
        }

        public string ExportType
        {
            get => MetaDrawSettings.ExportType;
            set { MetaDrawSettings.ExportType = value; OnPropertyChanged(nameof(ExportType)); }
        }
        public double Dpi
        {
            get => MetaDrawSettings.CanvasPdfExportDpi;
            set { MetaDrawSettings.CanvasPdfExportDpi = value; OnPropertyChanged(nameof(Dpi)); }
        }
        public bool DisplayIonAnnotations
        {
            get => MetaDrawSettings.DisplayIonAnnotations;
            set { MetaDrawSettings.DisplayIonAnnotations = value; OnPropertyChanged(nameof(DisplayIonAnnotations)); }
        }
        public bool AnnotateMzValues
        {
            get => MetaDrawSettings.AnnotateMzValues;
            set { MetaDrawSettings.AnnotateMzValues = value; OnPropertyChanged(nameof(AnnotateMzValues)); }
        }
        public bool SuppressMessageBoxes
        {
            get => MetaDrawSettings.SuppressMessageBoxes;
            set 
            {
                MessageBoxHelper.SuppressMessageBoxes = value;
                MetaDrawSettings.SuppressMessageBoxes = value; 
                OnPropertyChanged(nameof(AnnotateMzValues)); 
            }
        }
        public bool AnnotateCharges
        {
            get => MetaDrawSettings.AnnotateCharges;
            set { MetaDrawSettings.AnnotateCharges = value; OnPropertyChanged(nameof(AnnotateCharges)); }
        }
        public bool DisplayInternalIonAnnotations
        {
            get => MetaDrawSettings.DisplayInternalIonAnnotations;
            set { MetaDrawSettings.DisplayInternalIonAnnotations = value; OnPropertyChanged(nameof(DisplayInternalIonAnnotations)); }
        }
        public bool DisplayInternalIons
        {
            get => MetaDrawSettings.DisplayInternalIons;
            set { MetaDrawSettings.DisplayInternalIons = value; OnPropertyChanged(nameof(DisplayInternalIons)); }
        }
        public bool AnnotationBold
        {
            get => MetaDrawSettings.AnnotationBold;
            set { MetaDrawSettings.AnnotationBold = value; OnPropertyChanged(nameof(AnnotationBold)); }
        }
        public bool SubAndSuperScriptIons
        {
            get => MetaDrawSettings.SubAndSuperScriptIons;
            set { MetaDrawSettings.SubAndSuperScriptIons = value; OnPropertyChanged(nameof(SubAndSuperScriptIons)); }
        }

        public int AnnotatedFontSize
        {
            get => MetaDrawSettings.AnnotatedFontSize;
            set { MetaDrawSettings.AnnotatedFontSize = value; OnPropertyChanged(nameof(AnnotatedFontSize)); }
        }
        public int AxisLabelTextSize
        {
            get => MetaDrawSettings.AxisLabelTextSize;
            set { MetaDrawSettings.AxisLabelTextSize = value; OnPropertyChanged(nameof(AxisLabelTextSize)); }
        }
        public int AxisTitleTextSize
        {
            get => MetaDrawSettings.AxisTitleTextSize;
            set { MetaDrawSettings.AxisTitleTextSize = value; OnPropertyChanged(nameof(AxisTitleTextSize)); }
        }
        public double StrokeThicknessAnnotated
        {
            get => MetaDrawSettings.StrokeThicknessAnnotated;
            set { MetaDrawSettings.StrokeThicknessAnnotated = value; OnPropertyChanged(nameof(StrokeThicknessAnnotated)); }
        }
        public double StrokeThicknessUnannotated
        {
            get => MetaDrawSettings.StrokeThicknessUnannotated;
            set { MetaDrawSettings.StrokeThicknessUnannotated = value; OnPropertyChanged(nameof(StrokeThicknessUnannotated)); }
        }

        // Chimera Settings
        public bool DisplayChimeraLegend
        {
            get => MetaDrawSettings.DisplayChimeraLegend;
            set { MetaDrawSettings.DisplayChimeraLegend = value; OnPropertyChanged(nameof(DisplayChimeraLegend)); }
        }
        public bool ChimeraLegendTakeFirstIfAmbiguous
        {
            get => MetaDrawSettings.ChimeraLegendTakeFirstIfAmbiguous;
            set { MetaDrawSettings.ChimeraLegendTakeFirstIfAmbiguous = value; OnPropertyChanged(nameof(ChimeraLegendTakeFirstIfAmbiguous)); }
        }
        public double ChimeraLegendMaxWidth
        {
            get => MetaDrawSettings.ChimeraLegendMaxWidth;
            set { MetaDrawSettings.ChimeraLegendMaxWidth = value; OnPropertyChanged(nameof(ChimeraLegendMaxWidth)); }
        }
        public LegendDisplayProperty ChimeraLegendMainTextType
        {
            get => MetaDrawSettings.ChimeraLegendMainTextType;
            set { MetaDrawSettings.ChimeraLegendMainTextType = value; OnPropertyChanged(nameof(ChimeraLegendMainTextType)); }
        }
        public LegendDisplayProperty ChimeraLegendSubTextType
        {
            get => MetaDrawSettings.ChimeraLegendSubTextType;
            set { MetaDrawSettings.ChimeraLegendSubTextType = value; OnPropertyChanged(nameof(ChimeraLegendSubTextType)); }
        }

        // Data Visualization Settings
        public bool DisplayFilteredOnly
        {
            get => MetaDrawSettings.DisplayFilteredOnly;
            set { MetaDrawSettings.DisplayFilteredOnly = value; OnPropertyChanged(nameof(DisplayFilteredOnly)); }
        }

        public bool NormalizeHistogramToFile
        {
            get => MetaDrawSettings.NormalizeHistogramToFile;
            set { MetaDrawSettings.NormalizeHistogramToFile = value; OnPropertyChanged(nameof(NormalizeHistogramToFile)); }
        }

        #endregion

        #region Constructor

        /// <summary>
        /// Constructs the instance asynchronously
        /// </summary>
        /// <param name="loadAsync"></param>
        internal MetaDrawSettingsViewModel(bool loadAsync = true)
        {
            SpectrumDescriptors = [.. MetaDrawSettings.SpectrumDescription.Select(p => new SpectrumDescriptorViewModel(p.Key))];

            if (loadAsync)
                Initialization = InitializeAsync();
            else
            {
                if (HasDefaultSaved)
                    LoadSettings();

                PossibleColors = new ObservableCollection<string>(MetaDrawSettings.PossibleColors.Values.ToList());
                AddSpaces(PossibleColors);

                LoadPTMs();
                LoadIonTypes();
                LoadSequenceCoverage();
                Initialization = Task.CompletedTask;
            }

            // This defaults to classic decon, and we set the charge to ensure it will work for top-down and bottom-up.
            // This is not the best approach, in the future we could try to locate the search toml when loading in a psm file and use those decon params. 
            DeconHostViewModel = new();
             // Ensure it will work for top-down and bottom-up.
            DeconHostViewModel.SetAllPrecursorMaxChargeState(60); 
        }

        private async Task InitializeAsync()
        {
            if (HasDefaultSaved)
                LoadSettings();

            PossibleColors = new ObservableCollection<string>(MetaDrawSettings.PossibleColors.Values.ToList());
            AddSpaces(PossibleColors);

            LoadPTMs();
            LoadIonTypes();
            LoadSequenceCoverage();
            await Task.Delay(100);
        }

        #endregion

        #region Methods

        /// <summary>
        /// Method to be executed when the Save Command is fired
        /// Takes the dynamic data of the settings m
        /// </summary>
        public void Save()
        {
            // save ion colors if changed
            foreach (var group in IonGroups)
            {
                foreach (var ion in group)
                {
                    if (ion.HasChanged)
                    {
                        if (ion.IonName.Equals("Unannotated Peak"))
                            MetaDrawSettings.UnannotatedPeakColor = DrawnSequence.ParseOxyColorFromName(ion.SelectedColor.Replace(" ", ""));
                        else if (ion.IonName.Equals("Internal Ion"))
                            MetaDrawSettings.InternalIonColor = DrawnSequence.ParseOxyColorFromName(ion.SelectedColor.Replace(" ", ""));
                        else if (ion.IsBeta)
                            MetaDrawSettings.BetaProductTypeToColor[ion.IonType] = DrawnSequence.ParseOxyColorFromName(ion.SelectedColor.Replace(" ", ""));
                        else
                            MetaDrawSettings.ProductTypeToColor[ion.IonType] = DrawnSequence.ParseOxyColorFromName(ion.SelectedColor.Replace(" ", ""));
                    }
                }
            }

            // save modification colors if changed
            foreach (var group in Modifications)
            {
                foreach (var mod in group)
                {
                    if (mod.HasChanged)
                        MetaDrawSettings.ModificationTypeToColor[mod.ModName] = DrawnSequence.ParseOxyColorFromName(mod.SelectedColor.Replace(" ", ""));
                }
            }

            // save sequence coverage colors if changed
            foreach (var color in CoverageColors)
            {
                if (color.HasChanged)
                {
                    MetaDrawSettings.CoverageTypeToColor[color.Name] = DrawnSequence.ParseOxyColorFromName(color.SelectedColor.Replace(" ", ""));
                }
            }
        }

        /// <summary>
        /// Method to be executed when the SaveAsDefault Command is fired
        /// </summary>
        public void SaveAsDefault()
        {
            Save();
            MetaDrawSettingsSnapshot settings = MetaDrawSettings.MakeSnapShot();
            string directoryPath = Path.GetDirectoryName(SettingsPath);
            if (!Directory.Exists(directoryPath))
            {
                Directory.CreateDirectory(directoryPath);
            }

            XmlReaderWriter.WriteToXmlFile<MetaDrawSettingsSnapshot>(SettingsPath, settings);
        }

        /// <summary>
        /// Loads default settings from their location if called
        /// </summary>
        public void LoadSettings()
        {
            MetaDrawSettingsSnapshot settings = null;
            settings = XmlReaderWriter.ReadFromXmlFile<MetaDrawSettingsSnapshot>(SettingsPath);
            MetaDrawSettings.LoadSettings(settings, out bool flaggedErrorOnRead);

            if (flaggedErrorOnRead)
                SaveAsDefault();
        }

      
        public void LoadIonTypes()
        {
            var ions = ((ProductType[])Enum.GetValues(typeof(ProductType)));
            var common = ions.Where(p => p.ToString().Equals("a") || p.ToString().Equals("b") || p.ToString().Equals("c")
                                          || p.ToString().Equals("x") || p.ToString().Equals("y") || p.ToString().Equals("zDot"));
            var lessCommon = ions.Where(p => !common.Any(m => m == p));
            _IonGroups.Add(new IonTypeForTreeViewModel("Common Ions", common, false));
            _IonGroups.Add(new IonTypeForTreeViewModel("Less Common Ions", lessCommon, false));
            _IonGroups.Add(new IonTypeForTreeViewModel("Cross Linked Beta Peptide", ions, true));
        }

        public void LoadPTMs()
        {
            var modGroups = GlobalVariables.AllModsKnown.GroupBy(b => b.ModificationType);
            foreach (var group in modGroups)
            {
                var theModType = new ModTypeForTreeViewModel(group.Key, false);
                _Modifications.Add(theModType);
                foreach (var mod in group)
                {
                    theModType.Children.Add(new ModForTreeViewModel(mod.ToString(), false, mod.IdWithMotif, false, theModType));
                }
            }
        }

        public void LoadSequenceCoverage()
        {
            _CoverageColors.Add(new CoverageTypeForTreeViewModel("N-Terminal Color"));
            _CoverageColors.Add(new CoverageTypeForTreeViewModel("C-Terminal Color"));
            _CoverageColors.Add(new CoverageTypeForTreeViewModel("Internal Color"));
        }

        #endregion
    }


    /// <summary>
    /// Marks a type as requiring asynchronous initialization and provides the result of that initialization.
    /// </summary>
    public interface IAsyncInitialization
    {
        /// <summary>
        /// The result of the asynchronous initialization of this instance.
        /// </summary>
        Task Initialization { get; }
    }

    /// <summary>
    /// Class to represent each component of the Spectrum Description
    /// </summary>
    /// <param name="key"></param>
    public class SpectrumDescriptorViewModel(string key) : BaseViewModel
    {
        private string _displayName;
        public string DisplayName => _displayName ??= key.Split(":")[0];

        public bool IsSelected
        {
            get => MetaDrawSettings.SpectrumDescription[key];
            set
            {
                MetaDrawSettings.SpectrumDescription[key] = value;
                OnPropertyChanged(nameof(IsSelected));
            }
        }
    }
}
