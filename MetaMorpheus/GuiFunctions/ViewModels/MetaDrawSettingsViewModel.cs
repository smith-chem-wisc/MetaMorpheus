using EngineLayer;
using GuiFunctions;
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

        public ObservableCollection<string> PossibleColors { get; set; }
        public bool HasDefaultSaved { get { return File.Exists(SettingsPath); } }
        public bool CanOpen { get { return (_LoadedIons && _LoadedPTMs && _LoadedSequenceCoverage); } }
        public Task Initialization { get; private set; }
        public static string SettingsPath = Path.Combine(GlobalVariables.DataDir, "DefaultParameters", @"MetaDrawSettingsDefault.xml");

        #endregion

        #region Constructor

        /// <summary>
        /// Constructs the instance asynchronously
        /// </summary>
        /// <param name="loadAsync"></param>
        public MetaDrawSettingsViewModel(bool loadAsync = true)
        {
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
}
