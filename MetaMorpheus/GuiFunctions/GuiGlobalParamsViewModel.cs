#nullable enable
using System;
using System.Collections.ObjectModel;
using System.IO;
using System.Windows;
using Easy.Common.Extensions;
using EngineLayer;
using GuiFunctions.Util;
using Nett;
using TaskLayer;

namespace GuiFunctions;

public class GuiGlobalParamsViewModel : BaseViewModel
{
    private static GuiGlobalParamsViewModel _instance;
    private GuiGlobalParams _current;
    private GuiGlobalParams _loaded;

    public static GuiGlobalParamsViewModel Instance
    {
        get
        {
            if (_instance == null)
            {
                _instance = new GuiGlobalParamsViewModel();
                _instance.Load();
            }
            return _instance;
        }
    }

    // private to ensure only one instance ever exists. 
    private GuiGlobalParamsViewModel() { }

    public string MainWindowTitle => IsRnaMode
        ? $"MetaMorpheus RNA: {GlobalVariables.MetaMorpheusVersion}"
        : $"MetaMorpheus Protein: {GlobalVariables.MetaMorpheusVersion}";

    #region Parameters

    public bool AskAboutUpdating
    {
        get => _current.AskAboutUpdating;
        set { _current.AskAboutUpdating = value; OnPropertyChanged(nameof(AskAboutUpdating));}
    }

    public bool AskBeforeExitingMetaMorpheus
    {
        get => _current.AskBeforeExitingMetaMorpheus;
        set { _current.AskBeforeExitingMetaMorpheus = value; OnPropertyChanged(nameof(AskBeforeExitingMetaMorpheus)); }
    }

    public string ProteomeDirectory
    {
        get
        {
            if (_current.ProteomeDirectory.IsNullOrEmpty())
                _current.ProteomeDirectory = DefaultProteomePath;
            return _current.ProteomeDirectory;
        }
        set
        {
            if (value.IsNullOrEmpty())
                return;
            if (Directory.Exists(value)) 
                _current.ProteomeDirectory = value;
            OnPropertyChanged(nameof(ProteomeDirectory));
        }
    }

    public bool AskAboutTopDownParams
    {
        get => _current.AskAboutTopDownParams;
        set { _current.AskAboutTopDownParams = value; OnPropertyChanged(nameof(AskAboutTopDownParams)); }
    }

    public bool AskAboutChymotrypsinParams
    {
        get => _current.AskAboutChymotrypsinParams;
        set { _current.AskAboutChymotrypsinParams = value; OnPropertyChanged(nameof(AskAboutChymotrypsinParams)); }
    }

    public bool AskAboutElastaseParams
    {
        get => _current.AskAboutElastaseParams;
        set { _current.AskAboutElastaseParams = value; OnPropertyChanged(nameof(AskAboutElastaseParams)); }
    }

    public bool AskAboutNonSpecificParams
    {
        get => _current.AskAboutNonSpecificParams;
        set { _current.AskAboutNonSpecificParams = value; OnPropertyChanged(nameof(AskAboutNonSpecificParams)); }
    }

    public bool AskAboutSemiTrypsinParams
    {
        get => _current.AskAboutSemiTrypsinParams;
        set { _current.AskAboutSemiTrypsinParams = value; OnPropertyChanged(nameof(AskAboutSemiTrypsinParams)); }
    }

    public bool AskAboutArgCParams
    {
        get => _current.AskAboutArgCParams;
        set { _current.AskAboutArgCParams = value; OnPropertyChanged(nameof(AskAboutArgCParams)); }
    }

    public bool AskAboutModeSwitch
    {
        get => _current.AskAboutModeSwitch;
        set { _current.AskAboutModeSwitch = value; OnPropertyChanged(nameof(AskAboutModeSwitch)); }
    }

    public bool AskAboutOverwritingOutputDirectory
    {
        get => _current.AskAboutOverwritingOutputDirectory;
        set { _current.AskAboutOverwritingOutputDirectory = value; OnPropertyChanged(nameof(AskAboutOverwritingOutputDirectory)); }
    }

    public bool UseTopDownParams
    {
        get => _current.UseTopDownParams;
        set { _current.UseTopDownParams = value; OnPropertyChanged(nameof(UseTopDownParams)); }
    }

    public bool UseChymotrypsinParams
    {
        get => _current.UseChymotrypsinParams;
        set { _current.UseChymotrypsinParams = value; OnPropertyChanged(nameof(UseChymotrypsinParams)); }
    }

    public bool UseElastaseParams
    {
        get => _current.UseElastaseParams;
        set { _current.UseElastaseParams = value; OnPropertyChanged(nameof(UseElastaseParams)); }
    }

    public bool UseNonSpecificParams
    {
        get => _current.UseNonSpecificParams;
        set { _current.UseNonSpecificParams = value; OnPropertyChanged(nameof(UseNonSpecificParams)); }
    }

    public bool UseSemiTrypsinParams
    {
        get => _current.UseSemiTrypsinParams;
        set { _current.UseSemiTrypsinParams = value; OnPropertyChanged(nameof(UseSemiTrypsinParams)); }
    }

    public bool UseArgCParams
    {
        get => _current.UseArgCParams;
        set { _current.UseArgCParams = value; OnPropertyChanged(nameof(UseArgCParams)); }
    }

    public bool OverwriteOutputDirectory
    {
        get => _current.OverwriteOutputDirectory;
        set { _current.OverwriteOutputDirectory = value; OnPropertyChanged(nameof(OverwriteOutputDirectory)); }
    }


    private AnalyteType? _previous;
    public bool IsRnaMode
    {
        get => _current.IsRnaMode;
        set
        {   
            // Invoke the event to check if the user wants to switch modes
            var args = new ModeSwitchRequestEventArgs();

            // Ask the GUI how to move forward
            // - If we have a default saved and are told not to ask, it will skip the pop-up
            // - if no files are loaded it will tell us to switch, otherwise it will trigger a pop-up
            RequestModeSwitchConfirmation?.Invoke(this, args);

            if (args.RememberMyDecision)
            {
                AskAboutModeSwitch = false;
                CachedModeSwitchResult = args.Result;
            }

            // Do not switch - Force UI to refresh to original state
            if (args.Result == ModeSwitchResult.Cancel)
            {
                OnPropertyChanged(nameof(IsRnaMode));
                return;
            }
            
            _current.IsRnaMode = value;
            if (_current.IsRnaMode)
                GlobalVariables.AnalyteType = AnalyteType.Oligo;
            else
                GlobalVariables.AnalyteType = _previous ??= AnalyteType.Peptide;

            OnPropertyChanged(nameof(IsRnaMode));
            OnPropertyChanged(nameof(MainWindowTitle));
        }
    }

    public ModeSwitchResult CachedModeSwitchResult
    {
        get => _current.CachedModeSwitchResult;
        set {_current.CachedModeSwitchResult = value; OnPropertyChanged(nameof(CachedModeSwitchResult)); }
    }

    public ObservableCollection<ModeSwitchResult> AllModeSwitchValues { get; } = new(Enum.GetValues<ModeSwitchResult>());
    public static EventHandler<ModeSwitchRequestEventArgs>? RequestModeSwitchConfirmation { get; set; }

    #endregion

    #region IO

    // Lazy-initialize these to avoid null reference during designer mode
    private static string _settingsPath;
    private static string _defaultProteomePath;

    private static string SettingsPath
    {
        get
        {
            if (string.IsNullOrEmpty(_settingsPath))
            {
                _settingsPath = Path.Combine(GlobalVariables.DataDir, @"GUIsettings.toml");
            }
            return _settingsPath;
        }
    }

    private static string DefaultProteomePath
    {
        get
        {
            if (string.IsNullOrEmpty(_defaultProteomePath))
            {
                _defaultProteomePath = Path.Combine(GlobalVariables.DataDir, @"Proteomes");
            }
            return _defaultProteomePath;
        }
    }

    // Load from disk
    public void Load()
    {
        if (File.Exists(SettingsPath))
        {
            try
            {
                _current = Toml.ReadFile<GuiGlobalParams>(SettingsPath);
            }
            catch
            {
                _current = new GuiGlobalParams
                {
                    ProteomeDirectory = DefaultProteomePath
                };
                Toml.WriteFile(_current, SettingsPath, MetaMorpheusTask.tomlConfig);
            }
        }
        else
        {
            _current = new GuiGlobalParams
            {
                ProteomeDirectory = DefaultProteomePath
            };
            Toml.WriteFile(_current, SettingsPath, MetaMorpheusTask.tomlConfig);
        }        
        
        GlobalVariables.DecoyIdentifier = _current.DecoyIdentifier ??= "DECOY";
        IsRnaMode = _current.IsRnaMode;
        _loaded = _current.Clone();
    }

    // Save to disk
    public void Save()
    {
        Toml.WriteFile(_current, SettingsPath, MetaMorpheusTask.tomlConfig);
        _loaded = _current.Clone();
    }

    public bool IsDirty()
    {
        return !_current.Equals(_loaded);
    }

    public static bool SettingsFileExists() => File.Exists(SettingsPath);

    #endregion

}