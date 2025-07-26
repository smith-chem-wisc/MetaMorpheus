using System.IO;
using Easy.Common.Extensions;
using EngineLayer;
using GuiFunctions;
using Nett;
using TaskLayer;

namespace MetaMorpheusGUI;

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

    public bool AskAboutSpectralRecoveryParams
    {
        get => _current.AskAboutSpectralRecoveryParams;
        set { _current.AskAboutSpectralRecoveryParams = value; OnPropertyChanged(nameof(AskAboutSpectralRecoveryParams)); }
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

    public bool UseSpectralRecoveryParams
    {
        get => _current.UseSpectralRecoveryParams;
        set { _current.UseSpectralRecoveryParams = value; OnPropertyChanged(nameof(UseSpectralRecoveryParams)); }
    }

    #endregion

    #region IO
        
    private static readonly string SettingsPath = Path.Combine(GlobalVariables.DataDir, @"GUIsettings.toml");
    private static readonly string DefaultProteomePath = Path.Combine(GlobalVariables.DataDir, @"Proteomes");

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
        _loaded = Clone(_current);
    }

    // Save to disk
    public void Save()
    {
        Toml.WriteFile(_current, SettingsPath, MetaMorpheusTask.tomlConfig);
        _loaded = Clone(_current);
    }

    public bool IsDirty()
    {
        return !Equals(_current, _loaded);
    }

    public static bool SettingsFileExists() => File.Exists(SettingsPath);

    #endregion

    // Deep equality check (can be improved for more complex types)
    public override bool Equals(object obj)
    {
        if (obj is not GuiGlobalParams other) return false;
        return
            AskAboutUpdating == other.AskAboutUpdating &&
            AskBeforeExitingMetaMorpheus == other.AskBeforeExitingMetaMorpheus &&
            ProteomeDirectory == other.ProteomeDirectory &&
            AskAboutTopDownParams == other.AskAboutTopDownParams &&
            AskAboutChymotrypsinParams == other.AskAboutChymotrypsinParams &&
            AskAboutElastaseParams == other.AskAboutElastaseParams &&
            AskAboutNonSpecificParams == other.AskAboutNonSpecificParams &&
            AskAboutSemiTrypsinParams == other.AskAboutSemiTrypsinParams &&
            AskAboutArgCParams == other.AskAboutArgCParams &&
            AskAboutSpectralRecoveryParams == other.AskAboutSpectralRecoveryParams &&
            UseTopDownParams == other.UseTopDownParams &&
            UseChymotrypsinParams == other.UseChymotrypsinParams &&
            UseElastaseParams == other.UseElastaseParams &&
            UseNonSpecificParams == other.UseNonSpecificParams &&
            UseSemiTrypsinParams == other.UseSemiTrypsinParams &&
            UseArgCParams == other.UseArgCParams &&
            UseSpectralRecoveryParams == other.UseSpectralRecoveryParams;
    }

    public override int GetHashCode() => _instance.GetHashCode();

    // Helper for deep copy
    private static GuiGlobalParams Clone(GuiGlobalParams src)
    {
        return new GuiGlobalParams
        {
            AskAboutUpdating = src.AskAboutUpdating,
            AskBeforeExitingMetaMorpheus = src.AskBeforeExitingMetaMorpheus,
            ProteomeDirectory = src.ProteomeDirectory,
            AskAboutTopDownParams = src.AskAboutTopDownParams,
            AskAboutChymotrypsinParams = src.AskAboutChymotrypsinParams,
            AskAboutElastaseParams = src.AskAboutElastaseParams,
            AskAboutNonSpecificParams = src.AskAboutNonSpecificParams,
            AskAboutSemiTrypsinParams = src.AskAboutSemiTrypsinParams,
            AskAboutArgCParams = src.AskAboutArgCParams,
            AskAboutSpectralRecoveryParams = src.AskAboutSpectralRecoveryParams,
            UseTopDownParams = src.UseTopDownParams,
            UseChymotrypsinParams = src.UseChymotrypsinParams,
            UseElastaseParams = src.UseElastaseParams,
            UseNonSpecificParams = src.UseNonSpecificParams,
            UseSemiTrypsinParams = src.UseSemiTrypsinParams,
            UseArgCParams = src.UseArgCParams,
            UseSpectralRecoveryParams = src.UseSpectralRecoveryParams
        };
    }
}