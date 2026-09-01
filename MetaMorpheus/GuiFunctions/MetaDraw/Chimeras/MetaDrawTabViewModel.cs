using System.IO;

namespace GuiFunctions.MetaDraw;

public abstract class MetaDrawTabViewModel : BaseViewModel
{
    protected MetaDrawTabViewModel(bool isTabEnabled = true)
    {
        IsTabEnabled = isTabEnabled;
    }

    protected object ThreadLocker { get; } = new();
    public abstract string TabHeader { get; init; }

    // Used when tab is loading externally to disable the tab (from MetaDrawDataLoader)
    private bool _isTabEnabled;
    public bool IsTabEnabled
    {
        get => _isTabEnabled;
        set { _isTabEnabled = value; OnPropertyChanged(nameof(IsTabEnabled)); }
    }

    private string _exportDirectory;
    public string ExportDirectory
    {
        get
        {
            if (!Directory.Exists(_exportDirectory))
                Directory.CreateDirectory(_exportDirectory);
            return _exportDirectory;
        }
        set
        {
            _exportDirectory = value;
            OnPropertyChanged(nameof(ExportDirectory));
        }
    }

    // Used when the tab is loading internally (adding a database or loading a file)
    private bool _isLoading;
    public bool IsLoading
    {
        get => _isLoading;
        set
        {
            _isLoading = value;
            OnPropertyChanged(nameof(IsLoading));
        }
    }
}