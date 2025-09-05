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
}