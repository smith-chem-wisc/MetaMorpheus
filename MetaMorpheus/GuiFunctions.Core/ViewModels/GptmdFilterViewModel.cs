using EngineLayer;

namespace GuiFunctions;

public class GptmdFilterViewModel : BaseViewModel
{
    private bool _isSelected = false;
    public IGptmdFilter Filter { get; }
    public string Name { get; }
    public string Summary { get; }

    public bool IsSelected
    {
        get => _isSelected;
        set { _isSelected = value; OnPropertyChanged(nameof(IsSelected)); }
    }

    public GptmdFilterViewModel(IGptmdFilter filter, string name, string summary, bool isSelected = true)
    {
        Filter = filter;
        Name = name;
        Summary = summary;
        _isSelected = isSelected;
    }
}