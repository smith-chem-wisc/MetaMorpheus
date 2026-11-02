using System;
using System.Collections.ObjectModel;

namespace GuiFunctions;

public class LoadingProgressViewModel : BaseViewModel
{
    private static readonly Lazy<LoadingProgressViewModel> _instance = new(() => new LoadingProgressViewModel());
    public static LoadingProgressViewModel Instance => _instance.Value;

    public ObservableCollection<LoadingStepViewModel> Steps { get; } = new();

    private bool _isVisible;
    public bool IsVisible
    {
        get => _isVisible;
        set
        {
            if (_isVisible == value) return;
            _isVisible = value;
            OnPropertyChanged(nameof(IsVisible));
        }
    }

    private LoadingProgressViewModel() { }
}


public class LoadingStepViewModel : BaseViewModel
{
    private int _current;
    public int Current
    {
        get => _current;
        set { _current = value; OnPropertyChanged(nameof(Current)); }
    }

    private int _total;
    public int Total
    {
        get => _total;
        set 
        {
            if (value == _total) return;
            _total = value; 
            OnPropertyChanged(nameof(Total)); 
        }
    }

    private string _stepName;
    public string StepName
    {
        get => _stepName;
        set 
        {
            if (_stepName == value) return;
            _stepName = value; 
            OnPropertyChanged(nameof(StepName));
        }
    }
}
