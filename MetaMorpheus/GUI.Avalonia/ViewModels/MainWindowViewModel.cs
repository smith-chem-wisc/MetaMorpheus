using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.IO;
using System.Linq;
using System.Threading.Tasks;
using Avalonia.Threading;
using CommunityToolkit.Mvvm.ComponentModel;
using CommunityToolkit.Mvvm.Input;
using EngineLayer;
using EngineLayer.DatabaseLoading;
using TaskLayer;

namespace MetaMorpheus.Avalonia.ViewModels;

/// <summary>
/// The main window's state and behaviour, with no UI types.
///
/// The WPF MainWindow.xaml.cs is 2,284 lines of collections, event wiring and click handlers in one
/// file, none of it reachable without a Window. Keeping it here means the run flow and the progress
/// wiring can be tested headlessly - see GUI.Avalonia.Tests.
/// </summary>
internal sealed partial class MainWindowViewModel : ObservableObject, IDisposable
{
    // Held so they can be detached again; static events otherwise keep this instance alive forever.
    private EventHandler<StringEventArgs> _onLog;
    private EventHandler<StringEventArgs> _onWarn;
    private EventHandler<StringEventArgs> _onStatus;
    private EventHandler<ProgressEventArgs> _onProgress;
    private EventHandler<SingleTaskEventArgs> _onTaskStarting;
    private EventHandler<SingleTaskEventArgs> _onTaskFinished;
    private EventHandler<StringEventArgs> _onRunnerWarn;
    private EventHandler<StringEventArgs> _onAllFinished;

    public ObservableCollection<SpectraFileForDisplay> SpectraFiles { get; } = new();
    public ObservableCollection<DatabaseForDisplay> Databases { get; } = new();
    public ObservableCollection<TaskForDisplay> Tasks { get; } = new();

    [ObservableProperty] private string _outputFolder = string.Empty;
    [ObservableProperty] private string _log = string.Empty;
    [ObservableProperty] private string _status = "Ready.";
    [ObservableProperty] private double _progress;
    [ObservableProperty] private bool _isRunning;

    public string Version => GlobalVariables.MetaMorpheusVersion;

    public MainWindowViewModel()
    {
        // Every entry point must do this before anything else: it resolves the data directory and
        // loads modifications, proteases, crosslinkers and glycans. Verified to work on arm64.
        GlobalVariables.SetUpGlobalVariables();
        Append($"MetaMorpheus {GlobalVariables.MetaMorpheusVersion}");
        Append($"{GlobalVariables.AllModsKnown.Count()} modifications and "
            + $"{GlobalVariables.Crosslinkers.Count()} crosslinkers loaded.");
        SubscribeToEngineEvents();
    }

    /// <summary>
    /// The engines report through static events rather than return values, so a front end has to
    /// subscribe. Two things this has to get right, both of which cost real bugs when they were
    /// missing:
    ///
    ///   * The events are static, so a subscription outlives the view model unless it is removed. The
    ///     handlers are kept in fields and detached in Dispose, otherwise a discarded view model keeps
    ///     receiving engine callbacks - which shows up as duplicated log lines, and as one test's run
    ///     writing into another's view model.
    ///   * They are raised on the engines' worker threads, so every property update is marshalled to
    ///     the UI thread. The WPF window does the same thing with Dispatcher.Invoke; without it,
    ///     Avalonia throws "the calling thread cannot access this object".
    /// </summary>
    private void SubscribeToEngineEvents()
    {
        _onLog = (_, e) => OnUiThread(() => Append(e.S));
        _onWarn = (_, e) => OnUiThread(() => Append($"WARNING: {e.S}"));
        _onStatus = (_, e) => OnUiThread(() => Status = e.S);
        _onProgress = (_, e) => OnUiThread(() => Progress = e.NewProgress);
        _onTaskStarting = (_, e) => OnUiThread(() => Append($"Starting task: {e.DisplayName}"));
        _onTaskFinished = (_, e) => OnUiThread(() => Append($"Finished task: {e.DisplayName}"));
        _onRunnerWarn = (_, e) => OnUiThread(() => Append($"WARNING: {e.S}"));
        _onAllFinished = (_, _) => OnUiThread(() => Append("All tasks finished."));

        MetaMorpheusTask.LogHandler += _onLog;
        MetaMorpheusTask.WarnHandler += _onWarn;
        MetaMorpheusTask.OutLabelStatusHandler += _onStatus;
        MetaMorpheusTask.OutProgressHandler += _onProgress;
        MetaMorpheusTask.StartingSingleTaskHander += _onTaskStarting;
        MetaMorpheusTask.FinishedSingleTaskHandler += _onTaskFinished;
        EverythingRunnerEngine.WarnHandler += _onRunnerWarn;
        EverythingRunnerEngine.FinishedAllTasksEngineHandler += _onAllFinished;
    }

    /// <summary>
    /// Runs an update on the UI thread. Falls back to running inline when there is no dispatcher,
    /// which is the case for tests that exercise the engines without an Avalonia application.
    /// </summary>
    private static void OnUiThread(Action update)
    {
        if (Dispatcher.UIThread is null || Dispatcher.UIThread.CheckAccess())
        {
            update();
            return;
        }
        Dispatcher.UIThread.Post(update);
    }

    /// <summary>Detaches from the static engine events. Not doing this leaks the view model.</summary>
    public void Dispose()
    {
        MetaMorpheusTask.LogHandler -= _onLog;
        MetaMorpheusTask.WarnHandler -= _onWarn;
        MetaMorpheusTask.OutLabelStatusHandler -= _onStatus;
        MetaMorpheusTask.OutProgressHandler -= _onProgress;
        MetaMorpheusTask.StartingSingleTaskHander -= _onTaskStarting;
        MetaMorpheusTask.FinishedSingleTaskHandler -= _onTaskFinished;
        EverythingRunnerEngine.WarnHandler -= _onRunnerWarn;
        EverythingRunnerEngine.FinishedAllTasksEngineHandler -= _onAllFinished;
    }

    /// <summary>Adds spectra paths chosen by the view, which owns the file dialog.</summary>
    public void AddSpectraFiles(IEnumerable<string> paths)
    {
        foreach (string path in paths.Where(IsSupportedSpectraFile))
        {
            if (SpectraFiles.All(f => f.FilePath != path))
            {
                SpectraFiles.Add(new SpectraFileForDisplay(path));
            }
        }
        if (string.IsNullOrEmpty(OutputFolder) && SpectraFiles.Count > 0)
        {
            OutputFolder = Path.Combine(Path.GetDirectoryName(SpectraFiles[0].FilePath) ?? ".", "$DATETIME");
        }
    }

    /// <summary>
    /// Thermo .raw reading works on Windows and Linux only, so on macOS it is rejected up front with
    /// a reason rather than failing part-way into a run.
    /// </summary>
    internal static bool IsSupportedSpectraFile(string path)
    {
        string extension = Path.GetExtension(path).ToLowerInvariant();
        return extension is ".mzml" or ".mgf" or ".raw" or ".d";
    }

    internal static bool IsThermoRawUnsupportedHere(string path) =>
        Path.GetExtension(path).Equals(".raw", StringComparison.OrdinalIgnoreCase)
        && OperatingSystem.IsMacOS();

    public void AddDatabases(IEnumerable<string> paths)
    {
        foreach (string path in paths)
        {
            if (Databases.All(d => d.FilePath != path))
            {
                Databases.Add(new DatabaseForDisplay(path));
            }
        }
    }

    [RelayCommand]
    private void AddSearchTask() => AddTask("Search", new SearchTask());

    [RelayCommand]
    private void AddCalibrationTask() => AddTask("Calibrate", new CalibrationTask());

    [RelayCommand]
    private void AddGptmdTask() => AddTask("GPTMD", new GptmdTask());

    private void AddTask(string kind, MetaMorpheusTask task)
    {
        Tasks.Add(new TaskForDisplay($"Task{Tasks.Count + 1}-{kind}", kind, task));
        Append($"Added {kind} task. Select it and choose Edit settings to change its parameters.");
    }

    /// <summary>Lets the view append to the log without reaching into private members.</summary>
    public void Note(string line) => Append(line);

    /// <summary>
    /// Builds the editable settings for a task the view has selected. Taking object keeps the view
    /// from needing to know about TaskForDisplay, which is internal.
    /// </summary>
    public TaskSettingsViewModel CreateSettingsFor(object selected) =>
        selected is TaskForDisplay display
            ? new TaskSettingsViewModel(display.Task, display.Kind)
            : null;

    [RelayCommand]
    private void ClearTasks()
    {
        Tasks.Clear();
        Progress = 0;
        Status = "Ready.";
    }

    public bool CanRun => !IsRunning && Tasks.Count > 0 && SpectraFiles.Count > 0 && Databases.Count > 0;

    [RelayCommand]
    private async Task RunAsync()
    {
        if (Tasks.Count == 0 || SpectraFiles.Count == 0 || Databases.Count == 0)
        {
            Append("Add at least one spectra file, one database and one task before running.");
            return;
        }

        var unsupported = SpectraFiles.Where(f => IsThermoRawUnsupportedHere(f.FilePath)).ToList();
        if (unsupported.Count > 0)
        {
            Append($"Thermo .raw reading is not available on macOS. Convert to mzML first: "
                + string.Join(", ", unsupported.Select(f => f.FileName)));
            return;
        }

        IsRunning = true;
        Progress = 0;
        try
        {
            var runner = new EverythingRunnerEngine(
                Tasks.Select(t => (t.DisplayName, t.Task)).ToList(),
                SpectraFiles.Select(f => f.FilePath).ToList(),
                Databases.Select(d => new DbForTask(d.FilePath, d.IsContaminant)).ToList(),
                OutputFolder);

            Append($"Running {Tasks.Count} task(s) on {SpectraFiles.Count} file(s)...");
            await Task.Run(runner.Run);
        }
        catch (Exception exception)
        {
            Append($"ERROR: {exception.GetType().Name}: {exception.Message}");
        }
        finally
        {
            IsRunning = false;
            Status = "Ready.";
        }
    }

    [RelayCommand]
    private void Cancel()
    {
        // The engines poll this cooperatively; there is no forceful abort.
        GlobalVariables.StopLoops = true;
        Append("Cancellation requested; tasks stop at their next checkpoint.");
    }

    private void Append(string line) => Log += line + Environment.NewLine;
}

internal sealed class SpectraFileForDisplay
{
    public SpectraFileForDisplay(string filePath)
    {
        FilePath = filePath;
        FileName = Path.GetFileName(filePath);
    }

    public string FileName { get; }
    public string FilePath { get; }
}

internal sealed class DatabaseForDisplay
{
    public DatabaseForDisplay(string filePath)
    {
        FilePath = filePath;
        FileName = Path.GetFileName(filePath);
        // Matches the WPF window's convention for flagging contaminant databases.
        IsContaminant = filePath.Contains("contaminant", StringComparison.OrdinalIgnoreCase)
            || filePath.Contains("crap", StringComparison.OrdinalIgnoreCase);
    }

    public string FileName { get; }
    public string FilePath { get; }
    public bool IsContaminant { get; set; }
}

internal sealed class TaskForDisplay
{
    public TaskForDisplay(string displayName, string kind, MetaMorpheusTask task)
    {
        DisplayName = displayName;
        Kind = kind;
        Task = task;
    }

    public string DisplayName { get; }
    public string Kind { get; }
    public MetaMorpheusTask Task { get; }
}
