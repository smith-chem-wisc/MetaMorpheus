using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading;
using System.Threading.Tasks;
using CommunityToolkit.Mvvm.ComponentModel;
using CommunityToolkit.Mvvm.Input;
using EngineLayer;
using EngineLayer.DatabaseLoading;
using TaskLayer;

namespace MetaMorpheusAvalonia.ViewModels;

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
    private EventHandler<StringEventArgs> _onEngineStatus;
    private EventHandler<ProgressEventArgs> _onEngineProgress;
    private EventHandler<StringEventArgs> _onEngineWarn;
    private EventHandler<StringEventArgs> _onFileManagerWarn;
    private EventHandler _onAllStarting;
    private EventHandler<StringEventArgs> _onAllResultsWritten;
    private EventHandler<XmlForTaskListEventArgs> _onNewDbs;
    private EventHandler<StringListEventArgs> _onNewSpectra;

    // Whatever thread constructed this view model. See OnUiThread.
    private readonly SynchronizationContext _uiContext = SynchronizationContext.Current;

    public ObservableCollection<SpectraFileForDisplay> SpectraFiles { get; } = new();
    public ObservableCollection<DatabaseForDisplay> Databases { get; } = new();
    public ObservableCollection<TaskForDisplay> Tasks { get; } = new();

    [ObservableProperty] private string _outputFolder = string.Empty;
    [ObservableProperty] private string _status = "Ready.";
    [ObservableProperty] private double _progress;
    [ObservableProperty] private bool _isRunning;

    public string Version => GlobalVariables.MetaMorpheusVersion;

    // Set when the data directory could not be read. Nothing can run, but the window still opens so
    // the user is told why rather than watching the process vanish.
    private readonly bool _startupFailed;

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

    private MainWindowViewModel(Exception startupFailure)
    {
        _startupFailed = true;
        Append("MetaMorpheus could not read its data directory, so no search can run.");
        Append($"ERROR: {startupFailure.GetType().Name}: {startupFailure.Message}");
        Append($"Expected the Mods, Data and Contaminants folders under {AppDomain.CurrentDomain.BaseDirectory}");
        Status = "Failed to start.";
    }

    /// <summary>A view model that only reports why start-up failed. See App.OnFrameworkInitializationCompleted.</summary>
    internal static MainWindowViewModel ThatFailedToStart(Exception failure) => new(failure);

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

        // The inner engines report separately from the tasks that own them. Without these three, the
        // progress bar sat still and engine warnings never reached the log for the whole of a search.
        _onEngineStatus = (_, e) => OnUiThread(() => Status = e.S);
        _onEngineProgress = (_, e) => OnUiThread(() => Progress = e.NewProgress);
        _onEngineWarn = (_, e) => OnUiThread(() => Append($"WARNING: {e.S}"));

        MetaMorpheusTask.LogHandler += _onLog;
        MetaMorpheusTask.WarnHandler += _onWarn;
        MetaMorpheusTask.OutLabelStatusHandler += _onStatus;
        MetaMorpheusTask.OutProgressHandler += _onProgress;
        MetaMorpheusTask.StartingSingleTaskHander += _onTaskStarting;
        MetaMorpheusTask.FinishedSingleTaskHandler += _onTaskFinished;
        EverythingRunnerEngine.WarnHandler += _onRunnerWarn;
        EverythingRunnerEngine.FinishedAllTasksEngineHandler += _onAllFinished;
        // MyFileManager warns about the spectra files themselves - unreadable or unexpected content -
        // which is exactly the class of problem a new platform surfaces first.
        _onFileManagerWarn = (_, e) => OnUiThread(() => Append($"WARNING: {e.S}"));

        // EverythingRunnerEngine feeds each task's output into the next, so a Calibrate -> GPTMD ->
        // Search run silently swaps the database and spectra lists partway through. Reporting that is
        // how a user can tell the chaining happened.
        _onAllStarting = (_, _) => OnUiThread(() => Append("Starting all tasks..."));
        _onAllResultsWritten = (_, e) => OnUiThread(() => Append($"Wrote results: {e.S}"));
        _onNewDbs = (_, e) => OnUiThread(() => Append(
            $"Later tasks will use {e.NewDatabases.Count} database(s) written by the task that just finished."));
        _onNewSpectra = (_, e) => OnUiThread(() => Append(
            $"Later tasks will use {e.StringList.Count()} calibrated spectra file(s)."));

        MetaMorpheusEngine.OutLabelStatusHandler += _onEngineStatus;
        MetaMorpheusEngine.OutProgressHandler += _onEngineProgress;
        MetaMorpheusEngine.WarnHandler += _onEngineWarn;
        MyFileManager.WarnHandler += _onFileManagerWarn;
        EverythingRunnerEngine.StartingAllTasksEngineHandler += _onAllStarting;
        EverythingRunnerEngine.FinishedWritingAllResultsFileHandler += _onAllResultsWritten;
        EverythingRunnerEngine.NewDbsHandler += _onNewDbs;
        EverythingRunnerEngine.NewSpectrasHandler += _onNewSpectra;
    }

    /// <summary>
    /// Gate C. The static event surface a front end has to cover, split into what this window handles
    /// and what it deliberately does not, so a reflection test can require every event to appear in
    /// one list or the other. A new engine event then fails the build rather than going unnoticed.
    /// </summary>
    internal static IReadOnlyList<string> SubscribedEventNames { get; } = new[]
    {
        "MetaMorpheusTask.LogHandler",
        "MetaMorpheusTask.WarnHandler",
        "MetaMorpheusTask.OutLabelStatusHandler",
        "MetaMorpheusTask.OutProgressHandler",
        "MetaMorpheusTask.StartingSingleTaskHander",
        "MetaMorpheusTask.FinishedSingleTaskHandler",
        "EverythingRunnerEngine.WarnHandler",
        "EverythingRunnerEngine.FinishedAllTasksEngineHandler",
        "EverythingRunnerEngine.StartingAllTasksEngineHandler",
        "EverythingRunnerEngine.FinishedWritingAllResultsFileHandler",
        "EverythingRunnerEngine.NewDbsHandler",
        "EverythingRunnerEngine.NewSpectrasHandler",
        "MetaMorpheusEngine.OutLabelStatusHandler",
        "MetaMorpheusEngine.OutProgressHandler",
        "MetaMorpheusEngine.WarnHandler",
        "MyFileManager.WarnHandler",
    };

    /// <summary>
    /// Events left alone on purpose, each because it exists to populate a piece of WPF-specific UI this
    /// window does not have. Every one of these is informational; none carries a warning or a result.
    /// </summary>
    internal static IReadOnlyDictionary<string, string> DeliberatelyUnsubscribedEvents { get; } =
        new Dictionary<string, string>
        {
            ["MetaMorpheusTask.NewCollectionHandler"] = "adds a branch to the WPF results tree; no such tree here",
            ["MetaMorpheusTask.FinishedWritingFileHandler"] = "populates the WPF results-file tree",
            ["MetaMorpheusTask.StartingDataFileHandler"] = "per-file progress rows; this window has no per-file grid",
            ["MetaMorpheusTask.FinishedDataFileHandler"] = "per-file progress rows; this window has no per-file grid",
            ["EverythingRunnerEngine.NewFileSpecificTomlHandler"] = "file-specific parameters are not implemented yet",
            ["MetaMorpheusEngine.StartingSingleEngineHander"] = "per-engine start; OutLabelStatusHandler already reports it, and the WPF window ignores this too",
            ["MetaMorpheusEngine.FinishedSingleEngineHandler"] = "carries engine results for the WPF tree; the WPF window ignores this too",
        };

    /// <summary>
    /// Runs an update on the thread that owns this view model, or inline when nothing owns it.
    ///
    /// This has to be right before any engine event is subscribed, because engines raise theirs from
    /// inside Parallel.For - arbitrary thread-pool threads, not one predictable one. The earlier
    /// version tested Dispatcher.UIThread for null, which is never null in Avalonia 12, and then
    /// Post()ed to it. With no Avalonia platform registered, Dispatcher.UIThread is built over
    /// ManagedDispatcherImpl and binds to whichever thread first touches it; a Post from any other
    /// thread then queues onto a dispatcher nobody pumps and is discarded silently. Log lines went
    /// missing that way, and it would have looked like the events were never raised.
    ///
    /// Capturing the SynchronizationContext at construction settles it: Avalonia installs its own on
    /// the UI thread, so this marshals properly in the running app and under [AvaloniaTest], and is
    /// null in tests that drive the engines with no application, where inline is correct.
    /// </summary>
    private void OnUiThread(Action update)
    {
        if (_uiContext is null || ReferenceEquals(_uiContext, SynchronizationContext.Current))
        {
            update();
            return;
        }
        _uiContext.Post(_ => update(), null);
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
        MetaMorpheusEngine.OutLabelStatusHandler -= _onEngineStatus;
        MetaMorpheusEngine.OutProgressHandler -= _onEngineProgress;
        MetaMorpheusEngine.WarnHandler -= _onEngineWarn;
        MyFileManager.WarnHandler -= _onFileManagerWarn;
        EverythingRunnerEngine.StartingAllTasksEngineHandler -= _onAllStarting;
        EverythingRunnerEngine.FinishedWritingAllResultsFileHandler -= _onAllResultsWritten;
        EverythingRunnerEngine.NewDbsHandler -= _onNewDbs;
        EverythingRunnerEngine.NewSpectrasHandler -= _onNewSpectra;
    }

    /// <summary>Adds spectra paths chosen by the view, which owns the file dialog.</summary>
    // Lives on the view model, not on MainWindow, on purpose. MainWindow carries
    // [ExcludeFromCodeCoverage] because a file picker cannot be unit tested, and these two
    // properties are the testable part of that feature - the same reason GuiFunctions exists
    // separately from the WPF GUI project. Keeping them here means the two platform invariants
    // below are pinned by a test instead of only by a comment.
    /// <summary>
    /// Built from GlobalVariables rather than written out, so the picker cannot drift from what the
    /// application accepts. Two things this gets right that a hand-written list did not:
    ///
    ///   * The entries are lowercase. Avalonia's FreeDesktop backend passes these to the XDG portal as
    ///     GlobStyle globs, and portal matching is case-sensitive, so "*.mzML" hides sample.mzml on
    ///     Linux - the platform this exists for. Windows and macOS match case-insensitively, so it
    ///     would never reproduce locally on either.
    ///   * ".d" is excluded. Bruker data is a directory and OpenFilePickerAsync cannot select one;
    ///     users pick the .tdf inside it and AddSpectraFiles maps that back to the folder.
    /// </summary>
    public static string[] SpectraPatterns => GlobalVariables.AcceptedSpectraFormats
        .Where(extension => extension != ".d")
        .Select(extension => "*" + extension)
        .ToArray();

    public static string[] DatabasePatterns => GlobalVariables.AcceptedDatabaseFormats
        .SelectMany(extension => new[] { "*" + extension, "*" + extension + ".gz" })
        .ToArray();

    public void AddSpectraFiles(IEnumerable<string> paths)
    {
        foreach (string path in paths.Where(IsSupportedSpectraFile).Select(ToBrukerFolderIfInside))
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
    /// GlobalVariables owns the canonical list, so a format added there is accepted here without this
    /// having to be updated. GetFileExtension rather than Path.GetExtension, for .xml.gz and friends.
    /// </summary>
    internal static bool IsSupportedSpectraFile(string path) =>
        GlobalVariables.AcceptedSpectraFormats.Contains(
            GlobalVariables.GetFileExtension(path).ToLowerInvariant());

    /// <summary>
    /// Bruker data is a .d folder. Selecting a file inside one has to be read as selecting the folder,
    /// which is what CMD and the WPF windows both do; otherwise the choice is silently dropped.
    /// </summary>
    internal static string ToBrukerFolderIfInside(string path)
    {
        string extension = GlobalVariables.GetFileExtension(path).ToLowerInvariant();
        if (extension is not (".tdf" or ".tdf_bin"))
        {
            return path;
        }

        string parent = Path.GetDirectoryName(path);
        return parent is not null && Path.GetExtension(parent).Equals(".d", StringComparison.OrdinalIgnoreCase)
            ? parent
            : path;
    }

    internal static bool IsThermoRawUnsupportedHere(string path) =>
        Path.GetExtension(path).Equals(".raw", StringComparison.OrdinalIgnoreCase)
        && OperatingSystem.IsMacOS();

    /// <summary>Databases are filtered against GlobalVariables too, so .msp and .msl libraries load.</summary>
    public void AddDatabases(IEnumerable<string> paths)
    {
        foreach (string path in paths.Where(IsSupportedDatabaseFile))
        {
            if (Databases.All(d => d.FilePath != path))
            {
                Databases.Add(new DatabaseForDisplay(path));
            }
        }
    }

    internal static bool IsSupportedDatabaseFile(string path) =>
        GlobalVariables.AcceptedDatabaseFormats.Contains(
            GlobalVariables.GetFileExtension(path).ToLowerInvariant());

    [RelayCommand]
    private void AddSearchTask() => AddTask("Search", new SearchTask());

    [RelayCommand]
    private void AddCalibrationTask() => AddTask("Calibrate", new CalibrationTask());

    [RelayCommand]
    private void AddGptmdTask() => AddTask("GPTMD", new GptmdTask());

    [RelayCommand]
    private void AddXlSearchTask() => AddTask("XLSearch", new XLSearchTask());

    [RelayCommand]
    private void AddGlycoSearchTask() => AddTask("GlycoSearch", new GlycoSearchTask());

    [RelayCommand]
    private void AddAveragingTask() => AddTask("Average", new SpectralAveragingTask());

    /// <summary>
    /// Loads the contaminant databases shipped alongside the application, as the WPF window's
    /// "+ADD DEFAULT CONTAMINANTS" button does.
    /// </summary>
    [RelayCommand]
    private void AddDefaultContaminants()
    {
        string folder = Path.Combine(GlobalVariables.DataDir, "Contaminants");
        if (!Directory.Exists(folder))
        {
            Append($"No contaminants folder at {folder}");
            return;
        }

        string[] found = Directory.GetFiles(folder).Where(IsSupportedDatabaseFile).ToArray();
        if (found.Length == 0)
        {
            Append($"No contaminant databases found in {folder}");
            return;
        }

        AddDatabases(found);
        Append($"Added {found.Length} contaminant database(s).");
    }

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

    public bool CanRun => !_startupFailed && !IsRunning
        && Tasks.Count > 0 && SpectraFiles.Count > 0 && Databases.Count > 0;

    [RelayCommand]
    private async Task RunAsync()
    {
        if (_startupFailed)
        {
            Append("Cannot run: MetaMorpheus did not finish starting up. See the errors above.");
            return;
        }

        // Clear a previous Cancel() before validating, as the WPF GUI does. Nothing else resets this,
        // so without it one cancellation ends every later run in the process.
        GlobalVariables.StopLoops = false;

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

    /// <summary>
    /// The run log. Held as a StringBuilder rather than a string that gets concatenated: a search
    /// emits thousands of lines, and `Log += line` reallocates the whole document every time, which is
    /// quadratic in the number of lines. The lock is because Append runs on engine worker threads
    /// whenever there is no UI thread to marshal to, making the old read-modify-write a genuine race.
    /// </summary>
    public string Log
    {
        get { lock (_logLines) { return _logLines.ToString(); } }
    }

    private readonly StringBuilder _logLines = new();

    private void Append(string line)
    {
        lock (_logLines)
        {
            _logLines.AppendLine(line);
        }
        OnPropertyChanged(nameof(Log));
    }
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

/// <summary>
/// ObservableObject rather than a plain class because the contaminant column is editable: the filename
/// heuristic is a guess and the user has to be able to overrule it, which needs a two-way binding and
/// therefore change notification.
/// </summary>
internal sealed partial class DatabaseForDisplay : ObservableObject
{
    [ObservableProperty] private bool _isContaminant;

    public DatabaseForDisplay(string filePath)
    {
        FilePath = filePath;
        FileName = Path.GetFileName(filePath);
        IsContaminant = DbForTask.LooksLikeContaminant(filePath);
    }

    public string FileName { get; }
    public string FilePath { get; }
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
