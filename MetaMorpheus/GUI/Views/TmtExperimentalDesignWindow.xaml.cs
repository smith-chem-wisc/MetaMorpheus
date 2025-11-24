using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.ComponentModel;
using System.IO;
using System.Linq;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Input;
using System.Windows.Media;

namespace MetaMorpheusGUI
{
    public partial class TmtExperimentalDesignWindow : Window
    {
        private readonly ObservableCollection<TmtDesignRow> _rows = new();
        private readonly ObservableCollection<RawDataForDataGrid> _spectraFiles;
        public record TmtDesignResult(string FilePath, int Fraction, int TechnicalReplicate, string Plex);
        private Point _dragStart;
        private TmtDesignRow _dragSource;
        private readonly Dictionary<string, List<PlexAnnotation>> _plexAnnotations = new(StringComparer.OrdinalIgnoreCase);

        // Persist design across dialog opens during app lifetime
        private static readonly Dictionary<string, string> s_fileToPlex = new(StringComparer.OrdinalIgnoreCase);
        private static readonly Dictionary<string, List<PlexAnnotation>> s_plexAnnotations = new(StringComparer.OrdinalIgnoreCase);

        // NEW: Per-file full state (fraction, tech rep, plex)
        private record FileDesignState(int Fraction, int TechnicalReplicate, string Plex);
        private static readonly Dictionary<string, FileDesignState> s_fileState = new(StringComparer.OrdinalIgnoreCase);

        public TmtExperimentalDesignWindow(ObservableCollection<RawDataForDataGrid> spectraFiles)
        {
            InitializeComponent();
            _spectraFiles = spectraFiles;
            DgTmt.ItemsSource = _rows;

            // seed rows with existing SpectraFiles.Use == true
            foreach (var file in _spectraFiles.Where(p => p.Use).Select(p => p.FilePath))
            {
                var row = new TmtDesignRow(file);
                // prefer full saved state if available
                if (s_fileState.TryGetValue(file, out var st))
                {
                    row.Fraction = st.Fraction;
                    row.TechnicalReplicate = st.TechnicalReplicate;
                    row.Plex = st.Plex ?? "";
                }
                else if (s_fileToPlex.TryGetValue(file, out var savedPlex)) // legacy support if only Plex was saved
                {
                    row.Plex = savedPlex;
                }
                _rows.Add(row);
            }

            // seed annotations cache for this window instance
            foreach (var kv in s_plexAnnotations)
                _plexAnnotations[kv.Key] = new List<PlexAnnotation>(kv.Value);
        }

        public IReadOnlyList<TmtDesignResult> GetResults()
        {
            return _rows
                .Select(r => new TmtDesignResult(r.FilePath, r.Fraction, r.TechnicalReplicate, r.Plex ?? string.Empty))
                .ToList();
        }

        #region Drag/drop add files
        private void Window_Drop(object sender, DragEventArgs e)
        {
            if (!e.Data.GetDataPresent(DataFormats.FileDrop))
                return;

            CommitPendingEdits();

            var files = (string[])e.Data.GetData(DataFormats.FileDrop);
            if (files == null) return;

            foreach (var path in files)
            {
                if (Directory.Exists(path))
                {
                    foreach (var f in Directory.EnumerateFiles(path))
                        AddFileIfNotExists(f);
                }
                else if (File.Exists(path))
                    AddFileIfNotExists(path);
            }

            DgTmt.Items.Refresh();
        }

        private void AddFileIfNotExists(string path)
        {
            if (_rows.Any(r => string.Equals(r.FilePath, path, StringComparison.OrdinalIgnoreCase)))
                return;

            var row = new TmtDesignRow(path);
            if (s_fileState.TryGetValue(path, out var st))
            {
                row.Fraction = st.Fraction;
                row.TechnicalReplicate = st.TechnicalReplicate;
                row.Plex = st.Plex ?? "";
            }
            else
            {
                // initialize new: fraction increments from current max; techRep = 1; plex empty
                int nextFraction = _rows.Any() ? _rows.Max(r => r.Fraction) : 0; // first becomes 1
                row.Fraction = nextFraction + 1;
                row.TechnicalReplicate = 1;
            }
            _rows.Add(row);
        }
        #endregion

        #region Validation
        private string ValidateDesign()
        {
            if (!_rows.Any())
                return "No files defined.";

            CommitPendingEdits();

            var distinctFractions = _rows.Select(r => r.Fraction).Distinct().OrderBy(i => i).ToList();
            if (distinctFractions.First() != 1)
                return "Fractions must start at 1.";
            int maxFraction = distinctFractions.Last();
            for (int i = 1; i <= maxFraction; i++)
                if (!distinctFractions.Contains(i))
                    return $"Missing fraction number {i} in distinct set.";

            if (_rows.Any(r => r.TechnicalReplicate < 1))
                return "Technical Replicate values must be >= 1.";

            foreach (var grp in _rows.GroupBy(r => r.Fraction))
            {
                var techs = grp.Select(r => r.TechnicalReplicate).Distinct().OrderBy(t => t).ToList();
                if (techs.First() != 1)
                    return $"Fraction {grp.Key}: technical replicates must start at 1.";
                int maxTech = techs.Last();
                for (int t = 1; t <= maxTech; t++)
                    if (!techs.Contains(t))
                        return $"Fraction {grp.Key}: missing technical replicate {t}.";
            }

            var duplicatePair = _rows.GroupBy(r => (r.Fraction, r.TechnicalReplicate))
                                     .FirstOrDefault(g => g.Count() > 1);
            if (duplicatePair != null)
                return $"Duplicate Fraction/Technical Replicate combination: Fraction {duplicatePair.Key.Fraction}, Technical Replicate {duplicatePair.Key.TechnicalReplicate}.";

            return null;
        }
        #endregion

        #region Row drag reorder
        private void Row_PreviewMouseLeftButtonDown(object sender, MouseButtonEventArgs e)
        {
            _dragStart = e.GetPosition(null);
            if (e.OriginalSource is DependencyObject dep)
            {
                var row = FindAncestor<DataGridRow>(dep);
                _dragSource = row?.Item as TmtDesignRow;
            }
        }

        private void Row_MouseMove(object sender, MouseEventArgs e)
        {
            if (e.LeftButton != MouseButtonState.Pressed || _dragSource == null)
                return;

            var diff = e.GetPosition(null) - _dragStart;
            if (Math.Abs(diff.X) < SystemParameters.MinimumHorizontalDragDistance &&
                Math.Abs(diff.Y) < SystemParameters.MinimumVerticalDragDistance)
                return;

            DragDrop.DoDragDrop(DgTmt, _dragSource, DragDropEffects.Move);
        }

        private void Row_DragOver(object sender, DragEventArgs e)
        {
            if (!e.Data.GetDataPresent(typeof(TmtDesignRow)))
            {
                e.Effects = DragDropEffects.None;
                e.Handled = true;
                return;
            }
            e.Effects = DragDropEffects.Move;
            e.Handled = true;
        }

        private void Row_Drop(object sender, DragEventArgs e)
        {
            if (!e.Data.GetDataPresent(typeof(TmtDesignRow)))
                return;

            var source = (TmtDesignRow)e.Data.GetData(typeof(TmtDesignRow));
            if (source == null) return;

            var target = ResolveTargetRow(e.OriginalSource);
            if (target == null || target == source) return;

            int sourceIndex = _rows.IndexOf(source);
            int targetIndex = _rows.IndexOf(target);
            if (sourceIndex < 0 || targetIndex < 0) return;

            _rows.RemoveAt(sourceIndex);
            _rows.Insert(targetIndex, source);

            DgTmt.SelectedItem = source;
            DgTmt.Items.Refresh();
        }

        private TmtDesignRow ResolveTargetRow(object originalSource)
        {
            if (originalSource is DependencyObject dep)
            {
                var row = FindAncestor<DataGridRow>(dep);
                return row?.Item as TmtDesignRow;
            }
            return null;
        }

        private static T FindAncestor<T>(DependencyObject current) where T : DependencyObject
        {
            while (current != null)
            {
                if (current is T t) return t;
                current = VisualTreeHelper.GetParent(current);
            }
            return null;
        }
        #endregion

        #region Editing
        private void DgTmt_CellEditEnding(object sender, DataGridCellEditEndingEventArgs e)
        {
            if (e.EditAction != DataGridEditAction.Commit)
                return;

            if (e.Row.Item is TmtDesignRow row)
            {
                var colHeader = e.Column.Header?.ToString();
                var tb = e.EditingElement as TextBox;
                if (tb == null) return;

                if (colHeader == "Fraction")
                {
                    if (!int.TryParse(tb.Text, out var val) || val < 1)
                    {
                        MessageBox.Show("Fraction must be integer >= 1");
                        e.Cancel = true;
                        return;
                    }
                    row.Fraction = val;
                }
                else if (colHeader == "Technical Replicate")
                {
                    if (!int.TryParse(tb.Text, out var val) || val < 1)
                    {
                        MessageBox.Show("Technical Replicate must be integer >= 1");
                        e.Cancel = true;
                        return;
                    }
                    row.TechnicalReplicate = val;
                }
            }
        }
        #endregion

        #region Annotate Plex
        private void AnnotatePlexButton_Click(object sender, RoutedEventArgs e)
        {
            CommitPendingEdits();

            var plexNames = _rows
                .Select(r => (r.Plex ?? string.Empty).Trim())
                .Where(s => !string.IsNullOrEmpty(s))
                .Distinct(StringComparer.OrdinalIgnoreCase)
                .OrderBy(s => s)
                .ToList();

            if (!plexNames.Any())
            {
                MessageBox.Show("No Plex values have been set.");
                return;
            }

            var plexFileMap = _rows
                .Where(r => !string.IsNullOrWhiteSpace(r.Plex))
                .GroupBy(r => r.Plex.Trim())
                .ToDictionary(g => g.Key,
                    g => g
                        .OrderBy(x => x.Fraction)
                        .ThenBy(x => x.TechnicalReplicate)
                        .Select(x =>
                        {
                            _plexAnnotations.TryGetValue(x.Plex.Trim(), out var anns);
                            return new PlexFileEntry(x.FilePath, x.Fraction, x.Plex.Trim(), x.TechnicalReplicate, anns);
                        })
                        .ToList());

            var dialog = new AnnotatePlexWindow(plexNames,
                _plexAnnotations.ToDictionary(k => k.Key, v => v.Value.ToList()),
                plexFileMap)
            {
                Owner = this
            };

            if (dialog.ShowDialog() == true &&
                !string.IsNullOrEmpty(dialog.SavedPlexName) &&
                dialog.SavedAnnotations != null)
            {
                _plexAnnotations[dialog.SavedPlexName] = dialog.SavedAnnotations;
                s_plexAnnotations[dialog.SavedPlexName] = new List<PlexAnnotation>(dialog.SavedAnnotations);
            }
        }
        #endregion

        #region Buttons
        private void ValidateButton_Click(object sender, RoutedEventArgs e)
        {
            CommitPendingEdits();
            var err = ValidateDesign();
            MessageBox.Show(err ?? "Validation passed.",
                err == null ? "OK" : "Validation Error");
        }

        private void SaveButton_Click(object sender, RoutedEventArgs e)
        {
            CommitPendingEdits();
            var err = ValidateDesign();
            if (err != null)
            {
                MessageBox.Show(err, "Cannot Save");
                return;
            }

            // persist per-file full state and plex choices for next open
            foreach (var r in _rows)
            {
                s_fileState[r.FilePath] = new FileDesignState(r.Fraction, r.TechnicalReplicate, r.Plex?.Trim() ?? "");
                if (!string.IsNullOrWhiteSpace(r.Plex))
                    s_fileToPlex[r.FilePath] = r.Plex.Trim();
            }

            DialogResult = true;
            Close();
        }

        private void CancelButton_Click(object sender, RoutedEventArgs e)
        {
            DialogResult = false;
            Close();
        }
        #endregion

        #region Row model
        private class TmtDesignRow : INotifyPropertyChanged
        {
            private int _fraction;
            private int _technicalReplicate;
            private string _plex = "";

            public TmtDesignRow(string filePath)
            {
                FilePath = filePath;
                _fraction = 1;
                _technicalReplicate = 1;
            }

            public string FilePath { get; }

            public int Fraction
            {
                get => _fraction;
                set
                {
                    if (_fraction == value) return;
                    _fraction = value;
                    OnPropertyChanged(nameof(Fraction));
                }
            }

            public int TechnicalReplicate
            {
                get => _technicalReplicate;
                set
                {
                    if (_technicalReplicate == value) return;
                    _technicalReplicate = value;
                    OnPropertyChanged(nameof(TechnicalReplicate));
                }
            }

            public string Plex
            {
                get => _plex;
                set
                {
                    if (_plex == value) return;
                    _plex = value;
                    OnPropertyChanged(nameof(Plex));
                }
            }

            public event PropertyChangedEventHandler PropertyChanged;
            private void OnPropertyChanged(string name) =>
                PropertyChanged?.Invoke(this, new PropertyChangedEventArgs(name));
        }
        #endregion

        #region Helpers
        private void CommitPendingEdits()
        {
            DgTmt.CommitEdit(DataGridEditingUnit.Cell, true);
            DgTmt.CommitEdit(DataGridEditingUnit.Row, true);
        }

        // Load TMT design .txt files in same folders as the provided raw files and seed caches
        public static void SeedFromDesignFiles(IEnumerable<string> rawFilePaths)
        {
            if (rawFilePaths == null) return;

            var rawSet = new HashSet<string>(
                rawFilePaths.Where(p => !string.IsNullOrWhiteSpace(p))
                            .Select(p => Path.GetFullPath(p)),
                StringComparer.OrdinalIgnoreCase);

            var dirs = rawSet.Select(Path.GetDirectoryName)
                             .Where(d => !string.IsNullOrEmpty(d))
                             .Distinct(StringComparer.OrdinalIgnoreCase);

            foreach (var dir in dirs)
            {
                var designPath = Path.Combine(dir, "TmtDesign.txt");
                if (!File.Exists(designPath))
                    continue;

                TryParseDesignFile(designPath, rawSet);
            }
        }

        private static void TryParseDesignFile(string path, HashSet<string> interestedFiles)
        {
            try
            {
                using var sr = new StreamReader(path);
                var header = sr.ReadLine();
                if (header == null)
                    return;

                var headers = header.Split('\t');
                int idxFile = Array.FindIndex(headers, h => string.Equals(h.Trim(), "File", StringComparison.OrdinalIgnoreCase));
                int idxPlex = Array.FindIndex(headers, h => string.Equals(h.Trim(), "Plex", StringComparison.OrdinalIgnoreCase));
                int idxSample = Array.FindIndex(headers, h => string.Equals(h.Trim(), "Sample Name", StringComparison.OrdinalIgnoreCase));
                int idxChannel = Array.FindIndex(headers, h => string.Equals(h.Trim(), "TMT Channel", StringComparison.OrdinalIgnoreCase));
                int idxCond = Array.FindIndex(headers, h => string.Equals(h.Trim(), "Condition", StringComparison.OrdinalIgnoreCase));
                int idxBio = Array.FindIndex(headers, h => string.Equals(h.Trim(), "Biological Replicate", StringComparison.OrdinalIgnoreCase));
                int idxFrac = Array.FindIndex(headers, h => string.Equals(h.Trim(), "Fraction", StringComparison.OrdinalIgnoreCase));
                int idxTech = Array.FindIndex(headers, h => string.Equals(h.Trim(), "Technical Replicate", StringComparison.OrdinalIgnoreCase));

                if (new[] { idxFile, idxPlex, idxSample, idxChannel, idxCond, idxBio, idxFrac, idxTech }.Any(i => i < 0))
                    return; // not a valid design file

                var fileState = new Dictionary<string, (int frac, int tech, string plex)>(StringComparer.OrdinalIgnoreCase);
                var plexToAnnotations = new Dictionary<string, Dictionary<string, PlexAnnotation>>(StringComparer.OrdinalIgnoreCase);
                var plexOrder = new Dictionary<string, List<string>>(StringComparer.OrdinalIgnoreCase);

                string? line;
                while ((line = sr.ReadLine()) != null)
                {
                    if (string.IsNullOrWhiteSpace(line)) continue;
                    var cols = line.Split('\t');
                    if (cols.Length < headers.Length) continue;

                    var file = cols[idxFile].Trim();
                    if (string.IsNullOrEmpty(file)) continue;

                    var full = Path.GetFullPath(file);
                    if (!interestedFiles.Contains(full)) continue;

                    var plex = cols[idxPlex].Trim();
                    var sample = cols[idxSample].Trim();
                    var chan = cols[idxChannel].Trim();
                    var cond = cols[idxCond].Trim();
                    int.TryParse(cols[idxBio].Trim(), out var bio);
                    int.TryParse(cols[idxFrac].Trim(), out var frac);
                    int.TryParse(cols[idxTech].Trim(), out var tech);

                    fileState[full] = (frac, tech, plex);

                    if (!plexToAnnotations.TryGetValue(plex, out var byChannel))
                    {
                        byChannel = new Dictionary<string, PlexAnnotation>(StringComparer.OrdinalIgnoreCase);
                        plexToAnnotations[plex] = byChannel;
                        plexOrder[plex] = new List<string>();
                    }

                    if (!byChannel.ContainsKey(chan))
                    {
                        byChannel[chan] = new PlexAnnotation
                        {
                            Tag = chan,
                            SampleName = sample,
                            Condition = cond,
                            BiologicalReplicate = bio
                        };
                        plexOrder[plex].Add(chan);
                    }
                }

                // Persist parsed state to caches used by the TMT window
                foreach (var kv in fileState)
                {
                    var (frac, tech, plex) = kv.Value;
                    s_fileState[kv.Key] = new FileDesignState(frac, tech, plex ?? "");
                    if (!string.IsNullOrWhiteSpace(plex))
                        s_fileToPlex[kv.Key] = plex;
                }

                foreach (var kv in plexToAnnotations)
                {
                    var plex = kv.Key;
                    var order = plexOrder.TryGetValue(plex, out var o) ? o : kv.Value.Keys.ToList();
                    var list = order.Where(ch => kv.Value.ContainsKey(ch)).Select(ch => kv.Value[ch]).ToList();
                    s_plexAnnotations[plex] = list;
                }
            }
            catch
            {
                // ignore malformed files
            }
        }
        #endregion
    }
}