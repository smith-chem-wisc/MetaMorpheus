using System;
using System.Collections.ObjectModel;
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

        private Point _dragStart;
        private TmtDesignRow _dragSource;

        private readonly System.Collections.Generic.Dictionary<string, System.Collections.Generic.List<PlexAnnotation>> _plexAnnotations
            = new();

        public TmtExperimentalDesignWindow(ObservableCollection<RawDataForDataGrid> spectraFiles)
        {
            InitializeComponent();
            _spectraFiles = spectraFiles;
            DgTmt.ItemsSource = _rows;

            foreach (var file in _spectraFiles.Where(p => p.Use).Select(p => p.FilePath))
            {
                AddFileIfNotExists(file);
            }
        }

        #region Drag/drop add files
        private void Window_Drop(object sender, DragEventArgs e)
        {
            if (!e.Data.GetDataPresent(DataFormats.FileDrop))
                return;

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

            int nextFraction = _rows.Any() ? _rows.Max(r => r.Fraction) : 0;
            // Start at 1; user can later duplicate or edit.
            _rows.Add(new TmtDesignRow(path) { Fraction = nextFraction + 1 });
        }
        #endregion

        #region Validation
        private string ValidateDesign()
        {
            if (!_rows.Any())
                return "No files defined.";

            // Fractions must start at 1 (if any rows exist) and cover 1..MaxFraction with no gaps when considering the DISTINCT set.
            var distinctFractions = _rows.Select(r => r.Fraction).Distinct().OrderBy(i => i).ToList();
            if (distinctFractions.First() < 1)
                return "Fraction numbers must be >= 1.";
            int maxFraction = distinctFractions.Last();
            for (int i = 1; i <= maxFraction; i++)
                if (!distinctFractions.Contains(i))
                    return $"Missing fraction number {i} in distinct set.";

            // Technical replicates must be >=1
            if (_rows.Any(r => r.TechnicalReplicate < 1))
                return "Technical Replicate values must be >= 1.";

            // For each fraction: technical replicate set must start at 1 and have no gaps up to max replicate for that fraction
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

            // No duplicate (Fraction, TechnicalReplicate) pairs
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

            // DO NOT renumber fractions on reorder (user-managed)
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
                if (colHeader == "Fraction")
                {
                    if (!int.TryParse(((TextBox)e.EditingElement).Text, out var val) || val < 1)
                    {
                        MessageBox.Show("Fraction must be integer >= 1");
                        e.Cancel = true;
                        return;
                    }
                    row.Fraction = val;
                }
                else if (colHeader == "Technical Replicate")
                {
                    if (!int.TryParse(((TextBox)e.EditingElement).Text, out var val) || val < 1)
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
                    g => g.OrderBy(x => x.Fraction)
                          .Select(x => new PlexFileEntry(x.FilePath, x.Fraction, x.Plex.Trim(), x.TechnicalReplicate))
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
            }
        }
        #endregion

        #region Buttons
        private void ValidateButton_Click(object sender, RoutedEventArgs e)
        {
            var err = ValidateDesign();
            MessageBox.Show(err ?? "Validation passed.",
                err == null ? "OK" : "Validation Error");
        }

        private void SaveButton_Click(object sender, RoutedEventArgs e)
        {
            var err = ValidateDesign();
            if (err != null)
            {
                MessageBox.Show(err, "Cannot Save");
                return;
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
        private class TmtDesignRow
        {
            public TmtDesignRow(string filePath)
            {
                FilePath = filePath;
                Fraction = 1;           // set when added
                TechnicalReplicate = 1; // default
            }

            public string FilePath { get; }
            public int Fraction { get; set; }
            public int TechnicalReplicate { get; set; }
            public string Plex { get; set; } = "";
        }
        #endregion
    }
}