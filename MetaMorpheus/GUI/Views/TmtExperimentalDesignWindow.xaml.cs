using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.ComponentModel;
using System.IO;
using System.Linq;
using System.Windows;
using System.Windows.Controls;
using EngineLayer;

namespace MetaMorpheusGUI
{
    /// <summary>
    /// Authors TmtDesign.txt for the spectra files currently selected on the main window.
    /// </summary>
    /// <remarks>
    /// The shape follows <see cref="ExperimentalDesignWindow"/>, which owns ExperimentalDesign.tsv:
    /// the file on disk is the single source of truth, the dialog annotates a file set it does not
    /// own, and the whole design is written exactly once, on Save, only after it validates. This
    /// window used to diverge on all three counts -- it could add spectra files itself, it persisted
    /// into process-lifetime statics, and its Save wrote nothing at all while the annotate dialog
    /// wrote a design of its own -- which is why the grid and the file could disagree.
    /// </remarks>
    public partial class TmtExperimentalDesignWindow : Window
    {
        private readonly ObservableCollection<TmtDesignRow> _rows = new();
        private readonly ObservableCollection<RawDataForDataGrid> _spectraFiles;
        private readonly Dictionary<string, List<PlexAnnotation>> _plexAnnotations =
            new(StringComparer.OrdinalIgnoreCase);

        // Carried across dialog opens within one session so a part-authored design survives a
        // Cancel. The FILE is still the source of truth: SeedFromDesignFiles refreshes both of
        // these from disk whenever the spectra-file list changes.
        private record FileDesignState(int Fraction, int TechnicalReplicate, string Plex);
        private static readonly Dictionary<string, FileDesignState> s_fileState =
            new(StringComparer.OrdinalIgnoreCase);
        private static readonly Dictionary<string, List<PlexAnnotation>> s_plexAnnotations =
            new(StringComparer.OrdinalIgnoreCase);

        public TmtExperimentalDesignWindow(ObservableCollection<RawDataForDataGrid> spectraFiles)
        {
            InitializeComponent();
            _spectraFiles = spectraFiles;
            DgTmt.ItemsSource = _rows;

            if (s_seedErrors.Count > 0)
            {
                MessageBox.Show(
                    "A TMT design file was found next to the spectra files, but it could not be read in full. "
                    + "The grid below shows only what could be loaded." + Environment.NewLine + Environment.NewLine
                    + string.Join(Environment.NewLine, s_seedErrors),
                    "TMT experimental design", MessageBoxButton.OK, MessageBoxImage.Warning);

                // Cleared once shown, or the same warning reappears every time the window is opened
                // even though nothing was re-read.
                s_seedErrors.Clear();
            }

            // One row per checked spectra file. The row set is not editable here -- files are added
            // and removed on the main window, which is the only place that applies the extension
            // allowlist and the Thermo RawFileReader licence agreement.
            foreach (var file in _spectraFiles.Where(p => p.Use).Select(p => p.FilePath))
            {
                var row = new TmtDesignRow(file);
                if (s_fileState.TryGetValue(file, out var st))
                {
                    row.Fraction = st.Fraction;
                    row.TechnicalReplicate = st.TechnicalReplicate;
                    row.Plex = st.Plex ?? "";
                }
                _rows.Add(row);
            }

            foreach (var kv in s_plexAnnotations)
                _plexAnnotations[kv.Key] = new List<PlexAnnotation>(kv.Value);
        }

        #region Validation
        private string ValidateDesign()
        {
            if (!_rows.Any())
                return "No files defined.";

            CommitPendingEdits();

            // A file with no plex passes every check below -- they are all scoped to a plex -- and
            // then does not appear in the written design at all, because the design is built by
            // grouping on plex. Silently dropping a file the user can see in the grid is worse than
            // refusing to save.
            var unplexed = _rows.Where(r => string.IsNullOrWhiteSpace(r.Plex)).ToList();
            if (unplexed.Any())
            {
                return "Every file must be assigned a Plex. Missing for: "
                    + string.Join(", ", unplexed.Select(r => Path.GetFileName(r.FilePath)));
            }

            if (_rows.Any(r => r.TechnicalReplicate < 1))
                return "Technical Replicate values must be >= 1.";

            if (_rows.Any(r => r.Fraction < 1))
                return "Fraction values must be >= 1.";

            // Every check below is scoped to one plex, because _rows spans plexes and a plex is its
            // own labelling experiment. Unscoped, two files in different plexes both sitting at
            // (Fraction 1, Technical Replicate 1) -- which is the correct design for a two-plex run
            // with one fraction and one technical replicate, and the value the TmtDesignRow
            // constructor assigns -- read as a duplicate and Save was refused. It also broke the
            // round trip: SeedFromDesignFiles stores fraction, replicate and plex per file, so a
            // two-plex design the parser accepts loaded into this grid and then could not be saved
            // back out. TmtExperimentalDesign keys its own per-file state on plex the same way.
            foreach (var plexGroup in _rows.GroupBy(r => r.Plex.Trim(), StringComparer.OrdinalIgnoreCase))
            {
                string plexLabel = $"Plex {plexGroup.Key}";

                var distinctFractions = plexGroup.Select(r => r.Fraction).Distinct().OrderBy(i => i).ToList();
                if (distinctFractions.First() != 1)
                    return $"{plexLabel}: fractions must start at 1.";
                int maxFraction = distinctFractions.Last();
                for (int i = 1; i <= maxFraction; i++)
                    if (!distinctFractions.Contains(i))
                        return $"{plexLabel}: missing fraction number {i} in distinct set.";

                foreach (var grp in plexGroup.GroupBy(r => r.Fraction))
                {
                    var techs = grp.Select(r => r.TechnicalReplicate).Distinct().OrderBy(t => t).ToList();
                    if (techs.First() != 1)
                        return $"{plexLabel}, fraction {grp.Key}: technical replicates must start at 1.";
                    int maxTech = techs.Last();
                    for (int t = 1; t <= maxTech; t++)
                        if (!techs.Contains(t))
                            return $"{plexLabel}, fraction {grp.Key}: missing technical replicate {t}.";
                }

                var duplicatePair = plexGroup.GroupBy(r => (r.Fraction, r.TechnicalReplicate))
                                             .FirstOrDefault(g => g.Count() > 1);
                if (duplicatePair != null)
                    return $"{plexLabel}: duplicate Fraction/Technical Replicate combination: Fraction {duplicatePair.Key.Fraction}, Technical Replicate {duplicatePair.Key.TechnicalReplicate}.";
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

                // The two integer columns bind with the default UpdateSourceTrigger (LostFocus)
                // rather than PropertyChanged, so cancelling here genuinely rejects the edit. Under
                // PropertyChanged the source had already been written by the time this ran and
                // e.Cancel reverted nothing.
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

            // OrdinalIgnoreCase throughout. The grid deduped plex names case-insensitively while the
            // dictionary handed to the annotate window was rebuilt with the default comparer, so
            // "PLEX1" in the design file and "Plex1" in the grid opened a blank annotation grid and
            // then wrote it over the real annotations.
            var plexNames = _rows
                .Select(r => r.Plex.Trim())
                .Where(s => !string.IsNullOrEmpty(s))
                .Distinct(StringComparer.OrdinalIgnoreCase)
                .OrderBy(s => s, StringComparer.OrdinalIgnoreCase)
                .ToList();

            if (!plexNames.Any())
            {
                MessageBox.Show("No Plex values have been set.");
                return;
            }

            var dialog = new AnnotatePlexWindow(plexNames,
                _plexAnnotations.ToDictionary(k => k.Key, v => v.Value.ToList(), StringComparer.OrdinalIgnoreCase))
            {
                Owner = this
            };

            if (dialog.ShowDialog() == true)
            {
                foreach (var kv in dialog.SavedAnnotationsByPlex)
                    _plexAnnotations[kv.Key] = new List<PlexAnnotation>(kv.Value);
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

            var design = BuildDesign();

            // A run refuses to start when both design files sit in one folder, because nothing there
            // could say which of them the user meant (CMD.Program.ResolveExperimentalDesign, and the
            // main window's own pre-run check). Both writers target the first spectra file's
            // directory and the two authoring buttons sit side by side, so warn here rather than let
            // the run be the first place this is discovered.
            var designDirectory = Directory.GetParent(design[0].FullFilePathWithExtension)!.FullName;
            var classicDesignPath = Path.Combine(designDirectory, GlobalVariables.ExperimentalDesignFileName);
            if (File.Exists(classicDesignPath))
            {
                var proceed = MessageBox.Show(
                    $"{GlobalVariables.ExperimentalDesignFileName} is already in this folder:"
                    + Environment.NewLine + designDirectory + Environment.NewLine + Environment.NewLine
                    + "Only one experimental design file may be present when a run starts. With both, "
                    + "MetaMorpheus stops and asks for one of them to be removed." + Environment.NewLine + Environment.NewLine
                    + $"Save {GlobalVariables.TmtExperimentalDesignFileName} anyway?",
                    "Two experimental designs", MessageBoxButton.OKCancel, MessageBoxImage.Warning);

                if (proceed != MessageBoxResult.OK)
                    return;
            }

            string savePath;
            try
            {
                savePath = TmtExperimentalDesign.Write(design);
            }
            catch (Exception ex)
            {
                MessageBox.Show("Failed to write file: " + ex.Message, "Cannot Save");
                return;
            }

            // Read the design straight back through the same parser a command-line run uses, so a
            // file this window reports as saved is one a run can actually load. This window's own
            // validation and the parser's are not the same set of rules, and the parser is the one
            // that decides whether the design works.
            TmtExperimentalDesign.Read(savePath,
                design.Select(f => f.FullFilePathWithExtension).ToList(), out var readErrors);

            if (readErrors.Count > 0)
            {
                MessageBox.Show(
                    "The design was written, but reading it back reported:" + Environment.NewLine + Environment.NewLine
                    + string.Join(Environment.NewLine, readErrors),
                    "Saved design does not load", MessageBoxButton.OK, MessageBoxImage.Warning);
            }
            else
            {
                int plexCount = design.Select(f => f.Plex).Distinct(StringComparer.OrdinalIgnoreCase).Count();
                MessageBox.Show(
                    $"Saved TMT design ({design.Count} file(s), {plexCount} plex(es)) to:"
                    + Environment.NewLine + savePath,
                    "Saved");
            }

            foreach (var r in _rows)
                s_fileState[r.FilePath] = new FileDesignState(r.Fraction, r.TechnicalReplicate, r.Plex.Trim());
            foreach (var kv in _plexAnnotations)
                s_plexAnnotations[kv.Key] = new List<PlexAnnotation>(kv.Value);

            DialogResult = true;
            Close();
        }

        private void CancelButton_Click(object sender, RoutedEventArgs e)
        {
            DialogResult = false;
            Close();
        }
        #endregion

        #region Building the design
        /// <summary>
        /// The grid and the annotations as one <see cref="TmtFileInfo"/> list covering every plex, so
        /// that <see cref="TmtExperimentalDesign.Write"/> rewrites a complete design rather than a
        /// fragment. A plex with no annotations yet is still emitted, as the file-only placeholder row
        /// Read understands, so a part-authored design does not lose its files.
        /// </summary>
        private List<TmtFileInfo> BuildDesign()
        {
            var files = new List<TmtFileInfo>();

            foreach (var plexGroup in _rows
                         .GroupBy(r => r.Plex.Trim(), StringComparer.OrdinalIgnoreCase)
                         .OrderBy(g => g.Key, StringComparer.OrdinalIgnoreCase))
            {
                _plexAnnotations.TryGetValue(plexGroup.Key, out var annotationRows);
                var annotations = ToModelAnnotations(annotationRows);

                foreach (var row in plexGroup.OrderBy(r => r.Fraction).ThenBy(r => r.TechnicalReplicate))
                    files.Add(new TmtFileInfo(row.FilePath, plexGroup.Key, row.Fraction, row.TechnicalReplicate, annotations));
            }

            return files;
        }

        /// <summary>
        /// Grid rows to model annotations. Sample Type goes through
        /// <see cref="TmtExperimentalDesign.TryParseSampleType"/>, the single place that interprets the
        /// spelling, and the free-text cells keep the tab/newline scrubbing this window has always
        /// applied: the model's writer does not escape, so dropping it would let a pasted tab split a row.
        /// </summary>
        private static IReadOnlyList<TmtPlexAnnotation> ToModelAnnotations(List<PlexAnnotation> rows)
        {
            if (rows == null)
                return Array.Empty<TmtPlexAnnotation>();

            return rows.Select(a =>
            {
                TmtExperimentalDesign.TryParseSampleType(a.SampleType, out var sampleType);
                return new TmtPlexAnnotation
                {
                    Tag = a.Tag,
                    SampleName = Escape(a.SampleName),
                    Condition = Escape(a.Condition),
                    BiologicalReplicate = a.BiologicalReplicate,
                    SampleType = sampleType
                };
            }).ToList();
        }

        private static string Escape(string s) =>
            string.IsNullOrEmpty(s) ? "" :
            s.Replace('\t', ' ').Replace('\r', ' ').Replace('\n', ' ');
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
                    _plex = value ?? "";
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

        /// <summary>
        /// Every problem the last seeding run reported, in the order the design files were read.
        /// Surfaced when the window opens rather than at drag-drop time: seeding happens while the
        /// user is still adding files, and a message box per dropped folder would be noise.
        /// </summary>
        private static readonly List<string> s_seedErrors = new();

        /// <summary>
        /// Loads TmtDesign.txt from the folders holding the given spectra files and refreshes the
        /// per-file state and per-plex annotations from it.
        /// </summary>
        /// <param name="rawFilePaths">
        /// The CHECKED spectra files only. Passing every file in the grid produced a "design did not
        /// contain the file(s)" warning naming files the user had deliberately unchecked, and the
        /// window itself only ever lists the checked ones.
        /// </param>
        public static void SeedFromDesignFiles(IEnumerable<string> rawFilePaths)
        {
            s_seedErrors.Clear();
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
                var designPath = Path.Combine(dir!, GlobalVariables.TmtExperimentalDesignFileName);
                if (!File.Exists(designPath))
                    continue;

                var filesInDir = rawSet.Where(f => string.Equals(Path.GetDirectoryName(f), dir, StringComparison.OrdinalIgnoreCase))
                                       .ToList();

                // Read reports a design that is present but unusable -- a channel described twice with
                // disagreeing samples, a Biological Replicate below one, a row too short, an
                // unrecognised Sample Type, or no row naming a file in this run. Discarding the list
                // here opened a partial design, or an empty grid indistinguishable from having no
                // design at all, and said nothing.
                List<TmtFileInfo> tmtFiles;
                try
                {
                    tmtFiles = TmtExperimentalDesign.Read(designPath, filesInDir, out var readErrors);
                    foreach (var error in readErrors)
                        s_seedErrors.Add($"{designPath}: {error}");
                }
                catch (Exception ex)
                {
                    // Read returns the problems it can describe and throws the ones it cannot -- an
                    // unreadable or locked file. Both belong in the same list; swallowing the throw
                    // left the grid silently unseeded.
                    s_seedErrors.Add($"{designPath}: {ex.Message}");
                    continue;
                }

                foreach (var fi in tmtFiles)
                {
                    s_fileState[fi.FullFilePathWithExtension] =
                        new FileDesignState(fi.Fraction, fi.TechnicalReplicate, fi.Plex ?? "");

                    if (string.IsNullOrWhiteSpace(fi.Plex)) continue;

                    s_plexAnnotations[fi.Plex] = fi.Annotations.Select(a => new PlexAnnotation
                    {
                        Tag = a.Tag,
                        SampleName = a.SampleName,
                        Condition = a.Condition,
                        BiologicalReplicate = a.BiologicalReplicate,
                        SampleType = EngineLayer.TmtExperimentalDesign.ToDesignFileValue(a.SampleType)
                    }).ToList();
                }
            }
        }
        #endregion
    }
}
