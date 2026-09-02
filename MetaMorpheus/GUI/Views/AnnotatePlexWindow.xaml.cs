using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Input;
using System.Windows.Media;
using EngineLayer;

namespace MetaMorpheusGUI
{
    public partial class AnnotatePlexWindow : Window
    {
        private readonly List<string> _plexNames;
        private readonly Dictionary<string, List<PlexAnnotation>> _existingAnnotations;
        private readonly Dictionary<string, List<PlexFileEntry>> _plexFileMap;

        public string SavedPlexName { get; private set; } = "";
        public List<PlexAnnotation> SavedAnnotations { get; private set; }

        /// <summary>
        /// The choices offered by the Sample Type column. Taken straight from
        /// <see cref="EngineLayer.TmtExperimentalDesign.SampleTypeNames"/> so the drop-down cannot
        /// offer a value the design-file parser would reject.
        /// </summary>
        public static IReadOnlyList<string> SampleTypeChoices =>
            EngineLayer.TmtExperimentalDesign.SampleTypeNames;

        private List<PlexAnnotation> _currentRows = new();
        private Point _dragStartPoint;
        private PlexAnnotation _dragSourceItem;

        public AnnotatePlexWindow(List<string> plexNames,
            Dictionary<string, List<PlexAnnotation>> existingAnnotations,
            Dictionary<string, List<PlexFileEntry>> plexFileMap)
        {
            InitializeComponent();
            _plexNames = plexNames ?? new();
            _existingAnnotations = existingAnnotations ?? new();
            _plexFileMap = plexFileMap ?? new();

            PlexComboBox.ItemsSource = _plexNames;
            TagTypeComboBox.ItemsSource = Enum.GetValues(typeof(IsobaricMassTagType)).Cast<IsobaricMassTagType>();
        }

        private void PlexComboBox_SelectionChanged(object sender, SelectionChangedEventArgs e) => RegenerateRows();
        private void TagTypeComboBox_SelectionChanged(object sender, SelectionChangedEventArgs e) => RegenerateRows();

        private void RegenerateRows()
        {
            if (PlexComboBox.SelectedItem == null || TagTypeComboBox.SelectedItem == null)
            {
                AnnotationGrid.ItemsSource = null;
                return;
            }

            var plex = PlexComboBox.SelectedItem.ToString();
            var tagType = (IsobaricMassTagType)TagTypeComboBox.SelectedItem;
            var reporterLabels = GetReporterIonLabels(tagType);

            if (_existingAnnotations.TryGetValue(plex, out var existing) &&
                existing.Count == reporterLabels.Count)
            {
                _currentRows = existing
                    .Select(a => new PlexAnnotation
                    {
                        Tag = a.Tag,
                        SampleName = a.SampleName,
                        Condition = a.Condition,
                        BiologicalReplicate = a.BiologicalReplicate,
                        SampleType = a.SampleType
                    }).ToList();
            }
            else
            {
                _currentRows = reporterLabels
                    .Select(lbl => new PlexAnnotation
                    {
                        Tag = lbl,
                        SampleName = "",
                        Condition = "",
                        // 1, not 0: TmtExperimentalDesign.Read rejects a Biological Replicate below 1,
                        // so seeding 0 made the DEFAULT grid produce a design that would not load at all.
                        BiologicalReplicate = 1,
                        SampleType = EngineLayer.TmtExperimentalDesign.ToDesignFileValue(
                            EngineLayer.TmtSampleType.StudySample)
                    }).ToList();
            }

            AnnotationGrid.ItemsSource = _currentRows;
            AnnotationGrid.Items.Refresh();
        }

        /// <summary>
        /// The channel labels for a tag, from the single table in
        /// <see cref="IsobaricMassTag.GetReporterIonLabels(IsobaricMassTagType)"/>.
        /// </summary>
        /// <remarks>
        /// This used to be a second copy of that table, and the two had drifted: TMT10 ended "131"
        /// here against "131N" there, iTRAQ8 ended "121" against "120", and diLeu12 disagreed
        /// throughout. A label typed into the design file by this window then failed to match any
        /// channel of the tag, so the design could not be projected onto mzLib's
        /// <see cref="IExperimentalDesign"/> at all.
        ///
        /// The iTRAQ8 disagreement was the one where this window had it right: 8-plex has no 120
        /// channel -- it is skipped because the phenylalanine immonium ion sits at 120.081 -- and the
        /// eighth reagent is 121, which is what the reporter ion at that index has always been.
        /// IsobaricMassTag was corrected rather than this window changed to match it.
        /// </remarks>
        private static List<string> GetReporterIonLabels(IsobaricMassTagType type) =>
            IsobaricMassTag.GetReporterIonLabels(type) ?? new List<string>();

        private void Paste_CanExecute(object sender, CanExecuteRoutedEventArgs e)
        {
            e.CanExecute = _currentRows is { Count: > 0 };
        }

        private void Paste_Executed(object sender, ExecutedRoutedEventArgs e)
        {
            var raw = Clipboard.GetText();
            if (string.IsNullOrWhiteSpace(raw) || AnnotationGrid.ItemsSource == null)
                return;

            int targetColumn = -1;
            int startRow = 0;

            if (AnnotationGrid.CurrentCell != null)
            {
                targetColumn = AnnotationGrid.CurrentCell.Column?.DisplayIndex ?? -1;
                if (AnnotationGrid.CurrentCell.Item is PlexAnnotation pa)
                {
                    startRow = _currentRows.IndexOf(pa);
                    if (startRow < 0) startRow = 0;
                }
            }

            if (targetColumn < 0)
                targetColumn = 1; // Sample Name column (Tag=0, Sample=1, Condition=2, BioRep=3)

            var lines = raw
                .Split(new[] { '\r', '\n' }, StringSplitOptions.RemoveEmptyEntries)
                .Select(l => l.Trim())
                .Where(l => !string.IsNullOrEmpty(l))
                .ToList();

            if (!lines.Any())
                return;

            for (int i = 0; i < lines.Count && (startRow + i) < _currentRows.Count; i++)
            {
                var row = _currentRows[startRow + i];
                var cols = lines[i]
                    .Split(new[] { '\t', ',' }, StringSplitOptions.RemoveEmptyEntries)
                    .Select(c => c.Trim())
                    .ToArray();

                if (cols.Length == 0)
                    continue;

                switch (targetColumn)
                {
                    case 1: // Sample Name column
                        row.SampleName = cols[0];
                        if (cols.Length >= 2)
                            row.Condition = cols[1];
                        break;
                    case 2: // Condition column
                        row.Condition = cols[0];
                        break;
                    case 3: // Biological Replicate column (single integer)
                        if (int.TryParse(cols[0], out int bio) && bio >= 0)
                            row.BiologicalReplicate = bio;
                        break;
                    default:
                        row.SampleName = cols[0];
                        break;
                }
            }

            AnnotationGrid.Items.Refresh();
        }

        private void AutoFillReplicates_Click(object sender, RoutedEventArgs e)
        {
            // i + 1: replicates are 1-based in the design file, matching Fraction and Technical
            // Replicate. Filling from 0 silently dropped the first channel on read.
            for (int i = 0; i < _currentRows.Count; i++)
                _currentRows[i].BiologicalReplicate = i + 1;

            AnnotationGrid.Items.Refresh();
        }

        private void ValidateButton_Click(object sender, RoutedEventArgs e)
        {
            var error = ValidateCurrent();
            MessageBox.Show(error ?? "Validation passed.",
                error == null ? "OK" : "Validation Error");
        }

        private string ValidateCurrent()
        {
            if (!_currentRows.Any())
                return "No rows to validate.";

            foreach (var r in _currentRows)
            {
                if (r.BiologicalReplicate < 1)
                    return $"Biological replicate must be >= 1 (row tag {r.Tag}).";
            }
            return null;
        }

        private void SaveButton_Click(object sender, RoutedEventArgs e)
        {
            var error = ValidateCurrent();
            if (error != null)
            {
                MessageBox.Show(error, "Cannot Save");
                return;
            }

            var plex = PlexComboBox.SelectedItem?.ToString();
            if (string.IsNullOrEmpty(plex))
            {
                MessageBox.Show("Select a Plex before saving.");
                return;
            }

            if (!_plexFileMap.TryGetValue(plex, out var fileEntries) || fileEntries.Count == 0)
            {
                MessageBox.Show($"No files mapped to plex '{plex}'.");
                return;
            }

            SavedPlexName = plex;
            SavedAnnotations = _currentRows.Select(a => a).ToList();

            // Write through the model instead of a second hand-rolled writer. TmtExperimentalDesign.Write
            // owns the header -- including the Sample Type column this grid collects and the old writer
            // dropped -- and it takes the WHOLE design, so annotating one plex no longer truncates the file
            // and deletes the others. Those were the two defects behind "there is no path through this GUI
            // to a working two-plex experiment".
            var allFiles = BuildWholeDesign(plex);

            string savePath;
            try
            {
                savePath = TmtExperimentalDesign.Write(allFiles);
            }
            catch (Exception ex)
            {
                MessageBox.Show("Failed to write file: " + ex.Message);
                return;
            }

            int plexCount = allFiles.Select(f => f.Plex).Distinct(StringComparer.OrdinalIgnoreCase).Count();
            MessageBox.Show("Saved TMT design (" + allFiles.Count + " file(s), " + plexCount + " plex(es)) to:" + Environment.NewLine + savePath);
            DialogResult = true;
            Close();
        }

        /// <summary>
        /// Every plex's files, not only the one being annotated, so that
        /// <see cref="TmtExperimentalDesign.Write"/> rewrites a complete design rather than a fragment.
        /// The plex open in the grid takes its annotations from the grid; every other plex keeps what it
        /// was seeded with. A plex with no annotations yet is still emitted, as the file-only placeholder
        /// row Read understands, so a part-authored design does not lose its files.
        /// </summary>
        private List<TmtFileInfo> BuildWholeDesign(string currentPlex)
        {
            var files = new List<TmtFileInfo>();

            foreach (var kvp in _plexFileMap.OrderBy(k => k.Key, StringComparer.OrdinalIgnoreCase))
            {
                var rows = string.Equals(kvp.Key, currentPlex, StringComparison.OrdinalIgnoreCase)
                    ? _currentRows
                    : (_existingAnnotations.TryGetValue(kvp.Key, out var saved) ? saved : null);

                var annotations = ToModelAnnotations(rows);

                foreach (var fe in kvp.Value.OrderBy(f => f.Fraction).ThenBy(f => f.TechnicalReplicate))
                    files.Add(new TmtFileInfo(fe.FilePath, kvp.Key, fe.Fraction, fe.TechnicalReplicate, annotations));
            }

            return files;
        }

        /// <summary>
        /// Grid rows to model annotations. Sample Type goes through
        /// <see cref="TmtExperimentalDesign.TryParseSampleType"/>, the single place that interprets the
        /// spelling, and the free-text cells keep the tab/newline scrubbing this window has always applied:
        /// the model's writer does not escape, so dropping it here would let a pasted tab split a row.
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

        private void CancelButton_Click(object sender, RoutedEventArgs e)
        {
            DialogResult = false;
            Close();
        }

        private static string Escape(string s) =>
            string.IsNullOrEmpty(s) ? "" :
            s.Replace('\t', ' ').Replace('\r', ' ').Replace('\n', ' ');
    }
}