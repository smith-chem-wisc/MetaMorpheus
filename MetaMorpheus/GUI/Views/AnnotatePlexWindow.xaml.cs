using System;
using System.Collections.Generic;
using System.Linq;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Input;
using EngineLayer;

namespace MetaMorpheusGUI
{
    public partial class AnnotatePlexWindow : Window
    {
        private readonly List<string> _plexNames;
        private readonly Dictionary<string, List<PlexAnnotation>> _annotations;

        /// <summary>
        /// Every plex this window holds annotations for, including ones the user never opened.
        /// </summary>
        /// <remarks>
        /// This window no longer writes the design. <see cref="TmtExperimentalDesignWindow"/> owns the
        /// file and writes the whole of it once, on its own Save, the way
        /// <see cref="ExperimentalDesignWindow"/> owns ExperimentalDesign.tsv. Two writers was how a
        /// design came to be written from a grid nothing had reconciled it with.
        /// </remarks>
        public IReadOnlyDictionary<string, List<PlexAnnotation>> SavedAnnotationsByPlex => _annotations;

        /// <summary>
        /// The choices offered by the Sample Type column. Taken straight from
        /// <see cref="EngineLayer.TmtExperimentalDesign.SampleTypeNames"/> so the drop-down cannot
        /// offer a value the design-file parser would reject.
        /// </summary>
        public static IReadOnlyList<string> SampleTypeChoices =>
            EngineLayer.TmtExperimentalDesign.SampleTypeNames;

        private List<PlexAnnotation> _currentRows = new();
        private string _currentPlex;

        public AnnotatePlexWindow(List<string> plexNames,
            Dictionary<string, List<PlexAnnotation>> existingAnnotations)
        {
            InitializeComponent();
            _plexNames = plexNames ?? new();
            _annotations = new Dictionary<string, List<PlexAnnotation>>(
                existingAnnotations ?? new Dictionary<string, List<PlexAnnotation>>(),
                StringComparer.OrdinalIgnoreCase);

            PlexComboBox.ItemsSource = _plexNames;
            TagTypeComboBox.ItemsSource = Enum.GetValues(typeof(IsobaricMassTagType)).Cast<IsobaricMassTagType>();
        }

        private void PlexComboBox_SelectionChanged(object sender, SelectionChangedEventArgs e) => RegenerateRows();
        private void TagTypeComboBox_SelectionChanged(object sender, SelectionChangedEventArgs e) => RegenerateRows();

        /// <summary>
        /// Rebuilds the channel grid for the plex and tag type now selected, after storing whatever
        /// the user had typed for the previous one.
        /// </summary>
        /// <remarks>
        /// Nothing else wrote the grid back, and neither combo starts with a selection, so without the
        /// stash annotating plex 1 and then selecting plex 2 discarded plex 1's samples with no
        /// prompt -- and the design was then written with plex 1 as a file-only placeholder row.
        ///
        /// The Tag always comes from the tag type's own table rather than from the stored annotation.
        /// Channel counts across the eight tag types are 6/10/11/18/4/8/4/12, so 4 is the only count
        /// that occurs twice: iTRAQ4 is 114,115,116,117 and diLeu4 is 115,116,117,118. Reusing stored
        /// tags whenever the counts matched therefore carried iTRAQ4's labels into a diLeu4 plex, and
        /// <see cref="EngineLayer.TmtExperimentalDesign.ToMzLibDesign"/> rejects a tag that is not a
        /// channel of the selected type -- so that design would not load.
        /// </remarks>
        private void RegenerateRows()
        {
            CommitPendingEdits();
            StashCurrentRows();

            if (PlexComboBox.SelectedItem == null || TagTypeComboBox.SelectedItem == null)
            {
                AnnotationGrid.ItemsSource = null;
                _currentPlex = null;
                return;
            }

            var plex = PlexComboBox.SelectedItem.ToString();
            var tagType = (IsobaricMassTagType)TagTypeComboBox.SelectedItem;
            var reporterLabels = GetReporterIonLabels(tagType);

            _annotations.TryGetValue(plex, out var existing);
            bool canCarryOver = existing != null && existing.Count == reporterLabels.Count;

            if (existing != null && !canCarryOver && existing.Any(a => !string.IsNullOrWhiteSpace(a.SampleName)))
            {
                MessageBox.Show(
                    $"Plex '{plex}' was annotated with {existing.Count} channels and {tagType} has "
                    + $"{reporterLabels.Count}. The sample names already entered for it cannot be carried "
                    + "across and have been cleared.",
                    "Channel count changed", MessageBoxButton.OK, MessageBoxImage.Warning);
            }

            _currentRows = reporterLabels.Select((label, i) =>
            {
                var previous = canCarryOver ? existing[i] : null;
                return new PlexAnnotation
                {
                    Tag = label,
                    SampleName = previous?.SampleName ?? "",
                    Condition = previous?.Condition ?? "",
                    // 1, not 0: TmtExperimentalDesign.Read rejects a Biological Replicate below 1,
                    // so seeding 0 made the DEFAULT grid produce a design that would not load at all.
                    BiologicalReplicate = previous != null && previous.BiologicalReplicate >= 1
                        ? previous.BiologicalReplicate
                        : 1,
                    SampleType = previous?.SampleType
                        ?? EngineLayer.TmtExperimentalDesign.ToDesignFileValue(EngineLayer.TmtSampleType.StudySample)
                };
            }).ToList();

            _currentPlex = plex;
            AnnotationGrid.ItemsSource = _currentRows;
            AnnotationGrid.Items.Refresh();
        }

        private void StashCurrentRows()
        {
            if (_currentPlex == null || _currentRows.Count == 0)
                return;

            _annotations[_currentPlex] = _currentRows
                .Select(a => new PlexAnnotation
                {
                    Tag = a.Tag,
                    SampleName = a.SampleName,
                    Condition = a.Condition,
                    BiologicalReplicate = a.BiologicalReplicate,
                    SampleType = a.SampleType
                })
                .ToList();
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
            if (string.IsNullOrEmpty(raw) || AnnotationGrid.ItemsSource == null)
                return;

            // Column order: Tag=0, Sample Name=1, Condition=2, Biological Replicate=3, Sample Type=4.
            int targetColumn = 1;
            int startRow = 0;

            // DataGridCellInfo is a struct, so the old `CurrentCell != null` was always true (CS8073).
            var cell = AnnotationGrid.CurrentCell;
            if (cell.IsValid)
            {
                targetColumn = cell.Column?.DisplayIndex ?? 1;
                if (cell.Item is PlexAnnotation pa)
                {
                    startRow = _currentRows.IndexOf(pa);
                    if (startRow < 0) startRow = 0;
                }
            }

            // Row i of the paste must land on row i of the grid, because rows ARE reporter channels
            // and an unused channel part-way through a plex is an ordinary TMT design. Dropping blank
            // lines moved every sample below the gap onto the wrong channel, silently, in a design
            // that then looked entirely plausible.
            var lines = raw.Replace("\r\n", "\n").Replace('\r', '\n').Split('\n').ToList();
            while (lines.Count > 0 && string.IsNullOrWhiteSpace(lines[lines.Count - 1]))
                lines.RemoveAt(lines.Count - 1);

            if (lines.Count == 0)
                return;

            var rejected = new List<string>();

            for (int i = 0; i < lines.Count && (startRow + i) < _currentRows.Count; i++)
            {
                var row = _currentRows[startRow + i];

                // Tab only, and no RemoveEmptyEntries. This grid collects free text: a sample name of
                // "HeLa, 2h" is ordinary, and splitting on comma turned it into a name and a
                // condition, while a leading empty cell shifted every column left.
                var cols = lines[i].Split('\t').Select(c => c.Trim()).ToArray();

                switch (targetColumn)
                {
                    case 2: // Condition
                        row.Condition = cols[0];
                        break;

                    case 3: // Biological Replicate
                        if (int.TryParse(cols[0], out int bio) && bio >= 1)
                            row.BiologicalReplicate = bio;
                        else if (cols[0].Length > 0)
                            rejected.Add($"{row.Tag}: \"{cols[0]}\" is not a Biological Replicate of 1 or more.");
                        break;

                    case 4: // Sample Type
                        if (EngineLayer.TmtExperimentalDesign.TryParseSampleType(cols[0], out var sampleType))
                            row.SampleType = EngineLayer.TmtExperimentalDesign.ToDesignFileValue(sampleType);
                        else if (cols[0].Length > 0)
                            rejected.Add($"{row.Tag}: \"{cols[0]}\" is not one of {string.Join(", ", SampleTypeChoices)}.");
                        break;

                    default: // Tag or Sample Name -- paste lands in Sample Name
                        row.SampleName = cols[0];
                        if (cols.Length >= 2)
                            row.Condition = cols[1];
                        break;
                }
            }

            AnnotationGrid.Items.Refresh();

            if (rejected.Count > 0)
            {
                MessageBox.Show(string.Join(Environment.NewLine, rejected),
                    "Some pasted values were not applied", MessageBoxButton.OK, MessageBoxImage.Warning);
            }
        }

        private void AutoFillReplicates_Click(object sender, RoutedEventArgs e)
        {
            CommitPendingEdits();

            // i + 1: replicates are 1-based in the design file, matching Fraction and Technical
            // Replicate. Filling from 0 silently dropped the first channel on read.
            for (int i = 0; i < _currentRows.Count; i++)
                _currentRows[i].BiologicalReplicate = i + 1;

            AnnotationGrid.Items.Refresh();
        }

        private void ValidateButton_Click(object sender, RoutedEventArgs e)
        {
            CommitPendingEdits();
            StashCurrentRows();

            var error = ValidateAll();
            MessageBox.Show(error ?? "Validation passed.",
                error == null ? "OK" : "Validation Error");
        }

        /// <summary>
        /// Every plex that has been annotated, not only the one on screen -- the parent window writes
        /// all of them, so all of them have to be loadable.
        /// </summary>
        private string ValidateAll()
        {
            if (_annotations.Count == 0)
                return "No plex has been annotated yet.";

            foreach (var plex in _annotations.OrderBy(k => k.Key, StringComparer.OrdinalIgnoreCase))
            {
                foreach (var r in plex.Value)
                {
                    if (r.BiologicalReplicate < 1)
                        return $"Plex {plex.Key}, channel {r.Tag}: Biological Replicate must be 1 or more.";

                    if (!EngineLayer.TmtExperimentalDesign.TryParseSampleType(r.SampleType, out _))
                        return $"Plex {plex.Key}, channel {r.Tag}: \"{r.SampleType}\" is not a recognised Sample Type.";
                }
            }

            return null;
        }

        private void SaveButton_Click(object sender, RoutedEventArgs e)
        {
            CommitPendingEdits();
            StashCurrentRows();

            var error = ValidateAll();
            if (error != null)
            {
                MessageBox.Show(error, "Cannot Save");
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

        /// <summary>
        /// Pushes the cell being edited into the bound object before the grid is read. Without it a
        /// value typed into a cell that still had focus was simply not there on Save, unlike in
        /// <see cref="TmtExperimentalDesignWindow"/>, which has always committed before every read.
        /// </summary>
        private void CommitPendingEdits()
        {
            AnnotationGrid.CommitEdit(DataGridEditingUnit.Cell, true);
            AnnotationGrid.CommitEdit(DataGridEditingUnit.Row, true);
        }
    }
}
