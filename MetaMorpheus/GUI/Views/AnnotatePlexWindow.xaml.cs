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
                        BiologicalReplicate = a.BiologicalReplicate
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
                        BiologicalReplicate = 0
                    }).ToList();
            }

            AnnotationGrid.ItemsSource = _currentRows;
            AnnotationGrid.Items.Refresh();
        }

        private List<string> GetReporterIonLabels(IsobaricMassTagType type) => type switch
        {
            IsobaricMassTagType.TMT6 => new() { "126", "127", "128", "129", "130", "131" },
            IsobaricMassTagType.TMT10 => new() { "126", "127N", "127C", "128N", "128C", "129N", "129C", "130N", "130C", "131" },
            IsobaricMassTagType.TMT11 => new() { "126", "127N", "127C", "128N", "128C", "129N", "129C", "130N", "130C", "131N", "131C" },
            IsobaricMassTagType.TMT18 => new() { "126", "127N", "127C", "128N", "128C", "129N", "129C", "130N", "130C", "131N", "131C", "132N", "132C", "133N", "133C", "134N", "134C", "135N" },
            IsobaricMassTagType.iTRAQ4 => new() { "114", "115", "116", "117" },
            IsobaricMassTagType.iTRAQ8 => new() { "113", "114", "115", "116", "117", "118", "119", "121" },
            IsobaricMassTagType.diLeu4 => new() { "115", "116", "117", "118" },
            IsobaricMassTagType.diLeu12 => new() { "115", "116", "117", "118", "119", "120", "121", "122", "123", "124", "125", "126" },
            _ => new() { "126", "127", "128", "129", "130", "131" }
        };

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
            for (int i = 0; i < _currentRows.Count; i++)
                _currentRows[i].BiologicalReplicate = i;

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
                if (r.BiologicalReplicate < 0)
                    return $"Biological replicate must be >= 0 (row tag {r.Tag}).";
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

            var defaultDir = Path.GetDirectoryName(fileEntries.First().FilePath) ?? Environment.CurrentDirectory;
            var defaultName = $"TmtDesign_{plex}.txt";
            var savePath = Path.Combine(defaultDir, defaultName);

            try
            {
                using var writer = new StreamWriter(savePath, false);
                writer.WriteLine("File\tPlex\tSample Name\tTMT Channel\tCondition\tBiological Replicate\tFraction\tTechnical Replicate");

                foreach (var fe in fileEntries)
                {
                    int fractionValue = fe.Fraction; // USE USER-DEFINED FRACTION
                    int techRepValue = fe.TechnicalReplicate;

                    foreach (var ann in _currentRows)
                    {
                        writer.WriteLine(
                            $"{fe.FilePath}\t{plex}\t{Escape(ann.SampleName)}\t{ann.Tag}\t{Escape(ann.Condition)}\t{ann.BiologicalReplicate}\t{fractionValue}\t{techRepValue}");
                    }
                }
            }
            catch (Exception ex)
            {
                MessageBox.Show("Failed to write file: " + ex.Message);
                return;
            }

            MessageBox.Show($"Saved TMT design (files x channels = {fileEntries.Count} x {_currentRows.Count}) to:\n{savePath}");
            DialogResult = true;
            Close();
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