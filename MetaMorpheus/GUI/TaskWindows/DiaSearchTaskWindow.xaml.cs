// Copyright 2026 MetaMorpheus Contributors
// Licensed under the MIT License

using EngineLayer.DiaSearch;
using System;
using System.Globalization;
using System.Windows;
using System.Windows.Input;
using TaskLayer;

namespace MetaMorpheusGUI
{
    /// <summary>
    /// Interaction logic for DiaSearchTaskWindow.xaml
    ///
    /// Mirrors the pattern of SearchTaskWindow / XLSearchTaskWindow:
    ///   - Constructor populates controls from an existing task (or creates a default one)
    ///   - SaveButton_Click validates inputs, writes back to TheTask, sets DialogResult = true
    ///   - Caller reads TheTask after ShowDialog() == true
    /// </summary>
    public partial class DiaSearchTaskWindow : Window
    {
        // ── Public result ─────────────────────────────────────────────────────

        /// <summary>
        /// The configured task. Populated when the user clicks "Add the DIA Search Task".
        /// </summary>
        public DiaSearchTask TheTask { get; private set; }

        // ── Construction ──────────────────────────────────────────────────────

        /// <summary>
        /// Opens the window pre-populated from an existing task (edit mode).
        /// Pass <c>null</c> to open with defaults (new task mode).
        /// </summary>
        public DiaSearchTaskWindow(DiaSearchTask task = null)
        {
            InitializeComponent();
            TheTask = task ?? new DiaSearchTask();
            PopulateFromTask(TheTask);
        }

        // ── Populate UI from task ─────────────────────────────────────────────

        private void PopulateFromTask(DiaSearchTask task)
        {
            var p = task.DiaSearchParameters;

            // Task name
            OutputFileNameTextBox.Text = task.CommonParameters.TaskDescriptor;

            // Extraction
            PpmToleranceTextBox.Text   = p.PpmTolerance.ToString("G", CultureInfo.InvariantCulture);
            RtWindowTextBox.Text       = p.RtToleranceMinutes.ToString("G", CultureInfo.InvariantCulture);
            MinFragmentsTextBox.Text   = p.MinFragmentsRequired.ToString(CultureInfo.InvariantCulture);

            // Calibration
            UseIrtCalibrationCheckBox.IsChecked = p.UseIrtCalibration;

            // Performance
            MaxThreadsTextBox.Text  = p.MaxThreads.ToString(CultureInfo.InvariantCulture);
            PreferGpuCheckBox.IsChecked = p.PreferGpu;

            // Output
            WritePsmTsvCheckBox.IsChecked       = p.WritePsmTsv;
            WriteDiagnosticsCheckBox.IsChecked   = p.WriteDiagnostics;
            WriteDecoyResultsCheckBox.IsChecked  = p.WriteDecoyResults;
        }

        // ── Save / validate ───────────────────────────────────────────────────

        private void SaveButton_Click(object sender, RoutedEventArgs e)
        {
            if (!ValidateAndApply())
                return;

            DialogResult = true;
        }

        /// <summary>
        /// Validates all text-box inputs. On success, writes back to TheTask and returns true.
        /// On failure, sets ErrorTextBox.Text and returns false.
        /// </summary>
        private bool ValidateAndApply()
        {
            ErrorTextBox.Text = string.Empty;

            // ── ppm tolerance ──────────────────────────────────────────────
            if (!float.TryParse(PpmToleranceTextBox.Text, NumberStyles.Float,
                    CultureInfo.InvariantCulture, out float ppm) || ppm <= 0)
            {
                ErrorTextBox.Text = "Fragment m/z tolerance must be a positive number.";
                return false;
            }

            // ── RT window ──────────────────────────────────────────────────
            if (!float.TryParse(RtWindowTextBox.Text, NumberStyles.Float,
                    CultureInfo.InvariantCulture, out float rtWindow) || rtWindow <= 0)
            {
                ErrorTextBox.Text = "RT window must be a positive number.";
                return false;
            }

            // ── min fragments ──────────────────────────────────────────────
            if (!int.TryParse(MinFragmentsTextBox.Text, out int minFrags) || minFrags < 1)
            {
                ErrorTextBox.Text = "Minimum fragments required must be an integer ≥ 1.";
                return false;
            }

            // ── max threads ───────────────────────────────────────────────
            if (!int.TryParse(MaxThreadsTextBox.Text, out int maxThreads) || maxThreads < -1 || maxThreads == 0)
            {
                ErrorTextBox.Text = "Max threads must be -1 (all cores) or a positive integer.";
                return false;
            }

            // ── Write back ────────────────────────────────────────────────

            // Task name: if the user left it blank, the default descriptor is used
            string taskName = OutputFileNameTextBox.Text.Trim();
            if (!string.IsNullOrEmpty(taskName))
            {
                TheTask.CommonParameters =
                    new EngineLayer.CommonParameters(taskDescriptor: taskName);
            }

            var p = TheTask.DiaSearchParameters;
            p.PpmTolerance          = ppm;
            p.RtToleranceMinutes    = rtWindow;
            p.MinFragmentsRequired  = minFrags;
            p.MaxThreads            = maxThreads;
            p.PreferGpu             = PreferGpuCheckBox.IsChecked == true;
            p.UseIrtCalibration     = UseIrtCalibrationCheckBox.IsChecked == true;
            p.WritePsmTsv           = WritePsmTsvCheckBox.IsChecked == true;
            p.WriteDiagnostics      = WriteDiagnosticsCheckBox.IsChecked == true;
            p.WriteDecoyResults     = WriteDecoyResultsCheckBox.IsChecked == true;

            return true;
        }

        // ── Keyboard shortcut ─────────────────────────────────────────────────

        private void KeyPressed(object sender, KeyEventArgs e)
        {
            if (e.Key == Key.Return)
                SaveButton_Click(sender, e);
            else if (e.Key == Key.Escape)
                DialogResult = false;
        }
    }
}
