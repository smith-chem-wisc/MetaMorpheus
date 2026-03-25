using EngineLayer;
using GuiFunctions;
using MassSpectrometry;
using MzLibUtil;
using Omics.Modifications;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.Globalization;
using System.Linq;
using System.Windows;
using System.Windows.Controls;
using TaskLayer;

namespace MetaMorpheusGUI
{
    /// <summary>
    /// Code-behind for CircularSearchTaskWindow.
    /// Modification tree follows the identical pattern used by SearchTaskWindow:
    ///   • PopulateChoices() builds ModTypeForTreeViewModel groups from GlobalVariables.AllModsKnown
    ///   • UpdateFieldsFromTask() ticks the boxes that were saved in the task's CommonParameters
    ///   • AddTaskButton_Click collects checked boxes back into ListOfModsFixed / ListOfModsVariable
    ///   • SearchFixedMod / SearchVariableMod TextBoxes filter the trees via SearchModifications
    /// </summary>
    public partial class CircularSearchTaskWindow : Window
    {
        // ── Observable collections for the two mod trees ─────────────────────
        private readonly ObservableCollection<ModTypeForTreeViewModel>
            FixedModTypeForTreeViewObservableCollection = new();

        private readonly ObservableCollection<ModTypeForTreeViewModel>
            VariableModTypeForTreeViewObservableCollection = new();

        // ── The task being created / edited ──────────────────────────────────
        public CircularSearchTask TheTask { get; private set; }

        // ── Constructor ───────────────────────────────────────────────────────
        public CircularSearchTaskWindow(CircularSearchTask task = null)
        {
            InitializeComponent();   // ← MUST be first so named controls exist

            TheTask = task ?? new CircularSearchTask();

            PopulateChoices();
            UpdateFieldsFromTask(TheTask);

            // Wire the search-filter timer (same as SearchTaskWindow)
            SearchModifications.Timer.Tick += TextChangeTimerHandler;
            base.Closing += OnWindowClosing;
        }

        // ── PopulateChoices ───────────────────────────────────────────────────
        private void PopulateChoices()
        {
            // Dissociation types
            foreach (string dt in GlobalVariables.AllSupportedDissociationTypes.Keys)
                DissociationTypeComboBox.Items.Add(dt);
            DissociationTypeComboBox.SelectedItem = "HCD";

            // Precursor tolerance units
            PrecursorMassToleranceComboBox.Items.Add("Da");
            PrecursorMassToleranceComboBox.Items.Add("ppm");
            PrecursorMassToleranceComboBox.SelectedIndex = 1; // default ppm

            // Mass-diff acceptor types
            foreach (MassDiffAcceptorType t in Enum.GetValues(typeof(MassDiffAcceptorType)))
                MassDiffAcceptorComboBox.Items.Add(t);
            MassDiffAcceptorComboBox.SelectedItem = MassDiffAcceptorType.OneMM;

            // Proteases
            foreach (Protease p in ProteaseDictionary.Dictionary.Values)
                ProteaseComboBox.Items.Add(p);
            ProteaseComboBox.SelectedItem = ProteaseDictionary.Dictionary["trypsin"];

            // ── Build mod trees — identical logic to SearchTaskWindow ─────────
            var modsToUse = GlobalVariables.AllModsKnown.ToList();

            // Fixed
            foreach (var group in modsToUse
                         .Where(b => b.ValidModification == true)
                         .GroupBy(b => b.ModificationType))
            {
                var typeNode = new ModTypeForTreeViewModel(group.Key, false);
                FixedModTypeForTreeViewObservableCollection.Add(typeNode);
                foreach (var uah in group)
                    typeNode.Children.Add(
                        new ModForTreeViewModel(uah.ToString(), false, uah.IdWithMotif, false, typeNode));
            }
            FixedModsTreeView.DataContext = FixedModTypeForTreeViewObservableCollection;

            // Variable
            foreach (var group in modsToUse
                         .Where(b => b.ValidModification == true)
                         .GroupBy(b => b.ModificationType))
            {
                var typeNode = new ModTypeForTreeViewModel(group.Key, false);
                VariableModTypeForTreeViewObservableCollection.Add(typeNode);
                foreach (var uah in group)
                    typeNode.Children.Add(
                        new ModForTreeViewModel(uah.ToString(), false, uah.IdWithMotif, false, typeNode));
            }
            VariableModsTreeView.DataContext = VariableModTypeForTreeViewObservableCollection;
        }

        // ── UpdateFieldsFromTask ──────────────────────────────────────────────
        private void UpdateFieldsFromTask(CircularSearchTask task)
        {
            // Fragmentation
            DissociationTypeComboBox.SelectedItem =
                task.CommonParameters.DissociationType.ToString();

            MinInternalFragmentLengthTextBox.Text =
                task.CircularSearchParameters.MinInternalFragmentLength
                    .ToString(CultureInfo.InvariantCulture);

            // Precursor tolerance
            PrecursorMassToleranceTextBox.Text =
                task.CommonParameters.PrecursorMassTolerance.Value
                    .ToString(CultureInfo.InvariantCulture);
            PrecursorMassToleranceComboBox.SelectedIndex =
                task.CommonParameters.PrecursorMassTolerance is AbsoluteTolerance ? 0 : 1;

            // Mass-diff acceptor
            MassDiffAcceptorComboBox.SelectedItem =
                task.CircularSearchParameters.MassDiffAcceptorType;

            // Digestion
            if (task.CommonParameters.DigestionParams is DigestionParams dp
                && ProteaseDictionary.Dictionary.TryGetValue(
                    dp.SpecificProtease.Name, out var savedProtease))
                ProteaseComboBox.SelectedItem = savedProtease;

            MissedCleavagesTextBox.Text =
                task.CommonParameters.DigestionParams.MaxMissedCleavages == int.MaxValue
                    ? ""
                    : task.CommonParameters.DigestionParams.MaxMissedCleavages
                          .ToString(CultureInfo.InvariantCulture);

            MaxModsTextBox.Text =
                task.CommonParameters.DigestionParams.MaxMods
                    .ToString(CultureInfo.InvariantCulture);

            // ── Restore fixed mod check states ────────────────────────────────
            foreach (var mod in task.CommonParameters.ListOfModsFixed)
            {
                var typeNode = FixedModTypeForTreeViewObservableCollection
                    .FirstOrDefault(b => b.DisplayName == mod.Item1);

                if (typeNode != null)
                {
                    var modNode = typeNode.Children.FirstOrDefault(b => b.ModName == mod.Item2);
                    if (modNode != null)
                        modNode.Use = true;
                    else
                        typeNode.Children.Add(
                            new ModForTreeViewModel("UNKNOWN MODIFICATION!", true,
                                mod.Item2, true, typeNode));
                }
                else
                {
                    typeNode = new ModTypeForTreeViewModel(mod.Item1, true);
                    FixedModTypeForTreeViewObservableCollection.Add(typeNode);
                    typeNode.Children.Add(
                        new ModForTreeViewModel("UNKNOWN MODIFICATION!", true,
                            mod.Item2, true, typeNode));
                }
            }

            // ── Restore variable mod check states ─────────────────────────────
            foreach (var mod in task.CommonParameters.ListOfModsVariable)
            {
                var typeNode = VariableModTypeForTreeViewObservableCollection
                    .FirstOrDefault(b => b.DisplayName == mod.Item1);

                if (typeNode != null)
                {
                    var modNode = typeNode.Children.FirstOrDefault(b => b.ModName == mod.Item2);
                    if (modNode != null)
                        modNode.Use = true;
                    else
                        typeNode.Children.Add(
                            new ModForTreeViewModel("UNKNOWN MODIFICATION!", true,
                                mod.Item2, true, typeNode));
                }
                else
                {
                    typeNode = new ModTypeForTreeViewModel(mod.Item1, true);
                    VariableModTypeForTreeViewObservableCollection.Add(typeNode);
                    typeNode.Children.Add(
                        new ModForTreeViewModel("UNKNOWN MODIFICATION!", true,
                            mod.Item2, true, typeNode));
                }
            }

            // Sync parent tri-state checkboxes
            foreach (var node in FixedModTypeForTreeViewObservableCollection)
                node.VerifyCheckState();
            foreach (var node in VariableModTypeForTreeViewObservableCollection)
                node.VerifyCheckState();

            // Output options
            WriteDecoyCheckBox.IsChecked = task.CircularSearchParameters.WriteDecoys;
            WriteContaminantCheckBox.IsChecked = task.CircularSearchParameters.WriteContaminants;
            WriteHighQValueCheckBox.IsChecked = task.CircularSearchParameters.WriteHighQValuePsms;
        }

        // ── Add Task (Save) ───────────────────────────────────────────────────
        private void AddTaskButton_Click(object sender, RoutedEventArgs e)
        {
            // Validate numeric fields
            if (!double.TryParse(PrecursorMassToleranceTextBox.Text,
                    NumberStyles.Any, CultureInfo.InvariantCulture, out double precTolValue)
                || precTolValue <= 0)
            {
                MessageBox.Show("Please enter a valid positive precursor mass tolerance.");
                return;
            }

            if (!int.TryParse(MinInternalFragmentLengthTextBox.Text,
                    out int minFragLen) || minFragLen < 1)
            {
                MessageBox.Show("Min internal fragment length must be a positive integer.");
                return;
            }

            if (!int.TryParse(MissedCleavagesTextBox.Text, out int missedCleavages)
                && MissedCleavagesTextBox.Text != "")
            {
                MessageBox.Show("Please enter a valid number of missed cleavages.");
                return;
            }

            if (!int.TryParse(MaxModsTextBox.Text, out int maxMods) || maxMods < 0)
            {
                MessageBox.Show("Max mods per peptide must be a non-negative integer.");
                return;
            }

            // Collect fixed mods
            var listOfModsFixed = new List<(string, string)>();
            foreach (var typeNode in FixedModTypeForTreeViewObservableCollection)
                listOfModsFixed.AddRange(
                    typeNode.Children
                            .Where(c => c.Use)
                            .Select(c => (c.Parent.DisplayName, c.ModName)));

            // Collect variable mods
            var listOfModsVariable = new List<(string, string)>();
            foreach (var typeNode in VariableModTypeForTreeViewObservableCollection)
                listOfModsVariable.AddRange(
                    typeNode.Children
                            .Where(c => c.Use)
                            .Select(c => (c.Parent.DisplayName, c.ModName)));

            // Build tolerance
            Tolerance precursorTolerance =
                PrecursorMassToleranceComboBox.SelectedIndex == 0
                    ? (Tolerance)new AbsoluteTolerance(precTolValue)
                    : new PpmTolerance(precTolValue);

            // Dissociation type
            DissociationType dissociationType =
                GlobalVariables.AllSupportedDissociationTypes[
                    DissociationTypeComboBox.SelectedItem.ToString()];

            // Protease
            int missedCleavagesValue = MissedCleavagesTextBox.Text == ""
                ? int.MaxValue
                : missedCleavages;

            var digestionParams = new DigestionParams(
                protease: ((Protease)ProteaseComboBox.SelectedItem).Name,
                maxMissedCleavages: missedCleavagesValue,
                maxModsForPeptides: maxMods);

            // Build CommonParameters
            TheTask.CommonParameters = new CommonParameters(
                dissociationType: dissociationType,
                precursorMassTolerance: precursorTolerance,
                digestionParams: digestionParams,
                listOfModsFixed: listOfModsFixed,
                listOfModsVariable: listOfModsVariable,
                taskDescriptor: "CircularSearch");

            // Circular-specific parameters
            TheTask.CircularSearchParameters.MinInternalFragmentLength = minFragLen;
            TheTask.CircularSearchParameters.MassDiffAcceptorType =
                (MassDiffAcceptorType)MassDiffAcceptorComboBox.SelectedItem;
            TheTask.CircularSearchParameters.WriteDecoys =
                WriteDecoyCheckBox.IsChecked == true;
            TheTask.CircularSearchParameters.WriteContaminants =
                WriteContaminantCheckBox.IsChecked == true;
            TheTask.CircularSearchParameters.WriteHighQValuePsms =
                WriteHighQValueCheckBox.IsChecked == true;

            DialogResult = true;
        }

        private void CancelButton_Click(object sender, RoutedEventArgs e)
        {
            DialogResult = false;
        }

        // ── Mod search-filter handlers (identical pattern to SearchTaskWindow) ─
        private void SearchFixedMod_TextChanged(object sender, TextChangedEventArgs e)
        {
            SearchModifications.SetTimer();
            SearchModifications.FixedSearch = true;
        }

        private void SearchVariableMod_TextChanged(object sender, TextChangedEventArgs e)
        {
            SearchModifications.SetTimer();
            SearchModifications.VariableSearch = true;
        }

        private void TextChangeTimerHandler(object sender, EventArgs e)
        {
            if (SearchModifications.FixedSearch)
            {
                SearchModifications.FilterTree(
                    SearchFixedMod,
                    FixedModsTreeView,
                    FixedModTypeForTreeViewObservableCollection);
                SearchModifications.FixedSearch = false;
            }

            if (SearchModifications.VariableSearch)
            {
                SearchModifications.FilterTree(
                    SearchVariableMod,
                    VariableModsTreeView,
                    VariableModTypeForTreeViewObservableCollection);
                SearchModifications.VariableSearch = false;
            }
        }

        private void OnWindowClosing(object sender, System.ComponentModel.CancelEventArgs e)
        {
            SearchModifications.Timer.Tick -= TextChangeTimerHandler;
        }
    }
}
