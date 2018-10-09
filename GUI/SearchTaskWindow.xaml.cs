﻿using EngineLayer;
using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.ComponentModel;
using System.Globalization;
using System.Linq;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Input;
using TaskLayer;
using UsefulProteomicsDatabases;
using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using MassSpectrometry;

namespace MetaMorpheusGUI
{
    /// <summary>
    /// Interaction logic for SearchTaskWindow.xaml
    /// </summary>
    public partial class SearchTaskWindow : Window
    {
        private readonly DataContextForSearchTaskWindow DataContextForSearchTaskWindow;
        private readonly ObservableCollection<SearchModeForDataGrid> SearchModesForThisTask = new ObservableCollection<SearchModeForDataGrid>();
        private readonly ObservableCollection<ModTypeForTreeView> FixedModTypeForTreeViewObservableCollection = new ObservableCollection<ModTypeForTreeView>();
        private readonly ObservableCollection<ModTypeForTreeView> VariableModTypeForTreeViewObservableCollection = new ObservableCollection<ModTypeForTreeView>();
        private readonly ObservableCollection<ModTypeForLoc> LocalizeModTypeForTreeViewObservableCollection = new ObservableCollection<ModTypeForLoc>();
        private readonly ObservableCollection<ModTypeForGrid> ModSelectionGridItems = new ObservableCollection<ModTypeForGrid>();

        public SearchTaskWindow()
        {
            InitializeComponent();
            TheTask = new SearchTask();

            PopulateChoices();

            UpdateFieldsFromTask(TheTask);

            this.saveButton.Content = "Add the Search Task";

            DataContextForSearchTaskWindow = new DataContextForSearchTaskWindow
            {
                ExpanderTitle = string.Join(", ", SearchModesForThisTask.Where(b => b.Use).Select(b => b.Name)),
                AnalysisExpanderTitle = "Some analysis properties...",
                SearchModeExpanderTitle = "Some search properties..."
            };
            this.DataContext = DataContextForSearchTaskWindow;
        }

        public SearchTaskWindow(SearchTask task)
        {
            InitializeComponent();
            TheTask = task;

            PopulateChoices();

            UpdateFieldsFromTask(TheTask);

            DataContextForSearchTaskWindow = new DataContextForSearchTaskWindow
            {
                ExpanderTitle = string.Join(", ", SearchModesForThisTask.Where(b => b.Use).Select(b => b.Name)),
                AnalysisExpanderTitle = "Some analysis properties...",
                SearchModeExpanderTitle = "Some search properties..."
            };
            this.DataContext = DataContextForSearchTaskWindow;
        }

        internal SearchTask TheTask { get; private set; }

        private void CheckIfNumber(object sender, TextCompositionEventArgs e)
        {
            e.Handled = !GlobalGuiSettings.CheckIsNumber(e.Text);
        }

        private void Row_DoubleClick(object sender, MouseButtonEventArgs e)
        {
            var ye = sender as DataGridCell;
            if (ye.Content is TextBlock hm && !string.IsNullOrEmpty(hm.Text))
            {
                System.Diagnostics.Process.Start(hm.Text);
            }
        }

        private void PopulateChoices()
        {
            foreach (Protease protease in ProteaseDictionary.Dictionary.Values)
            {
                proteaseComboBox.Items.Add(protease);
            }
            proteaseComboBox.SelectedIndex = 12;

            foreach (string initiatior_methionine_behavior in Enum.GetNames(typeof(InitiatorMethionineBehavior)))
            {
                initiatorMethionineBehaviorComboBox.Items.Add(initiatior_methionine_behavior);
            }

            foreach (string dissassociationType in GlobalVariables.AllSupportedDissociationTypes.Keys)
            {
                dissociationTypeComboBox.Items.Add(dissassociationType);
            }

            productMassToleranceComboBox.Items.Add("Da");
            productMassToleranceComboBox.Items.Add("ppm");

            precursorMassToleranceComboBox.Items.Add("Da");
            precursorMassToleranceComboBox.Items.Add("ppm");

            foreach (var hm in GlobalVariables.AllModsKnown.Where(b => b.ValidModification == true).GroupBy(b => b.ModificationType))
            {
                var theModType = new ModTypeForGrid(hm.Key);
                ModSelectionGridItems.Add(theModType);
            }
            ModSelectionGrid.ItemsSource = ModSelectionGridItems;

            foreach (var hm in GlobalVariables.AllModsKnown.Where(b => b.ValidModification == true).GroupBy(b => b.ModificationType))
            {
                var theModType = new ModTypeForTreeView(hm.Key, false);
                FixedModTypeForTreeViewObservableCollection.Add(theModType);
                foreach (var uah in hm)
                {
                    theModType.Children.Add(new ModForTreeView(uah.ToString(), false, uah.IdWithMotif, false, theModType));
                }
            }
            fixedModsTreeView.DataContext = FixedModTypeForTreeViewObservableCollection;

            foreach (var hm in GlobalVariables.AllModsKnown.Where(b => b.ValidModification == true).GroupBy(b => b.ModificationType))
            {
                var theModType = new ModTypeForTreeView(hm.Key, false);
                VariableModTypeForTreeViewObservableCollection.Add(theModType);
                foreach (var uah in hm)
                {
                    theModType.Children.Add(new ModForTreeView(uah.ToString(), false, uah.IdWithMotif, false, theModType));
                }
            }
            variableModsTreeView.DataContext = VariableModTypeForTreeViewObservableCollection;

            foreach (var hm in GlobalVariables.AllModsKnown.Where(b => b.ValidModification == true).GroupBy(b => b.ModificationType))
            {
                LocalizeModTypeForTreeViewObservableCollection.Add(new ModTypeForLoc(hm.Key));
            }
        }

        private void UpdateFieldsFromTask(SearchTask task)
        {
            classicSearchRadioButton.IsChecked = task.SearchParameters.SearchType == SearchType.Classic;
            modernSearchRadioButton.IsChecked = task.SearchParameters.SearchType == SearchType.Modern;
            nonSpecificSearchRadioButton.IsChecked = task.SearchParameters.SearchType == SearchType.NonSpecific && task.CommonParameters.DigestionParams.SearchModeType == CleavageSpecificity.None;
            semiSpecificSearchRadioButton.IsChecked = task.SearchParameters.SearchType == SearchType.NonSpecific && task.CommonParameters.DigestionParams.SearchModeType != CleavageSpecificity.None;
            MaxFragmentMassTextBox.Text = task.SearchParameters.MaxFragmentSize.ToString(CultureInfo.InvariantCulture);
            checkBoxParsimony.IsChecked = task.SearchParameters.DoParsimony;
            checkBoxNoOneHitWonders.IsChecked = task.SearchParameters.NoOneHitWonders;
            checkBoxQuantification.IsChecked = task.SearchParameters.DoQuantification;
            peakFindingToleranceTextBox.Text = task.SearchParameters.QuantifyPpmTol.ToString(CultureInfo.InvariantCulture);
            checkBoxMatchBetweenRuns.IsChecked = task.SearchParameters.MatchBetweenRuns;
            checkBoxNormalize.IsChecked = task.SearchParameters.Normalize;
            modPepsAreUnique.IsChecked = task.SearchParameters.ModPeptidesAreDifferent;
            checkBoxHistogramAnalysis.IsChecked = task.SearchParameters.DoHistogramAnalysis;
            histogramBinWidthTextBox.Text = task.SearchParameters.HistogramBinTolInDaltons.ToString(CultureInfo.InvariantCulture);
            checkBoxTarget.IsChecked = task.SearchParameters.SearchTarget;
            checkBoxDecoy.IsChecked = task.SearchParameters.DecoyType != DecoyType.None;
            radioButtonReverseDecoy.IsChecked = task.SearchParameters.DecoyType == DecoyType.Reverse;
            radioButtonSlideDecoy.IsChecked = task.SearchParameters.DecoyType == DecoyType.Slide;
            missedCleavagesTextBox.Text = task.CommonParameters.DigestionParams.MaxMissedCleavages == int.MaxValue ? "" : task.CommonParameters.DigestionParams.MaxMissedCleavages.ToString(CultureInfo.InvariantCulture);
            MinPeptideLengthTextBox.Text = task.CommonParameters.DigestionParams.MinPeptideLength.ToString(CultureInfo.InvariantCulture);
            MaxPeptideLengthTextBox.Text = task.CommonParameters.DigestionParams.MaxPeptideLength == int.MaxValue ? "" : task.CommonParameters.DigestionParams.MaxPeptideLength.ToString(CultureInfo.InvariantCulture);
            proteaseComboBox.SelectedItem = task.CommonParameters.DigestionParams.SpecificProtease;
            maxModificationIsoformsTextBox.Text = task.CommonParameters.DigestionParams.MaxModificationIsoforms.ToString(CultureInfo.InvariantCulture);
            MaxModNumTextBox.Text = task.CommonParameters.DigestionParams.MaxModsForPeptide.ToString(CultureInfo.InvariantCulture);
            initiatorMethionineBehaviorComboBox.SelectedIndex = (int)task.CommonParameters.DigestionParams.InitiatorMethionineBehavior;
            dissociationTypeComboBox.SelectedItem = task.CommonParameters.DissociationType.ToString();
            nTerminalIons.IsChecked = task.CommonParameters.DigestionParams.FragmentationTerminus == FragmentationTerminus.Both || task.CommonParameters.DigestionParams.FragmentationTerminus == FragmentationTerminus.N;
            cTerminalIons.IsChecked = task.CommonParameters.DigestionParams.FragmentationTerminus == FragmentationTerminus.Both || task.CommonParameters.DigestionParams.FragmentationTerminus == FragmentationTerminus.C;
            productMassToleranceTextBox.Text = task.CommonParameters.ProductMassTolerance.Value.ToString(CultureInfo.InvariantCulture);
            productMassToleranceComboBox.SelectedIndex = task.CommonParameters.ProductMassTolerance is AbsoluteTolerance ? 0 : 1;
            precursorMassToleranceTextBox.Text = task.CommonParameters.PrecursorMassTolerance.Value.ToString(CultureInfo.InvariantCulture);
            precursorMassToleranceComboBox.SelectedIndex = task.CommonParameters.PrecursorMassTolerance is AbsoluteTolerance ? 0 : 1;
            addCompIonCheckBox.IsChecked = task.CommonParameters.AddCompIons;
            numberOfDatabaseSearchesTextBox.Text = task.CommonParameters.TotalPartitions.ToString(CultureInfo.InvariantCulture);
            deconvolutePrecursors.IsChecked = task.CommonParameters.DoPrecursorDeconvolution;
            useProvidedPrecursor.IsChecked = task.CommonParameters.UseProvidedPrecursorInfo;
            allAmbiguity.IsChecked = task.CommonParameters.ReportAllAmbiguity;
            DeconvolutionMaxAssumedChargeStateTextBox.Text = task.CommonParameters.DeconvolutionMaxAssumedChargeState.ToString();
            minScoreAllowed.Text = task.CommonParameters.ScoreCutoff.ToString(CultureInfo.InvariantCulture);
            eValueCheckBox.IsChecked = task.CommonParameters.CalculateEValue;
            deltaScoreCheckBox.IsChecked = task.CommonParameters.UseDeltaScore;
            trimMs1.IsChecked = task.CommonParameters.TrimMs1Peaks;
            trimMsMs.IsChecked = task.CommonParameters.TrimMsMsPeaks;
            TopNPeaksTextBox.Text = task.CommonParameters.TopNpeaks == int.MaxValue ? "" : task.CommonParameters.TopNpeaks.ToString(CultureInfo.InvariantCulture);
            MinRatioTextBox.Text = task.CommonParameters.MinRatio.ToString(CultureInfo.InvariantCulture);
            maxThreadsTextBox.Text = task.CommonParameters.MaxThreadsToUsePerFile.ToString(CultureInfo.InvariantCulture);

            if (task.CommonParameters.QValueOutputFilter < 1)
            {
                QValueTextBox.Text = task.CommonParameters.QValueOutputFilter.ToString(CultureInfo.InvariantCulture);
                QValueCheckBox.IsChecked = true;
            }
            else
            {
                QValueTextBox.Text = "0.01";
                QValueCheckBox.IsChecked = false;
            }

            OutputFileNameTextBox.Text = task.CommonParameters.TaskDescriptor;
            //ckbPepXML.IsChecked = task.SearchParameters.OutPepXML;
            ckbMzId.IsChecked = task.SearchParameters.WriteMzId;
            writeDecoyCheckBox.IsChecked = task.SearchParameters.WriteDecoys;
            writeContaminantCheckBox.IsChecked = task.SearchParameters.WriteContaminants;

            foreach (var mod in task.CommonParameters.ListOfModsFixed)
            {
                var theModType = FixedModTypeForTreeViewObservableCollection.FirstOrDefault(b => b.DisplayName.Equals(mod.Item1));
                if (theModType != null)
                {
                    var theMod = theModType.Children.FirstOrDefault(b => b.DisplayName.Equals(mod.Item2));
                    if (theMod != null)
                    {
                        theMod.Use = true;
                    }
                    else
                    {
                        theModType.Children.Add(new ModForTreeView("UNKNOWN MODIFICATION!", true, mod.Item2, true, theModType));
                    }
                }
                else
                {
                    theModType = new ModTypeForTreeView(mod.Item1, true);
                    FixedModTypeForTreeViewObservableCollection.Add(theModType);
                    theModType.Children.Add(new ModForTreeView("UNKNOWN MODIFICATION!", true, mod.Item2, true, theModType));
                }
            }
            foreach (var mod in task.CommonParameters.ListOfModsVariable)
            {
                var theModType = VariableModTypeForTreeViewObservableCollection.FirstOrDefault(b => b.DisplayName.Equals(mod.Item1));
                if (theModType != null)
                {
                    var theMod = theModType.Children.FirstOrDefault(b => b.DisplayName.Equals(mod.Item2));
                    if (theMod != null)
                    {
                        theMod.Use = true;
                    }
                    else
                    {
                        theModType.Children.Add(new ModForTreeView("UNKNOWN MODIFICATION!", true, mod.Item2, true, theModType));
                    }
                }
                else
                {
                    theModType = new ModTypeForTreeView(mod.Item1, true);
                    VariableModTypeForTreeViewObservableCollection.Add(theModType);
                    theModType.Children.Add(new ModForTreeView("UNKNOWN MODIFICATION!", true, mod.Item2, true, theModType));
                }
            }

            foreach (var heh in LocalizeModTypeForTreeViewObservableCollection)
            {
                heh.Use = false;
            }

            foreach (var ye in VariableModTypeForTreeViewObservableCollection)
            {
                ye.VerifyCheckState();
            }
            foreach (var ye in FixedModTypeForTreeViewObservableCollection)
            {
                ye.VerifyCheckState();
            }

            massDiffAcceptExact.IsChecked = task.SearchParameters.MassDiffAcceptorType == MassDiffAcceptorType.Exact;
            massDiffAccept1mm.IsChecked = task.SearchParameters.MassDiffAcceptorType == MassDiffAcceptorType.OneMM;
            massDiffAccept2mm.IsChecked = task.SearchParameters.MassDiffAcceptorType == MassDiffAcceptorType.TwoMM;
            massDiffAccept3mm.IsChecked = task.SearchParameters.MassDiffAcceptorType == MassDiffAcceptorType.ThreeMM;
            massDiffAccept187.IsChecked = task.SearchParameters.MassDiffAcceptorType == MassDiffAcceptorType.ModOpen;
            massDiffAcceptOpen.IsChecked = task.SearchParameters.MassDiffAcceptorType == MassDiffAcceptorType.Open;
            massDiffAcceptCustom.IsChecked = task.SearchParameters.MassDiffAcceptorType == MassDiffAcceptorType.Custom;

            if (task.SearchParameters.MassDiffAcceptorType == MassDiffAcceptorType.Custom)
            {
                customkMdacTextBox.Text = task.SearchParameters.CustomMdac;
            }
            writePrunedDBCheckBox.IsChecked = task.SearchParameters.WritePrunedDatabase;
            UpdateModSelectionGrid();
        }

        private void CancelButton_Click(object sender, RoutedEventArgs e)
        {
            DialogResult = false;
        }

        private void SaveButton_Click(object sender, RoutedEventArgs e)
        {
            CleavageSpecificity searchModeType = CleavageSpecificity.Full; //classic and modern by default
            if (semiSpecificSearchRadioButton.IsChecked.Value) //semi
            {
                searchModeType = CleavageSpecificity.Semi;
            }
            else if (nonSpecificSearchRadioButton.IsChecked.Value) //non
            {
                searchModeType = CleavageSpecificity.None;
            }
            //else it's the default of full

            if (searchModeType != CleavageSpecificity.Full)
            {
                if (((Protease)proteaseComboBox.SelectedItem).Name.Contains("non-specific"))
                {
                    if (cTerminalIons.IsChecked.Value)
                    {
                        while (!((Protease)proteaseComboBox.SelectedItem).Name.Equals("singleC"))
                        {
                            proteaseComboBox.Items.MoveCurrentToNext();
                            proteaseComboBox.SelectedItem = proteaseComboBox.Items.CurrentItem;
                        }
                    }
                    else //we're not allowing no ion types. It must have C if it doesn't have N.
                    {
                        while (!((Protease)proteaseComboBox.SelectedItem).Name.Equals("singleN"))
                        {
                            proteaseComboBox.Items.MoveCurrentToNext();
                            proteaseComboBox.SelectedItem = proteaseComboBox.Items.CurrentItem;
                        }
                    }
                    searchModeType = CleavageSpecificity.Full; //we're going to change this to override the semi, or the singleN/C proteases will be treated as semi instead of full
                }
                if (!addCompIonCheckBox.IsChecked.Value)
                {
                    MessageBox.Show("Warning: Complementary ions are strongly recommended when using this algorithm.");
                }
                //only use N or C termini, not both
                if(cTerminalIons.IsChecked.Value)
                {
                    nTerminalIons.IsChecked = false;
                }
                else
                {
                    nTerminalIons.IsChecked = true;
                }
            }

            if (!GlobalGuiSettings.CheckTaskSettingsValidity(precursorMassToleranceTextBox.Text, productMassToleranceTextBox.Text, missedCleavagesTextBox.Text,
                maxModificationIsoformsTextBox.Text, MinPeptideLengthTextBox.Text, MaxPeptideLengthTextBox.Text, maxThreadsTextBox.Text, minScoreAllowed.Text,
                peakFindingToleranceTextBox.Text, histogramBinWidthTextBox.Text, DeconvolutionMaxAssumedChargeStateTextBox.Text, TopNPeaksTextBox.Text,
                MinRatioTextBox.Text, numberOfDatabaseSearchesTextBox.Text, MaxModNumTextBox.Text, MaxFragmentMassTextBox.Text, QValueTextBox.Text))
            {
                return;
            }

            Protease protease = (Protease)proteaseComboBox.SelectedItem;

            DissociationType dissociationType = GlobalVariables.AllSupportedDissociationTypes[dissociationTypeComboBox.SelectedItem.ToString()];
            FragmentationTerminus fragmentationTerminus = FragmentationTerminus.Both;
            if (nTerminalIons.IsChecked.Value && !cTerminalIons.IsChecked.Value)
            {
                fragmentationTerminus = FragmentationTerminus.N;
            }
            else if (!nTerminalIons.IsChecked.Value && cTerminalIons.IsChecked.Value)
            {
                fragmentationTerminus = FragmentationTerminus.C;
            }
            else if (!nTerminalIons.IsChecked.Value && !cTerminalIons.IsChecked.Value) //why would you want this
            {
                fragmentationTerminus = FragmentationTerminus.None;
                MessageBox.Show("Warning: No ion types were selected. MetaMorpheus will be unable to search MS/MS spectra.");
            }
            //else both

            int maxMissedCleavages = string.IsNullOrEmpty(missedCleavagesTextBox.Text) ? int.MaxValue : (int.Parse(missedCleavagesTextBox.Text, NumberStyles.Any, CultureInfo.InvariantCulture));
            int minPeptideLengthValue = (int.Parse(MinPeptideLengthTextBox.Text, NumberStyles.Any, CultureInfo.InvariantCulture));
            int maxPeptideLengthValue = string.IsNullOrEmpty(MaxPeptideLengthTextBox.Text) ? int.MaxValue : (int.Parse(MaxPeptideLengthTextBox.Text, NumberStyles.Any, CultureInfo.InvariantCulture));
            int maxModificationIsoformsValue = (int.Parse(maxModificationIsoformsTextBox.Text, CultureInfo.InvariantCulture));
            int maxModsForPeptideValue = (int.Parse(MaxModNumTextBox.Text, CultureInfo.InvariantCulture));
            InitiatorMethionineBehavior initiatorMethionineBehavior = ((InitiatorMethionineBehavior)initiatorMethionineBehaviorComboBox.SelectedIndex);

            DigestionParams digestionParamsToSave = new DigestionParams(
                protease: protease.Name,
                maxMissedCleavages: maxMissedCleavages,
                minPeptideLength: minPeptideLengthValue,
                maxPeptideLength: maxPeptideLengthValue,
                maxModificationIsoforms: maxModificationIsoformsValue,
                initiatorMethionineBehavior: initiatorMethionineBehavior,
                maxModsForPeptides: maxModsForPeptideValue,
                searchModeType: searchModeType,
                fragmentationTerminus: fragmentationTerminus);

            Tolerance ProductMassTolerance;
            if (productMassToleranceComboBox.SelectedIndex == 0)
            {
                ProductMassTolerance = new AbsoluteTolerance(double.Parse(productMassToleranceTextBox.Text, CultureInfo.InvariantCulture));
            }
            else
            {
                ProductMassTolerance = new PpmTolerance(double.Parse(productMassToleranceTextBox.Text, CultureInfo.InvariantCulture));
            }

            Tolerance PrecursorMassTolerance;
            if (precursorMassToleranceComboBox.SelectedIndex == 0)
            {
                PrecursorMassTolerance = new AbsoluteTolerance(double.Parse(precursorMassToleranceTextBox.Text, CultureInfo.InvariantCulture));
            }
            else
            {
                PrecursorMassTolerance = new PpmTolerance(double.Parse(precursorMassToleranceTextBox.Text, CultureInfo.InvariantCulture));
            }
            TheTask.SearchParameters.MaxFragmentSize = Double.Parse(MaxFragmentMassTextBox.Text, CultureInfo.InvariantCulture);

            var listOfModsVariable = new List<(string, string)>();
            foreach (var heh in VariableModTypeForTreeViewObservableCollection)
            {
                listOfModsVariable.AddRange(heh.Children.Where(b => b.Use).Select(b => (b.Parent.DisplayName, b.DisplayName)));
            }

            var listOfModsFixed = new List<(string, string)>();
            foreach (var heh in FixedModTypeForTreeViewObservableCollection)
            {
                listOfModsFixed.AddRange(heh.Children.Where(b => b.Use).Select(b => (b.Parent.DisplayName, b.DisplayName)));
            }

            if (!GlobalGuiSettings.VariableModCheck(listOfModsVariable))
            {
                return;
            }

            bool TrimMs1Peaks = trimMs1.IsChecked.Value;
            bool TrimMsMsPeaks = trimMsMs.IsChecked.Value;
            int TopNpeaks = int.Parse(TopNPeaksTextBox.Text);
            double MinRatio = double.Parse(MinRatioTextBox.Text);

            bool parseMaxThreadsPerFile = !maxThreadsTextBox.Text.Equals("") && (int.Parse(maxThreadsTextBox.Text) <= Environment.ProcessorCount && int.Parse(maxThreadsTextBox.Text) > 0);

            CommonParameters commonParamsToSave = new CommonParameters(
                taskDescriptor: OutputFileNameTextBox.Text != "" ? OutputFileNameTextBox.Text : "SearchTask",
                maxThreadsToUsePerFile: parseMaxThreadsPerFile ? int.Parse(maxThreadsTextBox.Text, CultureInfo.InvariantCulture) : new CommonParameters().MaxThreadsToUsePerFile,
                useDeltaScore: deltaScoreCheckBox.IsChecked.Value,
                reportAllAmbiguity: allAmbiguity.IsChecked.Value,
                deconvolutionMaxAssumedChargeState: int.Parse(DeconvolutionMaxAssumedChargeStateTextBox.Text, CultureInfo.InvariantCulture),
                totalPartitions: int.Parse(numberOfDatabaseSearchesTextBox.Text, CultureInfo.InvariantCulture),
                doPrecursorDeconvolution: deconvolutePrecursors.IsChecked.Value,
                useProvidedPrecursorInfo: useProvidedPrecursor.IsChecked.Value,
                scoreCutoff: double.Parse(minScoreAllowed.Text, CultureInfo.InvariantCulture),
                calculateEValue: eValueCheckBox.IsChecked.Value,
                listOfModsFixed: listOfModsFixed,
                listOfModsVariable: listOfModsVariable,
                dissociationType: dissociationType,
                precursorMassTolerance: PrecursorMassTolerance,
                productMassTolerance: ProductMassTolerance,
                digestionParams: digestionParamsToSave,
                trimMs1Peaks: TrimMs1Peaks,
                trimMsMsPeaks: TrimMsMsPeaks,
                topNpeaks: TopNpeaks,
                minRatio: MinRatio,
                addCompIons: addCompIonCheckBox.IsChecked.Value,
                qValueOutputFilter: QValueCheckBox.IsChecked.Value ? double.Parse(QValueTextBox.Text, CultureInfo.InvariantCulture) : 1.0);

            if (classicSearchRadioButton.IsChecked.Value)
            {
                TheTask.SearchParameters.SearchType = SearchType.Classic;
            }
            else if (modernSearchRadioButton.IsChecked.Value)
            {
                TheTask.SearchParameters.SearchType = SearchType.Modern;
            }
            else //both semi and nonspecific are termed "nonspecific", because they both contain at least one nonspecific cleavage and they share the same algorithm
            {
                TheTask.SearchParameters.SearchType = SearchType.NonSpecific;
            }

            TheTask.SearchParameters.DoParsimony = checkBoxParsimony.IsChecked.Value;
            TheTask.SearchParameters.NoOneHitWonders = checkBoxNoOneHitWonders.IsChecked.Value;
            TheTask.SearchParameters.DoQuantification = checkBoxQuantification.IsChecked.Value;
            TheTask.SearchParameters.Normalize = checkBoxNormalize.IsChecked.Value;
            TheTask.SearchParameters.MatchBetweenRuns = checkBoxMatchBetweenRuns.IsChecked.Value;
            TheTask.SearchParameters.ModPeptidesAreDifferent = modPepsAreUnique.IsChecked.Value;
            TheTask.SearchParameters.QuantifyPpmTol = double.Parse(peakFindingToleranceTextBox.Text, CultureInfo.InvariantCulture);
            TheTask.SearchParameters.SearchTarget = checkBoxTarget.IsChecked.Value;
            TheTask.SearchParameters.WriteMzId = ckbMzId.IsChecked.Value;
            TheTask.SearchParameters.WriteDecoys = writeDecoyCheckBox.IsChecked.Value;
            TheTask.SearchParameters.WriteContaminants = writeContaminantCheckBox.IsChecked.Value;
            //TheTask.SearchParameters.OutPepXML = ckbPepXML.IsChecked.Value;

            if (checkBoxDecoy.IsChecked.Value)
            {
                if (radioButtonReverseDecoy.IsChecked.Value)
                {
                    TheTask.SearchParameters.DecoyType = DecoyType.Reverse;
                }
                else //if (radioButtonSlideDecoy.IsChecked.Value)
                {
                    TheTask.SearchParameters.DecoyType = DecoyType.Slide;
                }
            }
            else
            {
                TheTask.SearchParameters.DecoyType = DecoyType.None;
            }

            if (massDiffAcceptExact.IsChecked.HasValue && massDiffAcceptExact.IsChecked.Value)
            {
                TheTask.SearchParameters.MassDiffAcceptorType = MassDiffAcceptorType.Exact;
            }
            if (massDiffAccept1mm.IsChecked.HasValue && massDiffAccept1mm.IsChecked.Value)
            {
                TheTask.SearchParameters.MassDiffAcceptorType = MassDiffAcceptorType.OneMM;
            }
            if (massDiffAccept2mm.IsChecked.HasValue && massDiffAccept2mm.IsChecked.Value)
            {
                TheTask.SearchParameters.MassDiffAcceptorType = MassDiffAcceptorType.TwoMM;
            }
            if (massDiffAccept3mm.IsChecked.HasValue && massDiffAccept3mm.IsChecked.Value)
            {
                TheTask.SearchParameters.MassDiffAcceptorType = MassDiffAcceptorType.ThreeMM;
            }
            if (massDiffAccept187.IsChecked.HasValue && massDiffAccept187.IsChecked.Value)
            {
                TheTask.SearchParameters.MassDiffAcceptorType = MassDiffAcceptorType.ModOpen;
            }
            if (massDiffAcceptOpen.IsChecked.HasValue && massDiffAcceptOpen.IsChecked.Value)
            {
                TheTask.SearchParameters.MassDiffAcceptorType = MassDiffAcceptorType.Open;
            }
            if (massDiffAcceptCustom.IsChecked.HasValue && massDiffAcceptCustom.IsChecked.Value)
            {
                TheTask.SearchParameters.MassDiffAcceptorType = MassDiffAcceptorType.Custom;
                TheTask.SearchParameters.CustomMdac = customkMdacTextBox.Text;
            }

            //determine if semi or nonspecific with a specific protease. Full is already added by default.
            if(searchModeType == CleavageSpecificity.Semi || protease.CleavageSpecificity==CleavageSpecificity.Semi)
            {
                TheTask.SearchParameters.LocalFdrCategories.Add(FdrCategory.SemiSpecific);
            }
            else if(searchModeType==CleavageSpecificity.None && protease.CleavageSpecificity!=CleavageSpecificity.None)
            {
                TheTask.SearchParameters.LocalFdrCategories.Add(FdrCategory.SemiSpecific);
                TheTask.SearchParameters.LocalFdrCategories.Add(FdrCategory.NonSpecific);
            }
            // else do nothing


            // displays warning if classic search is enabled with an open search mode
            if (TheTask.SearchParameters.SearchType == SearchType.Classic &&
                (TheTask.SearchParameters.MassDiffAcceptorType == MassDiffAcceptorType.ModOpen || TheTask.SearchParameters.MassDiffAcceptorType == MassDiffAcceptorType.Open))
            {
                MessageBoxResult result = MessageBox.Show("Modern Search mode is recommended when conducting open precursor mass searches to reduce search time.\n\n" +
                    "Continue anyway?", "Modern search recommended", MessageBoxButton.OKCancel);

                if (result == MessageBoxResult.Cancel)
                {
                    return;
                }
            }

            TheTask.SearchParameters.DoHistogramAnalysis = checkBoxHistogramAnalysis.IsChecked.Value;
            TheTask.SearchParameters.HistogramBinTolInDaltons = double.Parse(histogramBinWidthTextBox.Text, CultureInfo.InvariantCulture);

            TheTask.SearchParameters.WritePrunedDatabase = writePrunedDBCheckBox.IsChecked.Value;

            SetModSelectionForPrunedDB();

            TheTask.CommonParameters = commonParamsToSave;

            DialogResult = true;
        }

        private void ApmdExpander_Collapsed(object sender, RoutedEventArgs e)
        {
            DataContextForSearchTaskWindow.ExpanderTitle = string.Join(", ", SearchModesForThisTask.Where(b => b.Use).Select(b => b.Name));
            DataContextForSearchTaskWindow.AnalysisExpanderTitle = "Some analysis properties...";
            DataContextForSearchTaskWindow.SearchModeExpanderTitle = "Some search properties...";
        }

        private void SetModSelectionForPrunedDB()
        {
            TheTask.SearchParameters.ModsToWriteSelection = new Dictionary<string, int>();
            //checks the grid values for which button is checked then sets paramaters accordingly
            foreach (var modTypeInGrid in ModSelectionGridItems)
            {
                if (modTypeInGrid.Item3)
                {
                    TheTask.SearchParameters.ModsToWriteSelection[modTypeInGrid.ModName] = 1;
                    continue;
                }
                if (modTypeInGrid.Item4)
                {
                    TheTask.SearchParameters.ModsToWriteSelection[modTypeInGrid.ModName] = 2;
                    continue;
                }
                if (modTypeInGrid.Item5)
                {
                    TheTask.SearchParameters.ModsToWriteSelection[modTypeInGrid.ModName] = 3;
                }
            }
        }

        private void UpdateModSelectionGrid()
        {
            foreach (var modType in TheTask.SearchParameters.ModsToWriteSelection)
            {
                var huhb = ModSelectionGridItems.FirstOrDefault(b => b.ModName == modType.Key);
                if (huhb != null)
                {
                    switch (modType.Value)
                    {
                        case (0):
                            huhb.Item2 = true;
                            huhb.Item3 = false;
                            huhb.Item4 = false;
                            huhb.Item5 = false;
                            break;

                        case (1):
                            huhb.Item2 = false;
                            huhb.Item3 = true;
                            huhb.Item4 = false;
                            huhb.Item5 = false;
                            break;

                        case (2):
                            huhb.Item2 = false;
                            huhb.Item3 = false;
                            huhb.Item4 = true;
                            huhb.Item5 = false;
                            break;

                        case (3):
                            huhb.Item2 = false;
                            huhb.Item3 = false;
                            huhb.Item4 = false;
                            huhb.Item5 = true;
                            break;
                    }
                }
            }
        }

        private void NonSpecificUsingNonSpecific(object sender, RoutedEventArgs e)
        {
            if (nonSpecificSearchRadioButton.IsChecked.Value)
            {
                proteaseComboBox.Items.MoveCurrentToFirst();
                proteaseComboBox.SelectedItem = proteaseComboBox.Items.CurrentItem;
                while (!((Protease)proteaseComboBox.SelectedItem).Name.Contains("non-specific"))
                {
                    proteaseComboBox.Items.MoveCurrentToNext();
                    proteaseComboBox.SelectedItem = proteaseComboBox.Items.CurrentItem;
                }
                addCompIonCheckBox.IsChecked = true;
            }
            else
            {
                addCompIonCheckBox.IsChecked = false;
                nTerminalIons.IsChecked = true;
                cTerminalIons.IsChecked = true;
            }
        }

        private void NonSpecificUpdate(object sender, TextChangedEventArgs e)
        {
            if (((Protease)proteaseComboBox.SelectedItem).Name.Contains("non-specific") || nonSpecificSearchRadioButton.IsChecked.Value)
            {
                try
                {
                    System.Windows.Controls.TextBox textBox = (TextBox)sender;
                    if (textBox.Name.Equals("MaxPeptideLengthTextBox")) //if maxPeptideLength was modified
                    {
                        if (!missedCleavagesTextBox.Text.Equals((Convert.ToInt32(MaxPeptideLengthTextBox.Text) - 1).ToString())) //prevents infinite loops
                        {
                            missedCleavagesTextBox.Text = (Convert.ToInt32(MaxPeptideLengthTextBox.Text) - 1).ToString();
                        }
                    }
                    else //if missedCleavagesTextBox was modified
                    {
                        if (!MaxPeptideLengthTextBox.Text.Equals((Convert.ToInt32(missedCleavagesTextBox.Text) + 1).ToString())) //prevents infinite loops
                        {
                            MaxPeptideLengthTextBox.Text = (Convert.ToInt32(missedCleavagesTextBox.Text) + 1).ToString();
                        }
                    }
                }
                catch
                {
                    //if not an entry, don't update the other box.
                }
            }
        }

        private void NonSpecificUpdate(object sender, SelectionChangedEventArgs e)
        {
            const int maxLength = 25;
            if (((Protease)proteaseComboBox.SelectedItem).Name.Contains("non-specific"))
            {
                MaxPeptideLengthTextBox.Text = maxLength.ToString();
            }
        }

        private void SemiSpecificUpdate(object sender, RoutedEventArgs e)
        {
            addCompIonCheckBox.IsChecked = semiSpecificSearchRadioButton.IsChecked.Value;
            if (semiSpecificSearchRadioButton.IsChecked.Value)
            {
                missedCleavagesTextBox.Text = "2";
                MaxPeptideLengthTextBox.Text = "50";
            }
            else
            {
                nTerminalIons.IsChecked = true;
                cTerminalIons.IsChecked = true;
            }
        }

        private void KeyPressed(object sender, KeyEventArgs e)
        {
            if (e.Key == Key.Return)
            {
                SaveButton_Click(sender, e);
            }
            else if (e.Key == Key.Escape)
            {
                CancelButton_Click(sender, e);
            }
        }

        private void dissassociationTypeComboBox_SelectionChanged(object sender, SelectionChangedEventArgs e)
        {

        }
    }

    public class DataContextForSearchTaskWindow : INotifyPropertyChanged
    {
        private string _ExpanderTitle;
        private string _SearchModeExpanderTitle;
        private string _ModExpanderTitle;
        private string _AnalysisExpanderTitle;

        public event PropertyChangedEventHandler PropertyChanged;

        public string ExpanderTitle
        {
            get
            { return _ExpanderTitle; }
            set
            {
                _ExpanderTitle = value;
                RaisePropertyChanged("ExpanderTitle");
            }
        }

        public string AnalysisExpanderTitle
        {
            get { return _AnalysisExpanderTitle; }
            set
            {
                _AnalysisExpanderTitle = value;
                RaisePropertyChanged("AnalysisExpanderTitle");
            }
        }

        public string ModExpanderTitle
        {
            get { return _ModExpanderTitle; }
            set
            {
                _ModExpanderTitle = value;
                RaisePropertyChanged("ModExpanderTitle");
            }
        }

        public string SearchModeExpanderTitle
        {
            get { return _SearchModeExpanderTitle; }
            set
            {
                _SearchModeExpanderTitle = value;
                RaisePropertyChanged("SearchModeExpanderTitle");
            }
        }

        protected void RaisePropertyChanged(string name)
        {
            PropertyChanged?.Invoke(this, new PropertyChangedEventArgs(name));
        }
    }
}
