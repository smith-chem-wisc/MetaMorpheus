﻿using EngineLayer;
using MassSpectrometry;
using MzLibUtil;
using Nett;
using Omics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.ComponentModel;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Input;
using TaskLayer;
using UsefulProteomicsDatabases;
using GuiFunctions;
using Proteomics;
using System.Threading.Tasks;
using Omics.Digestion;
using Omics.Modifications;

namespace MetaMorpheusGUI
{
    /// <summary>
    /// Interaction logic for SearchTaskWindow.xaml
    /// </summary>
    public partial class SearchTaskWindow : Window
    {
        private readonly DataContextForSearchTaskWindow DataContextForSearchTaskWindow;
        private readonly ObservableCollection<SearchModeForDataGrid> SearchModesForThisTask = new ObservableCollection<SearchModeForDataGrid>();
        private readonly ObservableCollection<ModTypeForTreeViewModel> FixedModTypeForTreeViewObservableCollection = new ObservableCollection<ModTypeForTreeViewModel>();
        private readonly ObservableCollection<ModTypeForTreeViewModel> VariableModTypeForTreeViewObservableCollection = new ObservableCollection<ModTypeForTreeViewModel>();
        private readonly ObservableCollection<ModTypeForLoc> LocalizeModTypeForTreeViewObservableCollection = new ObservableCollection<ModTypeForLoc>();
        private readonly ObservableCollection<ModTypeForGrid> ModSelectionGridItems = new ObservableCollection<ModTypeForGrid>();
        private readonly ObservableCollection<SilacInfoForDataGrid> StaticSilacLabelsObservableCollection = new ObservableCollection<SilacInfoForDataGrid>();
        private bool AutomaticallyAskAndOrUpdateParametersBasedOnProtease = true;
        private CustomFragmentationWindow CustomFragmentationWindow;
        private MassDifferenceAcceptorSelectionViewModel _massDifferenceAcceptorViewModel;
        private string _defaultMultiplexType = "TMT10";
        private DeconHostViewModel DeconHostViewModel;

        internal SearchTask TheTask { get; private set; }

        public SearchTaskWindow(SearchTask task)
        {
            InitializeComponent();
            TheTask = task ?? new SearchTask();

            AutomaticallyAskAndOrUpdateParametersBasedOnProtease = false;
            PopulateChoices();
            UpdateFieldsFromTask(TheTask);
            AutomaticallyAskAndOrUpdateParametersBasedOnProtease = true;
            DeisotopingControl.DataContext = DeconHostViewModel;
            MassDifferenceAcceptorControl.DataContext = _massDifferenceAcceptorViewModel;

            if (task == null)
            {
                this.saveButton.Content = "Add the Search Task";
            }

            DataContextForSearchTaskWindow = new DataContextForSearchTaskWindow
            {
                ExpanderTitle = string.Join(", ", SearchModesForThisTask.Where(b => b.Use).Select(b => b.Name)),
                AnalysisExpanderTitle = "Some analysis properties...",
                SearchModeExpanderTitle = "Some search properties..."
            };
            this.DataContext = DataContextForSearchTaskWindow;
            SearchModifications.Timer.Tick += new EventHandler(TextChangeTimerHandler);
            dataGridSilacLabels.DataContext = StaticSilacLabelsObservableCollection;
            base.Closing += this.OnClosing;
        }

        private void Row_DoubleClick(object sender, MouseButtonEventArgs e)
        {
            var ye = sender as DataGridCell;
            if (ye.Content is TextBlock hm && !string.IsNullOrEmpty(hm.Text))
            {
                GlobalVariables.StartProcess(hm.Text);
            }
        }

        private void PopulateChoices()
        {
            foreach (Protease protease in ProteaseDictionary.Dictionary.Values)
            {
                ProteaseComboBox.Items.Add(protease);
            }
            Protease trypsin = ProteaseDictionary.Dictionary["trypsin"];
            ProteaseComboBox.SelectedItem = trypsin;

            foreach (string initiatior_methionine_behavior in Enum.GetNames(typeof(InitiatorMethionineBehavior)))
            {
                InitiatorMethionineBehaviorComboBox.Items.Add(initiatior_methionine_behavior);
            }

            foreach (string dissassociationType in GlobalVariables.AllSupportedDissociationTypes.Keys)
            {
                DissociationTypeComboBox.Items.Add(dissassociationType);
            }


            foreach (string separationType in GlobalVariables.SeparationTypes)
            {
                SeparationTypeComboBox.Items.Add(separationType);
            }
            SeparationTypeComboBox.SelectedItem = "HPLC";


            ProductMassToleranceComboBox.Items.Add("Da");
            ProductMassToleranceComboBox.Items.Add("ppm");

            PrecursorMassToleranceComboBox.Items.Add("Da");
            PrecursorMassToleranceComboBox.Items.Add("ppm");

            foreach (var hm in GlobalVariables.AllModsKnown.Where(b => b.ValidModification == true).GroupBy(b => b.ModificationType))
            {
                var theModType = new ModTypeForGrid(hm.Key);
                ModSelectionGridItems.Add(theModType);
            }
            ModSelectionGrid.ItemsSource = ModSelectionGridItems;

            foreach (var hm in GlobalVariables.AllModsKnown.Where(b => b.ValidModification == true).GroupBy(b => b.ModificationType))
            {
                var theModType = new ModTypeForTreeViewModel(hm.Key, false);
                FixedModTypeForTreeViewObservableCollection.Add(theModType);
                foreach (var uah in hm)
                {
                    theModType.Children.Add(new ModForTreeViewModel(uah.ToString(), false, uah.IdWithMotif, false, theModType));
                }
            }
            FixedModsTreeView.DataContext = FixedModTypeForTreeViewObservableCollection;

            foreach (var hm in GlobalVariables.AllModsKnown.Where(b => b.ValidModification == true).GroupBy(b => b.ModificationType))
            {
                var theModType = new ModTypeForTreeViewModel(hm.Key, false);
                VariableModTypeForTreeViewObservableCollection.Add(theModType);
                foreach (var uah in hm)
                {
                    theModType.Children.Add(new ModForTreeViewModel(uah.ToString(), false, uah.IdWithMotif, false, theModType));
                }
            }
            VariableModsTreeView.DataContext = VariableModTypeForTreeViewObservableCollection;

            foreach (var hm in GlobalVariables.AllModsKnown.Where(b => b.ValidModification == true).GroupBy(b => b.ModificationType))
            {
                LocalizeModTypeForTreeViewObservableCollection.Add(new ModTypeForLoc(hm.Key));
            }

            foreach (string labelModType in GlobalVariables.AllModsKnown.Where(m => m.ModificationType.Equals("Multiplex Label"))
                         .Select(m => m.OriginalId).Distinct())
            {
                MultiplexComboBox.Items.Add(labelModType);
            }
            MultiplexComboBox.SelectedItem = _defaultMultiplexType;
        }

        /// <summary>
        /// Initializes the fields in the search task window upon opening to the settings of the param Task
        /// </summary>
        /// <param name="task"></param>
        private void UpdateFieldsFromTask(SearchTask task)
        {
            MetaMorpheusEngine.DetermineAnalyteType(TheTask.CommonParameters);
            if (task.CommonParameters.DigestionParams is DigestionParams digestionParams)
            {
                ProteaseComboBox.SelectedItem = digestionParams.SpecificProtease; //needs to be first, so nonspecific can override if necessary
                //do these in if statements so as not to trigger the change
                if (task.SearchParameters.SearchType == SearchType.NonSpecific && digestionParams.SearchModeType == CleavageSpecificity.None)
                {
                    NonSpecificSearchRadioButton.IsChecked = true; //when this is changed it overrides the protease
                    if (digestionParams.SpecificProtease.Name.Equals("singleC") || digestionParams.SpecificProtease.Name.Equals("singleN"))
                    {
                        ProteaseComboBox.SelectedItem = ProteaseDictionary.Dictionary["non-specific"];
                    }
                    else
                    {
                        ProteaseComboBox.SelectedItem = digestionParams.SpecificProtease;
                    }
                }

                //If SILAC turnover
                if (task.SearchParameters.StartTurnoverLabel != null || task.SearchParameters.EndTurnoverLabel != null)
                {
                    task.SearchParameters.SilacLabels = null; //reset if between runs
                    CheckBoxSILAC.IsChecked = true;
                    var startLabel = task.SearchParameters.StartTurnoverLabel;
                    if (startLabel != null)
                    {
                        SilacInfoForDataGrid infoToAdd = new SilacInfoForDataGrid(startLabel, SilacModificationWindow.ExperimentType.Start);
                        if (startLabel.AdditionalLabels != null)
                        {
                            foreach (SilacLabel additionalLabel in startLabel.AdditionalLabels)
                            {
                                infoToAdd.AddAdditionalLabel(new SilacInfoForDataGrid(additionalLabel, SilacModificationWindow.ExperimentType.Start));
                            }
                        }
                        StaticSilacLabelsObservableCollection.Add(infoToAdd);
                    }
                    else //it's unlabeled for the start condition
                    {
                        StaticSilacLabelsObservableCollection.Add(new SilacInfoForDataGrid(SilacModificationWindow.ExperimentType.Start));
                    }
                    var endLabel = task.SearchParameters.EndTurnoverLabel;
                    if (endLabel != null)
                    {
                        SilacInfoForDataGrid infoToAdd = new SilacInfoForDataGrid(endLabel, SilacModificationWindow.ExperimentType.End);
                        if (endLabel.AdditionalLabels != null)
                        {
                            foreach (SilacLabel additionalLabel in endLabel.AdditionalLabels)
                            {
                                infoToAdd.AddAdditionalLabel(new SilacInfoForDataGrid(additionalLabel, SilacModificationWindow.ExperimentType.End));
                            }
                        }
                        StaticSilacLabelsObservableCollection.Add(infoToAdd);
                    }
                    else //it's unlabeled for the end condition
                    {
                        StaticSilacLabelsObservableCollection.Add(new SilacInfoForDataGrid(SilacModificationWindow.ExperimentType.End));
                    }
                }
                //else if SILAC multiplex
                else if (task.SearchParameters.SilacLabels != null && task.SearchParameters.SilacLabels.Count != 0)
                {
                    CheckBoxSILAC.IsChecked = true;
                    List<SilacLabel> labels = task.SearchParameters.SilacLabels;
                    foreach (SilacLabel label in labels)
                    {
                        SilacInfoForDataGrid infoToAdd = new SilacInfoForDataGrid(label, SilacModificationWindow.ExperimentType.Multiplex);
                        if (label.AdditionalLabels != null)
                        {
                            foreach (SilacLabel additionalLabel in label.AdditionalLabels)
                            {
                                infoToAdd.AddAdditionalLabel(new SilacInfoForDataGrid(additionalLabel, SilacModificationWindow.ExperimentType.Multiplex));
                            }
                        }
                        StaticSilacLabelsObservableCollection.Add(infoToAdd);
                    }
                    if (digestionParams.GeneratehUnlabeledProteinsForSilac)
                    {
                        StaticSilacLabelsObservableCollection.Add(new SilacInfoForDataGrid(SilacModificationWindow.ExperimentType.Multiplex));
                    }
                }

                CheckBoxQuantifyUnlabeledForSilac.IsChecked = digestionParams.GeneratehUnlabeledProteinsForSilac;
                InitiatorMethionineBehaviorComboBox.SelectedIndex = (int)digestionParams.InitiatorMethionineBehavior;
            }
            ClassicSearchRadioButton.IsChecked = task.SearchParameters.SearchType == SearchType.Classic;
            ModernSearchRadioButton.IsChecked = task.SearchParameters.SearchType == SearchType.Modern;
            
            if (task.SearchParameters.SearchType == SearchType.NonSpecific && task.CommonParameters.DigestionParams.SearchModeType != CleavageSpecificity.None)
            {
                SemiSpecificSearchRadioButton.IsChecked = true;
            }
            CheckBoxParsimony.IsChecked = task.SearchParameters.DoParsimony;
            CheckBoxNoOneHitWonders.IsChecked = task.SearchParameters.NoOneHitWonders;
            CheckBoxNoQuant.IsChecked = !task.SearchParameters.DoLabelFreeQuantification;
            CheckBoxLFQ.IsChecked = task.SearchParameters.DoLabelFreeQuantification;
            CheckBoxUseSharedPeptidesForQuant.IsChecked = task.SearchParameters.UseSharedPeptidesForLFQ;
            CheckBoxMultiplex.IsChecked = task.SearchParameters.DoMultiplexQuantification;
            MultiplexComboBox.SelectedItem = task.SearchParameters.MultiplexModId ?? _defaultMultiplexType;
            // If Spectral Recovery is enabled
            if (task.SearchParameters.WriteSpectralLibrary & task.SearchParameters.MatchBetweenRuns)
            {
                CheckBoxLFQwSpectralRecovery.IsChecked = task.SearchParameters.DoSpectralRecovery;
            }
            
            PeakFindingToleranceTextBox.Text = task.SearchParameters.QuantifyPpmTol.ToString(CultureInfo.InvariantCulture);
            CheckBoxMatchBetweenRuns.IsChecked = task.SearchParameters.MatchBetweenRuns;
            MbrFdrThresholdTextBox.Text = task.SearchParameters.MbrFdrThreshold.ToString(CultureInfo.InvariantCulture);
            CheckBoxNormalize.IsChecked = task.SearchParameters.Normalize;
            ModPepsAreUnique.IsChecked = task.SearchParameters.ModPeptidesAreDifferent;
            CheckBoxHistogramAnalysis.IsChecked = task.SearchParameters.DoHistogramAnalysis;
            HistogramBinWidthTextBox.Text = task.SearchParameters.HistogramBinTolInDaltons.ToString(CultureInfo.InvariantCulture);
            CheckBoxTarget.IsChecked = task.SearchParameters.SearchTarget;
            CheckBoxDecoy.IsChecked = task.SearchParameters.DecoyType != DecoyType.None;
            RadioButtonReverseDecoy.IsChecked = task.SearchParameters.DecoyType == DecoyType.Reverse;
            RadioButtonSlideDecoy.IsChecked = task.SearchParameters.DecoyType == DecoyType.Slide;
            MissedCleavagesTextBox.Text = task.CommonParameters.DigestionParams.MaxMissedCleavages == int.MaxValue ? "" : task.CommonParameters.DigestionParams.MaxMissedCleavages.ToString(CultureInfo.InvariantCulture);
            MinPeptideLengthTextBox.Text = task.CommonParameters.DigestionParams.MinLength.ToString(CultureInfo.InvariantCulture);
            MaxPeptideLengthTextBox.Text = task.CommonParameters.DigestionParams.MaxLength == int.MaxValue ? "" : task.CommonParameters.DigestionParams.MaxLength.ToString(CultureInfo.InvariantCulture);
            MaxFragmentMassTextBox.Text = task.SearchParameters.MaxFragmentSize.ToString(CultureInfo.InvariantCulture); //put after max peptide length to allow for override of auto
            maxModificationIsoformsTextBox.Text = task.CommonParameters.DigestionParams.MaxModificationIsoforms.ToString(CultureInfo.InvariantCulture);
            MaxModNumTextBox.Text = task.CommonParameters.DigestionParams.MaxMods.ToString(CultureInfo.InvariantCulture);
            DissociationTypeComboBox.SelectedItem = task.CommonParameters.DissociationType.ToString();
            SeparationTypeComboBox.SelectedItem = task.CommonParameters.SeparationType.ToString();
            NTerminalIons.IsChecked = task.CommonParameters.DigestionParams.FragmentationTerminus == FragmentationTerminus.Both || task.CommonParameters.DigestionParams.FragmentationTerminus == FragmentationTerminus.N;
            CTerminalIons.IsChecked = task.CommonParameters.DigestionParams.FragmentationTerminus == FragmentationTerminus.Both || task.CommonParameters.DigestionParams.FragmentationTerminus == FragmentationTerminus.C;
            InternalIonsCheckBox.IsChecked = task.SearchParameters.MinAllowedInternalFragmentLength != 0;
            MinInternalFragmentLengthTextBox.Text = task.SearchParameters.MinAllowedInternalFragmentLength.ToString();
            ProductMassToleranceTextBox.Text = task.CommonParameters.ProductMassTolerance.Value.ToString(CultureInfo.InvariantCulture);
            ProductMassToleranceComboBox.SelectedIndex = task.CommonParameters.ProductMassTolerance is AbsoluteTolerance ? 0 : 1;
            PrecursorMassToleranceTextBox.Text = task.CommonParameters.PrecursorMassTolerance.Value.ToString(CultureInfo.InvariantCulture);
            PrecursorMassToleranceComboBox.SelectedIndex = task.CommonParameters.PrecursorMassTolerance is AbsoluteTolerance ? 0 : 1;
            AddCompIonCheckBox.IsChecked = task.CommonParameters.AddCompIons;
            NumberOfDatabaseSearchesTextBox.Text = task.CommonParameters.TotalPartitions.ToString(CultureInfo.InvariantCulture);
            RemoveContaminantRadioBox.IsChecked = task.SearchParameters.TCAmbiguity == TargetContaminantAmbiguity.RemoveContaminant;
            RemoveTargetRadioBox.IsChecked = task.SearchParameters.TCAmbiguity == TargetContaminantAmbiguity.RemoveTarget;
            RenameTCProteinsRadioBox.IsChecked = task.SearchParameters.TCAmbiguity == TargetContaminantAmbiguity.RenameProtein;
            AllAmbiguity.IsChecked = task.CommonParameters.ReportAllAmbiguity;
            MinScoreAllowed.Text = task.CommonParameters.ScoreCutoff.ToString(CultureInfo.InvariantCulture);
            TrimMs1.IsChecked = task.CommonParameters.TrimMs1Peaks;
            TrimMsMs.IsChecked = task.CommonParameters.TrimMsMsPeaks;
            AddTruncationsCheckBox.IsChecked = task.CommonParameters.AddTruncations;

            DeconHostViewModel = new DeconHostViewModel(TheTask.CommonParameters.PrecursorDeconvolutionParameters,
                TheTask.CommonParameters.ProductDeconvolutionParameters,
                TheTask.CommonParameters.UseProvidedPrecursorInfo, TheTask.CommonParameters.DoPrecursorDeconvolution);
            DeisotopingControl.DataContext = DeconHostViewModel;

            NumberOfPeaksToKeepPerWindowTextBox.Text = task.CommonParameters.NumberOfPeaksToKeepPerWindow == int.MaxValue || !task.CommonParameters.NumberOfPeaksToKeepPerWindow.HasValue ? "" : task.CommonParameters.NumberOfPeaksToKeepPerWindow.Value.ToString(CultureInfo.InvariantCulture);
            MinimumAllowedIntensityRatioToBasePeakTexBox.Text = task.CommonParameters.MinimumAllowedIntensityRatioToBasePeak == double.MaxValue || !task.CommonParameters.MinimumAllowedIntensityRatioToBasePeak.HasValue ? "" : task.CommonParameters.MinimumAllowedIntensityRatioToBasePeak.Value.ToString(CultureInfo.InvariantCulture);
            WindowWidthThomsonsTextBox.Text = task.CommonParameters.WindowWidthThomsons == double.MaxValue || !task.CommonParameters.WindowWidthThomsons.HasValue ? "" : task.CommonParameters.WindowWidthThomsons.Value.ToString(CultureInfo.InvariantCulture);
            NumberOfWindowsTextBox.Text = task.CommonParameters.NumberOfWindows == int.MaxValue || !task.CommonParameters.NumberOfWindows.HasValue ? "" : task.CommonParameters.NumberOfWindows.Value.ToString(CultureInfo.InvariantCulture);
            NormalizePeaksInWindowCheckBox.IsChecked = task.CommonParameters.NormalizePeaksAccrossAllWindows;

            MaxThreadsTextBox.Text = task.CommonParameters.MaxThreadsToUsePerFile.ToString(CultureInfo.InvariantCulture);
            MinVariantDepthTextBox.Text = task.CommonParameters.MinVariantDepth.ToString(CultureInfo.InvariantCulture);
            MaxHeterozygousVariantsTextBox.Text = task.CommonParameters.MaxHeterozygousVariants.ToString(CultureInfo.InvariantCulture);
            CustomFragmentationWindow = new CustomFragmentationWindow(task.CommonParameters.CustomIons);

            if (task.CommonParameters.QValueThreshold < 1)
            {
                QValueThresholdTextBox.Text = task.CommonParameters.QValueThreshold.ToString(CultureInfo.InvariantCulture);
            }
            else
            {
                QValueThresholdTextBox.Text = "0.01";
            }

            if (task.CommonParameters.PepQValueThreshold < 1)
            {
                PepQValueThresholdTextBox.Text = task.CommonParameters.PepQValueThreshold.ToString(CultureInfo.InvariantCulture);
                PepQValueThresholdCheckbox.IsChecked = true;
            }
            else
            {
                PepQValueThresholdTextBox.Text = "0.01";
            }

            OutputFileNameTextBox.Text = task.CommonParameters.TaskDescriptor;
            CkbMzId.IsChecked = task.SearchParameters.WriteMzId;
            WriteDigestCountCheckBox.IsChecked = task.SearchParameters.WriteDigestionProductCountFile;
            WriteHighQPsmsCheckBox.IsChecked = task.SearchParameters.WriteHighQValuePsms;
            WriteDecoyCheckBox.IsChecked = task.SearchParameters.WriteDecoys;
            WriteContaminantCheckBox.IsChecked = task.SearchParameters.WriteContaminants;
            WriteIndividualResultsCheckBox.IsChecked = task.SearchParameters.WriteIndividualFiles;
            WriteSpectralLibraryCheckBox.IsChecked = task.SearchParameters.WriteSpectralLibrary;
            UpdateSpectralLibraryCheckBox.IsChecked = task.SearchParameters.UpdateSpectralLibrary;
            CompressIndividualResultsCheckBox.IsChecked = task.SearchParameters.CompressIndividualFiles;
            IncludeMotifInModNamesCheckBox.IsChecked = task.SearchParameters.IncludeModMotifInMzid;

            foreach (var mod in task.CommonParameters.ListOfModsFixed)
            {
                var theModType = FixedModTypeForTreeViewObservableCollection.FirstOrDefault(b => b.DisplayName.Equals(mod.Item1));
                if (theModType != null)
                {
                    var theMod = theModType.Children.FirstOrDefault(b => b.ModName.Equals(mod.Item2));
                    if (theMod != null)
                    {
                        theMod.Use = true;
                    }
                    else
                    {
                        theModType.Children.Add(new ModForTreeViewModel("UNKNOWN MODIFICATION!", true, mod.Item2, true, theModType));
                    }
                }
                else
                {
                    theModType = new ModTypeForTreeViewModel(mod.Item1, true);
                    FixedModTypeForTreeViewObservableCollection.Add(theModType);
                    theModType.Children.Add(new ModForTreeViewModel("UNKNOWN MODIFICATION!", true, mod.Item2, true, theModType));
                }
            }
            foreach (var mod in task.CommonParameters.ListOfModsVariable)
            {
                var theModType = VariableModTypeForTreeViewObservableCollection.FirstOrDefault(b => b.DisplayName.Equals(mod.Item1));
                if (theModType != null)
                {
                    var theMod = theModType.Children.FirstOrDefault(b => b.ModName.Equals(mod.Item2));
                    if (theMod != null)
                    {
                        theMod.Use = true;
                    }
                    else
                    {
                        theModType.Children.Add(new ModForTreeViewModel("UNKNOWN MODIFICATION!", true, mod.Item2, true, theModType));
                    }
                }
                else
                {
                    theModType = new ModTypeForTreeViewModel(mod.Item1, true);
                    VariableModTypeForTreeViewObservableCollection.Add(theModType);
                    theModType.Children.Add(new ModForTreeViewModel("UNKNOWN MODIFICATION!", true, mod.Item2, true, theModType));
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

            _massDifferenceAcceptorViewModel = new(task.SearchParameters.MassDiffAcceptorType, task.SearchParameters.CustomMdac, task.CommonParameters.PrecursorMassTolerance.Value);
            WritePrunedDBCheckBox.IsChecked = task.SearchParameters.WritePrunedDatabase;
            UpdateModSelectionGrid();
        }

        private void CancelButton_Click(object sender, RoutedEventArgs e)
        {
            DialogResult = false;
            CustomFragmentationWindow.Close();
        }

        private void SaveButton_Click(object sender, RoutedEventArgs e)
        {
            CleavageSpecificity searchModeType = GetSearchModeType(); //change search type to semi or non if selected
            SnesUpdates(searchModeType); //decide on singleN/C, make comp ion changes

            if (!GlobalGuiSettings.CheckTaskSettingsValidity(
                PrecursorMassToleranceTextBox.Text, 
                ProductMassToleranceTextBox.Text, 
                MissedCleavagesTextBox.Text,
                maxModificationIsoformsTextBox.Text, 
                MinPeptideLengthTextBox.Text, 
                MaxPeptideLengthTextBox.Text,
                MaxThreadsTextBox.Text, 
                MinScoreAllowed.Text,
                PeakFindingToleranceTextBox.Text, 
                MbrFdrThresholdTextBox.Text,
                HistogramBinWidthTextBox.Text, 
                DeconHostViewModel.PrecursorDeconvolutionParameters.MaxAssumedChargeState.ToString(), 
                NumberOfPeaksToKeepPerWindowTextBox.Text,
                MinimumAllowedIntensityRatioToBasePeakTexBox.Text, 
                WindowWidthThomsonsTextBox.Text, 
                NumberOfWindowsTextBox.Text, 
                NumberOfDatabaseSearchesTextBox.Text, 
                MaxModNumTextBox.Text, 
                MaxFragmentMassTextBox.Text, 
                QValueThresholdTextBox.Text, 
                PepQValueThresholdTextBox.Text, 
                InternalIonsCheckBox.IsChecked.Value ? MinInternalFragmentLengthTextBox.Text : null))
            {
                return;
            }

            Protease protease = (Protease)ProteaseComboBox.SelectedItem;

            string separationType = SeparationTypeComboBox.SelectedItem.ToString();

            DissociationType dissociationType = GlobalVariables.AllSupportedDissociationTypes[DissociationTypeComboBox.SelectedItem.ToString()];
            CustomFragmentationWindow.Close();

            FragmentationTerminus fragmentationTerminus = GetFragmentationTerminus();

            if (CheckBoxSILAC.IsChecked.Value)
            {
                SilacUpdates(out string silacError);
                if (silacError.Length != 0)
                {
                    MessageBox.Show(silacError);
                    return;
                }
            }
            else
            {
                TheTask.SearchParameters.SilacLabels = null;
                TheTask.SearchParameters.StartTurnoverLabel = null;
                TheTask.SearchParameters.EndTurnoverLabel = null;
            }

            int maxMissedCleavages = string.IsNullOrEmpty(MissedCleavagesTextBox.Text) ? int.MaxValue : (int.Parse(MissedCleavagesTextBox.Text, NumberStyles.Any, CultureInfo.InvariantCulture));
            int minPeptideLengthValue = (int.Parse(MinPeptideLengthTextBox.Text, NumberStyles.Any, CultureInfo.InvariantCulture));
            int maxPeptideLengthValue = string.IsNullOrEmpty(MaxPeptideLengthTextBox.Text) ? int.MaxValue : (int.Parse(MaxPeptideLengthTextBox.Text, NumberStyles.Any, CultureInfo.InvariantCulture));
            int MinVariantDepth = int.Parse(MinVariantDepthTextBox.Text, NumberStyles.Any, CultureInfo.InvariantCulture);
            int MaxHeterozygousVariants = int.Parse(MaxHeterozygousVariantsTextBox.Text, NumberStyles.Any, CultureInfo.InvariantCulture);
            int maxModificationIsoformsValue = (int.Parse(maxModificationIsoformsTextBox.Text, CultureInfo.InvariantCulture));
            int maxModsForPeptideValue = (int.Parse(MaxModNumTextBox.Text, CultureInfo.InvariantCulture));
            InitiatorMethionineBehavior initiatorMethionineBehavior = ((InitiatorMethionineBehavior)InitiatorMethionineBehaviorComboBox.SelectedIndex);

            DigestionParams digestionParamsToSave = new DigestionParams(
                protease: protease.Name,
                maxMissedCleavages: maxMissedCleavages,
                minPeptideLength: minPeptideLengthValue,
                maxPeptideLength: maxPeptideLengthValue,
                maxModificationIsoforms: maxModificationIsoformsValue,
                initiatorMethionineBehavior: initiatorMethionineBehavior,
                maxModsForPeptides: maxModsForPeptideValue,
                searchModeType: searchModeType,
                fragmentationTerminus: fragmentationTerminus,
                generateUnlabeledProteinsForSilac: CheckBoxQuantifyUnlabeledForSilac.IsChecked.Value);

            Tolerance ProductMassTolerance;
            if (ProductMassToleranceComboBox.SelectedIndex == 0)
            {
                ProductMassTolerance = new AbsoluteTolerance(double.Parse(ProductMassToleranceTextBox.Text, CultureInfo.InvariantCulture));
            }
            else
            {
                ProductMassTolerance = new PpmTolerance(double.Parse(ProductMassToleranceTextBox.Text, CultureInfo.InvariantCulture));
            }

            Tolerance PrecursorMassTolerance;
            if (PrecursorMassToleranceComboBox.SelectedIndex == 0)
            {
                PrecursorMassTolerance = new AbsoluteTolerance(double.Parse(PrecursorMassToleranceTextBox.Text, CultureInfo.InvariantCulture));
            }
            else
            {
                PrecursorMassTolerance = new PpmTolerance(double.Parse(PrecursorMassToleranceTextBox.Text, CultureInfo.InvariantCulture));
            }
            TheTask.SearchParameters.MaxFragmentSize = double.Parse(MaxFragmentMassTextBox.Text, CultureInfo.InvariantCulture);

            var listOfModsVariable = new List<(string, string)>();
            foreach (var heh in VariableModTypeForTreeViewObservableCollection)
            {
                listOfModsVariable.AddRange(heh.Children.Where(b => b.Use).Select(b => (b.Parent.DisplayName, b.ModName)));
            }

            var listOfModsFixed = new List<(string, string)>();
            foreach (var heh in FixedModTypeForTreeViewObservableCollection)
            {
                listOfModsFixed.AddRange(heh.Children.Where(b => b.Use).Select(b => (b.Parent.DisplayName, b.ModName)));
            }

            if (!GlobalGuiSettings.VariableModCheck(listOfModsVariable))
            {
                return;
            }

            bool TrimMs1Peaks = TrimMs1.IsChecked.Value;
            bool TrimMsMsPeaks = TrimMsMs.IsChecked.Value;
            bool AddTruncations = AddTruncationsCheckBox.IsChecked.Value;

            int? numPeaksToKeep = null;
            if (int.TryParse(NumberOfPeaksToKeepPerWindowTextBox.Text, out int numberOfPeaksToKeeep))
            {
                numPeaksToKeep = numberOfPeaksToKeeep;
            }

            double? minimumAllowedIntensityRatioToBasePeak = null;
            if (double.TryParse(MinimumAllowedIntensityRatioToBasePeakTexBox.Text, NumberStyles.Any, CultureInfo.InvariantCulture, out double minimumAllowedIntensityRatio))
            {
                minimumAllowedIntensityRatioToBasePeak = minimumAllowedIntensityRatio;
            }

            double? windowWidthThompsons = null;
            if (double.TryParse(WindowWidthThomsonsTextBox.Text, NumberStyles.Any, CultureInfo.InvariantCulture, out double windowWidth))
            {
                windowWidthThompsons = windowWidth;
            }

            int? numberOfWindows = null;
            if (int.TryParse(NumberOfWindowsTextBox.Text, out int numWindows))
            {
                numberOfWindows = numWindows;
            }

            bool normalizePeaksAccrossAllWindows = NormalizePeaksInWindowCheckBox.IsChecked.Value;

            bool parseMaxThreadsPerFile = !MaxThreadsTextBox.Text.Equals("") && (int.Parse(MaxThreadsTextBox.Text) <= Environment.ProcessorCount && int.Parse(MaxThreadsTextBox.Text) > 0);


            DeconvolutionParameters precursorDeconvolutionParameters = DeconHostViewModel.PrecursorDeconvolutionParameters.Parameters;
            DeconvolutionParameters productDeconvolutionParameters = DeconHostViewModel.ProductDeconvolutionParameters.Parameters;
            bool useProvidedPrecursorInfo = DeconHostViewModel.UseProvidedPrecursors;
            bool doPrecursorDeconvolution = DeconHostViewModel.DoPrecursorDeconvolution;
            
            CommonParameters commonParamsToSave = new CommonParameters(
                taskDescriptor: OutputFileNameTextBox.Text != "" ? OutputFileNameTextBox.Text : "SearchTask",
                maxThreadsToUsePerFile: parseMaxThreadsPerFile ? int.Parse(MaxThreadsTextBox.Text, CultureInfo.InvariantCulture) : new CommonParameters().MaxThreadsToUsePerFile,
                reportAllAmbiguity: AllAmbiguity.IsChecked.Value,
                totalPartitions: int.Parse(NumberOfDatabaseSearchesTextBox.Text, CultureInfo.InvariantCulture),
                doPrecursorDeconvolution: doPrecursorDeconvolution,
                useProvidedPrecursorInfo: useProvidedPrecursorInfo,
                qValueThreshold: !PepQValueThresholdCheckbox.IsChecked.Value ? double.Parse(QValueThresholdTextBox.Text, CultureInfo.InvariantCulture) : 1.0,
                pepQValueThreshold: PepQValueThresholdCheckbox.IsChecked.Value ? double.Parse(PepQValueThresholdTextBox.Text, CultureInfo.InvariantCulture) : 1.0,
                scoreCutoff: double.Parse(MinScoreAllowed.Text, CultureInfo.InvariantCulture),
                listOfModsFixed: listOfModsFixed,
                listOfModsVariable: listOfModsVariable,
                dissociationType: dissociationType,
                precursorMassTolerance: PrecursorMassTolerance,
                productMassTolerance: ProductMassTolerance,
                digestionParams: digestionParamsToSave,
                separationType: separationType,
                trimMs1Peaks: TrimMs1Peaks,
                trimMsMsPeaks: TrimMsMsPeaks,
                addTruncations: AddTruncations,
                numberOfPeaksToKeepPerWindow: numPeaksToKeep,
                minimumAllowedIntensityRatioToBasePeak: minimumAllowedIntensityRatioToBasePeak,
                windowWidthThomsons: windowWidthThompsons,
                numberOfWindows: numberOfWindows,//maybe change this some day
                normalizePeaksAccrossAllWindows: normalizePeaksAccrossAllWindows,//maybe change this some day
                addCompIons: AddCompIonCheckBox.IsChecked.Value,
                assumeOrphanPeaksAreZ1Fragments: protease.Name != "top-down",
                minVariantDepth: MinVariantDepth,
                maxHeterozygousVariants: MaxHeterozygousVariants,
                precursorDeconParams: precursorDeconvolutionParameters,
                productDeconParams: productDeconvolutionParameters);

            if (ClassicSearchRadioButton.IsChecked.Value)
            {
                TheTask.SearchParameters.SearchType = SearchType.Classic;
            }
            else if (ModernSearchRadioButton.IsChecked.Value)
            {
                TheTask.SearchParameters.SearchType = SearchType.Modern;
            }
            else //both semi and nonspecific are termed "nonspecific", because they both contain at least one nonspecific cleavage and they share the same algorithm
            {
                TheTask.SearchParameters.SearchType = SearchType.NonSpecific;
            }

            if (TheTask.SearchParameters.SearchType != SearchType.Classic && dissociationType == DissociationType.Autodetect)
            {
                MessageBox.Show("Autodetection of dissociation type from scan headers is only available for classic search. Please choose a different dissociation type or search mode");
                return;
            }

            TheTask.SearchParameters.MinAllowedInternalFragmentLength = InternalIonsCheckBox.IsChecked.Value ? Convert.ToInt32(MinInternalFragmentLengthTextBox.Text) : 0;
            TheTask.SearchParameters.DoParsimony = CheckBoxParsimony.IsChecked.Value;
            TheTask.SearchParameters.NoOneHitWonders = CheckBoxNoOneHitWonders.IsChecked.Value;
            TheTask.SearchParameters.DoLabelFreeQuantification = !CheckBoxNoQuant.IsChecked.Value;
            TheTask.SearchParameters.UseSharedPeptidesForLFQ = CheckBoxUseSharedPeptidesForQuant.IsChecked.Value;
            TheTask.SearchParameters.DoSpectralRecovery = CheckBoxLFQwSpectralRecovery.IsChecked.Value;
            TheTask.SearchParameters.DoMultiplexQuantification = CheckBoxMultiplex.IsChecked.Value;
            TheTask.SearchParameters.MultiplexModId = (string)MultiplexComboBox.SelectedItem;
            TheTask.SearchParameters.Normalize = CheckBoxNormalize.IsChecked.Value;
            TheTask.SearchParameters.MatchBetweenRuns = CheckBoxMatchBetweenRuns.IsChecked.Value;
            TheTask.SearchParameters.MbrFdrThreshold = double.Parse(MbrFdrThresholdTextBox.Text, CultureInfo.InvariantCulture);
            TheTask.SearchParameters.ModPeptidesAreDifferent = ModPepsAreUnique.IsChecked.Value;
            TheTask.SearchParameters.QuantifyPpmTol = double.Parse(PeakFindingToleranceTextBox.Text, CultureInfo.InvariantCulture);
            TheTask.SearchParameters.SearchTarget = CheckBoxTarget.IsChecked.Value;
            TheTask.SearchParameters.WriteMzId = CkbMzId.IsChecked.Value;
            TheTask.SearchParameters.WriteHighQValuePsms = WriteHighQPsmsCheckBox.IsChecked.Value;
            TheTask.SearchParameters.WriteDecoys = WriteDecoyCheckBox.IsChecked.Value;
            TheTask.SearchParameters.WriteContaminants = WriteContaminantCheckBox.IsChecked.Value;
            TheTask.SearchParameters.WriteIndividualFiles = WriteIndividualResultsCheckBox.IsChecked.Value;
            TheTask.SearchParameters.WriteSpectralLibrary = WriteSpectralLibraryCheckBox.IsChecked.Value;
            TheTask.SearchParameters.UpdateSpectralLibrary = UpdateSpectralLibraryCheckBox.IsChecked.Value;
            TheTask.SearchParameters.CompressIndividualFiles = CompressIndividualResultsCheckBox.IsChecked.Value;
            TheTask.SearchParameters.IncludeModMotifInMzid = IncludeMotifInModNamesCheckBox.IsChecked.Value;
            TheTask.SearchParameters.WriteDigestionProductCountFile = WriteDigestCountCheckBox.IsChecked.Value;

            if (RemoveContaminantRadioBox.IsChecked.Value)
            {
                TheTask.SearchParameters.TCAmbiguity = TargetContaminantAmbiguity.RemoveContaminant;
            }
            else if (RemoveTargetRadioBox.IsChecked.Value)
            {
                TheTask.SearchParameters.TCAmbiguity = TargetContaminantAmbiguity.RemoveTarget;
            }
            else //RenameTCProteinsRadioBox.IsChecked.Value
            {
                TheTask.SearchParameters.TCAmbiguity = TargetContaminantAmbiguity.RenameProtein;
            }

            //TheTask.SearchParameters.OutPepXML = ckbPepXML.IsChecked.Value;

            if (CheckBoxDecoy.IsChecked.Value)
            {
                if (RadioButtonReverseDecoy.IsChecked.Value)
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

            // Custom Mdac will be "" for all non-custom types, so no need to check for those.
            if (_massDifferenceAcceptorViewModel.CustomMdac != string.Empty)
            {
                try
                {
                    MassDiffAcceptor customMassDiffAcceptor =
                        SearchTask.GetMassDiffAcceptor(null, MassDiffAcceptorType.Custom, _massDifferenceAcceptorViewModel.CustomMdac);
                }
                catch (Exception ex)
                {
                    MessageBox.Show("Could not parse custom mass difference acceptor: " + ex.Message, "Error",
                        MessageBoxButton.OK, MessageBoxImage.Error);
                    return;
                }
            }

            TheTask.SearchParameters.MassDiffAcceptorType = _massDifferenceAcceptorViewModel.SelectedType.Type;
            TheTask.SearchParameters.CustomMdac = _massDifferenceAcceptorViewModel.CustomMdac;

            //determine if semi or nonspecific with a specific protease.
            if (searchModeType == CleavageSpecificity.Semi || protease.CleavageSpecificity == CleavageSpecificity.Semi)
            {
                TheTask.SearchParameters.LocalFdrCategories = new List<FdrCategory> { FdrCategory.FullySpecific, FdrCategory.SemiSpecific };
            }
            else if (searchModeType == CleavageSpecificity.None && protease.CleavageSpecificity != CleavageSpecificity.None)
            {
                TheTask.SearchParameters.LocalFdrCategories = new List<FdrCategory> { FdrCategory.FullySpecific, FdrCategory.SemiSpecific, FdrCategory.NonSpecific };
            }
            else
            {
                TheTask.SearchParameters.LocalFdrCategories = new List<FdrCategory> { FdrCategory.FullySpecific };
            }

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

            TheTask.SearchParameters.DoHistogramAnalysis = CheckBoxHistogramAnalysis.IsChecked.Value;
            TheTask.SearchParameters.HistogramBinTolInDaltons = double.Parse(HistogramBinWidthTextBox.Text, CultureInfo.InvariantCulture);

            TheTask.SearchParameters.WritePrunedDatabase = WritePrunedDBCheckBox.IsChecked.Value;

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
            if (NonSpecificSearchRadioButton.IsChecked.Value)
            {
                ProteaseComboBox.SelectedItem = ProteaseDictionary.Dictionary["non-specific"];
                AddCompIonCheckBox.IsChecked = true;
            }
            else
            {
                AddCompIonCheckBox.IsChecked = false;
                NTerminalIons.IsChecked = true;
                CTerminalIons.IsChecked = true;
            }
        }

        /// <summary>
        /// When a new multiplex label is selected from the drop-down menu (MultiplexComboBox),
        /// the fixed mod tree view is updated so that only those mods are selected. 
        /// </summary>
        private void MultiplexUpdate(object sender, RoutedEventArgs routedEventArgs)
        {
            ModTypeForTreeViewModel multiplexModType = FixedModTypeForTreeViewObservableCollection
                .FirstOrDefault(b => b.DisplayName.Equals("Multiplex Label"));
            string selectedModId = CheckBoxMultiplex.IsChecked.Value
                ? (string)MultiplexComboBox.SelectedItem
                : "NoMatchingMod";

            foreach (ModForTreeViewModel mod in multiplexModType.Children)
            {
                mod.Use = mod.DisplayName.Contains(selectedModId);
            }
            UpdateModSelectionGrid();
        }

        //this one is used by the GUI
        private void ProteaseSpecificUpdate(object sender, SelectionChangedEventArgs e)
        {
            string proteaseName = ((Protease)ProteaseComboBox.SelectedItem).Name;
            MissedCleavagesTextBox.IsEnabled = !proteaseName.Equals("top-down");

            if (AutomaticallyAskAndOrUpdateParametersBasedOnProtease)
            {
                switch (proteaseName)
                {
                    case "non-specific":
                        if (UpdateGUISettings.UseNonSpecificRecommendedSettings())
                        {
                            MaxPeptideLengthTextBox.Text = "25";
                            ClassicSearchRadioButton.IsChecked = false;
                            NonSpecificSearchRadioButton.IsChecked = true;
                        }
                        break;
                    case "top-down":
                        if (UpdateGUISettings.UseTopDownRecommendedSettings())
                        {
                            DeconHostViewModel.DoPrecursorDeconvolution = true;
                            DeconHostViewModel.UseProvidedPrecursors = false;
                            DeconHostViewModel.SetAllPrecursorMaxChargeState(60);
                            DeconHostViewModel.SetAllProductMaxChargeState(20);
                            TrimMsMs.IsChecked = false;
                            CheckBoxNoQuant.IsChecked = true; 
                            _massDifferenceAcceptorViewModel.SelectedType =
                                _massDifferenceAcceptorViewModel.MassDiffAcceptorTypes.First(p => p.Type == MassDiffAcceptorType.PlusOrMinusThreeMM);
                            maxModificationIsoformsTextBox.Text = "4096";
                            InternalIonsCheckBox.IsChecked = true;
                            MinInternalFragmentLengthTextBox.Text = "10";
                            //uncheck all variable mods
                            foreach (var mod in VariableModTypeForTreeViewObservableCollection)
                            {
                                mod.Use = false;
                            }
                        }
                        break;
                    case "Arg-C":
                        if (UpdateGUISettings.UseArgCRecommendedSettings())
                        {
                            ProteaseComboBox.SelectedItem = ProteaseDictionary.Dictionary["trypsin"];
                        }
                        break;
                    case "chymotrypsin (don't cleave before proline)":
                    case "chymotrypsin (cleave before proline)":
                        {
                            if (UpdateGUISettings.UseChymotrypsinRecommendedSettings())
                            {
                                MissedCleavagesTextBox.Text = "3";
                                SemiSpecificSearchRadioButton.IsChecked = true;
                            }
                        }
                        break;
                    case "elastase":
                        {
                            if (UpdateGUISettings.UseElastaseRecommendedSettings())
                            {
                                MissedCleavagesTextBox.Text = "16";
                                SemiSpecificSearchRadioButton.IsChecked = true;
                            }
                        }
                        break;
                    case "semi-trypsin":
                        {
                            if (UpdateGUISettings.UseSemiTrypsinRecommendedSettings())
                            {
                                ProteaseComboBox.SelectedItem = ProteaseDictionary.Dictionary["trypsin"];
                                SemiSpecificSearchRadioButton.IsChecked = true;
                            }
                        }
                        break;
                    default:
                        break;
                }
            }
        }

        //this is an overload for co-dependent parameters
        private void ProteaseSpecificUpdate(object sender, TextChangedEventArgs e)
        {
            if (((Protease)ProteaseComboBox.SelectedItem).Name.Contains("non-specific"))
            {
                try
                {
                    System.Windows.Controls.TextBox textBox = (TextBox)sender;
                    if (textBox.Name.Equals("MaxPeptideLengthTextBox")) //if maxPeptideLength was modified
                    {
                        if (!MissedCleavagesTextBox.Text.Equals((Convert.ToInt32(MaxPeptideLengthTextBox.Text) - 1).ToString())) //prevents infinite loops
                        {
                            MissedCleavagesTextBox.Text = (Convert.ToInt32(MaxPeptideLengthTextBox.Text) - 1).ToString();
                        }
                    }
                    else //if missedCleavagesTextBox was modified
                    {
                        if (!MaxPeptideLengthTextBox.Text.Equals((Convert.ToInt32(MissedCleavagesTextBox.Text) + 1).ToString())) //prevents infinite loops
                        {
                            MaxPeptideLengthTextBox.Text = (Convert.ToInt32(MissedCleavagesTextBox.Text) + 1).ToString();
                        }
                    }
                }
                catch
                {
                    //if not an entry, don't update the other box.
                }
            }

            if (!ClassicSearchRadioButton.IsChecked.Value && MaxPeptideLengthTextBox.Text.Length != 0)
            {
                //trim maxFragmentMass to reduce index size ( reduces index size on disc and makes reading/writing faster)
                int maxLength = Convert.ToInt32(MaxPeptideLengthTextBox.Text);
                if (maxLength > 0 && maxLength < 100) //default is 30000; 30000/300=100
                {
                    MaxFragmentMassTextBox.Text = (maxLength * 300).ToString(); //assume the average residue doesn't have a mass over 300 Da (largest is W @ 204, but mods exist)
                }
            }
        }

        private void SemiSpecificUpdate(object sender, RoutedEventArgs e)
        {
            AddCompIonCheckBox.IsChecked = SemiSpecificSearchRadioButton.IsChecked.Value;
            if (SemiSpecificSearchRadioButton.IsChecked.Value)
            {
                MissedCleavagesTextBox.Text = "2";
                MaxPeptideLengthTextBox.Text = "50";
            }
            else
            {
                NTerminalIons.IsChecked = true;
                CTerminalIons.IsChecked = true;
            }
        }

        private void SpectralRecoveryUpdate(object sender, RoutedEventArgs e)
        {
            if (CheckBoxLFQwSpectralRecovery.IsChecked.Value)
            {
                if (UpdateGUISettings.UseSpectralRecoveryMandatorySettings())
                {
                    CheckBoxMatchBetweenRuns.IsChecked = true;
                    WriteSpectralLibraryCheckBox.IsChecked = true;
                }
            }
        }

        private void KeyPressed(object sender, KeyEventArgs e)
        {
            if (e.Key == Key.Return && e.OriginalSource is not TextBox)
            {
                SaveButton_Click(sender, e);
            }
            else if (e.Key == Key.Escape)
            {
                CancelButton_Click(sender, e);
            }
        }

        private void TextChanged_Fixed(object sender, TextChangedEventArgs args)
        {
            SearchModifications.SetTimer();
            SearchModifications.FixedSearch = true;
        }

        private void TextChanged_Var(object sender, TextChangedEventArgs args)
        {
            SearchModifications.SetTimer();
            SearchModifications.VariableSearch = true;
        }

        private void TextChangeTimerHandler(object sender, EventArgs e)
        {
            if (SearchModifications.FixedSearch)
            {
                SearchModifications.FilterTree(SearchFixMod, FixedModsTreeView, FixedModTypeForTreeViewObservableCollection);
                SearchModifications.FixedSearch = false;
            }

            if (SearchModifications.VariableSearch)
            {
                SearchModifications.FilterTree(SearchVarMod, VariableModsTreeView, VariableModTypeForTreeViewObservableCollection);
                SearchModifications.VariableSearch = false;
            }
        }

        private void CustomFragmentationHandler(object sender, EventArgs e)
        {
            if (DissociationTypeComboBox.SelectedItem.ToString().Equals(DissociationType.Custom.ToString()))
            {
                CustomFragmentationWindow.Show();
            }
        }

        private void AddSilac_Click(object sender, RoutedEventArgs e)
        {
            var dialog = new SilacModificationWindow();
            dialog.MultiplexRadioButton.IsChecked = true; //set default
            if (dialog.ShowDialog() == true)
            {
                StaticSilacLabelsObservableCollection.Add(dialog.SilacLabel);
                //refresh the unlabeled option to update to the choices made by the added label (e.g. make the unlabeled the end condition if a start condition was added)
                CheckBoxQuantifyUnlabeledForSilac_Checked(sender, e);
                dataGridSilacLabels.Items.Refresh();
            }
        }

        private void DeleteSilac_Click(object sender, RoutedEventArgs e)
        {
            var selectedTask = (SilacInfoForDataGrid)dataGridSilacLabels.SelectedItem;
            if (selectedTask != null)
            {
                if (selectedTask.SilacLabel == null) //if unlabeled, uncheck the box
                {
                    CheckBoxQuantifyUnlabeledForSilac.IsChecked = false;
                }
                StaticSilacLabelsObservableCollection.Remove(selectedTask);
                dataGridSilacLabels.Items.Refresh();
            }
        }

        private void ClearSilac_Click(object sender, RoutedEventArgs e)
        {
            StaticSilacLabelsObservableCollection.Clear();
            CheckBoxQuantifyUnlabeledForSilac.IsChecked = false; //remove the unlabeled condition
            dataGridSilacLabels.Items.Refresh();
        }

        private void OnClosing(object sender, CancelEventArgs e)
        {
            SearchModifications.Timer.Tick -= new EventHandler(TextChangeTimerHandler);
            // remove event handler from timer
            // keeping it will trigger an exception because the closed window stops existing

            CustomFragmentationWindow.Close();
        }

        private static SilacLabel ConvertSilacDataGridInfoToSilacLabel(SilacInfoForDataGrid info)
        {
            if (info == null)
            {
                return null;
            }
            else
            {
                SilacLabel label = info.SilacLabel[0];
                //This is needed to prevent double adding of additional labels. 
                //A quick test is to create a silac condition with two labels, save, reopen the task, save, and reopen again. 
                //Without this line, the second label will be doubled. Example: (K+8)&(R+10)&(R+10)
                if (label.AdditionalLabels != null)
                {
                    label.AdditionalLabels.Clear();
                }

                for (int infoIndex = 1; infoIndex < info.SilacLabel.Count; infoIndex++)
                {
                    label.AddAdditionalSilacLabel(info.SilacLabel[infoIndex]);
                }
                return label;
            }
        }

        private CleavageSpecificity GetSearchModeType()
        {
            if (SemiSpecificSearchRadioButton.IsChecked.Value) //semi
            {
                return CleavageSpecificity.Semi;
            }
            else if (NonSpecificSearchRadioButton.IsChecked.Value) //non
            {
                return CleavageSpecificity.None;
            }
            else
            {
                return CleavageSpecificity.Full;
            }
        }

        private void SnesUpdates(CleavageSpecificity searchModeType)
        {
            if (searchModeType != CleavageSpecificity.Full)
            {
                if (((Protease)ProteaseComboBox.SelectedItem).Name.Contains("non-specific"))
                {
                    searchModeType = CleavageSpecificity.None; //prevents an accidental semi attempt of a non-specific protease

                    if (CTerminalIons.IsChecked.Value)
                    {
                        Protease singleC = ProteaseDictionary.Dictionary["singleC"];
                        ProteaseComboBox.SelectedItem = singleC;
                    }
                    else //we're not allowing no ion types. It must have N if it doesn't have C.
                    {
                        Protease singleN = ProteaseDictionary.Dictionary["singleN"];
                        ProteaseComboBox.SelectedItem = singleN;
                    }
                }
                if (!AddCompIonCheckBox.IsChecked.Value)
                {
                    MessageBox.Show("Warning: Complementary ions are strongly recommended when using this algorithm.");
                }
                //only use N or C termini, not both
                if (CTerminalIons.IsChecked.Value)
                {
                    NTerminalIons.IsChecked = false;
                }
                else
                {
                    NTerminalIons.IsChecked = true;
                }
            }
        }

        private FragmentationTerminus GetFragmentationTerminus()
        {
            if (NTerminalIons.IsChecked.Value && !CTerminalIons.IsChecked.Value)
            {
                return FragmentationTerminus.N;
            }
            else if (!NTerminalIons.IsChecked.Value && CTerminalIons.IsChecked.Value)
            {
                return FragmentationTerminus.C;
            }
            else if (!NTerminalIons.IsChecked.Value && !CTerminalIons.IsChecked.Value) //why would you want this
            {
                MessageBox.Show("Warning: No ion types were selected. MetaMorpheus will be unable to search MS/MS spectra.");
                return FragmentationTerminus.None;
            }
            else
            {
                return FragmentationTerminus.Both;
            }
        }

        //string out is for error messages
        private void SilacUpdates(out string error)
        {
            error = "";
            if (StaticSilacLabelsObservableCollection.Count != 0)
            {
                //remove the unlabeled
                for (int i = 0; i < StaticSilacLabelsObservableCollection.Count; i++)
                {
                    if (StaticSilacLabelsObservableCollection[i].SilacLabel == null)
                    {
                        StaticSilacLabelsObservableCollection.RemoveAt(i);
                        break;
                    }
                }

                //Validate it really quick to determine if they're all multiplex (normal), or if it's a turnover experiment (requires one start label and one end label, either of which may be unlabeled)
                //if they're all multiplex
                if (StaticSilacLabelsObservableCollection.All(x => x.LabelType == SilacModificationWindow.ExperimentType.Multiplex))
                {
                    List<SilacLabel> labelsToSave = new List<SilacLabel>();
                    foreach (SilacInfoForDataGrid info in StaticSilacLabelsObservableCollection)
                    {
                        labelsToSave.Add(ConvertSilacDataGridInfoToSilacLabel(info));
                    }
                    TheTask.SearchParameters.SilacLabels = labelsToSave;
                    //Clear existing start and end labels, if any
                    TheTask.SearchParameters.StartTurnoverLabel = null;
                    TheTask.SearchParameters.EndTurnoverLabel = null;
                }
                //if it's a turnover experiment, there's a start and/or end and no others
                else if (StaticSilacLabelsObservableCollection.Count <= 2)
                {
                    TheTask.SearchParameters.SilacLabels = null; //clear it if it's populated from a previous run. Test is run turnover, stop it, open the file, check if it's still turnover after closing and reopening the window
                    SilacInfoForDataGrid startLabel = StaticSilacLabelsObservableCollection.Where(x => x.LabelType == SilacModificationWindow.ExperimentType.Start).FirstOrDefault();
                    SilacInfoForDataGrid endLabel = StaticSilacLabelsObservableCollection.Where(x => x.LabelType == SilacModificationWindow.ExperimentType.End).FirstOrDefault();

                    //check that two labels weren't set as the same thing
                    if ((startLabel == null && StaticSilacLabelsObservableCollection.Count == 2) || (endLabel == null && StaticSilacLabelsObservableCollection.Count == 2))
                    {
                        string missingLabel = startLabel == null ? "Start" : "End";
                        error = "A SILAC label could not be found with the '" + missingLabel + "' LabelType." +
                            "\nThis error occurs when one of the conditions was mislabeled (i.e. one of the LabelTypes is 'Multiplex' or a duplicate). " +
                            "\nPlease check your 'LabelTypes' and try again.";
                    }
                    else //we're good!
                    {
                        TheTask.SearchParameters.StartTurnoverLabel = ConvertSilacDataGridInfoToSilacLabel(startLabel);
                        TheTask.SearchParameters.EndTurnoverLabel = ConvertSilacDataGridInfoToSilacLabel(endLabel);
                        //set the unlabeled bool for consistency
                        CheckBoxQuantifyUnlabeledForSilac.IsChecked = (startLabel == null || endLabel == null);
                    }
                }
                else //there are too many labels for a turnover experiment, but it's not a multiplex?
                {
                    error = StaticSilacLabelsObservableCollection.Count.ToString() + " SILAC labeling conditions were specified, but only two (a start and an end) were expected for a turnover experiment." +
                        "\nIf you are not analyzing a turnover experiment, please maintain the default setting 'Multiplex' when specifying your labels.";
                }
            }
        }

        private void CheckBoxQuantifyUnlabeledForSilac_Checked(object sender, RoutedEventArgs e)
        {
            if (CheckBoxSILAC.IsChecked.Value) //if we're doing SILAC (check needed to start a new task)
            {
                //if checked, add a proxy label for the unlabeled condition
                if (CheckBoxQuantifyUnlabeledForSilac.IsChecked.Value)
                {
                    //remove the unlabeled if it previously existed
                    for (int i = 0; i < StaticSilacLabelsObservableCollection.Count; i++)
                    {
                        if (StaticSilacLabelsObservableCollection[i].SilacLabel == null)
                        {
                            StaticSilacLabelsObservableCollection.RemoveAt(i);
                        }
                    }

                    bool startConditionFound = false;
                    bool endConditionFound = false;
                    foreach (SilacInfoForDataGrid item in StaticSilacLabelsObservableCollection)
                    {
                        if (item.LabelType == SilacModificationWindow.ExperimentType.Start)
                        {
                            startConditionFound = true;
                        }
                        else if (item.LabelType == SilacModificationWindow.ExperimentType.End)
                        {
                            endConditionFound = true;
                        }
                        //else do nothing
                    }
                    if (!startConditionFound && !endConditionFound)
                    {
                        StaticSilacLabelsObservableCollection.Add(new SilacInfoForDataGrid(SilacModificationWindow.ExperimentType.Multiplex));
                    }
                    else if (startConditionFound && !endConditionFound)
                    {
                        StaticSilacLabelsObservableCollection.Add(new SilacInfoForDataGrid(SilacModificationWindow.ExperimentType.End));
                    }
                    else if (!startConditionFound && endConditionFound)
                    {
                        StaticSilacLabelsObservableCollection.Add(new SilacInfoForDataGrid(SilacModificationWindow.ExperimentType.Start));
                    }
                    else
                    {
                        CheckBoxQuantifyUnlabeledForSilac.IsChecked = false;
                        MessageBox.Show("Unable to add unlabeled condition. Two turnover conditions have already been specified.");
                    }
                }
                else //remove the unlabeled condition
                {
                    for (int i = 0; i < StaticSilacLabelsObservableCollection.Count; i++)
                    {
                        if (StaticSilacLabelsObservableCollection[i].SilacLabel == null)
                        {
                            StaticSilacLabelsObservableCollection.RemoveAt(i);
                        }
                    }
                }
                dataGridSilacLabels.Items.Refresh();
            }
        }

        //displays the unlabeled sequence when SILAC is selected
        private void CheckBoxSILAC_Checked(object sender, RoutedEventArgs e)
        {
            CheckBoxQuantifyUnlabeledForSilac_Checked(sender, e);
        }

        private void CheckBoxMultiplex_Checked(object sender, RoutedEventArgs e)
        {

        }

        private void SaveAsDefault_Click(object sender, RoutedEventArgs e)
        {
            SaveButton_Click(sender, e);
            Toml.WriteFile(TheTask, Path.Combine(GlobalVariables.DataDir, "DefaultParameters", @"SearchTaskDefault.toml"), MetaMorpheusTask.tomlConfig);
        }

        /// <summary>
        /// Retained/Lost Methionine is best handled through truncation search when truncation search is selected.
        /// </summary>
        /// <param name="sender"></param>
        /// <param name="e"></param>
        private void AddTruncationsCheckBox_Checked(object sender, RoutedEventArgs e)
        {
            if (AddTruncationsCheckBox.IsChecked.Value)
            {
                InitiatorMethionineBehaviorComboBox.SelectedIndex = (int)InitiatorMethionineBehavior.Retain;
            }
        }

        /// <summary>
        /// Sets the value of the Internal Ions TextBox upon being checked
        /// </summary>
        /// <param name="sender"></param>
        /// <param name="e"></param>
        private void InternalIonsCheckBox_Checked(object sender, RoutedEventArgs e)
        {
            MinInternalFragmentLengthTextBox.Text = "4";
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