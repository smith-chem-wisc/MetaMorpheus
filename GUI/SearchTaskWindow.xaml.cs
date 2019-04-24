using EngineLayer;
using MassSpectrometry;
using MzLibUtil;
using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;
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
        private readonly ObservableCollection<SilacInfoForDataGrid> StaticSilacLabelsObservableCollection = new ObservableCollection<SilacInfoForDataGrid>();

        private CustomFragmentationWindow CustomFragmentationWindow;

        internal SearchTask TheTask { get; private set; }

        public SearchTaskWindow() : this(null)
        {
        }

        public SearchTaskWindow(SearchTask task)
        {
            InitializeComponent();
            TheTask = task ?? new SearchTask();
            PopulateChoices();
            UpdateFieldsFromTask(TheTask);

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
            Protease trypsin = ProteaseDictionary.Dictionary["trypsin"];
            proteaseComboBox.SelectedItem = trypsin;

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
            proteaseComboBox.SelectedItem = task.CommonParameters.DigestionParams.SpecificProtease; //needs to be first, so nonspecific can override if necessary
            classicSearchRadioButton.IsChecked = task.SearchParameters.SearchType == SearchType.Classic;
            modernSearchRadioButton.IsChecked = task.SearchParameters.SearchType == SearchType.Modern;
            //do these in if statements so as not to trigger the change
            if (task.SearchParameters.SearchType == SearchType.NonSpecific && task.CommonParameters.DigestionParams.SearchModeType == CleavageSpecificity.None)
            {
                nonSpecificSearchRadioButton.IsChecked = true; //when this is changed it overrides the protease
                if (task.CommonParameters.DigestionParams.SpecificProtease.Name.Equals("singleC") || task.CommonParameters.DigestionParams.SpecificProtease.Name.Equals("singleN"))
                {
                    Protease nonspecific = ProteaseDictionary.Dictionary["non-specific"];
                    proteaseComboBox.SelectedItem = nonspecific;
                }
                else
                {
                    proteaseComboBox.SelectedItem = task.CommonParameters.DigestionParams.SpecificProtease;
                }
            }
            if (task.SearchParameters.SearchType == SearchType.NonSpecific && task.CommonParameters.DigestionParams.SearchModeType != CleavageSpecificity.None)
            {
                semiSpecificSearchRadioButton.IsChecked = true;
            }
            MaxFragmentMassTextBox.Text = task.SearchParameters.MaxFragmentSize.ToString(CultureInfo.InvariantCulture);
            checkBoxParsimony.IsChecked = task.SearchParameters.DoParsimony;
            checkBoxNoOneHitWonders.IsChecked = task.SearchParameters.NoOneHitWonders;
            checkBoxNoQuant.IsChecked = !task.SearchParameters.DoQuantification;
            checkBoxLFQ.IsChecked = task.SearchParameters.DoQuantification;
            if (task.SearchParameters.SilacLabels != null && task.SearchParameters.SilacLabels.Count != 0)
            {
                checkBoxSILAC.IsChecked = true;
                List<Proteomics.SilacLabel> labels = task.SearchParameters.SilacLabels;
                foreach (Proteomics.SilacLabel label in labels)
                {
                    SilacInfoForDataGrid infoToAdd = new SilacInfoForDataGrid(label);
                    if (label.AdditionalLabels != null)
                    {
                        foreach (Proteomics.SilacLabel additionalLabel in label.AdditionalLabels)
                        {
                            infoToAdd.AddAdditionalLabel(new SilacInfoForDataGrid(additionalLabel));
                        }
                    }
                    StaticSilacLabelsObservableCollection.Add(infoToAdd);
                }
            }
            CheckBoxQuantifyUnlabeledForSilac.IsChecked = task.CommonParameters.DigestionParams.GeneratehUnlabeledProteinsForSilac;
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

            NumberOfPeaksToKeepPerWindowTextBox.Text = task.CommonParameters.NumberOfPeaksToKeepPerWindow == int.MaxValue || !task.CommonParameters.NumberOfPeaksToKeepPerWindow.HasValue ? "" : task.CommonParameters.NumberOfPeaksToKeepPerWindow.Value.ToString(CultureInfo.InvariantCulture);
            MinimumAllowedIntensityRatioToBasePeakTexBox.Text = task.CommonParameters.MinimumAllowedIntensityRatioToBasePeak == double.MaxValue || !task.CommonParameters.MinimumAllowedIntensityRatioToBasePeak.HasValue ? "" : task.CommonParameters.MinimumAllowedIntensityRatioToBasePeak.Value.ToString(CultureInfo.InvariantCulture);
            WindowWidthThomsonsTextBox.Text = task.CommonParameters.WindowWidthThomsons == double.MaxValue || !task.CommonParameters.WindowWidthThomsons.HasValue ? "" : task.CommonParameters.WindowWidthThomsons.Value.ToString(CultureInfo.InvariantCulture);
            NumberOfWindowsTextBox.Text = task.CommonParameters.NumberOfWindows == int.MaxValue || !task.CommonParameters.NumberOfWindows.HasValue ? "" : task.CommonParameters.NumberOfWindows.Value.ToString(CultureInfo.InvariantCulture);
            normalizePeaksInWindowCheckBox.IsChecked = task.CommonParameters.NormalizePeaksAccrossAllWindows;

            maxThreadsTextBox.Text = task.CommonParameters.MaxThreadsToUsePerFile.ToString(CultureInfo.InvariantCulture);
            MinVariantDepthTextBox.Text = task.CommonParameters.MinVariantDepth.ToString(CultureInfo.InvariantCulture);
            MaxHeterozygousVariantsTextBox.Text = task.CommonParameters.MaxHeterozygousVariants.ToString(CultureInfo.InvariantCulture);
            CustomFragmentationWindow = new CustomFragmentationWindow(task.CommonParameters.CustomIons);

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
            ckbMzId.IsChecked = task.SearchParameters.WriteMzId;
            writeDecoyCheckBox.IsChecked = task.SearchParameters.WriteDecoys;
            writeContaminantCheckBox.IsChecked = task.SearchParameters.WriteContaminants;

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
                    var theMod = theModType.Children.FirstOrDefault(b => b.ModName.Equals(mod.Item2));
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
            CustomFragmentationWindow.Close();
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
                    searchModeType = CleavageSpecificity.None; //prevents an accidental semi attempt of a non-specific protease

                    if (cTerminalIons.IsChecked.Value)
                    {
                        Protease singleC = ProteaseDictionary.Dictionary["singleC"];
                        proteaseComboBox.SelectedItem = singleC;
                    }
                    else //we're not allowing no ion types. It must have N if it doesn't have C.
                    {
                        Protease singleN = ProteaseDictionary.Dictionary["singleN"];
                        proteaseComboBox.SelectedItem = singleN;
                    }
                }
                if (!addCompIonCheckBox.IsChecked.Value)
                {
                    MessageBox.Show("Warning: Complementary ions are strongly recommended when using this algorithm.");
                }
                //only use N or C termini, not both
                if (cTerminalIons.IsChecked.Value)
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
                peakFindingToleranceTextBox.Text, histogramBinWidthTextBox.Text, DeconvolutionMaxAssumedChargeStateTextBox.Text, NumberOfPeaksToKeepPerWindowTextBox.Text,
                MinimumAllowedIntensityRatioToBasePeakTexBox.Text, WindowWidthThomsonsTextBox.Text, NumberOfWindowsTextBox.Text, numberOfDatabaseSearchesTextBox.Text, MaxModNumTextBox.Text, MaxFragmentMassTextBox.Text, QValueTextBox.Text))
            {
                return;
            }

            Protease protease = (Protease)proteaseComboBox.SelectedItem;

            DissociationType dissociationType = GlobalVariables.AllSupportedDissociationTypes[dissociationTypeComboBox.SelectedItem.ToString()];
            CustomFragmentationWindow.Close();

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
            int MinVariantDepth = int.Parse(MinVariantDepthTextBox.Text, NumberStyles.Any, CultureInfo.InvariantCulture);
            int MaxHeterozygousVariants = int.Parse(MaxHeterozygousVariantsTextBox.Text, NumberStyles.Any, CultureInfo.InvariantCulture);
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
                fragmentationTerminus: fragmentationTerminus,
                generateUnlabeledProteinsForSilac: CheckBoxQuantifyUnlabeledForSilac.IsChecked.Value);

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

            bool TrimMs1Peaks = trimMs1.IsChecked.Value;
            bool TrimMsMsPeaks = trimMsMs.IsChecked.Value;

            int? numPeaksToKeep = null;
            if (int.TryParse(NumberOfPeaksToKeepPerWindowTextBox.Text, out int numberOfPeaksToKeeep))
            {
                numPeaksToKeep = numberOfPeaksToKeeep;
            }

            double? minimumAllowedIntensityRatioToBasePeak = null;
            if (double.TryParse(MinimumAllowedIntensityRatioToBasePeakTexBox.Text, out double minimumAllowedIntensityRatio))
            {
                minimumAllowedIntensityRatioToBasePeak = minimumAllowedIntensityRatio;
            }

            double? windowWidthThompsons = null;
            if (double.TryParse(WindowWidthThomsonsTextBox.Text, out double windowWidth))
            {
                windowWidthThompsons = windowWidth;
            }

            int? numberOfWindows = null;
            if (int.TryParse(NumberOfWindowsTextBox.Text, out int numWindows))
            {
                numberOfWindows = numWindows;
            }

            bool normalizePeaksAccrossAllWindows = normalizePeaksInWindowCheckBox.IsChecked.Value;

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
                numberOfPeaksToKeepPerWindow: numPeaksToKeep,
                minimumAllowedIntensityRatioToBasePeak: minimumAllowedIntensityRatioToBasePeak,
                windowWidthThomsons: windowWidthThompsons,
                numberOfWindows: numberOfWindows,//maybe change this some day
                normalizePeaksAccrossAllWindows: normalizePeaksAccrossAllWindows,//maybe change this some day
                addCompIons: addCompIonCheckBox.IsChecked.Value,
                qValueOutputFilter: QValueCheckBox.IsChecked.Value ? double.Parse(QValueTextBox.Text, CultureInfo.InvariantCulture) : 1.0,
                assumeOrphanPeaksAreZ1Fragments: protease.Name != "top-down",
                minVariantDepth: MinVariantDepth,
                maxHeterozygousVariants: MaxHeterozygousVariants);

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
            TheTask.SearchParameters.DoQuantification = !checkBoxNoQuant.IsChecked.Value;

            //SilacLabel deconvolution
            {
                if (StaticSilacLabelsObservableCollection.Count == 0)
                {
                    TheTask.SearchParameters.SilacLabels = null;
                }
                else
                {
                    List<Proteomics.SilacLabel> labelsToSave = new List<Proteomics.SilacLabel>();
                    foreach (SilacInfoForDataGrid info in StaticSilacLabelsObservableCollection)
                    {
                        Proteomics.SilacLabel labelToAdd = info.SilacLabel[0];

                        //This is needed to prevent double adding of additional labels. 
                        //A quick test is to create a silac condition with two labels, save, reopen the task, save, and reopen again. 
                        //Without this line, the second label will be doubled (K+8)&(R+10)&(R+10)
                        if (labelToAdd.AdditionalLabels != null)
                        {
                            labelToAdd.AdditionalLabels.Clear();
                        }

                        for (int infoIndex = 1; infoIndex < info.SilacLabel.Count; infoIndex++)
                        {
                            labelToAdd.AddAdditionalSilacLabel(info.SilacLabel[infoIndex]);
                        }
                        labelsToSave.Add(labelToAdd);
                    }
                    TheTask.SearchParameters.SilacLabels = labelsToSave;
                }
            }

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
                try
                {
                    MassDiffAcceptor customMassDiffAcceptor = SearchTask.GetMassDiffAcceptor(null, MassDiffAcceptorType.Custom, customkMdacTextBox.Text);
                }
                catch (Exception ex)
                {
                    MessageBox.Show("Could not parse custom mass difference acceptor: " + ex.Message, "Error", MessageBoxButton.OK, MessageBoxImage.Error);
                    return;
                }

                TheTask.SearchParameters.MassDiffAcceptorType = MassDiffAcceptorType.Custom;
                TheTask.SearchParameters.CustomMdac = customkMdacTextBox.Text;
            }

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
                Protease nonspecific = ProteaseDictionary.Dictionary["non-specific"];
                proteaseComboBox.SelectedItem = nonspecific;
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
                SearchModifications.FilterTree(SearchFixMod, fixedModsTreeView, FixedModTypeForTreeViewObservableCollection);
                SearchModifications.FixedSearch = false;
            }

            if (SearchModifications.VariableSearch)
            {
                SearchModifications.FilterTree(SearchVarMod, variableModsTreeView, VariableModTypeForTreeViewObservableCollection);
                SearchModifications.VariableSearch = false;
            }
        }

        private void CustomFragmentationHandler(object sender, EventArgs e)
        {
            if (dissociationTypeComboBox.SelectedItem.ToString().Equals(DissociationType.Custom.ToString()))
            {
                CustomFragmentationWindow.Show();
            }
        }

        private void AddSilac_Click(object sender, RoutedEventArgs e)
        {
            var dialog = new SilacModificationWindow();
            if (dialog.ShowDialog() == true)
            {
                if (GetNumberOfSilacMods() + dialog.SilacLabel.SilacLabel.Count <= 26)
                {
                    StaticSilacLabelsObservableCollection.Add(dialog.SilacLabel);
                    dataGridSilacLabels.Items.Refresh();
                }
                else
                {
                    MessageBox.Show("More than 26 total SILAC labels have been specified, which is the maximum for this implementation.\t" +
                        "The most recent label was not added.");
                }
            }
        }

        private void DeleteSilac_Click(object sender, RoutedEventArgs e)
        {
            var selectedTask = (SilacInfoForDataGrid)dataGridSilacLabels.SelectedItem;
            if (selectedTask != null)
            {
                StaticSilacLabelsObservableCollection.Remove(selectedTask);
                dataGridSilacLabels.Items.Refresh();
            }
        }

        private void ClearSilac_Click(object sender, RoutedEventArgs e)
        {
            StaticSilacLabelsObservableCollection.Clear();
            dataGridSilacLabels.Items.Refresh();
        }

        private void OnClosing(object sender, CancelEventArgs e)
        {
            CustomFragmentationWindow.Close();
        }

        public int GetNumberOfSilacMods()
        {
            return StaticSilacLabelsObservableCollection.Sum(x => x.SilacLabel.Count);
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