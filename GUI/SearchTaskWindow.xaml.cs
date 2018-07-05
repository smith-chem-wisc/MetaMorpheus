using EngineLayer;
using MzLibUtil;
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

        public SearchTaskWindow()
        {
            InitializeComponent();
            TheTask = new SearchTask();

            PopulateChoices();
            UpdateFieldsFromTask(TheTask);
            SaveButton.Content = "Add the Search Task";

            DataContextForSearchTaskWindow = new DataContextForSearchTaskWindow
            {
                ExpanderTitle = string.Join(", ", SearchModesForThisTask.Where(b => b.Use).Select(b => b.Name)),
                AnalysisExpanderTitle = "Some analysis properties...",
                SearchModeExpanderTitle = "Some search properties..."
            };
            DataContext = DataContextForSearchTaskWindow;
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
            DataContext = DataContextForSearchTaskWindow;
        }

        internal SearchTask TheTask { get; private set; }

        private void PreviewIfInt(object sender, TextCompositionEventArgs e)
        {
            e.Handled = !TextBoxIntAllowed(e.Text);
        }

        private static bool TextBoxIntAllowed(String Text2)
        {
            return Array.TrueForAll(Text2.ToCharArray(),
                delegate (Char c) { return Char.IsDigit(c) || Char.IsControl(c); });
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
                ProteaseComboBox.Items.Add(protease);
            }
            ProteaseComboBox.SelectedIndex = 12;

            foreach (string initiatior_methionine_behavior in Enum.GetNames(typeof(InitiatorMethionineBehavior)))
            {
                InitiatorMethionineBehaviorComboBox.Items.Add(initiatior_methionine_behavior);
            }

            ProductMassToleranceComboBox.Items.Add("Da");
            ProductMassToleranceComboBox.Items.Add("ppm");

            PrecursorMassToleranceComboBox.Items.Add("Da");
            PrecursorMassToleranceComboBox.Items.Add("ppm");

            foreach (var hm in GlobalVariables.AllModsKnown.GroupBy(b => b.modificationType))
            {
                var theModType = new ModTypeForGrid(hm.Key);
                ModSelectionGridItems.Add(theModType);
            }
            ModSelectionGrid.ItemsSource = ModSelectionGridItems;

            foreach (var hm in GlobalVariables.AllModsKnown.GroupBy(b => b.modificationType))
            {
                var theModType = new ModTypeForTreeView(hm.Key, false);
                FixedModTypeForTreeViewObservableCollection.Add(theModType);
                foreach (var uah in hm)
                {
                    theModType.Children.Add(new ModForTreeView(uah.ToString(), false, uah.id, false, theModType));
                }
            }
            fixedModsTreeView.DataContext = FixedModTypeForTreeViewObservableCollection;

            foreach (var hm in GlobalVariables.AllModsKnown.GroupBy(b => b.modificationType))
            {
                var theModType = new ModTypeForTreeView(hm.Key, false);
                VariableModTypeForTreeViewObservableCollection.Add(theModType);
                foreach (var uah in hm)
                {
                    theModType.Children.Add(new ModForTreeView(uah.ToString(), false, uah.id, false, theModType));
                }
            }
            variableModsTreeView.DataContext = VariableModTypeForTreeViewObservableCollection;

            foreach (var hm in GlobalVariables.AllModsKnown.GroupBy(b => b.modificationType))
            {
                LocalizeModTypeForTreeViewObservableCollection.Add(new ModTypeForLoc(hm.Key));
            }
        }

        private void UpdateFieldsFromTask(SearchTask task)
        {
            ClassicSearchRadioButton.IsChecked = task.SearchParameters.SearchType == SearchType.Classic;
            ModernSearchRadioButton.IsChecked = task.SearchParameters.SearchType == SearchType.Modern;
            NonSpecificSearchRadioButton.IsChecked = task.SearchParameters.SearchType == SearchType.NonSpecific && task.CommonParameters.DigestionParams.Protease.Name.Contains("non-specific");
            SemiSpecificSearchRadioButton.IsChecked = task.SearchParameters.SearchType == SearchType.NonSpecific && !task.CommonParameters.DigestionParams.Protease.Name.Contains("non-specific");
            TxtMaxFragmentSize.Text = task.SearchParameters.MaxFragmentSize.ToString(CultureInfo.InvariantCulture);
            CheckBoxParsimony.IsChecked = task.SearchParameters.DoParsimony;
            CheckBoxNoOneHitWonders.IsChecked = task.SearchParameters.NoOneHitWonders;
            CheckBoxQuantification.IsChecked = task.SearchParameters.DoQuantification;
            QuantPpmTolerance.Text = task.SearchParameters.QuantifyPpmTol.ToString(CultureInfo.InvariantCulture);
            CheckBoxMatchBetweenRuns.IsChecked = task.SearchParameters.MatchBetweenRuns;
            CheckBoxNormalize.IsChecked = task.SearchParameters.Normalize;
            ModPepsAreUnique.IsChecked = task.SearchParameters.ModPeptidesAreDifferent;
            CheckBoxHistogramAnalysis.IsChecked = task.SearchParameters.DoHistogramAnalysis;
            HistogramBinWidthTextBox.Text = task.SearchParameters.HistogramBinTolInDaltons.ToString(CultureInfo.InvariantCulture);
            CheckBoxTarget.IsChecked = task.SearchParameters.SearchTarget;
            CheckBoxDecoy.IsChecked = task.SearchParameters.DecoyType != DecoyType.None;
            RadioButtonReverseDecoy.IsChecked = task.SearchParameters.DecoyType == DecoyType.Reverse;
            RadioButtonSlideDecoy.IsChecked = task.SearchParameters.DecoyType == DecoyType.Slide;
            MissedCleavagesTextBox.Text = task.CommonParameters.DigestionParams.MaxMissedCleavages == int.MaxValue ? "" : task.CommonParameters.DigestionParams.MaxMissedCleavages.ToString(CultureInfo.InvariantCulture);
            TxtMinPeptideLength.Text = task.CommonParameters.DigestionParams.MinPeptideLength.ToString(CultureInfo.InvariantCulture);
            TxtMaxPeptideLength.Text = task.CommonParameters.DigestionParams.MaxPeptideLength == int.MaxValue ? "" : task.CommonParameters.DigestionParams.MaxPeptideLength.ToString(CultureInfo.InvariantCulture);
            ProteaseComboBox.SelectedItem = task.CommonParameters.DigestionParams.Protease;
            MaxModificationIsoformsTextBox.Text = task.CommonParameters.DigestionParams.MaxModificationIsoforms.ToString(CultureInfo.InvariantCulture);
            TxtMaxModNum.Text = task.CommonParameters.DigestionParams.MaxModsForPeptide.ToString(CultureInfo.InvariantCulture);
            InitiatorMethionineBehaviorComboBox.SelectedIndex = (int)task.CommonParameters.DigestionParams.InitiatorMethionineBehavior;
            ProductMassToleranceTextBox.Text = task.CommonParameters.ProductMassTolerance.Value.ToString(CultureInfo.InvariantCulture);
            ProductMassToleranceComboBox.SelectedIndex = task.CommonParameters.ProductMassTolerance is AbsoluteTolerance ? 0 : 1;
            PrecursorMassToleranceTextBox.Text = task.CommonParameters.PrecursorMassTolerance.Value.ToString(CultureInfo.InvariantCulture);
            PrecursorMassToleranceComboBox.SelectedIndex = task.CommonParameters.PrecursorMassTolerance is AbsoluteTolerance ? 0 : 1;
            AddCompIonCheckBox.IsChecked = task.CommonParameters.AddCompIons;
            BCheckBox.IsChecked = task.CommonParameters.BIons;
            YCheckBox.IsChecked = task.CommonParameters.YIons;
            CCheckBox.IsChecked = task.CommonParameters.CIons;
            ZdotCheckBox.IsChecked = task.CommonParameters.ZdotIons;
            NumberOfDatabaseSearchesTextBox.Text = task.CommonParameters.TotalPartitions.ToString(CultureInfo.InvariantCulture);
            DeconvolutePrecursors.IsChecked = task.CommonParameters.DoPrecursorDeconvolution;
            UseProvidedPrecursor.IsChecked = task.CommonParameters.UseProvidedPrecursorInfo;
            AllAmbiguity.IsChecked = task.CommonParameters.ReportAllAmbiguity;
            DeconvolutionMaxAssumedChargeStateTextBox.Text = task.CommonParameters.DeconvolutionMaxAssumedChargeState.ToString();
            MinScoreAllowed.Text = task.CommonParameters.ScoreCutoff.ToString(CultureInfo.InvariantCulture);
            EValueCheckBox.IsChecked = task.CommonParameters.CalculateEValue;
            DeltaScoreCheckBox.IsChecked = task.CommonParameters.UseDeltaScore;
            TrimMs1.IsChecked = task.CommonParameters.TrimMs1Peaks;
            TrimMsMs.IsChecked = task.CommonParameters.TrimMsMsPeaks;
            TopNPeaksTextBox.Text = task.CommonParameters.TopNpeaks == int.MaxValue ? "" : task.CommonParameters.TopNpeaks.ToString(CultureInfo.InvariantCulture);
            MinRatioTextBox.Text = task.CommonParameters.MinRatio.ToString(CultureInfo.InvariantCulture);
            maxThreadsTextBox.Text = task.CommonParameters.MaxThreadsToUsePerFile.ToString(CultureInfo.InvariantCulture);

            OutputFileNameTextBox.Text = task.CommonParameters.TaskDescriptor;
            //ckbPepXML.IsChecked = task.SearchParameters.OutPepXML;
            CkbMzId.IsChecked = task.SearchParameters.OutMzId;
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
            // Check Task Validity
            if (NonSpecificSearchRadioButton.IsChecked.Value || SemiSpecificSearchRadioButton.IsChecked.Value)
            {
                if ((BCheckBox.IsChecked.Value || CCheckBox.IsChecked.Value) && (YCheckBox.IsChecked.Value || ZdotCheckBox.IsChecked.Value))
                {
                    //MessageBox.Show("Only ion types from a single terminus are allowed for this search algorithm. \ne.g. b- and/or c-ions OR y- and/or zdot-ions. \nC-terminal ions (y and/or zdot) will be chosen by default.");
                    BCheckBox.IsChecked = false;
                    CCheckBox.IsChecked = false;
                }
                if (((Protease)ProteaseComboBox.SelectedItem).Name.Contains("non-specific"))
                {
                    ProteaseComboBox.Items.MoveCurrentToFirst();
                    ProteaseComboBox.SelectedItem = ProteaseComboBox.Items.CurrentItem;
                    if ((BCheckBox.IsChecked.Value || CCheckBox.IsChecked.Value))
                    {
                        while (!((Protease)ProteaseComboBox.SelectedItem).Name.Equals("singleN"))
                        {
                            ProteaseComboBox.Items.MoveCurrentToNext();
                            ProteaseComboBox.SelectedItem = ProteaseComboBox.Items.CurrentItem;
                        }
                    }
                    else
                    {
                        while (!((Protease)ProteaseComboBox.SelectedItem).Name.Equals("singleC"))
                        {
                            ProteaseComboBox.Items.MoveCurrentToNext();
                            ProteaseComboBox.SelectedItem = ProteaseComboBox.Items.CurrentItem;
                        }
                    }
                }
                if (((Protease)ProteaseComboBox.SelectedItem).Name.Contains("semi-trypsin"))
                {
                    ProteaseComboBox.Items.MoveCurrentToFirst();
                    ProteaseComboBox.SelectedItem = ProteaseComboBox.Items.CurrentItem;
                    while (!((Protease)ProteaseComboBox.SelectedItem).Name.Equals("trypsin"))
                    {
                        ProteaseComboBox.Items.MoveCurrentToNext();
                        ProteaseComboBox.SelectedItem = ProteaseComboBox.Items.CurrentItem;
                    }
                }
                if (!AddCompIonCheckBox.IsChecked.Value)
                    MessageBox.Show("Warning: Complementary ions are strongly recommended when using this algorithm.");
            }
            if (!int.TryParse(DeconvolutionMaxAssumedChargeStateTextBox.Text, out int deconMaxAssumedCharge) || deconMaxAssumedCharge < 1)
            {
                MessageBox.Show("The maximum assumed charge state for deconvolution contains unrecognized characters or was zero. \n You entered " + '"' + DeconvolutionMaxAssumedChargeStateTextBox.Text + '"' + "\n Please enter a positive number.");
                return;
            }
            if (TopNPeaksTextBox.Text.Length == 0)
            {
                TopNPeaksTextBox.Text = int.MaxValue.ToString();
            }
            if (!int.TryParse(TopNPeaksTextBox.Text, out int numPeaks) || numPeaks < 1)
            {
                MessageBox.Show("The Top N Peaks to be retained must be greater than zero. \n You entered " + '"' + TopNPeaksTextBox.Text + '"' + "\n Please enter a positive number.");
                return;
            }
            if (!double.TryParse(MinRatioTextBox.Text, out double minRatio) || minRatio < 0 || minRatio > 1)
            {
                MessageBox.Show("The minimum ratio was not set to a number between zero and one. \n You entered " + '"' + MinRatioTextBox.Text + '"');
                return;
            }
            if (!int.TryParse(NumberOfDatabaseSearchesTextBox.Text, out int numberOfDatabaseSearches) || numberOfDatabaseSearches == 0)
            {
                MessageBox.Show("The number of database partitions was set to zero. At least one database is required for searching.");
                return;
            }
            if (string.IsNullOrEmpty(MissedCleavagesTextBox.Text))
            {
                MissedCleavagesTextBox.Text = int.MaxValue.ToString();
            }
            if (!int.TryParse(TxtMinPeptideLength.Text, out int minPeptideLength) || minPeptideLength < 1)
            {
                MessageBox.Show("The minimum peptide length must be a positive integer");
                return;
            }
            if (string.IsNullOrEmpty(TxtMaxPeptideLength.Text))
            {
                TxtMaxPeptideLength.Text = int.MaxValue.ToString();
            }
            if (!int.TryParse(TxtMaxPeptideLength.Text, out int maxPeptideLength) || maxPeptideLength < 1)
            {
                MessageBox.Show("The minimum peptide length must be a positive integer");
                return;
            }
            if (maxPeptideLength < minPeptideLength)
            {
                MessageBox.Show("The maximum peptide length must be greater than or equal to the minimum peptide length.");
                return;
            }
            if (!double.TryParse(PrecursorMassToleranceTextBox.Text, out double precursorMassTolerance) || precursorMassTolerance <= 0)
            {
                MessageBox.Show("The precursor mass tolerance contains unrecognized characters. \n You entered " + '"' + PrecursorMassToleranceTextBox.Text + '"' + "\n Please enter a positive number.");
                return;
            }
            if (!double.TryParse(ProductMassToleranceTextBox.Text, out double productMassTolerance) || productMassTolerance <= 0)
            {
                MessageBox.Show("The product mass tolerance contains unrecognized characters. \n You entered " + '"' + ProductMassToleranceTextBox.Text + '"' + "\n Please enter a positive number.");
                return;
            }
            if (!double.TryParse(MinScoreAllowed.Text, out double minScore) || minScore < 1)
            {
                MessageBox.Show("The minimum score allowed contains unrecognized characters. \n You entered " + '"' + MinScoreAllowed.Text + '"' + "\n Please enter a positive, non-zero number.");
                return;
            }
            if (!int.TryParse(MaxModificationIsoformsTextBox.Text, out int maxModificationIsoforms) || maxModificationIsoforms < 1)
            {
                MessageBox.Show("The maximum number of modification isoforms contains unrecognized characters. \n You entered " + '"' + MaxModificationIsoformsTextBox.Text + '"' + "\n Please enter a positive, non-zero number.");
                return;
            }
            if (!int.TryParse(maxThreadsTextBox.Text, out int maxThreads) || maxThreads > Environment.ProcessorCount || maxThreads <= 0)
            {
                MessageBox.Show("Your current device has " + Environment.ProcessorCount + " processors. \n Please select a positive value less than or equal to this number.");
                return;
            }

            // Save Parameters
            Protease protease = (Protease)ProteaseComboBox.SelectedItem;
            bool semiProteaseDigestion = (SemiSpecificSearchRadioButton.IsChecked.Value && ((Protease)ProteaseComboBox.SelectedItem).CleavageSpecificity != CleavageSpecificity.SingleN && ((Protease)ProteaseComboBox.SelectedItem).CleavageSpecificity != CleavageSpecificity.SingleC);
            TerminusType terminusTypeSemiProtease = (BCheckBox.IsChecked.Value || CCheckBox.IsChecked.Value ? TerminusType.N : TerminusType.C);
            int maxMissedCleavages = (int.Parse(MissedCleavagesTextBox.Text, CultureInfo.InvariantCulture));
            int minPeptideLengthValue = (int.Parse(TxtMinPeptideLength.Text, NumberStyles.Any, CultureInfo.InvariantCulture));
            int maxPeptideLengthValue = (int.Parse(TxtMaxPeptideLength.Text, NumberStyles.Any, CultureInfo.InvariantCulture));
            int maxModificationIsoformsValue = (int.Parse(MaxModificationIsoformsTextBox.Text, CultureInfo.InvariantCulture));
            int maxModsForPeptideValue = (int.Parse(TxtMaxModNum.Text, CultureInfo.InvariantCulture));
            InitiatorMethionineBehavior initiatorMethionineBehavior = ((InitiatorMethionineBehavior)InitiatorMethionineBehaviorComboBox.SelectedIndex);
            DigestionParams digestionParamsToSave = new DigestionParams(
                protease: protease.Name,
                semiProteaseDigestion: semiProteaseDigestion,
                terminusTypeSemiProtease: terminusTypeSemiProtease,
                maxMissedCleavages: maxMissedCleavages,
                minPeptideLength: minPeptideLengthValue,
                maxPeptideLength: maxPeptideLengthValue,
                maxModificationIsoforms: maxModificationIsoforms,
                initiatorMethionineBehavior: initiatorMethionineBehavior,
                maxModsForPeptides: maxModsForPeptideValue);

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
            TheTask.SearchParameters.MaxFragmentSize = Double.Parse(TxtMaxFragmentSize.Text, CultureInfo.InvariantCulture);

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

            bool TrimMs1Peaks = TrimMs1.IsChecked.Value;
            bool TrimMsMsPeaks = TrimMsMs.IsChecked.Value;
            int TopNpeaks = int.Parse(TopNPeaksTextBox.Text);
            double MinRatio = double.Parse(MinRatioTextBox.Text);

            CommonParameters CommonParamsToSave = new CommonParameters(
                useDeltaScore: DeltaScoreCheckBox.IsChecked.Value,
                reportAllAmbiguity: AllAmbiguity.IsChecked.Value,
                deconvolutionMaxAssumedChargeState: int.Parse(DeconvolutionMaxAssumedChargeStateTextBox.Text, CultureInfo.InvariantCulture),
                totalPartitions: int.Parse(NumberOfDatabaseSearchesTextBox.Text, CultureInfo.InvariantCulture),
                doPrecursorDeconvolution: DeconvolutePrecursors.IsChecked.Value,
                useProvidedPrecursorInfo: UseProvidedPrecursor.IsChecked.Value,
                scoreCutoff: double.Parse(MinScoreAllowed.Text, CultureInfo.InvariantCulture),
                calculateEValue: EValueCheckBox.IsChecked.Value,
                listOfModsFixed: listOfModsFixed,
                listOfModsVariable: listOfModsVariable,
                bIons: BCheckBox.IsChecked.Value,
                yIons: YCheckBox.IsChecked.Value,
                cIons: CCheckBox.IsChecked.Value,
                zDotIons: ZdotCheckBox.IsChecked.Value,
                precursorMassTolerance: PrecursorMassTolerance,
                productMassTolerance: ProductMassTolerance,
                digestionParams: digestionParamsToSave,
                trimMs1Peaks: TrimMs1Peaks,
                trimMsMsPeaks: TrimMsMsPeaks,
                topNpeaks: TopNpeaks,
                minRatio: MinRatio,
                addCompIons: AddCompIonCheckBox.IsChecked.Value);

            if (OutputFileNameTextBox.Text != "")
            {
                CommonParamsToSave.TaskDescriptor = OutputFileNameTextBox.Text;
            }
            else
            {
                CommonParamsToSave.TaskDescriptor = "SearchTask";
            }

            if (ClassicSearchRadioButton.IsChecked.Value)
            {
                TheTask.SearchParameters.SearchType = SearchType.Classic;
            }
            else if (ModernSearchRadioButton.IsChecked.Value)
            {
                TheTask.SearchParameters.SearchType = SearchType.Modern;
            }
            else //if (nonSpecificSearchRadioButton.IsChecked.Value)
            {
                TheTask.SearchParameters.SearchType = SearchType.NonSpecific;
            }

            TheTask.SearchParameters.DoParsimony = CheckBoxParsimony.IsChecked.Value;
            TheTask.SearchParameters.NoOneHitWonders = CheckBoxNoOneHitWonders.IsChecked.Value;
            TheTask.SearchParameters.DoQuantification = CheckBoxQuantification.IsChecked.Value;
            TheTask.SearchParameters.Normalize = CheckBoxNormalize.IsChecked.Value;
            TheTask.SearchParameters.MatchBetweenRuns = CheckBoxMatchBetweenRuns.IsChecked.Value;
            TheTask.SearchParameters.ModPeptidesAreDifferent = ModPepsAreUnique.IsChecked.Value;
            TheTask.SearchParameters.QuantifyPpmTol = double.Parse(QuantPpmTolerance.Text, CultureInfo.InvariantCulture);
            TheTask.SearchParameters.SearchTarget = CheckBoxTarget.IsChecked.Value;
            TheTask.SearchParameters.OutMzId = CkbMzId.IsChecked.Value;
            //TheTask.SearchParameters.OutPepXML = ckbPepXML.IsChecked.Value;

            if (CheckBoxDecoy.IsChecked.Value)
            {
                if (RadioButtonReverseDecoy.IsChecked.Value)
                {
                    TheTask.SearchParameters.DecoyType = DecoyType.Reverse;
                }
                else if (RadioButtonSlideDecoy.IsChecked.Value)
                {
                    TheTask.SearchParameters.DecoyType = DecoyType.Slide;
                }
                else if (RadioButtonShuffleDecoy.IsChecked.Value)
                {
                    TheTask.SearchParameters.DecoyType = DecoyType.Reverse;
                }
                else
                {
                    throw new ArgumentException("Decoys should be generated, but decoy selected either doesn't exist or hasn't been implemented.");
                }
            }
            else
            {
                TheTask.SearchParameters.DecoyType = DecoyType.None;
            }

            if (!maxThreadsTextBox.Text.Equals("") && (int.Parse(maxThreadsTextBox.Text) <= Environment.ProcessorCount && int.Parse(maxThreadsTextBox.Text) > 0))
            {
                CommonParamsToSave.MaxThreadsToUsePerFile = int.Parse(maxThreadsTextBox.Text, CultureInfo.InvariantCulture);
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

            TheTask.SearchParameters.DoHistogramAnalysis = CheckBoxHistogramAnalysis.IsChecked.Value;
            TheTask.SearchParameters.HistogramBinTolInDaltons = double.Parse(HistogramBinWidthTextBox.Text, CultureInfo.InvariantCulture);

            TheTask.SearchParameters.WritePrunedDatabase = writePrunedDBCheckBox.IsChecked.Value;

            SetModSelectionForPrunedDB();

            TheTask.CommonParameters = CommonParamsToSave;

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
                ProteaseComboBox.Items.MoveCurrentToFirst();
                ProteaseComboBox.SelectedItem = ProteaseComboBox.Items.CurrentItem;
                while (!((Protease)ProteaseComboBox.SelectedItem).Name.Contains("non-specific"))
                {
                    ProteaseComboBox.Items.MoveCurrentToNext();
                    ProteaseComboBox.SelectedItem = ProteaseComboBox.Items.CurrentItem;
                }
                ProteaseComboBox.IsEnabled = false;
                AddCompIonCheckBox.IsChecked = true;
            }
            else
            {
                ProteaseComboBox.IsEnabled = true;
                AddCompIonCheckBox.IsChecked = false;
            }
        }

        private void NonSpecificUpdate(object sender, TextChangedEventArgs e)
        {
            if (((Protease)ProteaseComboBox.SelectedItem).Name.Contains("non-specific"))
            {
                try
                {
                    TextBox textBox = (TextBox)sender;
                    if (textBox.Name.Equals("txtMaxPeptideLength")) //if maxPeptideLength was modified
                    {
                        if (!MissedCleavagesTextBox.Text.Equals((Convert.ToInt32(TxtMaxPeptideLength.Text) - 1).ToString())) //prevents infinite loops
                        {
                            MissedCleavagesTextBox.Text = (Convert.ToInt32(TxtMaxPeptideLength.Text) - 1).ToString();
                        }
                    }
                    else //if missedCleavagesTextBox was modified
                    {
                        if (!TxtMaxPeptideLength.Text.Equals((Convert.ToInt32(MissedCleavagesTextBox.Text) + 1).ToString())) //prevents infinite loops
                        {
                            TxtMaxPeptideLength.Text = (Convert.ToInt32(MissedCleavagesTextBox.Text) + 1).ToString();
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
            if (((Protease)ProteaseComboBox.SelectedItem).Name.Contains("non-specific"))
            {
                TxtMaxPeptideLength.Text = maxLength.ToString();
            }
        }

        private void SemiSpecificUpdate(object sender, RoutedEventArgs e)
        {
            AddCompIonCheckBox.IsChecked = SemiSpecificSearchRadioButton.IsChecked.Value;
        }

        private void RadioButtonReverseDecoy_Checked(object sender, RoutedEventArgs e)
        {
            TheTask.SearchParameters.NumDecoyDatabases = 1;
            NumDecoyDatabasesTextBox.Text = 1.ToString(CultureInfo.InvariantCulture);
        }

        private void RadioButtonSlideDecoy_Checked(object sender, RoutedEventArgs e)
        {
            TheTask.SearchParameters.NumDecoyDatabases = 1;
            NumDecoyDatabasesTextBox.Text = 1.ToString(CultureInfo.InvariantCulture);
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