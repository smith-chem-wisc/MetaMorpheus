using EngineLayer;
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

namespace MetaMorpheusGUI
{
    /// <summary>
    /// Interaction logic for SearchTaskWindow.xaml
    /// </summary>
    public partial class SearchTaskWindow : Window
    {
        #region Private Fields

        private readonly DataContextForSearchTaskWindow dataContextForSearchTaskWindow;
        private readonly ObservableCollection<SearchModeForDataGrid> SearchModesForThisTask = new ObservableCollection<SearchModeForDataGrid>();
        private readonly ObservableCollection<ModTypeForTreeView> fixedModTypeForTreeViewObservableCollection = new ObservableCollection<ModTypeForTreeView>();
        private readonly ObservableCollection<ModTypeForTreeView> variableModTypeForTreeViewObservableCollection = new ObservableCollection<ModTypeForTreeView>();
        private readonly ObservableCollection<ModTypeForLoc> localizeModTypeForTreeViewObservableCollection = new ObservableCollection<ModTypeForLoc>();
        private readonly ObservableCollection<ModTypeForGrid> modSelectionGridItems = new ObservableCollection<ModTypeForGrid>();

        #endregion Private Fields

        #region Public Constructors

        public SearchTaskWindow()
        {
            InitializeComponent();
            TheTask = new SearchTask();

            PopulateChoices();

            UpdateFieldsFromTask(TheTask);

            this.saveButton.Content = "Add the Search Task";

            dataContextForSearchTaskWindow = new DataContextForSearchTaskWindow
            {
                ExpanderTitle = string.Join(", ", SearchModesForThisTask.Where(b => b.Use).Select(b => b.Name)),
                AnalysisExpanderTitle = "Some analysis properties...",
                SearchModeExpanderTitle = "Some search properties..."
            };
            this.DataContext = dataContextForSearchTaskWindow;
        }

        public SearchTaskWindow(SearchTask task)
        {
            InitializeComponent();
            TheTask = task;

            PopulateChoices();

            UpdateFieldsFromTask(TheTask);

            dataContextForSearchTaskWindow = new DataContextForSearchTaskWindow
            {
                ExpanderTitle = string.Join(", ", SearchModesForThisTask.Where(b => b.Use).Select(b => b.Name)),
                AnalysisExpanderTitle = "Some analysis properties...",
                SearchModeExpanderTitle = "Some search properties..."
            };
            this.DataContext = dataContextForSearchTaskWindow;
        }

        #endregion Public Constructors

        #region Internal Properties

        internal SearchTask TheTask { get; private set; }

        #endregion Internal Properties

        #region Private Methods

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
            foreach (Protease protease in GlobalVariables.ProteaseDictionary.Values)
            {
                proteaseComboBox.Items.Add(protease);
            }
            proteaseComboBox.SelectedIndex = 12;

            foreach (string initiatior_methionine_behavior in Enum.GetNames(typeof(InitiatorMethionineBehavior)))
            {
                initiatorMethionineBehaviorComboBox.Items.Add(initiatior_methionine_behavior);
            }

            productMassToleranceComboBox.Items.Add("Da");
            productMassToleranceComboBox.Items.Add("ppm");

            precursorMassToleranceComboBox.Items.Add("Da");
            precursorMassToleranceComboBox.Items.Add("ppm");

            foreach (var hm in GlobalVariables.AllModsKnown.GroupBy(b => b.modificationType))
            {
                var theModType = new ModTypeForGrid(hm.Key);
                modSelectionGridItems.Add(theModType);
            }
            ModSelectionGrid.ItemsSource = modSelectionGridItems;

            foreach (var hm in GlobalVariables.AllModsKnown.GroupBy(b => b.modificationType))
            {
                var theModType = new ModTypeForTreeView(hm.Key, false);
                fixedModTypeForTreeViewObservableCollection.Add(theModType);
                foreach (var uah in hm)
                {
                    theModType.Children.Add(new ModForTreeView(uah.ToString(), false, uah.id, false, theModType));
                }
            }
            fixedModsTreeView.DataContext = fixedModTypeForTreeViewObservableCollection;

            foreach (var hm in GlobalVariables.AllModsKnown.GroupBy(b => b.modificationType))
            {
                var theModType = new ModTypeForTreeView(hm.Key, false);
                variableModTypeForTreeViewObservableCollection.Add(theModType);
                foreach (var uah in hm)
                {
                    theModType.Children.Add(new ModForTreeView(uah.ToString(), false, uah.id, false, theModType));
                }
            }
            variableModsTreeView.DataContext = variableModTypeForTreeViewObservableCollection;

            foreach (var hm in GlobalVariables.AllModsKnown.GroupBy(b => b.modificationType))
            {
                localizeModTypeForTreeViewObservableCollection.Add(new ModTypeForLoc(hm.Key));
            }
        }

        private void UpdateFieldsFromTask(SearchTask task)
        {
            classicSearchRadioButton.IsChecked = task.SearchParameters.SearchType == SearchType.Classic;
            modernSearchRadioButton.IsChecked = task.SearchParameters.SearchType == SearchType.Modern;
            nonSpecificSearchRadioButton1.IsChecked = task.SearchParameters.SearchType == SearchType.NonSpecific && task.CommonParameters.DigestionParams.Protease.Name.Contains("non-specific");
            semiSpecificSearchRadioButton.IsChecked = task.SearchParameters.SearchType == SearchType.NonSpecific && !task.CommonParameters.DigestionParams.Protease.Name.Contains("non-specific");
            txtMaxFragmentSize.Text = task.SearchParameters.MaxFragmentSize.ToString(CultureInfo.InvariantCulture);
            checkBoxParsimony.IsChecked = task.SearchParameters.DoParsimony;
            checkBoxNoOneHitWonders.IsChecked = task.SearchParameters.NoOneHitWonders;
            checkBoxQuantification.IsChecked = task.SearchParameters.DoQuantification;
            quantPpmTolerance.Text = task.SearchParameters.QuantifyPpmTol.ToString(CultureInfo.InvariantCulture);
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
            txtMinPeptideLength.Text = task.CommonParameters.DigestionParams.MinPeptideLength.ToString(CultureInfo.InvariantCulture);
            txtMaxPeptideLength.Text = task.CommonParameters.DigestionParams.MaxPeptideLength == int.MaxValue ? "" : task.CommonParameters.DigestionParams.MaxPeptideLength.ToString(CultureInfo.InvariantCulture);
            proteaseComboBox.SelectedItem = task.CommonParameters.DigestionParams.Protease;
            maxModificationIsoformsTextBox.Text = task.CommonParameters.DigestionParams.MaxModificationIsoforms.ToString(CultureInfo.InvariantCulture);
            txtMaxModNum.Text = task.CommonParameters.DigestionParams.MaxModsForPeptide.ToString(CultureInfo.InvariantCulture);
            initiatorMethionineBehaviorComboBox.SelectedIndex = (int)task.CommonParameters.DigestionParams.InitiatorMethionineBehavior;
            productMassToleranceTextBox.Text = task.CommonParameters.ProductMassTolerance.Value.ToString(CultureInfo.InvariantCulture);
            productMassToleranceComboBox.SelectedIndex = task.CommonParameters.ProductMassTolerance is AbsoluteTolerance ? 0 : 1;
            precursorMassToleranceTextBox.Text = task.CommonParameters.PrecursorMassTolerance.Value.ToString(CultureInfo.InvariantCulture);
            precursorMassToleranceComboBox.SelectedIndex = task.CommonParameters.PrecursorMassTolerance is AbsoluteTolerance ? 0 : 1;
            addCompIonCheckBox.IsChecked = task.SearchParameters.AddCompIons;
            bCheckBox.IsChecked = task.CommonParameters.BIons;
            yCheckBox.IsChecked = task.CommonParameters.YIons;
            cCheckBox.IsChecked = task.CommonParameters.CIons;
            zdotCheckBox.IsChecked = task.CommonParameters.ZdotIons;
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

            OutputFileNameTextBox.Text = task.CommonParameters.TaskDescriptor;
            //ckbPepXML.IsChecked = task.SearchParameters.OutPepXML;
            ckbMzId.IsChecked = task.SearchParameters.OutMzId;
            foreach (var mod in task.CommonParameters.ListOfModsFixed)
            {
                var theModType = fixedModTypeForTreeViewObservableCollection.FirstOrDefault(b => b.DisplayName.Equals(mod.Item1));
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
                    fixedModTypeForTreeViewObservableCollection.Add(theModType);
                    theModType.Children.Add(new ModForTreeView("UNKNOWN MODIFICATION!", true, mod.Item2, true, theModType));
                }
            }
            foreach (var mod in task.CommonParameters.ListOfModsVariable)
            {
                var theModType = variableModTypeForTreeViewObservableCollection.FirstOrDefault(b => b.DisplayName.Equals(mod.Item1));
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
                    variableModTypeForTreeViewObservableCollection.Add(theModType);
                    theModType.Children.Add(new ModForTreeView("UNKNOWN MODIFICATION!", true, mod.Item2, true, theModType));
                }
            }
            
            foreach (var heh in localizeModTypeForTreeViewObservableCollection)
            {
                    heh.Use = false;
            }
           
            foreach (var ye in variableModTypeForTreeViewObservableCollection)
            {
                ye.VerifyCheckState();
            }
            foreach (var ye in fixedModTypeForTreeViewObservableCollection)
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
            #region Check Task Validity

            if (nonSpecificSearchRadioButton1.IsChecked.Value || semiSpecificSearchRadioButton.IsChecked.Value)
            {
                if ((bCheckBox.IsChecked.Value || cCheckBox.IsChecked.Value) && (yCheckBox.IsChecked.Value || zdotCheckBox.IsChecked.Value))
                {
                    //MessageBox.Show("Only ion types from a single terminus are allowed for this search algorithm. \ne.g. b- and/or c-ions OR y- and/or zdot-ions. \nC-terminal ions (y and/or zdot) will be chosen by default.");
                    bCheckBox.IsChecked = false;
                    cCheckBox.IsChecked = false;
                }
                if (((Protease)proteaseComboBox.SelectedItem).Name.Contains("non-specific"))
                {
                    proteaseComboBox.Items.MoveCurrentToFirst();
                    proteaseComboBox.SelectedItem = proteaseComboBox.Items.CurrentItem;
                    if ((bCheckBox.IsChecked.Value || cCheckBox.IsChecked.Value))
                    {
                        while (!((Protease)proteaseComboBox.SelectedItem).Name.Equals("singleN"))
                        {
                            proteaseComboBox.Items.MoveCurrentToNext();
                            proteaseComboBox.SelectedItem = proteaseComboBox.Items.CurrentItem;
                        }
                    }
                    else
                    {
                        while (!((Protease)proteaseComboBox.SelectedItem).Name.Equals("singleC"))
                        {
                            proteaseComboBox.Items.MoveCurrentToNext();
                            proteaseComboBox.SelectedItem = proteaseComboBox.Items.CurrentItem;
                        }
                    }
                }
                if (((Protease)proteaseComboBox.SelectedItem).Name.Contains("semi-trypsin"))
                {
                    proteaseComboBox.Items.MoveCurrentToFirst();
                    proteaseComboBox.SelectedItem = proteaseComboBox.Items.CurrentItem;
                    while (!((Protease)proteaseComboBox.SelectedItem).Name.Equals("trypsin"))
                    {
                        proteaseComboBox.Items.MoveCurrentToNext();
                        proteaseComboBox.SelectedItem = proteaseComboBox.Items.CurrentItem;
                    }
                }
                if (!addCompIonCheckBox.IsChecked.Value)
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
            if (!int.TryParse(numberOfDatabaseSearchesTextBox.Text, out int numberOfDatabaseSearches) || numberOfDatabaseSearches == 0)
            {
                MessageBox.Show("The number of database partitions was set to zero. At least one database is required for searching.");
                return;
            }
            if (string.IsNullOrEmpty(missedCleavagesTextBox.Text))
            {
                missedCleavagesTextBox.Text = int.MaxValue.ToString();
            }
            if (!int.TryParse(txtMinPeptideLength.Text, out int minPeptideLength) || minPeptideLength < 1)
            {
                MessageBox.Show("The minimum peptide length must be a positive integer");
                return;
            }
            if (string.IsNullOrEmpty(txtMaxPeptideLength.Text))
            {
                txtMaxPeptideLength.Text = int.MaxValue.ToString();
            }
            if (!int.TryParse(txtMaxPeptideLength.Text, out int maxPeptideLength) || maxPeptideLength < 1)
            {
                MessageBox.Show("The minimum peptide length must be a positive integer");
                return;
            }
            if (maxPeptideLength < minPeptideLength)
            {
                MessageBox.Show("The maximum peptide length must be greater than or equal to the minimum peptide length.");
                return;
            }
            if (!double.TryParse(precursorMassToleranceTextBox.Text, out double precursorMassTolerance) || precursorMassTolerance <= 0)
            {
                MessageBox.Show("The precursor mass tolerance contains unrecognized characters. \n You entered " + '"' + precursorMassToleranceTextBox.Text + '"' + "\n Please enter a positive number.");
                return;
            }
            if (!double.TryParse(productMassToleranceTextBox.Text, out double productMassTolerance) || productMassTolerance <= 0)
            {
                MessageBox.Show("The product mass tolerance contains unrecognized characters. \n You entered " + '"' + productMassToleranceTextBox.Text + '"' + "\n Please enter a positive number.");
                return;
            }
            if (!double.TryParse(minScoreAllowed.Text, out double minScore) || minScore < 1)
            {
                MessageBox.Show("The minimum score allowed contains unrecognized characters. \n You entered " + '"' + minScoreAllowed.Text + '"' + "\n Please enter a positive, non-zero number.");
                return;
            }
            if (!int.TryParse(maxModificationIsoformsTextBox.Text, out int maxModificationIsoforms) || maxModificationIsoforms < 1)
            {
                MessageBox.Show("The maximum number of modification isoforms contains unrecognized characters. \n You entered " + '"' + maxModificationIsoformsTextBox.Text + '"' + "\n Please enter a positive, non-zero number.");
                return;
            }

            #endregion Check Task Validity

            #region Save Parameters

            Protease protease = (Protease)proteaseComboBox.SelectedItem;
            bool semiProteaseDigestion = (semiSpecificSearchRadioButton.IsChecked.Value && ((Protease)proteaseComboBox.SelectedItem).CleavageSpecificity != CleavageSpecificity.SingleN && ((Protease)proteaseComboBox.SelectedItem).CleavageSpecificity != CleavageSpecificity.SingleC);
            TerminusType terminusTypeSemiProtease = (bCheckBox.IsChecked.Value || cCheckBox.IsChecked.Value ? TerminusType.N : TerminusType.C);
            int maxMissedCleavages = (int.Parse(missedCleavagesTextBox.Text, CultureInfo.InvariantCulture));
            int minPeptideLengthValue = (int.Parse(txtMinPeptideLength.Text, NumberStyles.Any, CultureInfo.InvariantCulture));
            int maxPeptideLengthValue = (int.Parse(txtMaxPeptideLength.Text, NumberStyles.Any, CultureInfo.InvariantCulture));
            int maxModificationIsoformsValue = (int.Parse(maxModificationIsoformsTextBox.Text, CultureInfo.InvariantCulture));
            int maxModsForPeptideValue = (int.Parse(txtMaxModNum.Text, CultureInfo.InvariantCulture));
            InitiatorMethionineBehavior initiatorMethionineBehavior = ((InitiatorMethionineBehavior)initiatorMethionineBehaviorComboBox.SelectedIndex);
            DigestionParams digestionParamsToSave = new DigestionParams(
                protease: protease.Name,
                SemiProteaseDigestion: semiProteaseDigestion,
                TerminusTypeSemiProtease: terminusTypeSemiProtease,
                MaxMissedCleavages: maxMissedCleavages,
                MinPeptideLength: minPeptideLengthValue,
                MaxPeptideLength: maxPeptideLengthValue,
                MaxModificationIsoforms: maxModificationIsoforms,
                InitiatorMethionineBehavior: initiatorMethionineBehavior,
                MaxModsForPeptides: maxModsForPeptideValue);

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
            TheTask.SearchParameters.MaxFragmentSize = Double.Parse(txtMaxFragmentSize.Text, CultureInfo.InvariantCulture);
            TheTask.SearchParameters.AddCompIons = addCompIonCheckBox.IsChecked.Value;
           
            var listOfModsVariable = new List<(string, string)>();
            foreach (var heh in variableModTypeForTreeViewObservableCollection)
            {
                listOfModsVariable.AddRange(heh.Children.Where(b => b.Use).Select(b => (b.Parent.DisplayName, b.DisplayName)));
            }
            
            var listOfModsFixed = new List<(string, string)>();
            foreach (var heh in fixedModTypeForTreeViewObservableCollection)
            {
                listOfModsFixed.AddRange(heh.Children.Where(b => b.Use).Select(b => (b.Parent.DisplayName, b.DisplayName)));
            }
           
            bool TrimMs1Peaks = trimMs1.IsChecked.Value;
            bool TrimMsMsPeaks = trimMsMs.IsChecked.Value;
            int TopNpeaks = int.Parse(TopNPeaksTextBox.Text);
            double MinRatio = double.Parse(MinRatioTextBox.Text);

            CommonParameters CommonParamsToSave = new CommonParameters(
                UseDeltaScore: deltaScoreCheckBox.IsChecked.Value,
                ReportAllAmbiguity: allAmbiguity.IsChecked.Value,
                DeconvolutionMaxAssumedChargeState: int.Parse(DeconvolutionMaxAssumedChargeStateTextBox.Text, CultureInfo.InvariantCulture), 
                TotalPartitions: int.Parse(numberOfDatabaseSearchesTextBox.Text, CultureInfo.InvariantCulture),
                DoPrecursorDeconvolution: deconvolutePrecursors.IsChecked.Value, 
                UseProvidedPrecursorInfo: useProvidedPrecursor.IsChecked.Value,
                ScoreCutoff: double.Parse(minScoreAllowed.Text, CultureInfo.InvariantCulture), 
                CalculateEValue: eValueCheckBox.IsChecked.Value, 
                ListOfModsFixed: listOfModsFixed, 
                ListOfModsVariable:  listOfModsVariable, 
                BIons: bCheckBox.IsChecked.Value, 
                YIons: yCheckBox.IsChecked.Value, 
                CIons: cCheckBox.IsChecked.Value, 
                ZdotIons: zdotCheckBox.IsChecked.Value,
                PrecursorMassTolerance: PrecursorMassTolerance,
                ProductMassTolerance:ProductMassTolerance,
                DigestionParams: digestionParamsToSave ,
                TrimMs1Peaks: TrimMs1Peaks, 
                TrimMsMsPeaks: TrimMsMsPeaks, 
                TopNpeaks: TopNpeaks, 
                MinRatio: MinRatio);

            if (OutputFileNameTextBox.Text != "")
            {
                CommonParamsToSave.TaskDescriptor = OutputFileNameTextBox.Text;
            }
            else
            {
                CommonParamsToSave.TaskDescriptor = "SearchTask";
            }

            if (classicSearchRadioButton.IsChecked.Value)
            {
                TheTask.SearchParameters.SearchType = SearchType.Classic;
            }
            else if (modernSearchRadioButton.IsChecked.Value)
            {
                TheTask.SearchParameters.SearchType = SearchType.Modern;
            }
            else //if (nonSpecificSearchRadioButton.IsChecked.Value)
            {
                TheTask.SearchParameters.SearchType = SearchType.NonSpecific;
            }

            TheTask.SearchParameters.DoParsimony = checkBoxParsimony.IsChecked.Value;
            TheTask.SearchParameters.NoOneHitWonders = checkBoxNoOneHitWonders.IsChecked.Value;
            TheTask.SearchParameters.DoQuantification = checkBoxQuantification.IsChecked.Value;
            TheTask.SearchParameters.Normalize = checkBoxNormalize.IsChecked.Value;
            TheTask.SearchParameters.MatchBetweenRuns = checkBoxMatchBetweenRuns.IsChecked.Value;
            TheTask.SearchParameters.ModPeptidesAreDifferent = modPepsAreUnique.IsChecked.Value;
            TheTask.SearchParameters.QuantifyPpmTol = double.Parse(quantPpmTolerance.Text, CultureInfo.InvariantCulture);
            TheTask.SearchParameters.SearchTarget = checkBoxTarget.IsChecked.Value;
            TheTask.SearchParameters.OutMzId = ckbMzId.IsChecked.Value;
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

            TheTask.SearchParameters.DoHistogramAnalysis = checkBoxHistogramAnalysis.IsChecked.Value;
            TheTask.SearchParameters.HistogramBinTolInDaltons = double.Parse(histogramBinWidthTextBox.Text, CultureInfo.InvariantCulture);

            TheTask.SearchParameters.WritePrunedDatabase = writePrunedDBCheckBox.IsChecked.Value;

            SetModSelectionForPrunedDB();

            TheTask.CommonParameters = CommonParamsToSave;

            #endregion Save Parameters

            DialogResult = true;
        }

        private void ApmdExpander_Collapsed(object sender, RoutedEventArgs e)
        {
            dataContextForSearchTaskWindow.ExpanderTitle = string.Join(", ", SearchModesForThisTask.Where(b => b.Use).Select(b => b.Name));
            dataContextForSearchTaskWindow.AnalysisExpanderTitle = "Some analysis properties...";
            dataContextForSearchTaskWindow.SearchModeExpanderTitle = "Some search properties...";
        }

        private void SetModSelectionForPrunedDB()
        {
            TheTask.SearchParameters.ModsToWriteSelection = new Dictionary<string, int>();
            //checks the grid values for which button is checked then sets paramaters accordingly
            foreach (var modTypeInGrid in modSelectionGridItems)
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
                var huhb = modSelectionGridItems.FirstOrDefault(b => b.ModName == modType.Key);
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
            if (nonSpecificSearchRadioButton1.IsChecked.Value)
            {
                proteaseComboBox.Items.MoveCurrentToFirst();
                proteaseComboBox.SelectedItem = proteaseComboBox.Items.CurrentItem;
                while (!((Protease)proteaseComboBox.SelectedItem).Name.Contains("non-specific"))
                {
                    proteaseComboBox.Items.MoveCurrentToNext();
                    proteaseComboBox.SelectedItem = proteaseComboBox.Items.CurrentItem;
                }
                proteaseComboBox.IsEnabled = false;
                addCompIonCheckBox.IsChecked = true;
            }
            else
            {
                proteaseComboBox.IsEnabled = true;
                addCompIonCheckBox.IsChecked = false;
            }
        }

        private void NonSpecificUpdate(object sender, TextChangedEventArgs e)
        {
            if(((Protease)proteaseComboBox.SelectedItem).Name.Contains("non-specific"))
            {
                try
                {
                    System.Windows.Controls.TextBox textBox = (TextBox)sender;
                    if (textBox.Name.Equals("txtMaxPeptideLength")) //if maxPeptideLength was modified
                    {
                        if (!missedCleavagesTextBox.Text.Equals((Convert.ToInt32(txtMaxPeptideLength.Text) - 1).ToString())) //prevents infinite loops
                        {
                            missedCleavagesTextBox.Text = (Convert.ToInt32(txtMaxPeptideLength.Text) - 1).ToString();
                        }
                    }
                    else //if missedCleavagesTextBox was modified
                    {
                        if (!txtMaxPeptideLength.Text.Equals((Convert.ToInt32(missedCleavagesTextBox.Text) + 1).ToString())) //prevents infinite loops
                        {
                            txtMaxPeptideLength.Text = (Convert.ToInt32(missedCleavagesTextBox.Text) + 1).ToString();
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
                txtMaxPeptideLength.Text = maxLength.ToString();
            }
        }

        private void SemiSpecificUpdate(object sender, RoutedEventArgs e)
        {
            addCompIonCheckBox.IsChecked = semiSpecificSearchRadioButton.IsChecked.Value;
        }
        
        #endregion Private Methods

    }

    public class DataContextForSearchTaskWindow : INotifyPropertyChanged
    {
        #region Private Fields

        private string expanderTitle;
        private string searchModeExpanderTitle;
        private string modExpanderTitle;
        private string analysisExpanderTitle;

        #endregion Private Fields

        #region Public Events

        public event PropertyChangedEventHandler PropertyChanged;

        #endregion Public Events

        #region Public Properties

        public string ExpanderTitle
        {
            get
            { return expanderTitle; }
            set
            {
                expanderTitle = value;
                RaisePropertyChanged("ExpanderTitle");
            }
        }

        public string AnalysisExpanderTitle
        {
            get { return analysisExpanderTitle; }
            set
            {
                analysisExpanderTitle = value;
                RaisePropertyChanged("AnalysisExpanderTitle");
            }
        }

        public string ModExpanderTitle
        {
            get { return modExpanderTitle; }
            set
            {
                modExpanderTitle = value;
                RaisePropertyChanged("ModExpanderTitle");
            }
        }

        public string SearchModeExpanderTitle
        {
            get { return searchModeExpanderTitle; }
            set
            {
                searchModeExpanderTitle = value;
                RaisePropertyChanged("SearchModeExpanderTitle");
            }
        }

        #endregion Public Properties

        #region Protected Methods

        protected void RaisePropertyChanged(string name)
        {
            PropertyChanged?.Invoke(this, new PropertyChangedEventArgs(name));
        }

        #endregion Protected Methods
    }
}