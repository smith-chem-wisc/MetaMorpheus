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
        private readonly ObservableCollection<ModTypeForTreeView> localizeModTypeForTreeViewObservableCollection = new ObservableCollection<ModTypeForTreeView>();

        #endregion Private Fields

        #region Public Constructors

        public SearchTaskWindow()
        {
            InitializeComponent();
            PopulateChoices();

            TheTask = new SearchTask();
            UpdateFieldsFromTask(TheTask);

            this.saveButton.Content = "Add the Search Task";

            dataContextForSearchTaskWindow = new DataContextForSearchTaskWindow
            {
                ExpanderTitle = string.Join(", ", SearchModesForThisTask.Where(b => b.Use).Select(b => b.Name)),
                //ModExpanderTitle =
                //"fixed: "
                //+ string.Join(",", ModFileListInWindow.Where(b => b.Fixed).Select(b => b.FileName))
                //+ " variable: "
                //+ string.Join(",", ModFileListInWindow.Where(b => b.Variable).Select(b => b.FileName))
                //+ " localize: "
                //+ string.Join(",", ModFileListInWindow.Where(b => b.Localize).Select(b => b.FileName)),
                AnalysisExpanderTitle = "Some analysis properties...",
                SearchModeExpanderTitle = "Some search properties..."
            };
            this.DataContext = dataContextForSearchTaskWindow;
        }

        public SearchTaskWindow(SearchTask task)
        {
            InitializeComponent();
            PopulateChoices();

            TheTask = task;
            UpdateFieldsFromTask(TheTask);

            dataContextForSearchTaskWindow = new DataContextForSearchTaskWindow
            {
                ExpanderTitle = string.Join(", ", SearchModesForThisTask.Where(b => b.Use).Select(b => b.Name)),
                //ModExpanderTitle =
                //"fixed: "
                //+ string.Join(",", ModFileListInWindow.Where(b => b.Fixed).Select(b => b.FileName))
                //+ " variable: "
                //+ string.Join(",", ModFileListInWindow.Where(b => b.Variable).Select(b => b.FileName))
                //+ " localize: "
                //+ string.Join(",", ModFileListInWindow.Where(b => b.Localize).Select(b => b.FileName)),
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

        private static Boolean TextBoxIntAllowed(String Text2)
        {
            return Array.TrueForAll<Char>(Text2.ToCharArray(),
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
            foreach (Protease protease in GlobalEngineLevelSettings.ProteaseDictionary.Values)
                proteaseComboBox.Items.Add(protease);
            proteaseComboBox.SelectedIndex = 12;

            foreach (string initiatior_methionine_behavior in Enum.GetNames(typeof(InitiatorMethionineBehavior)))
                initiatorMethionineBehaviorComboBox.Items.Add(initiatior_methionine_behavior);

            productMassToleranceComboBox.Items.Add("Absolute");
            productMassToleranceComboBox.Items.Add("Ppm");

            foreach (var hm in GlobalEngineLevelSettings.AllModsKnown.GroupBy(b => b.modificationType))
            {
                var theModType = new ModTypeForTreeView(hm.Key, false);
                fixedModTypeForTreeViewObservableCollection.Add(theModType);
                foreach (var uah in hm)
                    theModType.Children.Add(new ModForTreeView(uah.ToString(), false, uah.id, false, theModType));
            }
            fixedModsTreeView.DataContext = fixedModTypeForTreeViewObservableCollection;
            foreach (var hm in GlobalEngineLevelSettings.AllModsKnown.GroupBy(b => b.modificationType))
            {
                var theModType = new ModTypeForTreeView(hm.Key, false);
                variableModTypeForTreeViewObservableCollection.Add(theModType);
                foreach (var uah in hm)
                    theModType.Children.Add(new ModForTreeView(uah.ToString(), false, uah.id, false, theModType));
            }
            variableModsTreeView.DataContext = variableModTypeForTreeViewObservableCollection;
            foreach (var hm in GlobalEngineLevelSettings.AllModsKnown.GroupBy(b => b.modificationType))
            {
                var theModType = new ModTypeForTreeView(hm.Key, false);
                localizeModTypeForTreeViewObservableCollection.Add(theModType);
                foreach (var uah in hm)
                    theModType.Children.Add(new ModForTreeView(uah.ToString(), false, uah.id, false, theModType));
            }
            localizeModsTreeView.DataContext = localizeModTypeForTreeViewObservableCollection;

            foreach (var uu in GlobalEngineLevelSettings.SearchModesKnown)
                SearchModesForThisTask.Add(new SearchModeForDataGrid(uu));
            searchModesDataGrid.DataContext = SearchModesForThisTask;
        }

        private void UpdateFieldsFromTask(SearchTask task)
        {
            classicSearchRadioButton.IsChecked = task.SearchParameters.SearchType == SearchType.Classic;
            modernSearchRadioButton.IsChecked = task.SearchParameters.SearchType == SearchType.Modern;
            nonSpecificSearchRadioButton.IsChecked = task.SearchParameters.SearchType == SearchType.NonSpecific;
            checkBoxParsimony.IsChecked = task.SearchParameters.DoParsimony;
            checkBoxNoOneHitWonders.IsChecked = task.SearchParameters.NoOneHitWonders;
            checkBoxQuantification.IsChecked = task.SearchParameters.DoQuantification;
            quantPpmTolerance.Text = task.SearchParameters.QuantifyPpmTol.ToString(CultureInfo.InvariantCulture);
            checkBoxMatchBetweenRuns.IsChecked = task.SearchParameters.MatchBetweenRuns;
            modPepsAreUnique.IsChecked = task.SearchParameters.ModPeptidesAreUnique;
            checkBoxHistogramAnalysis.IsChecked = task.SearchParameters.DoHistogramAnalysis;
            checkBoxTarget.IsChecked = task.SearchParameters.SearchTarget;
            checkBoxDecoy.IsChecked = task.SearchParameters.SearchDecoy;
            missedCleavagesTextBox.Text = task.CommonParameters.DigestionParams.MaxMissedCleavages.ToString(CultureInfo.InvariantCulture);
            txtMinPeptideLength.Text = task.CommonParameters.DigestionParams.MinPeptideLength.HasValue ? task.CommonParameters.DigestionParams.MinPeptideLength.Value.ToString(CultureInfo.InvariantCulture) : "";
            txtMaxPeptideLength.Text = task.CommonParameters.DigestionParams.MaxPeptideLength.HasValue ? task.CommonParameters.DigestionParams.MaxPeptideLength.Value.ToString(CultureInfo.InvariantCulture) : "";
            proteaseComboBox.SelectedItem = task.CommonParameters.DigestionParams.Protease;
            maxModificationIsoformsTextBox.Text = task.CommonParameters.DigestionParams.MaxModificationIsoforms.ToString(CultureInfo.InvariantCulture);
            txtMaxModNum.Text = task.CommonParameters.DigestionParams.MaxModsForPeptide.ToString(CultureInfo.InvariantCulture);
            initiatorMethionineBehaviorComboBox.SelectedIndex = (int)task.CommonParameters.DigestionParams.InitiatorMethionineBehavior;
            productMassToleranceTextBox.Text = task.CommonParameters.ProductMassTolerance.Value.ToString(CultureInfo.InvariantCulture);
            productMassToleranceComboBox.SelectedIndex = task.CommonParameters.ProductMassTolerance is AbsoluteTolerance ? 0 : 1;
            addCompIonCheckBox.IsChecked = task.SearchParameters.AddCompIons;
            bCheckBox.IsChecked = task.CommonParameters.BIons;
            yCheckBox.IsChecked = task.CommonParameters.YIons;
            cCheckBox.IsChecked = task.CommonParameters.CIons;
            zdotCheckBox.IsChecked = task.CommonParameters.ZdotIons;
            conserveMemoryCheckBox.IsChecked = task.CommonParameters.ConserveMemory;
            numberOfDatabaseSearchesTextBox.Text = task.CommonParameters.TotalPartitions.ToString(CultureInfo.InvariantCulture);
            deconvolutePrecursors.IsChecked = task.CommonParameters.DoPrecursorDeconvolution;
            useProvidedPrecursor.IsChecked = task.CommonParameters.UseProvidedPrecursorInfo;
            maxDegreesOfParallelism.Text = task.CommonParameters.MaxDegreeOfParallelism.ToString();
            disposeOfFilesWhenDone.IsChecked = task.SearchParameters.DisposeOfFileWhenDone;
            allAmbiguity.IsChecked = task.CommonParameters.ReportAllAmbiguity;
            excelCompatible.IsChecked = task.CommonParameters.ExcelCompatible;
            DeconvolutionIntensityRatioTextBox.Text = task.CommonParameters.DeconvolutionIntensityRatio.ToString();
            DeconvolutionMaxAssumedChargeStateTextBox.Text = task.CommonParameters.DeconvolutionMaxAssumedChargeState.ToString();
            DeconvolutionMassToleranceInPpmTextBox.Text = task.CommonParameters.DeconvolutionMassTolerance.Value.ToString();
            minScoreAllowed.Text = task.CommonParameters.ScoreCutoff.ToString(CultureInfo.InvariantCulture);

            trimMs1.IsChecked = task.CommonParameters.TrimMs1Peaks;
            trimMsMs.IsChecked = task.CommonParameters.TrimMsMsPeaks;
            TopNPeaksCheckBox.Text = task.CommonParameters.TopNpeaks.HasValue ? task.CommonParameters.TopNpeaks.Value.ToString(CultureInfo.InvariantCulture) : "";
            MinRatioCheckBox.Text = task.CommonParameters.MinRatio.HasValue ? task.CommonParameters.MinRatio.Value.ToString(CultureInfo.InvariantCulture) : "";

            foreach (var mod in task.CommonParameters.ListOfModsFixed)
            {
                var theModType = fixedModTypeForTreeViewObservableCollection.FirstOrDefault(b => b.DisplayName.Equals(mod.Item1));
                if (theModType != null)
                {
                    var theMod = theModType.Children.FirstOrDefault(b => b.DisplayName.Equals(mod.Item2));
                    if (theMod != null)
                        theMod.Use = true;
                    else
                        theModType.Children.Add(new ModForTreeView("UNKNOWN MODIFICATION!", true, mod.Item2, true, theModType));
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
                        theMod.Use = true;
                    else
                        theModType.Children.Add(new ModForTreeView("UNKNOWN MODIFICATION!", true, mod.Item2, true, theModType));
                }
                else
                {
                    theModType = new ModTypeForTreeView(mod.Item1, true);
                    variableModTypeForTreeViewObservableCollection.Add(theModType);
                    theModType.Children.Add(new ModForTreeView("UNKNOWN MODIFICATION!", true, mod.Item2, true, theModType));
                }
            }

            localizeAllCheckBox.IsChecked = task.CommonParameters.LocalizeAll;
            if (task.CommonParameters.LocalizeAll)
            {
                foreach (var heh in localizeModTypeForTreeViewObservableCollection)
                    heh.Use = true;
            }
            else
            {
                foreach (var mod in task.CommonParameters.ListOfModsLocalize)
                {
                    var theModType = localizeModTypeForTreeViewObservableCollection.FirstOrDefault(b => b.DisplayName.Equals(mod.Item1));
                    if (theModType != null)
                    {
                        var theMod = theModType.Children.FirstOrDefault(b => b.DisplayName.Equals(mod.Item2));
                        if (theMod != null)
                            theMod.Use = true;
                        else
                            theModType.Children.Add(new ModForTreeView("UNKNOWN MODIFICATION!", true, mod.Item2, true, theModType));
                    }
                    else
                    {
                        theModType = new ModTypeForTreeView(mod.Item1, true);
                        localizeModTypeForTreeViewObservableCollection.Add(theModType);
                        theModType.Children.Add(new ModForTreeView("UNKNOWN MODIFICATION!", true, mod.Item2, true, theModType));
                    }
                }
            }

            foreach (var ye in variableModTypeForTreeViewObservableCollection)
                ye.VerifyCheckState();
            foreach (var ye in fixedModTypeForTreeViewObservableCollection)
                ye.VerifyCheckState();
            foreach (var ye in localizeModTypeForTreeViewObservableCollection)
                ye.VerifyCheckState();

            SearchModesForThisTask.First(b => b.searchMode.FileNameAddition.Equals(task.SearchParameters.MassDiffAcceptor.FileNameAddition)).Use = true;

            writePrunedDatabaseCheckBox.IsChecked = task.SearchParameters.WritePrunedDatabase;
            keepAllUniprotModsCheckBox.IsChecked = task.SearchParameters.KeepAllUniprotMods;

            searchModesDataGrid.Items.Refresh();
        }

        private void CancelButton_Click(object sender, RoutedEventArgs e)
        {
            DialogResult = false;
        }

        private void SaveButton_Click(object sender, RoutedEventArgs e)
        {
            #region Check Task Validity

            if (nonSpecificSearchRadioButton.IsChecked.Value)
            {
                if (((Protease)proteaseComboBox.SelectedItem).Name.Equals("singleC") && (bCheckBox.IsChecked.Value || cCheckBox.IsChecked.Value))
                    MessageBox.Show("Warning: N-terminal ions were chosen for the C-terminal protease 'singleC'");
                if (((Protease)proteaseComboBox.SelectedItem).Name.Equals("singleN") && (yCheckBox.IsChecked.Value || zdotCheckBox.IsChecked.Value))
                    MessageBox.Show("Warning: C-terminal ions were chosen for the N-terminal protease 'singleN'");
                if (((Protease)proteaseComboBox.SelectedItem).Name.Contains("non-specific"))
                {
                    MessageBox.Show("The non-specific protease is designed for classic/modern searches and should not be assigned for the non-specific search. \n Please use 'singleN' or 'singleC'.");
                    return;
                }
                if (((Protease)proteaseComboBox.SelectedItem).Name.Contains("semi-trypsin"))
                {
                    MessageBox.Show("The semi-trypsin protease is designed for classic/modern searches and should not be assigned for the non-specific search. \n Please use 'trypsin'.");
                    return;
                }
                if (!addCompIonCheckBox.IsChecked.Value)
                    MessageBox.Show("Warning: Complementary ions are recommended for non-specific searches");
            }
            if (int.Parse(numberOfDatabaseSearchesTextBox.Text, CultureInfo.InvariantCulture) == 0)
            {
                MessageBox.Show("The number of database partitions was set to zero. At least one database is required for searching.");
                return;
            }
            if (missedCleavagesTextBox.Text.Length == 0)
            {
                MessageBox.Show("The number of missed cleavages was left empty. For no missed cleavages, please enter zero.");
                return;
            }
            if (!double.TryParse(DeconvolutionIntensityRatioTextBox.Text, out double dir) || dir <= 0)
            {
                MessageBox.Show("The deconvolution intensity ratio contains unrecognized characters. \n You entered " + '"' + DeconvolutionIntensityRatioTextBox.Text + '"' + "\n Please enter a positive number.");
                return;
            }
            if (!double.TryParse(DeconvolutionMassToleranceInPpmTextBox.Text, out double dmtip) || dmtip <= 0)
            {
                MessageBox.Show("The deconvolution mass tolerance (in ppm) contains unrecognized characters. \n You entered " + '"' + DeconvolutionMassToleranceInPpmTextBox.Text + '"' + "\n Please enter a positive number.");
                return;
            }
            if (!double.TryParse(productMassToleranceTextBox.Text, out double pmt) || pmt <= 0)
            {
                MessageBox.Show("The product mass tolerance contains unrecognized characters. \n You entered " + '"' + productMassToleranceTextBox.Text + '"' + "\n Please enter a positive number.");
                return;
            }
            if (!double.TryParse(minScoreAllowed.Text, out double msa) || msa < 0)
            {
                MessageBox.Show("The minimum score allowed contains unrecognized characters. \n You entered " + '"' + minScoreAllowed.Text + '"' + "\n Please enter a positive number.");
                return;
            }

            #endregion Check Task Validity

            #region Save Parameters

            TheTask.CommonParameters.TrimMs1Peaks = trimMs1.IsChecked.Value;
            TheTask.CommonParameters.TrimMsMsPeaks = trimMsMs.IsChecked.Value;
            TheTask.CommonParameters.TopNpeaks = int.TryParse(TopNPeaksCheckBox.Text, out int TopNPeak) ? (int?)TopNPeak : null;
            TheTask.CommonParameters.MinRatio = double.TryParse(MinRatioCheckBox.Text, out double MinRatio) ? (double?)MinRatio : null;

            if (classicSearchRadioButton.IsChecked.Value)
                TheTask.SearchParameters.SearchType = SearchType.Classic;
            else if (modernSearchRadioButton.IsChecked.Value)
                TheTask.SearchParameters.SearchType = SearchType.Modern;
            else //if (nonSpecificSearchRadioButton.IsChecked.Value)
                TheTask.SearchParameters.SearchType = SearchType.NonSpecific;

            //Code for determining SemiSpecific
            TheTask.CommonParameters.DigestionParams.SemiProteaseDigestion = nonSpecificSearchRadioButton.IsChecked.Value && ((Protease)proteaseComboBox.SelectedItem).CleavageSpecificity != CleavageSpecificity.SingleN && ((Protease)proteaseComboBox.SelectedItem).CleavageSpecificity != CleavageSpecificity.SingleC;
            TheTask.CommonParameters.DigestionParams.TerminusTypeSemiProtease = bCheckBox.IsChecked.Value || cCheckBox.IsChecked.Value ? TerminusType.N : TerminusType.C;

            TheTask.SearchParameters.DoParsimony = checkBoxParsimony.IsChecked.Value;
            TheTask.SearchParameters.NoOneHitWonders = checkBoxNoOneHitWonders.IsChecked.Value;
            TheTask.SearchParameters.DoQuantification = checkBoxQuantification.IsChecked.Value;
            TheTask.SearchParameters.MatchBetweenRuns = checkBoxMatchBetweenRuns.IsChecked.Value;
            TheTask.SearchParameters.ModPeptidesAreUnique = modPepsAreUnique.IsChecked.Value;
            TheTask.SearchParameters.QuantifyPpmTol = double.Parse(quantPpmTolerance.Text, CultureInfo.InvariantCulture);
            TheTask.SearchParameters.SearchTarget = checkBoxTarget.IsChecked.Value;
            TheTask.SearchParameters.SearchDecoy = checkBoxDecoy.IsChecked.Value;
            TheTask.CommonParameters.DigestionParams.MaxMissedCleavages = int.Parse(missedCleavagesTextBox.Text, CultureInfo.InvariantCulture);
            TheTask.CommonParameters.DigestionParams.MinPeptideLength = int.TryParse(txtMinPeptideLength.Text, NumberStyles.Any, CultureInfo.InvariantCulture, out int temp) ? (int?)temp : null;
            TheTask.CommonParameters.DigestionParams.MaxPeptideLength = int.TryParse(txtMaxPeptideLength.Text, NumberStyles.Any, CultureInfo.InvariantCulture, out temp) ? (int?)temp : null;
            TheTask.CommonParameters.DigestionParams.Protease = (Protease)proteaseComboBox.SelectedItem;
            TheTask.CommonParameters.DigestionParams.MaxModificationIsoforms = int.Parse(maxModificationIsoformsTextBox.Text, CultureInfo.InvariantCulture);
            TheTask.CommonParameters.DigestionParams.MaxModsForPeptide = int.Parse(txtMaxModNum.Text, CultureInfo.InvariantCulture);
            TheTask.CommonParameters.DigestionParams.InitiatorMethionineBehavior = (InitiatorMethionineBehavior)initiatorMethionineBehaviorComboBox.SelectedIndex;
            if (productMassToleranceComboBox.SelectedIndex == 0)
                TheTask.CommonParameters.ProductMassTolerance = new AbsoluteTolerance(double.Parse(productMassToleranceTextBox.Text, CultureInfo.InvariantCulture));
            else
                TheTask.CommonParameters.ProductMassTolerance = new PpmTolerance(double.Parse(productMassToleranceTextBox.Text, CultureInfo.InvariantCulture));

            TheTask.SearchParameters.AddCompIons = addCompIonCheckBox.IsChecked.Value;
            TheTask.CommonParameters.BIons = bCheckBox.IsChecked.Value;
            TheTask.CommonParameters.YIons = yCheckBox.IsChecked.Value;
            TheTask.CommonParameters.CIons = cCheckBox.IsChecked.Value;
            TheTask.CommonParameters.ZdotIons = zdotCheckBox.IsChecked.Value;
            TheTask.CommonParameters.ConserveMemory = conserveMemoryCheckBox.IsChecked.Value;
            TheTask.CommonParameters.TotalPartitions = int.Parse(numberOfDatabaseSearchesTextBox.Text, CultureInfo.InvariantCulture);

            TheTask.CommonParameters.DoPrecursorDeconvolution = deconvolutePrecursors.IsChecked.Value;
            TheTask.CommonParameters.UseProvidedPrecursorInfo = useProvidedPrecursor.IsChecked.Value;

            TheTask.CommonParameters.ScoreCutoff = double.Parse(minScoreAllowed.Text, CultureInfo.InvariantCulture);

            TheTask.CommonParameters.ReportAllAmbiguity = allAmbiguity.IsChecked.Value;
            TheTask.CommonParameters.ExcelCompatible = excelCompatible.IsChecked.Value;

            TheTask.CommonParameters.DeconvolutionIntensityRatio = double.Parse(DeconvolutionIntensityRatioTextBox.Text, CultureInfo.InvariantCulture);
            TheTask.CommonParameters.DeconvolutionMaxAssumedChargeState = int.Parse(DeconvolutionMaxAssumedChargeStateTextBox.Text, CultureInfo.InvariantCulture);
            TheTask.CommonParameters.DeconvolutionMassTolerance = new PpmTolerance(double.Parse(DeconvolutionMassToleranceInPpmTextBox.Text, CultureInfo.InvariantCulture));
            TheTask.SearchParameters.DisposeOfFileWhenDone = disposeOfFilesWhenDone.IsChecked.Value;
            TheTask.CommonParameters.ListOfModsVariable = new List<Tuple<string, string>>();
            foreach (var heh in variableModTypeForTreeViewObservableCollection)
                TheTask.CommonParameters.ListOfModsVariable.AddRange(heh.Children.Where(b => b.Use).Select(b => new Tuple<string, string>(b.Parent.DisplayName, b.DisplayName)));
            TheTask.CommonParameters.ListOfModsFixed = new List<Tuple<string, string>>();
            foreach (var heh in fixedModTypeForTreeViewObservableCollection)
                TheTask.CommonParameters.ListOfModsFixed.AddRange(heh.Children.Where(b => b.Use).Select(b => new Tuple<string, string>(b.Parent.DisplayName, b.DisplayName)));

            if (localizeAllCheckBox.IsChecked.Value)
            {
                TheTask.CommonParameters.ListOfModsLocalize = null;
                TheTask.CommonParameters.LocalizeAll = true;
            }
            else
            {
                TheTask.CommonParameters.LocalizeAll = false;
                TheTask.CommonParameters.ListOfModsLocalize = new List<Tuple<string, string>>();
                foreach (var heh in localizeModTypeForTreeViewObservableCollection)
                    TheTask.CommonParameters.ListOfModsLocalize.AddRange(heh.Children.Where(b => b.Use).Select(b => new Tuple<string, string>(b.Parent.DisplayName, b.DisplayName)));
            }

            TheTask.SearchParameters.MassDiffAcceptor = SearchModesForThisTask.Where(b => b.Use).Select(b => b.searchMode).First();
            TheTask.SearchParameters.DoHistogramAnalysis = checkBoxHistogramAnalysis.IsChecked.Value;

            TheTask.SearchParameters.WritePrunedDatabase = writePrunedDatabaseCheckBox.IsChecked.Value;
            TheTask.SearchParameters.KeepAllUniprotMods = keepAllUniprotModsCheckBox.IsChecked.Value;
            if (int.TryParse(maxDegreesOfParallelism.Text, out int jsakdf))
                TheTask.CommonParameters.MaxDegreeOfParallelism = jsakdf;

            #endregion Save Parameters

            DialogResult = true;
        }

        private void AddNewAllowedPrecursorMassDiffsButton_Click(object sender, RoutedEventArgs e)
        {
            try
            {
                var ye = MetaMorpheusTask.ParseSearchMode(newAllowedPrecursorMassDiffsTextBox.Text);
                GlobalEngineLevelSettings.SearchModesKnown.Add(ye);
                SearchModesForThisTask.Add(new SearchModeForDataGrid(ye));
                searchModesDataGrid.Items.Refresh();
            }
            catch
            {
                MessageBox.Show("Examples:" + Environment.NewLine + "name dot 5 ppm 0,1.003,2.006" + Environment.NewLine + "name interval [-4;-3],[-0.5;0.5],[101;102]", "Error parsing search mode text box", MessageBoxButton.OK, MessageBoxImage.Warning);
            }
        }

        private void ApmdExpander_Collapsed(object sender, RoutedEventArgs e)
        {
            dataContextForSearchTaskWindow.ExpanderTitle = string.Join(", ", SearchModesForThisTask.Where(b => b.Use).Select(b => b.Name));
            //dataContextForSearchTaskWindow.ModExpanderTitle =
            //    "fixed: "
            //    + string.Join(",", ModFileListInWindow.Where(b => b.Fixed).Select(b => b.FileName))
            //    + " variable: "
            //    + string.Join(",", ModFileListInWindow.Where(b => b.Variable).Select(b => b.FileName))
            //    + " localize: "
            //    + string.Join(",", ModFileListInWindow.Where(b => b.Localize).Select(b => b.FileName));
            dataContextForSearchTaskWindow.AnalysisExpanderTitle = "Some analysis properties...";
            dataContextForSearchTaskWindow.SearchModeExpanderTitle = "Some search properties...";
        }

        private void PreviewIfInt(object sender, TextCompositionEventArgs e)
        {
            e.Handled = !TextBoxIntAllowed(e.Text);
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
            get { return expanderTitle; }
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