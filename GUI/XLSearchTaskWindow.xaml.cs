using EngineLayer;
using EngineLayer.CrosslinkSearch;
using MzLibUtil;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
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
    public partial class XLSearchTaskWindow : Window
    {
        private readonly DataContextForSearchTaskWindow DataContextForSearchTaskWindow;
        private readonly ObservableCollection<SearchModeForDataGrid> SearchModesForThisTask = new ObservableCollection<SearchModeForDataGrid>();
        private readonly ObservableCollection<ModTypeForTreeView> FixedModTypeForTreeViewObservableCollection = new ObservableCollection<ModTypeForTreeView>();
        private readonly ObservableCollection<ModTypeForTreeView> VariableModTypeForTreeViewObservableCollection = new ObservableCollection<ModTypeForTreeView>();
        private readonly ObservableCollection<ModTypeForLoc> LocalizeModTypeForTreeViewObservableCollection = new ObservableCollection<ModTypeForLoc>();

        public XLSearchTaskWindow()
        {
            InitializeComponent();
            PopulateChoices();

            TheTask = new XLSearchTask();
            UpdateFieldsFromTask(TheTask);

            saveButton.Content = "Add the XLSearch Task";

            DataContextForSearchTaskWindow = new DataContextForSearchTaskWindow()
            {
                ExpanderTitle = string.Join(", ", SearchModesForThisTask.Where(b => b.Use).Select(b => b.Name)),
                AnalysisExpanderTitle = "Some analysis properties...",
                SearchModeExpanderTitle = "Some search properties..."
            };
            DataContext = DataContextForSearchTaskWindow;
        }

        public XLSearchTaskWindow(XLSearchTask task)
        {
            InitializeComponent();
            PopulateChoices();

            TheTask = task;
            UpdateFieldsFromTask(TheTask);

            DataContextForSearchTaskWindow = new DataContextForSearchTaskWindow()
            {
                ExpanderTitle = string.Join(", ", SearchModesForThisTask.Where(b => b.Use).Select(b => b.Name)),
                AnalysisExpanderTitle = "Some analysis properties...",
                SearchModeExpanderTitle = "Some search properties..."
            };
            DataContext = DataContextForSearchTaskWindow;
        }

        internal XLSearchTask TheTask { get; private set; }

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
            foreach (string crosslinkerName in Enum.GetNames(typeof(CrosslinkerType)))
                cbCrosslinker.Items.Add(crosslinkerName);

            //foreach (string fragmentationType in Enum.GetNames(typeof(FragmentaionType)))
            //    cbFragmentation.Items.Add(fragmentationType);

            cbbXLprecusorMsTl.Items.Add("Da");
            cbbXLprecusorMsTl.Items.Add("ppm");

            //cbbXLBetaprecusorMsTl.Items.Add("Absolute");
            //cbbXLBetaprecusorMsTl.Items.Add("ppm");

            foreach (Protease protease in ProteaseDictionary.Dictionary.Values)
                proteaseComboBox.Items.Add(protease);
            proteaseComboBox.SelectedIndex = 12;

            foreach (string initiatior_methionine_behavior in Enum.GetNames(typeof(InitiatorMethionineBehavior)))
                initiatorMethionineBehaviorComboBox.Items.Add(initiatior_methionine_behavior);

            productMassToleranceComboBox.Items.Add("Da");
            productMassToleranceComboBox.Items.Add("ppm");

            //foreach (string toleranceUnit in Enum.GetNames(typeof(ToleranceUnit)))
            //    productMassToleranceComboBox.Items.Add(toleranceUnit);

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
            localizeModsTreeView.DataContext = LocalizeModTypeForTreeViewObservableCollection;
        }

        private void UpdateFieldsFromTask(XLSearchTask task)
        {
            //Crosslink search para
            cbCrosslinker.SelectedIndex = (int)task.XlSearchParameters.CrosslinkerType;
            //cbFragmentation.SelectedIndex = (int)task.XlSearchParameters.FragmentationType;
            ckbXLTopNum.IsChecked = task.XlSearchParameters.CrosslinkSearchTop;
            txtXLTopNum.Text = task.XlSearchParameters.CrosslinkSearchTopNum.ToString(CultureInfo.InvariantCulture);
            ckbQuenchH2O.IsChecked = task.XlSearchParameters.XlQuench_H2O;
            ckbQuenchNH2.IsChecked = task.XlSearchParameters.XlQuench_NH2;
            ckbQuenchTris.IsChecked = task.XlSearchParameters.XlQuench_Tris;
            //ckbSearchWithXLAllBeta.IsChecked = task.XlSearchParameters.CrosslinkSearchWithAllBeta;
            txtUdXLKerName.Text = task.XlSearchParameters.UdXLkerName;
            ckbUdXLkerCleavable.IsChecked = task.XlSearchParameters.UdXLkerCleavable;
            txtUdXLkerTotalMs.Text = task.XlSearchParameters.UdXLkerTotalMass.HasValue ? task.XlSearchParameters.UdXLkerTotalMass.Value.ToString(CultureInfo.InvariantCulture) : "";
            txtUdXLkerShortMass.Text = task.XlSearchParameters.UdXLkerShortMass.HasValue ? task.XlSearchParameters.UdXLkerShortMass.Value.ToString(CultureInfo.InvariantCulture) : "";
            txtUdXLkerLongMass.Text = task.XlSearchParameters.UdXLkerLongMass.HasValue ? task.XlSearchParameters.UdXLkerLongMass.Value.ToString(CultureInfo.InvariantCulture) : "";
            txtUdXLkerLoopMass.Text = task.XlSearchParameters.UdXLkerLoopMass.HasValue ? task.XlSearchParameters.UdXLkerLoopMass.Value.ToString(CultureInfo.InvariantCulture) : "";
            txtUdXLkerAminoAcids.Text = task.XlSearchParameters.UdXLkerResidues;
            txtUdXLkerAminoAcids2.Text = task.XlSearchParameters.UdXLkerResidues2;
            txtUdXLkerDeadendH2O.Text = task.XlSearchParameters.UdXLkerDeadendMassH2O.HasValue ? task.XlSearchParameters.UdXLkerDeadendMassH2O.Value.ToString(CultureInfo.InvariantCulture) : "";
            txtUdXLkerDeadendNH2.Text = task.XlSearchParameters.UdXLkerDeadendMassNH2.HasValue ? task.XlSearchParameters.UdXLkerDeadendMassNH2.Value.ToString(CultureInfo.InvariantCulture) : "";
            txtUdXLkerDeadendTris.Text = task.XlSearchParameters.UdXLkerDeadendMassTris.HasValue ? task.XlSearchParameters.UdXLkerDeadendMassTris.Value.ToString(CultureInfo.InvariantCulture) : "";
            cbbXLprecusorMsTl.SelectedIndex = task.XlSearchParameters.XlPrecusorMsTl is AbsoluteTolerance ? 0 : 1;
            txtXLPrecusorMsTl.Text = task.XlSearchParameters.XlPrecusorMsTl.Value.ToString(CultureInfo.InvariantCulture);
            //cbbXLBetaprecusorMsTl.SelectedIndex = task.XlSearchParameters.XlPrecusorMsTl is AbsoluteTolerance ? 0 : 1;
            //txtXLBetaPrecusorMsTl.Text = task.XlSearchParameters.XlPrecusorMsTl.Value.ToString(CultureInfo.InvariantCulture);
            trimMs1.IsChecked = task.CommonParameters.TrimMs1Peaks;
            trimMsMs.IsChecked = task.CommonParameters.TrimMsMsPeaks;
            TopNPeaksCheckBox.Text = task.CommonParameters.TopNpeaks.ToString(CultureInfo.InvariantCulture);
            MinRatioCheckBox.Text = task.CommonParameters.MinRatio.ToString(CultureInfo.InvariantCulture);

            ckbCharge_2_3.IsChecked = task.XlSearchParameters.XlCharge_2_3;
            ckbCharge_2_3_PrimeFragments.IsChecked = task.XlSearchParameters.XlCharge_2_3_PrimeFragment;

            checkBoxDecoy.IsChecked = task.XlSearchParameters.DecoyType != DecoyType.None;
            deconvolutePrecursors.IsChecked = task.CommonParameters.DoPrecursorDeconvolution;
            useProvidedPrecursor.IsChecked = task.CommonParameters.UseProvidedPrecursorInfo;
            missedCleavagesTextBox.Text = task.CommonParameters.DigestionParams.MaxMissedCleavages.ToString(CultureInfo.InvariantCulture);
            txtMinPeptideLength.Text = task.CommonParameters.DigestionParams.MinPeptideLength.ToString(CultureInfo.InvariantCulture);
            txtMaxPeptideLength.Text = task.CommonParameters.DigestionParams.MaxPeptideLength.ToString(CultureInfo.InvariantCulture);
            proteaseComboBox.SelectedItem = task.CommonParameters.DigestionParams.Protease;
            maxModificationIsoformsTextBox.Text = task.CommonParameters.DigestionParams.MaxModificationIsoforms.ToString(CultureInfo.InvariantCulture);
            initiatorMethionineBehaviorComboBox.SelectedIndex = (int)task.CommonParameters.DigestionParams.InitiatorMethionineBehavior;
            productMassToleranceTextBox.Text = task.CommonParameters.ProductMassTolerance.Value.ToString(CultureInfo.InvariantCulture);
            productMassToleranceComboBox.SelectedIndex = task.CommonParameters.ProductMassTolerance is AbsoluteTolerance ? 0 : 1;
            bCheckBox.IsChecked = task.CommonParameters.BIons;
            yCheckBox.IsChecked = task.CommonParameters.YIons;
            cCheckBox.IsChecked = task.CommonParameters.CIons;
            zdotCheckBox.IsChecked = task.CommonParameters.ZdotIons;
            minScoreAllowed.Text = task.CommonParameters.ScoreCutoff.ToString(CultureInfo.InvariantCulture);
            txtNumberOfDatabaseSearches.Text = task.CommonParameters.TotalPartitions.ToString(CultureInfo.InvariantCulture);
            maxThreadsTextBox.Text = task.CommonParameters.MaxThreadsToUsePerFile.ToString(CultureInfo.InvariantCulture);

            ckbAllResults.IsChecked = task.XlSearchParameters.XlOutAll;
            ckbPercolator.IsChecked = task.XlSearchParameters.XlOutPercolator;
            ckbCrosslink.IsChecked = task.XlSearchParameters.XlOutCrosslink;
            ckbPepXML.IsChecked = task.XlSearchParameters.XlOutPepXML;

            OutputFileNameTextBox.Text = task.CommonParameters.TaskDescriptor;

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
        }

        private void CancelButton_Click(object sender, RoutedEventArgs e)
        {
            DialogResult = false;
        }

        private void SaveButton_Click(object sender, RoutedEventArgs e)
        {
            TheTask.XlSearchParameters.CrosslinkSearchTop = ckbXLTopNum.IsChecked.Value;
            TheTask.XlSearchParameters.CrosslinkSearchTopNum = int.Parse(txtXLTopNum.Text, CultureInfo.InvariantCulture);
            //TheTask.XlSearchParameters.CrosslinkSearchWithAllBeta = ckbSearchWithXLAllBeta.IsChecked.Value;
            TheTask.XlSearchParameters.CrosslinkerType = (CrosslinkerType)cbCrosslinker.SelectedIndex;
            //TheTask.XlSearchParameters.FragmentationType = (FragmentaionType)cbFragmentation.SelectedIndex;
            if (cbbXLprecusorMsTl.SelectedIndex == 0)
            {
                TheTask.XlSearchParameters.XlPrecusorMsTl = new AbsoluteTolerance(double.Parse(txtXLPrecusorMsTl.Text, CultureInfo.InvariantCulture));
            }
            else
            {
                TheTask.XlSearchParameters.XlPrecusorMsTl = new PpmTolerance(double.Parse(txtXLPrecusorMsTl.Text, CultureInfo.InvariantCulture));
            }
            TheTask.XlSearchParameters.XlCharge_2_3 = ckbCharge_2_3.IsChecked.Value;
            TheTask.XlSearchParameters.XlCharge_2_3_PrimeFragment = ckbCharge_2_3_PrimeFragments.IsChecked.Value;
            TheTask.XlSearchParameters.XlQuench_H2O = ckbQuenchH2O.IsChecked.Value;
            TheTask.XlSearchParameters.XlQuench_NH2 = ckbQuenchNH2.IsChecked.Value;
            TheTask.XlSearchParameters.XlQuench_Tris = ckbQuenchTris.IsChecked.Value;

            if (TheTask.XlSearchParameters.CrosslinkerType == CrosslinkerType.UserDefined)
            {
                TheTask.XlSearchParameters.UdXLkerName = txtUdXLKerName.Text;
                TheTask.XlSearchParameters.UdXLkerCleavable = ckbUdXLkerCleavable.IsChecked.Value;
                TheTask.XlSearchParameters.UdXLkerResidues = txtUdXLkerAminoAcids.Text;
                TheTask.XlSearchParameters.UdXLkerResidues2 = txtUdXLkerAminoAcids2.Text;
                TheTask.XlSearchParameters.UdXLkerLongMass = string.IsNullOrEmpty(txtUdXLkerLongMass.Text) ? (double?)null : double.Parse(txtUdXLkerLongMass.Text, CultureInfo.InvariantCulture);
                TheTask.XlSearchParameters.UdXLkerShortMass = string.IsNullOrEmpty(txtUdXLkerShortMass.Text) ? (double?)null : double.Parse(txtUdXLkerShortMass.Text, CultureInfo.InvariantCulture);
                TheTask.XlSearchParameters.UdXLkerTotalMass = string.IsNullOrEmpty(txtUdXLkerTotalMs.Text) ? (double?)null : double.Parse(txtUdXLkerTotalMs.Text, CultureInfo.InvariantCulture);
                TheTask.XlSearchParameters.UdXLkerLoopMass = string.IsNullOrEmpty(txtUdXLkerLoopMass.Text) ? (double?)null : double.Parse(txtUdXLkerLoopMass.Text, CultureInfo.InvariantCulture);
                TheTask.XlSearchParameters.UdXLkerDeadendMassH2O = string.IsNullOrEmpty(txtUdXLkerDeadendH2O.Text) ? (double?)null : double.Parse(txtUdXLkerDeadendH2O.Text, CultureInfo.InvariantCulture);
                TheTask.XlSearchParameters.UdXLkerDeadendMassNH2 = string.IsNullOrEmpty(txtUdXLkerDeadendNH2.Text) ? (double?)null : double.Parse(txtUdXLkerDeadendNH2.Text, CultureInfo.InvariantCulture);
                TheTask.XlSearchParameters.UdXLkerDeadendMassTris = string.IsNullOrEmpty(txtUdXLkerDeadendTris.Text) ? (double?)null : double.Parse(txtUdXLkerDeadendTris.Text, CultureInfo.InvariantCulture);
            }

            // either of these does not require the protease to be entered
            // this may change if we want to use the aTDC method for xl tasks
            TheTask.XlSearchParameters.DecoyType = checkBoxDecoy.IsChecked.Value ? DecoyType.Reverse : DecoyType.None;

            Protease protease = (Protease)proteaseComboBox.SelectedItem;
            int MaxMissedCleavages = (int.Parse(missedCleavagesTextBox.Text, CultureInfo.InvariantCulture));
            int MinPeptideLength = (int.Parse(txtMinPeptideLength.Text, NumberStyles.Any, CultureInfo.InvariantCulture));
            int MaxPeptideLength = (int.Parse(txtMaxPeptideLength.Text, NumberStyles.Any, CultureInfo.InvariantCulture));
            int MaxModificationIsoforms = (int.Parse(maxModificationIsoformsTextBox.Text, CultureInfo.InvariantCulture));
            InitiatorMethionineBehavior InitiatorMethionineBehavior = ((InitiatorMethionineBehavior)initiatorMethionineBehaviorComboBox.SelectedIndex);
            DigestionParams digestionParamsToSave = new DigestionParams(
                protease: protease.Name,
                maxMissedCleavages: MaxMissedCleavages,
                minPeptideLength: MinPeptideLength,
                maxPeptideLength: MaxPeptideLength,
                maxModificationIsoforms: MaxModificationIsoforms,
                initiatorMethionineBehavior: InitiatorMethionineBehavior);

            Tolerance ProductMassTolerance;
            if (productMassToleranceComboBox.SelectedIndex == 0)
            {
                ProductMassTolerance = new AbsoluteTolerance(double.Parse(productMassToleranceTextBox.Text, CultureInfo.InvariantCulture));
            }
            else
            {
                ProductMassTolerance = new PpmTolerance(double.Parse(productMassToleranceTextBox.Text, CultureInfo.InvariantCulture));
            }

            TheTask.XlSearchParameters.XlOutPercolator = ckbPercolator.IsChecked.Value;
            TheTask.XlSearchParameters.XlOutPepXML = ckbPepXML.IsChecked.Value;
            TheTask.XlSearchParameters.XlOutAll = ckbAllResults.IsChecked.Value;
            TheTask.XlSearchParameters.XlOutCrosslink = ckbCrosslink.IsChecked.Value;
            //TheTask.UseProvidedPrecursorInfo = useProvidedPrecursor.IsChecked.Value;

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

            CommonParameters CommonParamsToSave = new CommonParameters(
                productMassTolerance: ProductMassTolerance,
                doPrecursorDeconvolution: deconvolutePrecursors.IsChecked.Value,
                useProvidedPrecursorInfo: useProvidedPrecursor.IsChecked.Value,
                digestionParams: digestionParamsToSave,
                trimMs1Peaks: trimMs1.IsChecked.Value,
                trimMsMsPeaks: trimMsMs.IsChecked.Value,
                topNpeaks: int.Parse(TopNPeaksCheckBox.Text),
                minRatio: double.Parse(MinRatioCheckBox.Text),
                bIons: bCheckBox.IsChecked.Value,
                yIons: yCheckBox.IsChecked.Value,
                cIons: cCheckBox.IsChecked.Value,
                zDotIons: zdotCheckBox.IsChecked.Value,
                scoreCutoff: double.Parse(minScoreAllowed.Text, CultureInfo.InvariantCulture),
                totalPartitions: int.Parse(txtNumberOfDatabaseSearches.Text, CultureInfo.InvariantCulture),
                listOfModsVariable: listOfModsVariable,
                listOfModsFixed: listOfModsFixed);
            if (OutputFileNameTextBox.Text != "")
            {
                CommonParamsToSave.TaskDescriptor = OutputFileNameTextBox.Text;
            }
            else
            {
                CommonParamsToSave.TaskDescriptor = "XLSearchTask";
            }
            TheTask.CommonParameters = CommonParamsToSave;

            DialogResult = true;
        }

        private void ApmdExpander_Collapsed(object sender, RoutedEventArgs e)
        {
            DataContextForSearchTaskWindow.ExpanderTitle = string.Join(", ", SearchModesForThisTask.Where(b => b.Use).Select(b => b.Name));
            DataContextForSearchTaskWindow.AnalysisExpanderTitle = "Some analysis properties...";
            DataContextForSearchTaskWindow.SearchModeExpanderTitle = "Some search properties...";
        }

        private void ModExpander_Expanded(object sender, RoutedEventArgs e)
        {
            //Envent Handler not used #1
        }

        private void ModificationsDataGrid_Loaded(object sender, RoutedEventArgs e)
        {
            //Envent Handler not used #2
        }

        private void ModificationsDataGrid_DataContextChanged(object sender, DependencyPropertyChangedEventArgs e)
        {
            //Envent Handler not used #3
        }

        private void ModificationsDataGrid_AutoGeneratedColumns(object sender, EventArgs e)
        {
            //if (!TheTask.WritePrunedDatabase)
            //    modificationsDataGrid.Columns[3].Visibility = Visibility.Collapsed;
        }
    }
}