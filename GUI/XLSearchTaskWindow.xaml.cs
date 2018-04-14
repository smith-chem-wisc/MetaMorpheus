﻿using EngineLayer;
using EngineLayer.CrosslinkSearch;
using MzLibUtil;
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
        #region Private Fields

        private readonly DataContextForSearchTaskWindow dataContextForSearchTaskWindow;

        private readonly ObservableCollection<SearchModeForDataGrid> SearchModesForThisTask = new ObservableCollection<SearchModeForDataGrid>();

        private readonly ObservableCollection<ModTypeForTreeView> fixedModTypeForTreeViewObservableCollection = new ObservableCollection<ModTypeForTreeView>();
        private readonly ObservableCollection<ModTypeForTreeView> variableModTypeForTreeViewObservableCollection = new ObservableCollection<ModTypeForTreeView>();
        private readonly ObservableCollection<ModTypeForLoc> localizeModTypeForTreeViewObservableCollection = new ObservableCollection<ModTypeForLoc>();

        #endregion Private Fields

        #region Public Constructors

        public XLSearchTaskWindow()
        {
            InitializeComponent();
            PopulateChoices();

            TheTask = new XLSearchTask();
            UpdateFieldsFromTask(TheTask);

            this.saveButton.Content = "Add the XLSearch Task";

            dataContextForSearchTaskWindow = new DataContextForSearchTaskWindow()
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

        public XLSearchTaskWindow(XLSearchTask task)
        {
            InitializeComponent();
            PopulateChoices();

            TheTask = task;
            UpdateFieldsFromTask(TheTask);

            dataContextForSearchTaskWindow = new DataContextForSearchTaskWindow()
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

        internal XLSearchTask TheTask { get; private set; }

        #endregion Internal Properties

        #region Private Methods

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

            foreach (Protease protease in GlobalVariables.ProteaseDictionary.Values)
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
                fixedModTypeForTreeViewObservableCollection.Add(theModType);
                foreach (var uah in hm)
                    theModType.Children.Add(new ModForTreeView(uah.ToString(), false, uah.id, false, theModType));
            }
            fixedModsTreeView.DataContext = fixedModTypeForTreeViewObservableCollection;
            foreach (var hm in GlobalVariables.AllModsKnown.GroupBy(b => b.modificationType))
            {
                var theModType = new ModTypeForTreeView(hm.Key, false);
                variableModTypeForTreeViewObservableCollection.Add(theModType);
                foreach (var uah in hm)
                    theModType.Children.Add(new ModForTreeView(uah.ToString(), false, uah.id, false, theModType));
            }
            variableModsTreeView.DataContext = variableModTypeForTreeViewObservableCollection;

            foreach (var hm in GlobalVariables.AllModsKnown.GroupBy(b => b.modificationType))
                localizeModTypeForTreeViewObservableCollection.Add(new ModTypeForLoc(hm.Key));
            localizeModsTreeView.DataContext = localizeModTypeForTreeViewObservableCollection;
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
            txtUdXLkerAminoAcid.Text = task.XlSearchParameters.UdXLkerResidue.ToString();
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

            ckbAllResults.IsChecked = task.XlSearchParameters.XlOutAll;
            ckbPercolator.IsChecked = task.XlSearchParameters.XlOutPercolator;
            ckbCrosslink.IsChecked = task.XlSearchParameters.XlOutCrosslink;
            ckbPepXML.IsChecked = task.XlSearchParameters.XlOutPepXML;

            OutputFileNameTextBox.Text = task.CommonParameters.TaskDescriptor;

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

            localizeAllCheckBox.IsChecked = task.CommonParameters.LocalizeAll;
            if (task.CommonParameters.LocalizeAll)
            {
                foreach (var heh in localizeModTypeForTreeViewObservableCollection)
                {
                    heh.Use = false;
                }
            }
            else
            {
                foreach (var heh in localizeModTypeForTreeViewObservableCollection)
                {
                    heh.Use = task.CommonParameters.ListOfModTypesLocalize.Contains(heh.DisplayName);
                }
            }

            foreach (var ye in variableModTypeForTreeViewObservableCollection)
            {
                ye.VerifyCheckState();
            }
            foreach (var ye in fixedModTypeForTreeViewObservableCollection)
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
            CommonParameters CommonParamsToSave = (TheTask.CommonParameters as CommonParameters).Clone();
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
            //if (cbbXLBetaprecusorMsTl.SelectedIndex == 0)
            //    TheTask.XlSearchParameters.XlBetaPrecusorMsTl = new AbsoluteTolerance(double.Parse(txtXLBetaPrecusorMsTl.Text, CultureInfo.InvariantCulture));
            //else
            //    TheTask.XlSearchParameters.XlBetaPrecusorMsTl = new PpmTolerance(double.Parse(txtXLBetaPrecusorMsTl.Text, CultureInfo.InvariantCulture));

            if (TheTask.XlSearchParameters.CrosslinkerType == CrosslinkerType.UserDefined)
            {
                TheTask.XlSearchParameters.UdXLkerName = txtUdXLKerName.Text;
                TheTask.XlSearchParameters.UdXLkerCleavable = ckbUdXLkerCleavable.IsChecked.Value;
                TheTask.XlSearchParameters.UdXLkerLongMass = string.IsNullOrEmpty(txtUdXLkerLongMass.Text) ? (double?)null : double.Parse(txtUdXLkerLongMass.Text, CultureInfo.InvariantCulture);
                TheTask.XlSearchParameters.UdXLkerShortMass = string.IsNullOrEmpty(txtUdXLkerShortMass.Text) ? (double?)null : double.Parse(txtUdXLkerShortMass.Text, CultureInfo.InvariantCulture);
                TheTask.XlSearchParameters.UdXLkerTotalMass = string.IsNullOrEmpty(txtUdXLkerTotalMs.Text) ? (double?)null : double.Parse(txtUdXLkerTotalMs.Text, CultureInfo.InvariantCulture);
                TheTask.XlSearchParameters.UdXLkerLoopMass = string.IsNullOrEmpty(txtUdXLkerLoopMass.Text) ? (double?)null : double.Parse(txtUdXLkerLoopMass.Text, CultureInfo.InvariantCulture);
                TheTask.XlSearchParameters.UdXLkerDeadendMassH2O = string.IsNullOrEmpty(txtUdXLkerDeadendH2O.Text) ? (double?)null : double.Parse(txtUdXLkerDeadendH2O.Text, CultureInfo.InvariantCulture);
                TheTask.XlSearchParameters.UdXLkerDeadendMassNH2 = string.IsNullOrEmpty(txtUdXLkerDeadendNH2.Text) ? (double?)null : double.Parse(txtUdXLkerDeadendNH2.Text, CultureInfo.InvariantCulture);
                TheTask.XlSearchParameters.UdXLkerDeadendMassTris = string.IsNullOrEmpty(txtUdXLkerDeadendTris.Text) ? (double?)null : double.Parse(txtUdXLkerDeadendTris.Text, CultureInfo.InvariantCulture);
            }
            CommonParamsToSave.TrimMs1Peaks = trimMs1.IsChecked.Value;
            CommonParamsToSave.TrimMsMsPeaks = trimMsMs.IsChecked.Value;
            CommonParamsToSave.TopNpeaks = int.Parse(TopNPeaksCheckBox.Text);
            CommonParamsToSave.MinRatio = double.Parse(MinRatioCheckBox.Text);

            CommonParamsToSave.DoPrecursorDeconvolution = deconvolutePrecursors.IsChecked.Value;
            CommonParamsToSave.UseProvidedPrecursorInfo = useProvidedPrecursor.IsChecked.Value;
            TheTask.XlSearchParameters.DecoyType = checkBoxDecoy.IsChecked.Value ? DecoyType.Reverse : DecoyType.None;

            DigestionParams digestionParamsToSave = new DigestionParams();
            digestionParamsToSave.MaxMissedCleavages = int.Parse(missedCleavagesTextBox.Text, CultureInfo.InvariantCulture);
            digestionParamsToSave.MinPeptideLength = int.Parse(txtMinPeptideLength.Text, NumberStyles.Any, CultureInfo.InvariantCulture);
            digestionParamsToSave.MaxPeptideLength = int.Parse(txtMaxPeptideLength.Text, NumberStyles.Any, CultureInfo.InvariantCulture);
            digestionParamsToSave.Protease = (Protease)proteaseComboBox.SelectedItem;
            digestionParamsToSave.MaxModificationIsoforms = int.Parse(maxModificationIsoformsTextBox.Text, CultureInfo.InvariantCulture);
            digestionParamsToSave.InitiatorMethionineBehavior = (InitiatorMethionineBehavior)initiatorMethionineBehaviorComboBox.SelectedIndex;
            CommonParamsToSave.DigestionParams = digestionParamsToSave;

            if (productMassToleranceComboBox.SelectedIndex == 0)
            {
                CommonParamsToSave.ProductMassTolerance = new AbsoluteTolerance(double.Parse(productMassToleranceTextBox.Text, CultureInfo.InvariantCulture));
            }
            else
            {
                CommonParamsToSave.ProductMassTolerance = new PpmTolerance(double.Parse(productMassToleranceTextBox.Text, CultureInfo.InvariantCulture));
            }
            CommonParamsToSave.BIons = bCheckBox.IsChecked.Value;
            CommonParamsToSave.YIons = yCheckBox.IsChecked.Value;
            CommonParamsToSave.CIons = cCheckBox.IsChecked.Value;
            CommonParamsToSave.ZdotIons = zdotCheckBox.IsChecked.Value;
            CommonParamsToSave.ScoreCutoff = double.Parse(minScoreAllowed.Text, CultureInfo.InvariantCulture);

            TheTask.XlSearchParameters.XlOutPercolator = ckbPercolator.IsChecked.Value;
            TheTask.XlSearchParameters.XlOutPepXML = ckbPepXML.IsChecked.Value;
            TheTask.XlSearchParameters.XlOutAll = ckbAllResults.IsChecked.Value;
            TheTask.XlSearchParameters.XlOutCrosslink = ckbCrosslink.IsChecked.Value;
            //TheTask.UseProvidedPrecursorInfo = useProvidedPrecursor.IsChecked.Value;

            if (OutputFileNameTextBox.Text != "")
            {
                CommonParamsToSave.TaskDescriptor = OutputFileNameTextBox.Text;
            }
            else
            {
                CommonParamsToSave.TaskDescriptor = "XLSearchTask";
            }

            var listOfModsVariable = new List<(string, string)>();
            foreach (var heh in variableModTypeForTreeViewObservableCollection)
            {
                listOfModsVariable.AddRange(heh.Children.Where(b => b.Use).Select(b => (b.Parent.DisplayName, b.DisplayName)));
            }
            CommonParamsToSave.ListOfModsVariable = listOfModsVariable;

            var listOfModsFixed = new List<(string, string)>();
            foreach (var heh in fixedModTypeForTreeViewObservableCollection)
            {
                listOfModsFixed.AddRange(heh.Children.Where(b => b.Use).Select(b => (b.Parent.DisplayName, b.DisplayName)));
            }
            CommonParamsToSave.ListOfModsFixed = listOfModsFixed;

            if (localizeAllCheckBox.IsChecked.Value)
            {
                CommonParamsToSave.ListOfModTypesLocalize = null;
                CommonParamsToSave.LocalizeAll = true;
            }
            else
            {
                CommonParamsToSave.LocalizeAll = false;
                CommonParamsToSave.ListOfModTypesLocalize = localizeModTypeForTreeViewObservableCollection.Where(b => b.Use.HasValue && b.Use.Value).Select(b => b.DisplayName).ToList();
            }

            TheTask.CommonParameters = CommonParamsToSave;

            DialogResult = true;
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

        #endregion Private Methods
    }
}