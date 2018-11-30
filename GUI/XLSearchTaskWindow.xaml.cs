using EngineLayer;
using EngineLayer.CrosslinkSearch;
using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.Globalization;
using System.Linq;
using System.Windows;
using System.Windows.Input;
using TaskLayer;
using UsefulProteomicsDatabases;
using Proteomics.ProteolyticDigestion;
using MassSpectrometry;
using System.Windows.Controls;

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
        private CustomFragmentationWindow CustomFragmentationWindow;

        public XLSearchTaskWindow() : this(null)
        {
        }

        public XLSearchTaskWindow(XLSearchTask task)
        {
            InitializeComponent();
            PopulateChoices();
            TheTask = task ?? new XLSearchTask();
            UpdateFieldsFromTask(TheTask);

            if (task == null)
            {
                this.saveButton.Content = "Add the XLSearch Task";
            }
            DataContextForSearchTaskWindow = new DataContextForSearchTaskWindow()
            {
                ExpanderTitle = string.Join(", ", SearchModesForThisTask.Where(b => b.Use).Select(b => b.Name)),
                AnalysisExpanderTitle = "Some analysis properties...",
                SearchModeExpanderTitle = "Some search properties..."
            };
            this.DataContext = DataContextForSearchTaskWindow;
            SearchModifications.Timer.Tick += new EventHandler(TextChangeTimerHandler);
        }

        internal XLSearchTask TheTask { get; private set; }

        private void CheckIfNumber(object sender, TextCompositionEventArgs e)
        {
            e.Handled = !GlobalGuiSettings.CheckIsNumber(e.Text);
        }

        private void PopulateChoices()
        {
            foreach (string crosslinkerName in Enum.GetNames(typeof(CrosslinkerType)))
                cbCrosslinker.Items.Add(crosslinkerName);

            foreach (string dissassociationType in GlobalVariables.AllSupportedDissociationTypes.Keys)
            {
                DissociationTypeComboBox.Items.Add(dissassociationType);
            }

            cbbXLprecusorMsTl.Items.Add("Da");
            cbbXLprecusorMsTl.Items.Add("ppm");

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

            productMassToleranceComboBox.Items.Add("Da");
            productMassToleranceComboBox.Items.Add("ppm");

            foreach (var hm in GlobalVariables.AllModsKnown.GroupBy(b => b.ModificationType))
            {
                var theModType = new ModTypeForTreeView(hm.Key, false);
                FixedModTypeForTreeViewObservableCollection.Add(theModType);
                foreach (var uah in hm)
                    theModType.Children.Add(new ModForTreeView(uah.ToString(), false, uah.IdWithMotif, false, theModType));
            }
            fixedModsTreeView.DataContext = FixedModTypeForTreeViewObservableCollection;
            foreach (var hm in GlobalVariables.AllModsKnown.GroupBy(b => b.ModificationType))
            {
                var theModType = new ModTypeForTreeView(hm.Key, false);
                VariableModTypeForTreeViewObservableCollection.Add(theModType);
                foreach (var uah in hm)
                {
                    theModType.Children.Add(new ModForTreeView(uah.ToString(), false, uah.IdWithMotif, false, theModType));
                }
            }
            variableModsTreeView.DataContext = VariableModTypeForTreeViewObservableCollection;
        }

        private void UpdateFieldsFromTask(XLSearchTask task)
        {
            //Crosslink search para
            //RbSearchCrosslink.IsChecked = !task.XlSearchParameters.SearchGlyco;
            //RbSearchGlyco.IsChecked = task.XlSearchParameters.SearchGlyco;
            //CkbSearchGlycoWithBgYgIndex.IsChecked = task.XlSearchParameters.SearchGlycoWithBgYgIndex;
            cbCrosslinker.SelectedIndex = (int)task.XlSearchParameters.CrosslinkerType;
            ckbXLTopNum.IsChecked = task.XlSearchParameters.RestrictToTopNHits;
            txtXLTopNum.Text = task.XlSearchParameters.CrosslinkSearchTopNum.ToString(CultureInfo.InvariantCulture);
            ckbQuenchH2O.IsChecked = task.XlSearchParameters.XlQuench_H2O;
            ckbQuenchNH2.IsChecked = task.XlSearchParameters.XlQuench_NH2;
            ckbQuenchTris.IsChecked = task.XlSearchParameters.XlQuench_Tris;
            txtUdXLKerName.Text = task.XlSearchParameters.CrosslinkerName;
            ckbUdXLkerCleavable.IsChecked = task.XlSearchParameters.IsCleavable;
            txtUdXLkerTotalMs.Text = task.XlSearchParameters.CrosslinkerTotalMass.HasValue ? 
                task.XlSearchParameters.CrosslinkerTotalMass.Value.ToString(CultureInfo.InvariantCulture) : "";
            txtUdXLkerShortMass.Text = task.XlSearchParameters.CrosslinkerShortMass.HasValue ? 
                task.XlSearchParameters.CrosslinkerShortMass.Value.ToString(CultureInfo.InvariantCulture) : "";
            txtUdXLkerLongMass.Text = task.XlSearchParameters.CrosslinkerLongMass.HasValue ? 
                task.XlSearchParameters.CrosslinkerLongMass.Value.ToString(CultureInfo.InvariantCulture) : "";
            txtH2OQuenchMass.Text = task.XlSearchParameters.CrosslinkerDeadEndMassH2O.HasValue ?
                task.XlSearchParameters.CrosslinkerDeadEndMassH2O.Value.ToString(CultureInfo.InvariantCulture) : "";
            txtNH2QuenchMass.Text = task.XlSearchParameters.CrosslinkerDeadEndMassNH2.HasValue ?
                task.XlSearchParameters.CrosslinkerDeadEndMassNH2.Value.ToString(CultureInfo.InvariantCulture) : "";
            txtTrisQuenchMass.Text = task.XlSearchParameters.CrosslinkerDeadEndMassTris.HasValue ?
                task.XlSearchParameters.CrosslinkerDeadEndMassTris.Value.ToString(CultureInfo.InvariantCulture) : "";

            txtUdXLkerAminoAcids.Text = task.XlSearchParameters.CrosslinkerResidues;
            txtUdXLkerAminoAcids2.Text = task.XlSearchParameters.CrosslinkerResidues2;
            cbbXLprecusorMsTl.SelectedIndex = task.CommonParameters.PrecursorMassTolerance is AbsoluteTolerance ? 0 : 1;
            XLPrecusorMsTlTextBox.Text = task.CommonParameters.PrecursorMassTolerance.Value.ToString(CultureInfo.InvariantCulture);
            trimMs1.IsChecked = task.CommonParameters.TrimMs1Peaks;
            trimMsMs.IsChecked = task.CommonParameters.TrimMsMsPeaks;
            TopNPeaksTextBox.Text = task.CommonParameters.TopNpeaks.ToString(CultureInfo.InvariantCulture);
            MinRatioTextBox.Text = task.CommonParameters.MinRatio.ToString(CultureInfo.InvariantCulture);
            DissociationTypeComboBox.SelectedItem = task.CommonParameters.DissociationType.ToString();

            //ckbCharge_2_3.IsChecked = task.XlSearchParameters.XlCharge_2_3;
            checkBoxDecoy.IsChecked = task.XlSearchParameters.DecoyType != DecoyType.None;
            deconvolutePrecursors.IsChecked = task.CommonParameters.DoPrecursorDeconvolution;
            useProvidedPrecursor.IsChecked = task.CommonParameters.UseProvidedPrecursorInfo;
            missedCleavagesTextBox.Text = task.CommonParameters.DigestionParams.MaxMissedCleavages.ToString(CultureInfo.InvariantCulture);
            MinPeptideLengthTextBox.Text = task.CommonParameters.DigestionParams.MinPeptideLength.ToString(CultureInfo.InvariantCulture);
            MaxPeptideLengthTextBox.Text = task.CommonParameters.DigestionParams.MaxPeptideLength == int.MaxValue ? "" : task.CommonParameters.DigestionParams.MaxPeptideLength.ToString(CultureInfo.InvariantCulture);
            proteaseComboBox.SelectedItem = task.CommonParameters.DigestionParams.Protease;
            maxModificationIsoformsTextBox.Text = task.CommonParameters.DigestionParams.MaxModificationIsoforms.ToString(CultureInfo.InvariantCulture);
            initiatorMethionineBehaviorComboBox.SelectedIndex = (int)task.CommonParameters.DigestionParams.InitiatorMethionineBehavior;
            productMassToleranceTextBox.Text = task.CommonParameters.ProductMassTolerance.Value.ToString(CultureInfo.InvariantCulture);
            productMassToleranceComboBox.SelectedIndex = task.CommonParameters.ProductMassTolerance is AbsoluteTolerance ? 0 : 1;
            DissociationTypeComboBox.SelectedItem = task.CommonParameters.DissociationType.ToString();
            minScoreAllowed.Text = task.CommonParameters.ScoreCutoff.ToString(CultureInfo.InvariantCulture);
            numberOfDatabaseSearchesTextBox.Text = task.CommonParameters.TotalPartitions.ToString(CultureInfo.InvariantCulture);
            maxThreadsTextBox.Text = task.CommonParameters.MaxThreadsToUsePerFile.ToString(CultureInfo.InvariantCulture);
            CustomFragmentationWindow = new CustomFragmentationWindow(task.CommonParameters.CustomIons);
            ckbPercolator.IsChecked = task.XlSearchParameters.WriteOutputForPercolator;
            ckbPepXML.IsChecked = task.XlSearchParameters.WritePepXml;

            OutputFileNameTextBox.Text = task.CommonParameters.TaskDescriptor;

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
            string fieldNotUsed = "1";

            if (!GlobalGuiSettings.CheckTaskSettingsValidity(XLPrecusorMsTlTextBox.Text, productMassToleranceTextBox.Text, missedCleavagesTextBox.Text,
                maxModificationIsoformsTextBox.Text, MinPeptideLengthTextBox.Text, MaxPeptideLengthTextBox.Text, maxThreadsTextBox.Text, minScoreAllowed.Text,
                fieldNotUsed, fieldNotUsed, fieldNotUsed, TopNPeaksTextBox.Text, MinRatioTextBox.Text, numberOfDatabaseSearchesTextBox.Text, fieldNotUsed, fieldNotUsed, fieldNotUsed))
            {
                return;
            }

            DissociationType dissociationType = GlobalVariables.AllSupportedDissociationTypes[DissociationTypeComboBox.SelectedItem.ToString()];
            CustomFragmentationWindow.Close();
            
            //TheTask.XlSearchParameters.SearchGlyco = RbSearchGlyco.IsChecked.Value;
            //TheTask.XlSearchParameters.SearchGlycoWithBgYgIndex = CkbSearchGlycoWithBgYgIndex.IsChecked.Value;
            TheTask.XlSearchParameters.RestrictToTopNHits = ckbXLTopNum.IsChecked.Value;
            TheTask.XlSearchParameters.CrosslinkSearchTopNum = int.Parse(txtXLTopNum.Text, CultureInfo.InvariantCulture);
            TheTask.XlSearchParameters.CrosslinkerType = (CrosslinkerType)cbCrosslinker.SelectedIndex;

            //TheTask.XlSearchParameters.XlCharge_2_3 = ckbCharge_2_3.IsChecked.Value;
            TheTask.XlSearchParameters.XlQuench_H2O = ckbQuenchH2O.IsChecked.Value;
            TheTask.XlSearchParameters.XlQuench_NH2 = ckbQuenchNH2.IsChecked.Value;
            TheTask.XlSearchParameters.XlQuench_Tris = ckbQuenchTris.IsChecked.Value;

            if (TheTask.XlSearchParameters.CrosslinkerType == CrosslinkerType.UserDefined)
            {
                TheTask.XlSearchParameters.CrosslinkerName = txtUdXLKerName.Text;
                TheTask.XlSearchParameters.IsCleavable = ckbUdXLkerCleavable.IsChecked.Value;
                TheTask.XlSearchParameters.CrosslinkerResidues = txtUdXLkerAminoAcids.Text;
                TheTask.XlSearchParameters.CrosslinkerResidues2 = txtUdXLkerAminoAcids2.Text;
                TheTask.XlSearchParameters.CrosslinkerLongMass = string.IsNullOrEmpty(txtUdXLkerLongMass.Text) ?
                    (double?)null : double.Parse(txtUdXLkerLongMass.Text, CultureInfo.InvariantCulture);

                TheTask.XlSearchParameters.CrosslinkerShortMass = string.IsNullOrEmpty(txtUdXLkerShortMass.Text) ?
                    (double?)null : double.Parse(txtUdXLkerShortMass.Text, CultureInfo.InvariantCulture);

                TheTask.XlSearchParameters.CrosslinkerTotalMass = string.IsNullOrEmpty(txtUdXLkerTotalMs.Text) ?
                    (double?)null : double.Parse(txtUdXLkerTotalMs.Text, CultureInfo.InvariantCulture);

                TheTask.XlSearchParameters.CrosslinkerDeadEndMassH2O = string.IsNullOrEmpty(txtH2OQuenchMass.Text) ?
                    (double?)null : double.Parse(txtH2OQuenchMass.Text, CultureInfo.InvariantCulture);

                TheTask.XlSearchParameters.CrosslinkerDeadEndMassNH2 = string.IsNullOrEmpty(txtNH2QuenchMass.Text) ?
                    (double?)null : double.Parse(txtNH2QuenchMass.Text, CultureInfo.InvariantCulture);

                TheTask.XlSearchParameters.CrosslinkerDeadEndMassTris = string.IsNullOrEmpty(txtTrisQuenchMass.Text) ?
                    (double?)null : double.Parse(txtTrisQuenchMass.Text, CultureInfo.InvariantCulture);
            }

            TheTask.XlSearchParameters.DecoyType = checkBoxDecoy.IsChecked.Value ? DecoyType.Reverse : DecoyType.None;

            Protease protease = (Protease)proteaseComboBox.SelectedItem;
            int MaxMissedCleavages = string.IsNullOrEmpty(missedCleavagesTextBox.Text) ? int.MaxValue : (int.Parse(missedCleavagesTextBox.Text, NumberStyles.Any, CultureInfo.InvariantCulture));
            int MinPeptideLength = (int.Parse(MinPeptideLengthTextBox.Text, NumberStyles.Any, CultureInfo.InvariantCulture));
            int MaxPeptideLength = string.IsNullOrEmpty(MaxPeptideLengthTextBox.Text) ? int.MaxValue : (int.Parse(MaxPeptideLengthTextBox.Text, NumberStyles.Any, CultureInfo.InvariantCulture));
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

            Tolerance PrecursorMassTolerance;
            if (cbbXLprecusorMsTl.SelectedIndex == 0)
            {
                PrecursorMassTolerance = new AbsoluteTolerance(double.Parse(XLPrecusorMsTlTextBox.Text, CultureInfo.InvariantCulture));
            }
            else
            {
                PrecursorMassTolerance = new PpmTolerance(double.Parse(XLPrecusorMsTlTextBox.Text, CultureInfo.InvariantCulture));
            }

            TheTask.XlSearchParameters.WriteOutputForPercolator = ckbPercolator.IsChecked.Value;
            TheTask.XlSearchParameters.WritePepXml = ckbPepXML.IsChecked.Value;
            //TheTask.UseProvidedPrecursorInfo = useProvidedPrecursor.IsChecked.Value;

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

            CommonParameters commonParamsToSave = new CommonParameters(
                precursorMassTolerance: PrecursorMassTolerance,
                taskDescriptor: OutputFileNameTextBox.Text != "" ? OutputFileNameTextBox.Text : "XLSearchTask",
                productMassTolerance: ProductMassTolerance,
                doPrecursorDeconvolution: deconvolutePrecursors.IsChecked.Value,
                useProvidedPrecursorInfo: useProvidedPrecursor.IsChecked.Value,
                digestionParams: digestionParamsToSave,
                trimMs1Peaks: trimMs1.IsChecked.Value,
                trimMsMsPeaks: trimMsMs.IsChecked.Value,
                topNpeaks: int.Parse(TopNPeaksTextBox.Text),
                minRatio: double.Parse(MinRatioTextBox.Text),
                dissociationType: dissociationType,
                scoreCutoff: double.Parse(minScoreAllowed.Text, CultureInfo.InvariantCulture),
                totalPartitions: int.Parse(numberOfDatabaseSearchesTextBox.Text, CultureInfo.InvariantCulture),
                listOfModsVariable: listOfModsVariable,
                listOfModsFixed: listOfModsFixed,
                assumeOrphanPeaksAreZ1Fragments: protease.Name != "top-down");

            TheTask.CommonParameters = commonParamsToSave;

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
            if (DissociationTypeComboBox.SelectedItem.ToString().Equals(DissociationType.Custom.ToString()))
            {
                CustomFragmentationWindow.Show();
            }
        }
    }
}