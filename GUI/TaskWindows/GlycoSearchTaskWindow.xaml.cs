﻿using EngineLayer;
using EngineLayer.CrosslinkSearch;
using MassSpectrometry;
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
using System.IO;


namespace MetaMorpheusGUI
{
    /// <summary>
    /// Interaction logic for SearchTaskWindow.xaml
    /// </summary>
    public partial class GlycoSearchTaskWindow : Window
    {
        private readonly DataContextForSearchTaskWindow DataContextForSearchTaskWindow;
        private readonly ObservableCollection<SearchModeForDataGrid> SearchModesForThisTask = new ObservableCollection<SearchModeForDataGrid>();
        private readonly ObservableCollection<ModTypeForTreeView> FixedModTypeForTreeViewObservableCollection = new ObservableCollection<ModTypeForTreeView>();
        private readonly ObservableCollection<ModTypeForTreeView> VariableModTypeForTreeViewObservableCollection = new ObservableCollection<ModTypeForTreeView>();
        private CustomFragmentationWindow CustomFragmentationWindow;

        public GlycoSearchTaskWindow() : this(null)
        {
        }

        public GlycoSearchTaskWindow(GlycoSearchTask task)
        {
            InitializeComponent();
            PopulateChoices();
            TheTask = task ?? new GlycoSearchTask();
            UpdateFieldsFromTask(TheTask);

            if (task == null)
            {
                this.saveButton.Content = "Add GlycoSearch Task";
            }
            DataContextForSearchTaskWindow = new DataContextForSearchTaskWindow()
            {
                ExpanderTitle = string.Join(", ", SearchModesForThisTask.Where(b => b.Use).Select(b => b.Name)),
                AnalysisExpanderTitle = "Some analysis properties...",
                SearchModeExpanderTitle = "Some search properties..."
            };
            this.DataContext = DataContextForSearchTaskWindow;
            SearchModifications.Timer.Tick += new EventHandler(TextChangeTimerHandler);
            base.Closing += this.OnClosing;
        }

        internal GlycoSearchTask TheTask { get; private set; }

        private void CheckIfNumber(object sender, TextCompositionEventArgs e)
        {
            e.Handled = !GlobalGuiSettings.CheckIsNumber(e.Text);
        }

        private void PopulateChoices()
        {
            ChildScanDissociationTypeComboBox.Items.Add("Null");
            foreach (string dissassociationType in GlobalVariables.AllSupportedDissociationTypes.Keys)
            {
                DissociationTypeComboBox.Items.Add(dissassociationType);
                ChildScanDissociationTypeComboBox.Items.Add(dissassociationType);
            }          

            cbbPrecusorMsTl.Items.Add("Da");
            cbbPrecusorMsTl.Items.Add("ppm");

            CmbOGlycanDatabase.ItemsSource = GlobalVariables.OGlycanLocations.Select(p=> Path.GetFileName(p));
            CmbNGlycanDatabase.ItemsSource = GlobalVariables.NGlycanLocations.Select(p => Path.GetFileName(p));

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
                {
                    theModType.Children.Add(new ModForTreeView(uah.ToString(), false, uah.IdWithMotif, false, theModType));
                }
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

        private void UpdateFieldsFromTask(GlycoSearchTask task)
        {
            RbtOGlycoSearch.IsChecked = task._glycoSearchParameters.GlycoSearchType == EngineLayer.GlycoSearch.GlycoSearchType.OGlycanSearch;
            RbtNGlycoSearch.IsChecked = task._glycoSearchParameters.GlycoSearchType == EngineLayer.GlycoSearch.GlycoSearchType.NGlycanSearch;
            Rbt_N_O_GlycoSearch.IsChecked = task._glycoSearchParameters.GlycoSearchType == EngineLayer.GlycoSearch.GlycoSearchType.N_O_GlycanSearch;
            TbMaxOGlycanNum.Text = task._glycoSearchParameters.MaximumOGlycanAllowed.ToString(CultureInfo.InvariantCulture);
            CkbOxoniumIonFilt.IsChecked = task._glycoSearchParameters.OxoniumIonFilt;

            txtTopNum.Text = task._glycoSearchParameters.GlycoSearchTopNum.ToString(CultureInfo.InvariantCulture);
            CmbOGlycanDatabase.SelectedItem = task._glycoSearchParameters.OGlycanDatabasefile;
            CmbNGlycanDatabase.SelectedItem = task._glycoSearchParameters.NGlycanDatabasefile;

            cbbPrecusorMsTl.SelectedIndex = task.CommonParameters.PrecursorMassTolerance is AbsoluteTolerance ? 0 : 1;
            PrecusorMsTlTextBox.Text = task.CommonParameters.PrecursorMassTolerance.Value.ToString(CultureInfo.InvariantCulture);
            trimMs1.IsChecked = task.CommonParameters.TrimMs1Peaks;
            trimMsMs.IsChecked = task.CommonParameters.TrimMsMsPeaks;

            TopNPeaksTextBox.Text = task.CommonParameters.NumberOfPeaksToKeepPerWindow == int.MaxValue || !task.CommonParameters.NumberOfPeaksToKeepPerWindow.HasValue ? "" : task.CommonParameters.NumberOfPeaksToKeepPerWindow.Value.ToString(CultureInfo.InvariantCulture);
            MinRatioTextBox.Text = task.CommonParameters.MinimumAllowedIntensityRatioToBasePeak == double.MaxValue || !task.CommonParameters.MinimumAllowedIntensityRatioToBasePeak.HasValue ? "" : task.CommonParameters.MinimumAllowedIntensityRatioToBasePeak.Value.ToString(CultureInfo.InvariantCulture);

            DissociationTypeComboBox.SelectedItem = task.CommonParameters.DissociationType.ToString();

            ChildScanDissociationTypeComboBox.SelectedItem = "Null";
            if (task.CommonParameters.MS2ChildScanDissociationType != DissociationType.Unknown)
            {
                ChildScanDissociationTypeComboBox.SelectedItem = task.CommonParameters.MS2ChildScanDissociationType.ToString();
            }

            CheckBoxDecoy.IsChecked = task._glycoSearchParameters.DecoyType != DecoyType.None;
            RadioButtonReverseDecoy.IsChecked = task._glycoSearchParameters.DecoyType == DecoyType.Reverse;
            RadioButtonSlideDecoy.IsChecked = task._glycoSearchParameters.DecoyType == DecoyType.Slide;
            deconvolutePrecursors.IsChecked = task.CommonParameters.DoPrecursorDeconvolution;
            useProvidedPrecursor.IsChecked = task.CommonParameters.UseProvidedPrecursorInfo;
            missedCleavagesTextBox.Text = task.CommonParameters.DigestionParams.MaxMissedCleavages.ToString(CultureInfo.InvariantCulture);
            MinPeptideLengthTextBox.Text = task.CommonParameters.DigestionParams.MinPeptideLength.ToString(CultureInfo.InvariantCulture);
            MaxPeptideLengthTextBox.Text = task.CommonParameters.DigestionParams.MaxPeptideLength == int.MaxValue ? "" : task.CommonParameters.DigestionParams.MaxPeptideLength.ToString(CultureInfo.InvariantCulture);
            proteaseComboBox.SelectedItem = task.CommonParameters.DigestionParams.Protease;
            maxModificationIsoformsTextBox.Text = task.CommonParameters.DigestionParams.MaxModificationIsoforms.ToString(CultureInfo.InvariantCulture);
            TxtBoxMaxModPerPep.Text = task.CommonParameters.DigestionParams.MaxModsForPeptide.ToString(CultureInfo.InvariantCulture);
            initiatorMethionineBehaviorComboBox.SelectedIndex = (int)task.CommonParameters.DigestionParams.InitiatorMethionineBehavior;
            productMassToleranceTextBox.Text = task.CommonParameters.ProductMassTolerance.Value.ToString(CultureInfo.InvariantCulture);
            productMassToleranceComboBox.SelectedIndex = task.CommonParameters.ProductMassTolerance is AbsoluteTolerance ? 0 : 1;
            minScoreAllowed.Text = task.CommonParameters.ScoreCutoff.ToString(CultureInfo.InvariantCulture);
            numberOfDatabaseSearchesTextBox.Text = task.CommonParameters.TotalPartitions.ToString(CultureInfo.InvariantCulture);
            maxThreadsTextBox.Text = task.CommonParameters.MaxThreadsToUsePerFile.ToString(CultureInfo.InvariantCulture);
            CustomFragmentationWindow = new CustomFragmentationWindow(task.CommonParameters.CustomIons);

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

            if (!GlobalGuiSettings.CheckTaskSettingsValidity(PrecusorMsTlTextBox.Text, productMassToleranceTextBox.Text, missedCleavagesTextBox.Text,
                maxModificationIsoformsTextBox.Text, MinPeptideLengthTextBox.Text, MaxPeptideLengthTextBox.Text, maxThreadsTextBox.Text, minScoreAllowed.Text,
                fieldNotUsed, fieldNotUsed, fieldNotUsed, TopNPeaksTextBox.Text, MinRatioTextBox.Text, null, null, numberOfDatabaseSearchesTextBox.Text, TxtBoxMaxModPerPep.Text, fieldNotUsed, fieldNotUsed))
            {
                return;
            }

            DissociationType dissociationType = GlobalVariables.AllSupportedDissociationTypes[DissociationTypeComboBox.SelectedItem.ToString()];

            DissociationType childDissociationType = DissociationType.Unknown;
            if (ChildScanDissociationTypeComboBox.SelectedItem.ToString() != "Null")
            {
                childDissociationType = GlobalVariables.AllSupportedDissociationTypes[ChildScanDissociationTypeComboBox.SelectedItem.ToString()];
            }
            CustomFragmentationWindow.Close();


            if (RbtOGlycoSearch.IsChecked.Value)
            {
                TheTask._glycoSearchParameters.GlycoSearchType = EngineLayer.GlycoSearch.GlycoSearchType.OGlycanSearch;
            }
            else if (RbtNGlycoSearch.IsChecked.Value)
            {
                TheTask._glycoSearchParameters.GlycoSearchType = EngineLayer.GlycoSearch.GlycoSearchType.NGlycanSearch;
            }
            else if (Rbt_N_O_GlycoSearch.IsChecked.Value)
            {
                TheTask._glycoSearchParameters.GlycoSearchType = EngineLayer.GlycoSearch.GlycoSearchType.N_O_GlycanSearch;
            }


            TheTask._glycoSearchParameters.OGlycanDatabasefile = CmbOGlycanDatabase.SelectedItem.ToString();
            TheTask._glycoSearchParameters.NGlycanDatabasefile = CmbNGlycanDatabase.SelectedItem.ToString();
            TheTask._glycoSearchParameters.GlycoSearchTopNum = int.Parse(txtTopNum.Text, CultureInfo.InvariantCulture);
            TheTask._glycoSearchParameters.MaximumOGlycanAllowed = int.Parse(TbMaxOGlycanNum.Text, CultureInfo.InvariantCulture);
            TheTask._glycoSearchParameters.OxoniumIonFilt = CkbOxoniumIonFilt.IsChecked.Value;

            if (CheckBoxDecoy.IsChecked.Value)
            {
                if (RadioButtonReverseDecoy.IsChecked.Value)
                {
                    TheTask._glycoSearchParameters.DecoyType = DecoyType.Reverse;
                }
                else //if (radioButtonSlideDecoy.IsChecked.Value)
                {
                    TheTask._glycoSearchParameters.DecoyType = DecoyType.Slide;
                }
            }
            else
            {
                TheTask._glycoSearchParameters.DecoyType = DecoyType.None;
            }

            Protease protease = (Protease)proteaseComboBox.SelectedItem;
            int MaxMissedCleavages = string.IsNullOrEmpty(missedCleavagesTextBox.Text) ? int.MaxValue : (int.Parse(missedCleavagesTextBox.Text, NumberStyles.Any, CultureInfo.InvariantCulture));
            int MinPeptideLength = (int.Parse(MinPeptideLengthTextBox.Text, NumberStyles.Any, CultureInfo.InvariantCulture));
            int MaxPeptideLength = string.IsNullOrEmpty(MaxPeptideLengthTextBox.Text) ? int.MaxValue : (int.Parse(MaxPeptideLengthTextBox.Text, NumberStyles.Any, CultureInfo.InvariantCulture));
            int MaxModificationIsoforms = (int.Parse(maxModificationIsoformsTextBox.Text, CultureInfo.InvariantCulture));
            int MaxModPerPep = (int.Parse(TxtBoxMaxModPerPep.Text, CultureInfo.InvariantCulture));
            InitiatorMethionineBehavior InitiatorMethionineBehavior = ((InitiatorMethionineBehavior)initiatorMethionineBehaviorComboBox.SelectedIndex);
            DigestionParams digestionParamsToSave = new DigestionParams(
                protease: protease.Name,
                maxMissedCleavages: MaxMissedCleavages,
                minPeptideLength: MinPeptideLength,
                maxPeptideLength: MaxPeptideLength,
                maxModificationIsoforms: MaxModificationIsoforms,
                maxModsForPeptides: MaxModPerPep,
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
            if (cbbPrecusorMsTl.SelectedIndex == 0)
            {
                PrecursorMassTolerance = new AbsoluteTolerance(double.Parse(PrecusorMsTlTextBox.Text, CultureInfo.InvariantCulture));
            }
            else
            {
                PrecursorMassTolerance = new PpmTolerance(double.Parse(PrecusorMsTlTextBox.Text, CultureInfo.InvariantCulture));
            }


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
                taskDescriptor: OutputFileNameTextBox.Text != "" ? OutputFileNameTextBox.Text : "GlycoSearchTask",
                productMassTolerance: ProductMassTolerance,
                doPrecursorDeconvolution: deconvolutePrecursors.IsChecked.Value,
                useProvidedPrecursorInfo: useProvidedPrecursor.IsChecked.Value,
                digestionParams: digestionParamsToSave,
                trimMs1Peaks: trimMs1.IsChecked.Value,
                trimMsMsPeaks: trimMsMs.IsChecked.Value,
                numberOfPeaksToKeepPerWindow: int.Parse(TopNPeaksTextBox.Text),
                minimumAllowedIntensityRatioToBasePeak: double.Parse(MinRatioTextBox.Text, CultureInfo.InvariantCulture),
                dissociationType: dissociationType,
                ms2childScanDissociationType: childDissociationType,
                scoreCutoff: double.Parse(minScoreAllowed.Text, CultureInfo.InvariantCulture),
                totalPartitions: int.Parse(numberOfDatabaseSearchesTextBox.Text, CultureInfo.InvariantCulture),
                maxThreadsToUsePerFile: int.Parse(maxThreadsTextBox.Text, CultureInfo.InvariantCulture),
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

        private void OnClosing(object sender, CancelEventArgs e)
        {
            SearchModifications.Timer.Tick -= new EventHandler(TextChangeTimerHandler);
            // remove event handler from timer
            // keeping it will trigger an exception because the closed window stops existing

            CustomFragmentationWindow.Close();
        }

        private void NonSpecificUpdate(object sender, SelectionChangedEventArgs e)
        {
            const int maxLength = 25;
            if (((Protease)proteaseComboBox.SelectedItem).Name.Contains("non-specific"))
            {
                MaxPeptideLengthTextBox.Text = maxLength.ToString();
            }
        }

        private void NonSpecificUpdate(object sender, TextChangedEventArgs e)
        {
            if (((Protease)proteaseComboBox.SelectedItem).Name.Contains("non-specific"))
            {
                try
                {
                    TextBox textBox = (TextBox)sender;
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
    }
}