﻿using EngineLayer;
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

namespace MetaMorpheusGUI
{
    /// <summary>
    /// Interaction logic for GptmdTaskWindow.xaml
    /// </summary>
    public partial class GptmdTaskWindow : Window
    {
        private readonly ObservableCollection<ModTypeForTreeView> fixedModTypeForTreeViewObservableCollection = new ObservableCollection<ModTypeForTreeView>();
        private readonly ObservableCollection<ModTypeForTreeView> variableModTypeForTreeViewObservableCollection = new ObservableCollection<ModTypeForTreeView>();
        private readonly ObservableCollection<ModTypeForLoc> localizeModTypeForTreeViewObservableCollection = new ObservableCollection<ModTypeForLoc>();
        private readonly ObservableCollection<ModTypeForTreeView> gptmdModTypeForTreeViewObservableCollection = new ObservableCollection<ModTypeForTreeView>();
        private CustomFragmentationWindow CustomFragmentationWindow;

        public GptmdTaskWindow() : this(null)
        {
        }

        public GptmdTaskWindow(GptmdTask myGPTMDtask)
        {
            InitializeComponent();
            PopulateChoices();

            TheTask = myGPTMDtask ?? new GptmdTask();
            UpdateFieldsFromTask(TheTask);

            if (myGPTMDtask == null)
            {
                this.saveButton.Content = "Add the GPTMD Task";
            }

            SearchModifications.Timer.Tick += new EventHandler(TextChangeTimerHandler);
            base.Closing += this.OnClosing;
        }

        internal GptmdTask TheTask { get; private set; }

        private void Row_DoubleClick(object sender, MouseButtonEventArgs e)
        {
            var ye = sender as DataGridCell;
            if (ye.Content is TextBlock hm && !string.IsNullOrEmpty(hm.Text))
            {
                System.Diagnostics.Process.Start(hm.Text);
            }
        }

        private void UpdateFieldsFromTask(GptmdTask task)
        {
            useProvidedPrecursor.IsChecked = task.CommonParameters.UseProvidedPrecursorInfo;
            deconvolutePrecursors.IsChecked = task.CommonParameters.DoPrecursorDeconvolution;
            DeconvolutionMaxAssumedChargeStateTextBox.Text = task.CommonParameters.DeconvolutionMaxAssumedChargeState.ToString();
            missedCleavagesTextBox.Text = task.CommonParameters.DigestionParams.MaxMissedCleavages == int.MaxValue ? "" : task.CommonParameters.DigestionParams.MaxMissedCleavages.ToString(CultureInfo.InvariantCulture);
            MinPeptideLengthTextBox.Text = task.CommonParameters.DigestionParams.MinPeptideLength.ToString(CultureInfo.InvariantCulture);
            MaxPeptideLengthTextBox.Text = task.CommonParameters.DigestionParams.MaxPeptideLength == int.MaxValue ? "" : task.CommonParameters.DigestionParams.MaxPeptideLength.ToString(CultureInfo.InvariantCulture);
            proteaseComboBox.SelectedItem = task.CommonParameters.DigestionParams.Protease;
            maxModificationIsoformsTextBox.Text = task.CommonParameters.DigestionParams.MaxModificationIsoforms.ToString(CultureInfo.InvariantCulture);
            initiatorMethionineBehaviorComboBox.SelectedIndex = (int)task.CommonParameters.DigestionParams.InitiatorMethionineBehavior;
            DissociationTypeComboBox.SelectedItem = task.CommonParameters.DissociationType.ToString();
            productMassToleranceTextBox.Text = task.CommonParameters.ProductMassTolerance.Value.ToString(CultureInfo.InvariantCulture);
            productMassToleranceComboBox.SelectedIndex = task.CommonParameters.ProductMassTolerance is AbsoluteTolerance ? 0 : 1;
            precursorMassToleranceTextBox.Text = task.CommonParameters.PrecursorMassTolerance.Value.ToString(CultureInfo.InvariantCulture);
            precursorMassToleranceComboBox.SelectedIndex = task.CommonParameters.PrecursorMassTolerance is AbsoluteTolerance ? 0 : 1;
            minScoreAllowed.Text = task.CommonParameters.ScoreCutoff.ToString(CultureInfo.InvariantCulture);
            maxThreadsTextBox.Text = task.CommonParameters.MaxThreadsToUsePerFile.ToString(CultureInfo.InvariantCulture);
            addCompIonCheckBox.IsChecked = task.CommonParameters.AddCompIons;
            MinVariantDepthTextBox.Text = task.CommonParameters.MinVariantDepth.ToString(CultureInfo.InvariantCulture);
            MaxHeterozygousVariantsTextBox.Text = task.CommonParameters.MaxHeterozygousVariants.ToString(CultureInfo.InvariantCulture);
            CustomFragmentationWindow = new CustomFragmentationWindow(task.CommonParameters.CustomIons);
            OutputFileNameTextBox.Text = task.CommonParameters.TaskDescriptor;

            trimMs1.IsChecked = task.CommonParameters.TrimMs1Peaks;
            trimMsMs.IsChecked = task.CommonParameters.TrimMsMsPeaks;
            NumberOfPeaksToKeepPerWindowTextBox.Text = task.CommonParameters.NumberOfPeaksToKeepPerWindow == int.MaxValue || !task.CommonParameters.NumberOfPeaksToKeepPerWindow.HasValue ? "" : task.CommonParameters.NumberOfPeaksToKeepPerWindow.Value.ToString(CultureInfo.InvariantCulture);
            MinimumAllowedIntensityRatioToBasePeakTexBox.Text = task.CommonParameters.MinimumAllowedIntensityRatioToBasePeak == double.MaxValue || !task.CommonParameters.MinimumAllowedIntensityRatioToBasePeak.HasValue ? "" : task.CommonParameters.MinimumAllowedIntensityRatioToBasePeak.Value.ToString(CultureInfo.InvariantCulture);

            foreach (var mod in task.CommonParameters.ListOfModsFixed)
            {
                var theModType = fixedModTypeForTreeViewObservableCollection.FirstOrDefault(b => b.DisplayName.Equals(mod.Item1));
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
                    fixedModTypeForTreeViewObservableCollection.Add(theModType);
                    theModType.Children.Add(new ModForTreeView("UNKNOWN MODIFICATION!", true, mod.Item2, true, theModType));
                }
            }
            foreach (var mod in task.CommonParameters.ListOfModsVariable)
            {
                var theModType = variableModTypeForTreeViewObservableCollection.FirstOrDefault(b => b.DisplayName.Equals(mod.Item1));
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
                    variableModTypeForTreeViewObservableCollection.Add(theModType);
                    theModType.Children.Add(new ModForTreeView("UNKNOWN MODIFICATION!", true, mod.Item2, true, theModType));
                }
            }

            foreach (var heh in localizeModTypeForTreeViewObservableCollection)
            {
                heh.Use = false;
            }

            foreach (var mod in task.GptmdParameters.ListOfModsGptmd)
            {
                var theModType = gptmdModTypeForTreeViewObservableCollection.FirstOrDefault(b => b.DisplayName.Equals(mod.Item1));
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
                    gptmdModTypeForTreeViewObservableCollection.Add(theModType);
                    theModType.Children.Add(new ModForTreeView("UNKNOWN MODIFICATION!", true, mod.Item2, true, theModType));
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

            foreach (var ye in gptmdModTypeForTreeViewObservableCollection)
            {
                ye.VerifyCheckState();
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
                DissociationTypeComboBox.Items.Add(dissassociationType);
            }

            productMassToleranceComboBox.Items.Add("Da");
            productMassToleranceComboBox.Items.Add("ppm");
            precursorMassToleranceComboBox.Items.Add("Da");
            precursorMassToleranceComboBox.Items.Add("ppm");

            foreach (var hm in GlobalVariables.AllModsKnown.GroupBy(b => b.ModificationType))
            {
                var theModType = new ModTypeForTreeView(hm.Key, false);
                fixedModTypeForTreeViewObservableCollection.Add(theModType);
                foreach (var uah in hm)
                {
                    theModType.Children.Add(new ModForTreeView(uah.ToString(), false, uah.IdWithMotif, false, theModType));
                }
            }
            fixedModsTreeView.DataContext = fixedModTypeForTreeViewObservableCollection;
            foreach (var hm in GlobalVariables.AllModsKnown.GroupBy(b => b.ModificationType))
            {
                var theModType = new ModTypeForTreeView(hm.Key, false);
                variableModTypeForTreeViewObservableCollection.Add(theModType);
                foreach (var uah in hm)
                {
                    theModType.Children.Add(new ModForTreeView(uah.ToString(), false, uah.IdWithMotif, false, theModType));
                }
            }
            variableModsTreeView.DataContext = variableModTypeForTreeViewObservableCollection;

            foreach (var hm in GlobalVariables.AllModsKnown.GroupBy(b => b.ModificationType))
            {
                var theModType = new ModTypeForTreeView(hm.Key, false);
                gptmdModTypeForTreeViewObservableCollection.Add(theModType);
                foreach (var uah in hm)
                {
                    theModType.Children.Add(new ModForTreeView(uah.ToString(), false, uah.IdWithMotif, false, theModType));
                }
            }
            gptmdModsTreeView.DataContext = gptmdModTypeForTreeViewObservableCollection;
        }

        private void CancelButton_Click(object sender, RoutedEventArgs e)
        {
            DialogResult = false;
            CustomFragmentationWindow.Close();
        }

        private void SaveButton_Click(object sender, RoutedEventArgs e)
        {
            string fieldNotUsed = "1";

            if (!GlobalGuiSettings.CheckTaskSettingsValidity(precursorMassToleranceTextBox.Text, productMassToleranceTextBox.Text, missedCleavagesTextBox.Text,
                 maxModificationIsoformsTextBox.Text, MinPeptideLengthTextBox.Text, MaxPeptideLengthTextBox.Text, maxThreadsTextBox.Text, minScoreAllowed.Text,
                fieldNotUsed, fieldNotUsed, DeconvolutionMaxAssumedChargeStateTextBox.Text, NumberOfPeaksToKeepPerWindowTextBox.Text, MinimumAllowedIntensityRatioToBasePeakTexBox.Text, null, null, fieldNotUsed, fieldNotUsed, fieldNotUsed, fieldNotUsed))

            {
                return;
            }

            Protease protease = (Protease)proteaseComboBox.SelectedItem;
            int MaxMissedCleavages = string.IsNullOrEmpty(missedCleavagesTextBox.Text) ? int.MaxValue : (int.Parse(missedCleavagesTextBox.Text, NumberStyles.Any, CultureInfo.InvariantCulture));
            int MinPeptideLength = int.Parse(MinPeptideLengthTextBox.Text, NumberStyles.Any, CultureInfo.InvariantCulture);
            int MaxPeptideLength = string.IsNullOrEmpty(MaxPeptideLengthTextBox.Text) ? int.MaxValue : (int.Parse(MaxPeptideLengthTextBox.Text, NumberStyles.Any, CultureInfo.InvariantCulture));
            int MinVariantDepth = int.Parse(MinVariantDepthTextBox.Text, NumberStyles.Any, CultureInfo.InvariantCulture);
            int MaxHeterozygousVariants = int.Parse(MaxHeterozygousVariantsTextBox.Text, NumberStyles.Any, CultureInfo.InvariantCulture);
            int MaxModificationIsoforms = int.Parse(maxModificationIsoformsTextBox.Text, CultureInfo.InvariantCulture);
            InitiatorMethionineBehavior InitiatorMethionineBehavior = (InitiatorMethionineBehavior)initiatorMethionineBehaviorComboBox.SelectedIndex;
            DissociationType dissociationType = GlobalVariables.AllSupportedDissociationTypes[DissociationTypeComboBox.SelectedItem.ToString()];
            CustomFragmentationWindow.Close();

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

            var listOfModsVariable = new List<(string, string)>();
            foreach (var heh in variableModTypeForTreeViewObservableCollection)
            {
                listOfModsVariable.AddRange(heh.Children.Where(b => b.Use).Select(b => (b.Parent.DisplayName, b.ModName)));
            }

            if (!GlobalGuiSettings.VariableModCheck(listOfModsVariable))
            {
                return;
            }

            bool TrimMs1Peaks = trimMs1.IsChecked.Value;
            bool TrimMsMsPeaks = trimMsMs.IsChecked.Value;

            int? numPeaksToKeep = null;
            if (!string.IsNullOrWhiteSpace(NumberOfPeaksToKeepPerWindowTextBox.Text))
            {
                if (int.TryParse(NumberOfPeaksToKeepPerWindowTextBox.Text, out int numberOfPeaksToKeeep))
                {
                    numPeaksToKeep = numberOfPeaksToKeeep;
                }
                else
                {
                    MessageBox.Show("The value that you entered for number of peaks to keep is not acceptable. Try again.");
                    return;
                }
            }

            double? minimumAllowedIntensityRatioToBasePeak = null;
            if (!string.IsNullOrWhiteSpace(MinimumAllowedIntensityRatioToBasePeakTexBox.Text))
            {
                if (double.TryParse(MinimumAllowedIntensityRatioToBasePeakTexBox.Text, out double minimumAllowedIntensityRatio))
                {
                    minimumAllowedIntensityRatioToBasePeak = minimumAllowedIntensityRatio;
                }
                else
                {
                    MessageBox.Show("The value that you entered for minimum allowed intensity ratio to keep is not acceptable. Try again.");
                    return;
                }
            }

            var listOfModsFixed = new List<(string, string)>();
            foreach (var heh in fixedModTypeForTreeViewObservableCollection)
            {
                listOfModsFixed.AddRange(heh.Children.Where(b => b.Use).Select(b => (b.Parent.DisplayName, b.ModName)));
            }
            bool parseMaxThreadsPerFile = int.Parse(maxThreadsTextBox.Text, CultureInfo.InvariantCulture) <= Environment.ProcessorCount && int.Parse(maxThreadsTextBox.Text, CultureInfo.InvariantCulture) > 0;

            CommonParameters commonParamsToSave = new CommonParameters(
                useProvidedPrecursorInfo: useProvidedPrecursor.IsChecked.Value,
                deconvolutionMaxAssumedChargeState: int.Parse(DeconvolutionMaxAssumedChargeStateTextBox.Text, CultureInfo.InvariantCulture),
                doPrecursorDeconvolution: deconvolutePrecursors.IsChecked.Value,
                taskDescriptor: OutputFileNameTextBox.Text != "" ? OutputFileNameTextBox.Text : "GPTMDTask",
                maxThreadsToUsePerFile: parseMaxThreadsPerFile ? int.Parse(maxThreadsTextBox.Text, CultureInfo.InvariantCulture) : new CommonParameters().MaxThreadsToUsePerFile,
                digestionParams: new DigestionParams(
                    protease: protease.Name,
                    maxMissedCleavages: MaxMissedCleavages,
                    minPeptideLength: MinPeptideLength,
                    maxPeptideLength: MaxPeptideLength,
                    maxModificationIsoforms: MaxModificationIsoforms,
                    initiatorMethionineBehavior: InitiatorMethionineBehavior),
                    dissociationType: dissociationType,
                    scoreCutoff: double.Parse(minScoreAllowed.Text, CultureInfo.InvariantCulture),
                    precursorMassTolerance: PrecursorMassTolerance,
                    productMassTolerance: ProductMassTolerance,                    
                    trimMs1Peaks: TrimMs1Peaks,
                    trimMsMsPeaks: TrimMsMsPeaks,
                    numberOfPeaksToKeepPerWindow: numPeaksToKeep,
                    minimumAllowedIntensityRatioToBasePeak: minimumAllowedIntensityRatioToBasePeak,
                    listOfModsFixed: listOfModsFixed,
                    listOfModsVariable: listOfModsVariable,
                    assumeOrphanPeaksAreZ1Fragments: protease.Name != "top-down",
                    addCompIons: addCompIonCheckBox.IsChecked.Value,
                    minVariantDepth: MinVariantDepth,
                    maxHeterozygousVariants: MaxHeterozygousVariants);

            TheTask.GptmdParameters.ListOfModsGptmd = new List<(string, string)>();
            foreach (var heh in gptmdModTypeForTreeViewObservableCollection)
            {
                TheTask.GptmdParameters.ListOfModsGptmd.AddRange(heh.Children.Where(b => b.Use).Select(b => (b.Parent.DisplayName, b.ModName)));
            }

            TheTask.CommonParameters = commonParamsToSave;

            DialogResult = true;
        }

        private void CheckIfNumber(object sender, TextCompositionEventArgs e)
        {
            e.Handled = !GlobalGuiSettings.CheckIsNumber(e.Text);
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

        private void TextChanged_GPTMD(object sender, TextChangedEventArgs args)
        {
            SearchModifications.SetTimer();
            SearchModifications.GptmdSearch = true;
        }

        private void TextChangeTimerHandler(object sender, EventArgs e)
        {
            if (SearchModifications.FixedSearch)
            {
                SearchModifications.FilterTree(SearchFixMod, fixedModsTreeView, fixedModTypeForTreeViewObservableCollection);
                SearchModifications.FixedSearch = false;
            }

            if (SearchModifications.VariableSearch)
            {
                SearchModifications.FilterTree(SearchVarMod, variableModsTreeView, variableModTypeForTreeViewObservableCollection);
                SearchModifications.VariableSearch = false;
            }

            if (SearchModifications.GptmdSearch)
            {
                SearchModifications.FilterTree(SearchGPTMD, gptmdModsTreeView, gptmdModTypeForTreeViewObservableCollection);
                SearchModifications.GptmdSearch = false;
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
            CustomFragmentationWindow.Close();
        }
    }
}