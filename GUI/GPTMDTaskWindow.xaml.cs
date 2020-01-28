using EngineLayer;
using MassSpectrometry;
using MzLibUtil;
using Nett;
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

namespace MetaMorpheusGUI
{
    /// <summary>
    /// Interaction logic for GptmdTaskWindow.xaml
    /// </summary>
    public partial class GptmdTaskWindow : Window
    {
        private readonly ObservableCollection<ModTypeForTreeView> FixedModTypeForTreeViewObservableCollection = new ObservableCollection<ModTypeForTreeView>();
        private readonly ObservableCollection<ModTypeForTreeView> VariableModTypeForTreeViewObservableCollection = new ObservableCollection<ModTypeForTreeView>();
        private readonly ObservableCollection<ModTypeForLoc> LocalizeModTypeForTreeViewObservableCollection = new ObservableCollection<ModTypeForLoc>();
        private readonly ObservableCollection<ModTypeForTreeView> GptmdModTypeForTreeViewObservableCollection = new ObservableCollection<ModTypeForTreeView>();
        private bool AutomaticallyAskAndOrUpdateParametersBasedOnProtease = true;
        private CustomFragmentationWindow CustomFragmentationWindow;

        public GptmdTaskWindow(GptmdTask myGPTMDtask)
        {
            InitializeComponent();
            TheTask = myGPTMDtask ?? new GptmdTask();

            AutomaticallyAskAndOrUpdateParametersBasedOnProtease = false;
            PopulateChoices();
            UpdateFieldsFromTask(TheTask);
            AutomaticallyAskAndOrUpdateParametersBasedOnProtease = true;

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
                GlobalVariables.StartProcess(hm.Text);
            }
        }

        private void UpdateFieldsFromTask(GptmdTask task)
        {
            ProteaseComboBox.SelectedItem = task.CommonParameters.DigestionParams.Protease; //protease needs to come first or recommended settings can overwrite the actual settings
            UseProvidedPrecursor.IsChecked = task.CommonParameters.UseProvidedPrecursorInfo;
            DeconvolutePrecursors.IsChecked = task.CommonParameters.DoPrecursorDeconvolution;
            DeconvolutionMaxAssumedChargeStateTextBox.Text = task.CommonParameters.DeconvolutionMaxAssumedChargeState.ToString();
            MissedCleavagesTextBox.Text = task.CommonParameters.DigestionParams.MaxMissedCleavages == int.MaxValue ? "" : task.CommonParameters.DigestionParams.MaxMissedCleavages.ToString(CultureInfo.InvariantCulture);
            MinPeptideLengthTextBox.Text = task.CommonParameters.DigestionParams.MinPeptideLength.ToString(CultureInfo.InvariantCulture);
            MaxPeptideLengthTextBox.Text = task.CommonParameters.DigestionParams.MaxPeptideLength == int.MaxValue ? "" : task.CommonParameters.DigestionParams.MaxPeptideLength.ToString(CultureInfo.InvariantCulture);

            MaxModificationIsoformsTextBox.Text = task.CommonParameters.DigestionParams.MaxModificationIsoforms.ToString(CultureInfo.InvariantCulture);
            MaxModsPerPeptideTextBox.Text = task.CommonParameters.DigestionParams.MaxModsForPeptide.ToString(CultureInfo.InvariantCulture);
            InitiatorMethionineBehaviorComboBox.SelectedIndex = (int)task.CommonParameters.DigestionParams.InitiatorMethionineBehavior;
            DissociationTypeComboBox.SelectedItem = task.CommonParameters.DissociationType.ToString();
            ProductMassToleranceTextBox.Text = task.CommonParameters.ProductMassTolerance.Value.ToString(CultureInfo.InvariantCulture);
            ProductMassToleranceComboBox.SelectedIndex = task.CommonParameters.ProductMassTolerance is AbsoluteTolerance ? 0 : 1;
            PrecursorMassToleranceTextBox.Text = task.CommonParameters.PrecursorMassTolerance.Value.ToString(CultureInfo.InvariantCulture);
            PrecursorMassToleranceComboBox.SelectedIndex = task.CommonParameters.PrecursorMassTolerance is AbsoluteTolerance ? 0 : 1;
            MinScoreAllowed.Text = task.CommonParameters.ScoreCutoff.ToString(CultureInfo.InvariantCulture);
            MaxThreadsTextBox.Text = task.CommonParameters.MaxThreadsToUsePerFile.ToString(CultureInfo.InvariantCulture);
            AddCompIonCheckBox.IsChecked = task.CommonParameters.AddCompIons;
            MinVariantDepthTextBox.Text = task.CommonParameters.MinVariantDepth.ToString(CultureInfo.InvariantCulture);
            MaxHeterozygousVariantsTextBox.Text = task.CommonParameters.MaxHeterozygousVariants.ToString(CultureInfo.InvariantCulture);
            CustomFragmentationWindow = new CustomFragmentationWindow(task.CommonParameters.CustomIons);
            OutputFileNameTextBox.Text = task.CommonParameters.TaskDescriptor;

            TrimMs1.IsChecked = task.CommonParameters.TrimMs1Peaks;
            TrimMsMs.IsChecked = task.CommonParameters.TrimMsMsPeaks;
            NumberOfPeaksToKeepPerWindowTextBox.Text = task.CommonParameters.NumberOfPeaksToKeepPerWindow == int.MaxValue || !task.CommonParameters.NumberOfPeaksToKeepPerWindow.HasValue ? "" : task.CommonParameters.NumberOfPeaksToKeepPerWindow.Value.ToString(CultureInfo.InvariantCulture);
            MinimumAllowedIntensityRatioToBasePeakTexBox.Text = task.CommonParameters.MinimumAllowedIntensityRatioToBasePeak == double.MaxValue || !task.CommonParameters.MinimumAllowedIntensityRatioToBasePeak.HasValue ? "" : task.CommonParameters.MinimumAllowedIntensityRatioToBasePeak.Value.ToString(CultureInfo.InvariantCulture);

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
                ModTypeForTreeView theModType = VariableModTypeForTreeViewObservableCollection.FirstOrDefault(b => b.DisplayName.Equals(mod.Item1));
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

            foreach (var mod in task.GptmdParameters.ListOfModsGptmd)
            {
                ModTypeForTreeView theModType = GptmdModTypeForTreeViewObservableCollection.FirstOrDefault(b => b.DisplayName.Equals(mod.Item1));
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
                    GptmdModTypeForTreeViewObservableCollection.Add(theModType);
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

            foreach (var ye in GptmdModTypeForTreeViewObservableCollection)
            {
                ye.VerifyCheckState();
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

            ProductMassToleranceComboBox.Items.Add("Da");
            ProductMassToleranceComboBox.Items.Add("ppm");
            PrecursorMassToleranceComboBox.Items.Add("Da");
            PrecursorMassToleranceComboBox.Items.Add("ppm");

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

            foreach (var hm in GlobalVariables.AllModsKnown.GroupBy(b => b.ModificationType))
            {
                var theModType = new ModTypeForTreeView(hm.Key, false);
                GptmdModTypeForTreeViewObservableCollection.Add(theModType);
                foreach (var uah in hm)
                {
                    theModType.Children.Add(new ModForTreeView(uah.ToString(), false, uah.IdWithMotif, false, theModType));
                }
            }
            gptmdModsTreeView.DataContext = GptmdModTypeForTreeViewObservableCollection;
        }

        private void CancelButton_Click(object sender, RoutedEventArgs e)
        {
            DialogResult = false;
            CustomFragmentationWindow.Close();
        }
        
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
                        }
                        break;
                    case "top-down":
                        if (UpdateGUISettings.UseTopDownRecommendedSettings())
                        {
                            UseProvidedPrecursor.IsChecked = false;
                            DeconvolutionMaxAssumedChargeStateTextBox.Text = "60";
                            TrimMsMs.IsChecked = false;
                            //uncheck all variable mods
                            foreach (var mod in VariableModTypeForTreeViewObservableCollection)
                            {
                                mod.Use = false;
                            }

                            //clear GPTMD mods and replace them with a subset
                            foreach (var mod in GptmdModTypeForTreeViewObservableCollection)
                            {
                                mod.Use = false;
                            }
                            //populate the recommended mods
                            foreach (var mod in UpdateGUISettings.TopDownModsForGPTMD)
                            {
                                ModTypeForTreeView theModType = GptmdModTypeForTreeViewObservableCollection.FirstOrDefault(b => b.DisplayName.Equals(mod.Item1));
                                if (theModType != null)
                                {
                                    var theMod = theModType.Children.FirstOrDefault(b => b.ModName.Equals(mod.Item2));
                                    if (theMod != null)
                                    {
                                        theMod.Use = true;
                                    }
                                }
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
                            }
                        }
                        break;
                    case "elastase":
                        {
                            if (UpdateGUISettings.UseElastaseRecommendedSettings())
                            {
                                MissedCleavagesTextBox.Text = "16";
                            }
                        }
                        break;
                    //nothing to change for semi-trypsin
                    default:
                        break;
                }
            }
        }

        private void ProteaseSpecificUpdate(object sender, TextChangedEventArgs e)
        {
            if (((Protease)ProteaseComboBox.SelectedItem).Name.Contains("non-specific"))
            {
                try
                {
                    TextBox textBox = (TextBox)sender;
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
        }

        private void SaveButton_Click(object sender, RoutedEventArgs e)
        {
            string fieldNotUsed = "1";

            if (!GlobalGuiSettings.CheckTaskSettingsValidity(PrecursorMassToleranceTextBox.Text, ProductMassToleranceTextBox.Text, MissedCleavagesTextBox.Text,
                 MaxModificationIsoformsTextBox.Text, MinPeptideLengthTextBox.Text, MaxPeptideLengthTextBox.Text, MaxThreadsTextBox.Text, MinScoreAllowed.Text,
                fieldNotUsed, fieldNotUsed, DeconvolutionMaxAssumedChargeStateTextBox.Text, NumberOfPeaksToKeepPerWindowTextBox.Text, MinimumAllowedIntensityRatioToBasePeakTexBox.Text, null, null, fieldNotUsed, fieldNotUsed, fieldNotUsed, fieldNotUsed))
            {
                return;
            }

            Protease protease = (Protease)ProteaseComboBox.SelectedItem;
            int maxMissedCleavages = string.IsNullOrEmpty(MissedCleavagesTextBox.Text) ? int.MaxValue : (int.Parse(MissedCleavagesTextBox.Text, NumberStyles.Any, CultureInfo.InvariantCulture));
            int minPeptideLength = int.Parse(MinPeptideLengthTextBox.Text, NumberStyles.Any, CultureInfo.InvariantCulture);
            int maxPeptideLength = string.IsNullOrEmpty(MaxPeptideLengthTextBox.Text) ? int.MaxValue : (int.Parse(MaxPeptideLengthTextBox.Text, NumberStyles.Any, CultureInfo.InvariantCulture));
            int minVariantDepth = int.Parse(MinVariantDepthTextBox.Text, NumberStyles.Any, CultureInfo.InvariantCulture);
            int maxHeterozygousVariants = int.Parse(MaxHeterozygousVariantsTextBox.Text, NumberStyles.Any, CultureInfo.InvariantCulture);
            int maxModificationIsoforms = int.Parse(MaxModificationIsoformsTextBox.Text, CultureInfo.InvariantCulture);
            int maxModsPerPeptide = int.Parse(MaxModsPerPeptideTextBox.Text, CultureInfo.InvariantCulture);
            InitiatorMethionineBehavior initiatorMethionineBehavior = (InitiatorMethionineBehavior)InitiatorMethionineBehaviorComboBox.SelectedIndex;
            DissociationType dissociationType = GlobalVariables.AllSupportedDissociationTypes[DissociationTypeComboBox.SelectedItem.ToString()];
            CustomFragmentationWindow.Close();

            Tolerance productMassTolerance;
            if (ProductMassToleranceComboBox.SelectedIndex == 0)
            {
                productMassTolerance = new AbsoluteTolerance(double.Parse(ProductMassToleranceTextBox.Text, CultureInfo.InvariantCulture));
            }
            else
            {
                productMassTolerance = new PpmTolerance(double.Parse(ProductMassToleranceTextBox.Text, CultureInfo.InvariantCulture));
            }

            Tolerance precursorMassTolerance;
            if (PrecursorMassToleranceComboBox.SelectedIndex == 0)
            {
                precursorMassTolerance = new AbsoluteTolerance(double.Parse(PrecursorMassToleranceTextBox.Text, CultureInfo.InvariantCulture));
            }
            else
            {
                precursorMassTolerance = new PpmTolerance(double.Parse(PrecursorMassToleranceTextBox.Text, CultureInfo.InvariantCulture));
            }

            List<(string, string)> listOfModsVariable = new List<(string, string)>();
            foreach (var heh in VariableModTypeForTreeViewObservableCollection)
            {
                listOfModsVariable.AddRange(heh.Children.Where(b => b.Use).Select(b => (b.Parent.DisplayName, b.ModName)));
            }

            if (!GlobalGuiSettings.VariableModCheck(listOfModsVariable))
            {
                return;
            }

            bool TrimMs1Peaks = TrimMs1.IsChecked.Value;
            bool TrimMsMsPeaks = TrimMsMs.IsChecked.Value;

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

            List<(string, string)> listOfModsFixed = new List<(string, string)>();
            foreach (var heh in FixedModTypeForTreeViewObservableCollection)
            {
                listOfModsFixed.AddRange(heh.Children.Where(b => b.Use).Select(b => (b.Parent.DisplayName, b.ModName)));
            }
            bool parseMaxThreadsPerFile = int.Parse(MaxThreadsTextBox.Text, CultureInfo.InvariantCulture) <= Environment.ProcessorCount && int.Parse(MaxThreadsTextBox.Text, CultureInfo.InvariantCulture) > 0;

            CommonParameters commonParamsToSave = new CommonParameters(
                useProvidedPrecursorInfo: UseProvidedPrecursor.IsChecked.Value,
                deconvolutionMaxAssumedChargeState: int.Parse(DeconvolutionMaxAssumedChargeStateTextBox.Text, CultureInfo.InvariantCulture),
                doPrecursorDeconvolution: DeconvolutePrecursors.IsChecked.Value,
                taskDescriptor: OutputFileNameTextBox.Text != "" ? OutputFileNameTextBox.Text : "GPTMDTask",
                maxThreadsToUsePerFile: parseMaxThreadsPerFile ? int.Parse(MaxThreadsTextBox.Text, CultureInfo.InvariantCulture) : new CommonParameters().MaxThreadsToUsePerFile,
                digestionParams: new DigestionParams(
                    protease: protease.Name,
                    maxMissedCleavages: maxMissedCleavages,
                    minPeptideLength: minPeptideLength,
                    maxPeptideLength: maxPeptideLength,
                    maxModificationIsoforms: maxModificationIsoforms,
                    maxModsForPeptides: maxModsPerPeptide,
                    initiatorMethionineBehavior: initiatorMethionineBehavior),
                    dissociationType: dissociationType,
                    scoreCutoff: double.Parse(MinScoreAllowed.Text, CultureInfo.InvariantCulture),
                    precursorMassTolerance: precursorMassTolerance,
                    productMassTolerance: productMassTolerance,                    
                    trimMs1Peaks: TrimMs1Peaks,
                    trimMsMsPeaks: TrimMsMsPeaks,
                    numberOfPeaksToKeepPerWindow: numPeaksToKeep,
                    minimumAllowedIntensityRatioToBasePeak: minimumAllowedIntensityRatioToBasePeak,
                    listOfModsFixed: listOfModsFixed,
                    listOfModsVariable: listOfModsVariable,
                    assumeOrphanPeaksAreZ1Fragments: protease.Name != "top-down",
                    addCompIons: AddCompIonCheckBox.IsChecked.Value,
                    minVariantDepth: minVariantDepth,
                    maxHeterozygousVariants: maxHeterozygousVariants);

            TheTask.GptmdParameters.ListOfModsGptmd = new List<(string, string)>();
            foreach (var heh in GptmdModTypeForTreeViewObservableCollection)
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
                SearchModifications.FilterTree(SearchFixMod, fixedModsTreeView, FixedModTypeForTreeViewObservableCollection);
                SearchModifications.FixedSearch = false;
            }

            if (SearchModifications.VariableSearch)
            {
                SearchModifications.FilterTree(SearchVarMod, variableModsTreeView, VariableModTypeForTreeViewObservableCollection);
                SearchModifications.VariableSearch = false;
            }

            if (SearchModifications.GptmdSearch)
            {
                SearchModifications.FilterTree(SearchGPTMD, gptmdModsTreeView, GptmdModTypeForTreeViewObservableCollection);
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
            SearchModifications.Timer.Tick -= new EventHandler(TextChangeTimerHandler);
            // remove event handler from timer
            // keeping it will trigger an exception because the closed window stops existing

            CustomFragmentationWindow.Close();
        }

        private void SaveAsDefault_Click(object sender, RoutedEventArgs e)
        {
            SaveButton_Click(sender, e);
            Toml.WriteFile(TheTask, Path.Combine(GlobalVariables.DataDir, "DefaultParameters", @"GptmdTaskDefault.toml"), MetaMorpheusTask.tomlConfig);
        }
    }
}