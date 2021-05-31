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
    /// Interaction logic for CalibrateTaskWindow.xaml
    /// </summary>
    public partial class CalibrateTaskWindow : Window
    {
        private readonly ObservableCollection<ModTypeForTreeView> FixedModTypeForTreeViewObservableCollection = new ObservableCollection<ModTypeForTreeView>();
        private readonly ObservableCollection<ModTypeForTreeView> VariableModTypeForTreeViewObservableCollection = new ObservableCollection<ModTypeForTreeView>();
        private readonly ObservableCollection<ModTypeForLoc> LocalizeModTypeForTreeViewObservableCollection = new ObservableCollection<ModTypeForLoc>();
        private bool AutomaticallyAskAndOrUpdateParametersBasedOnProtease = true;
        private CustomFragmentationWindow CustomFragmentationWindow;

        public CalibrateTaskWindow(CalibrationTask myCalibrateTask)
        {
            InitializeComponent();
            TheTask = myCalibrateTask ?? new CalibrationTask();

            AutomaticallyAskAndOrUpdateParametersBasedOnProtease = false;
            PopulateChoices();
            UpdateFieldsFromTask(TheTask);
            AutomaticallyAskAndOrUpdateParametersBasedOnProtease = true;

            if (myCalibrateTask == null)
            {
                SaveButton.Content = "Add the Calibration Task";
            }
            SearchModifications.Timer.Tick += new EventHandler(TextChangeTimerHandler);
            base.Closing += this.OnClosing;
        }

        internal CalibrationTask TheTask { get; private set; }

        private void UpdateFieldsFromTask(CalibrationTask task)
        {
            ProteaseComboBox.SelectedItem = task.CommonParameters.DigestionParams.Protease; //protease needs to come first or recommended settings can overwrite the actual settings
            MissedCleavagesTextBox.Text = task.CommonParameters.DigestionParams.MaxMissedCleavages == int.MaxValue ? "" : task.CommonParameters.DigestionParams.MaxMissedCleavages.ToString(CultureInfo.InvariantCulture);
            MinPeptideLengthTextBox.Text = task.CommonParameters.DigestionParams.MinPeptideLength.ToString(CultureInfo.InvariantCulture);
            MaxPeptideLengthTextBox.Text = task.CommonParameters.DigestionParams.MaxPeptideLength == int.MaxValue ? "" : task.CommonParameters.DigestionParams.MaxPeptideLength.ToString(CultureInfo.InvariantCulture);

            MaxModsPerPeptideTextBox.Text = task.CommonParameters.DigestionParams.MaxModsForPeptide.ToString(CultureInfo.InvariantCulture);
            MaxModificationIsoformsTextBox.Text = task.CommonParameters.DigestionParams.MaxModificationIsoforms.ToString(CultureInfo.InvariantCulture);
            InitiatorMethionineBehaviorComboBox.SelectedIndex = (int)task.CommonParameters.DigestionParams.InitiatorMethionineBehavior;
            DissociationTypeComboBox.SelectedItem = task.CommonParameters.DissociationType.ToString();
            MaxThreadsTextBox.Text = task.CommonParameters.MaxThreadsToUsePerFile.ToString(CultureInfo.InvariantCulture);
            MinVariantDepthTextBox.Text = task.CommonParameters.MinVariantDepth.ToString(CultureInfo.InvariantCulture);
            MaxHeterozygousVariantsTextBox.Text = task.CommonParameters.MaxHeterozygousVariants.ToString(CultureInfo.InvariantCulture);

            ProductMassToleranceTextBox.Text = task.CommonParameters.ProductMassTolerance.Value.ToString(CultureInfo.InvariantCulture);
            ProductMassToleranceComboBox.SelectedIndex = task.CommonParameters.ProductMassTolerance is AbsoluteTolerance ? 0 : 1;
            PrecursorMassToleranceTextBox.Text = task.CommonParameters.PrecursorMassTolerance.Value.ToString(CultureInfo.InvariantCulture);
            PrecursorMassToleranceComboBox.SelectedIndex = task.CommonParameters.PrecursorMassTolerance is AbsoluteTolerance ? 0 : 1;
            CustomFragmentationWindow = new CustomFragmentationWindow(task.CommonParameters.CustomIons);
            writeIndexMzmlCheckbox.IsChecked = task.CalibrationParameters.WriteIndexedMzml;

            //writeIntermediateFilesCheckBox.IsChecked = task.CalibrationParameters.WriteIntermediateFiles;

            MinScoreAllowed.Text = task.CommonParameters.ScoreCutoff.ToString(CultureInfo.InvariantCulture);

            OutputFileNameTextBox.Text = task.CommonParameters.TaskDescriptor;

            foreach (var mod in task.CommonParameters.ListOfModsFixed)
            {
                var theModType = FixedModTypeForTreeViewObservableCollection.FirstOrDefault(b => b.DisplayName.Equals(mod.Item1));
                if (theModType != null)
                {
                    var theMod = theModType.Children.FirstOrDefault(b => b.ModName.Equals(mod.Item2));
                    if (theMod != null)
                        theMod.Use = true;
                    else
                        theModType.Children.Add(new ModForTreeView("UNKNOWN MODIFICATION!", true, mod.Item2, true, theModType));
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
                        theMod.Use = true;
                    else
                        theModType.Children.Add(new ModForTreeView("UNKNOWN MODIFICATION!", true, mod.Item2, true, theModType));
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

            FixedModsTreeView.DataContext = FixedModTypeForTreeViewObservableCollection;

            foreach (var hm in GlobalVariables.AllModsKnown.GroupBy(b => b.ModificationType))
            {
                var theModType = new ModTypeForTreeView(hm.Key, false);
                VariableModTypeForTreeViewObservableCollection.Add(theModType);

                foreach (var uah in hm)
                {
                    theModType.Children.Add(new ModForTreeView(uah.ToString(), false, uah.IdWithMotif, false, theModType));
                }
            }

            VariableModsTreeView.DataContext = VariableModTypeForTreeViewObservableCollection;
        }

        private void CancelButton_Click(object sender, RoutedEventArgs e)
        {
            DialogResult = false;
            CustomFragmentationWindow.Close();
        }

        private void SaveButton_Click(object sender, RoutedEventArgs e)
        {
            string fieldNotUsed = "1";

            if (!GlobalGuiSettings.CheckTaskSettingsValidity(PrecursorMassToleranceTextBox.Text, ProductMassToleranceTextBox.Text, MissedCleavagesTextBox.Text,
                 MaxModificationIsoformsTextBox.Text, MinPeptideLengthTextBox.Text, MaxPeptideLengthTextBox.Text, MaxThreadsTextBox.Text, MinScoreAllowed.Text,
                 fieldNotUsed, fieldNotUsed, fieldNotUsed, fieldNotUsed, fieldNotUsed, null, null, fieldNotUsed, MaxModsPerPeptideTextBox.Text, fieldNotUsed, fieldNotUsed))
            {
                return;
            }

            Protease protease = (Protease)ProteaseComboBox.SelectedItem;
            int maxMissedCleavages = string.IsNullOrEmpty(MissedCleavagesTextBox.Text) ? int.MaxValue : (int.Parse(MissedCleavagesTextBox.Text, CultureInfo.InvariantCulture));
            int minPeptideLength = int.Parse(MinPeptideLengthTextBox.Text, NumberStyles.Any, CultureInfo.InvariantCulture);
            int maxPeptideLength = string.IsNullOrEmpty(MaxPeptideLengthTextBox.Text) ? int.MaxValue : (int.Parse(MaxPeptideLengthTextBox.Text, NumberStyles.Any, CultureInfo.InvariantCulture));
            int minVariantDepth = int.Parse(MinVariantDepthTextBox.Text, NumberStyles.Any, CultureInfo.InvariantCulture);
            int maxHeterozygousVariants = int.Parse(MaxHeterozygousVariantsTextBox.Text, NumberStyles.Any, CultureInfo.InvariantCulture);
            int maxModificationIsoforms = int.Parse(MaxModificationIsoformsTextBox.Text, CultureInfo.InvariantCulture);
            int maxModsPerPeptide = int.Parse(MaxModsPerPeptideTextBox.Text, CultureInfo.InvariantCulture);
            DissociationType dissociationType = GlobalVariables.AllSupportedDissociationTypes[DissociationTypeComboBox.SelectedItem.ToString()];
            CustomFragmentationWindow.Close();

            DigestionParams digestionParamsToSave = new DigestionParams(
                protease: protease.Name,
                maxMissedCleavages: maxMissedCleavages,
                minPeptideLength: minPeptideLength,
                maxPeptideLength: maxPeptideLength,
                maxModificationIsoforms: maxModificationIsoforms,
                maxModsForPeptides: maxModsPerPeptide);

            var listOfModsVariable = new List<(string, string)>();
            foreach (var heh in VariableModTypeForTreeViewObservableCollection)
            {
                listOfModsVariable.AddRange(heh.Children.Where(b => b.Use).Select(b => (b.Parent.DisplayName, b.ModName)));
            }

            if (!GlobalGuiSettings.VariableModCheck(listOfModsVariable))
            {
                return;
            }

            var listOfModsFixed = new List<(string, string)>();
            foreach (var heh in FixedModTypeForTreeViewObservableCollection)
            {
                listOfModsFixed.AddRange(heh.Children.Where(b => b.Use).Select(b => (b.Parent.DisplayName, b.ModName)));
            }
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

            bool parseMaxThreadsPerFile = int.Parse(MaxThreadsTextBox.Text, CultureInfo.InvariantCulture) <= Environment.ProcessorCount && int.Parse(MaxThreadsTextBox.Text, CultureInfo.InvariantCulture) > 0;

            //the below parameters are optimized for top-down but do not exist in the GUI as of Nov. 13, 2019
            if (((Protease)ProteaseComboBox.SelectedItem).Name.Contains("top-down"))
            {
                CommonParameters commonParamsToSave = new CommonParameters(
                    taskDescriptor: OutputFileNameTextBox.Text != "" ? OutputFileNameTextBox.Text : "CalibrateTask",
                    maxThreadsToUsePerFile: parseMaxThreadsPerFile ? int.Parse(MaxThreadsTextBox.Text, CultureInfo.InvariantCulture) : new CommonParameters().MaxThreadsToUsePerFile,
                    digestionParams: digestionParamsToSave,
                    dissociationType: dissociationType,
                    scoreCutoff: double.Parse(MinScoreAllowed.Text, CultureInfo.InvariantCulture),
                    listOfModsFixed: listOfModsFixed,
                    listOfModsVariable: listOfModsVariable,
                    productMassTolerance: productMassTolerance,
                    precursorMassTolerance: precursorMassTolerance,
                    assumeOrphanPeaksAreZ1Fragments: protease.Name != "top-down",
                    minVariantDepth: minVariantDepth,
                    maxHeterozygousVariants: maxHeterozygousVariants,
                    useProvidedPrecursorInfo: false, //Updated
                    deconvolutionMaxAssumedChargeState: 60, //Updated
                    trimMsMsPeaks: false); //Updated
                TheTask.CommonParameters = commonParamsToSave;
            }
            else //bottom-up
            {
                CommonParameters commonParamsToSave = new CommonParameters(
                    taskDescriptor: OutputFileNameTextBox.Text != "" ? OutputFileNameTextBox.Text : "CalibrateTask",
                    maxThreadsToUsePerFile: parseMaxThreadsPerFile ? int.Parse(MaxThreadsTextBox.Text, CultureInfo.InvariantCulture) : new CommonParameters().MaxThreadsToUsePerFile,
                    digestionParams: digestionParamsToSave,
                    dissociationType: dissociationType,
                    scoreCutoff: double.Parse(MinScoreAllowed.Text, CultureInfo.InvariantCulture),
                    listOfModsFixed: listOfModsFixed,
                    listOfModsVariable: listOfModsVariable,
                    productMassTolerance: productMassTolerance,
                    precursorMassTolerance: precursorMassTolerance,
                    assumeOrphanPeaksAreZ1Fragments: protease.Name != "top-down",
                    minVariantDepth: minVariantDepth,
                    maxHeterozygousVariants: maxHeterozygousVariants);
                TheTask.CommonParameters = commonParamsToSave;
            }

            TheTask.CalibrationParameters.WriteIndexedMzml = writeIndexMzmlCheckbox.IsChecked.Value;
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

        private void TextChangeTimerHandler(object sender, EventArgs e)
        {
            if (SearchModifications.FixedSearch)
            {
                SearchModifications.FilterTree(SearchFixMod, FixedModsTreeView, FixedModTypeForTreeViewObservableCollection);
                SearchModifications.FixedSearch = false;
            }

            if (SearchModifications.VariableSearch)
            {
                SearchModifications.FilterTree(SearchVarMod, VariableModsTreeView, VariableModTypeForTreeViewObservableCollection);
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

        private void SaveAsDefault_Click(object sender, RoutedEventArgs e)
        {
            SaveButton_Click(sender, e);
            Toml.WriteFile(TheTask, Path.Combine(GlobalVariables.DataDir, "DefaultParameters", @"CalibrationTaskDefault.toml"), MetaMorpheusTask.tomlConfig);
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
                            //many variables are not present in the calibrate task gui, but we modify them when saving
                            //uncheck all variable mods
                            foreach (var mod in VariableModTypeForTreeViewObservableCollection)
                            {
                                mod.Use = false;
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
    }
}