﻿using EngineLayer;
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

namespace MetaMorpheusGUI
{
    /// <summary>
    /// Interaction logic for GptmdTaskWindow.xaml
    /// </summary>
    public partial class GptmdTaskWindow : Window
    {
        #region Private Fields

        private readonly ObservableCollection<ModTypeForTreeView> fixedModTypeForTreeViewObservableCollection = new ObservableCollection<ModTypeForTreeView>();
        private readonly ObservableCollection<ModTypeForTreeView> variableModTypeForTreeViewObservableCollection = new ObservableCollection<ModTypeForTreeView>();
        private readonly ObservableCollection<ModTypeForLoc> localizeModTypeForTreeViewObservableCollection = new ObservableCollection<ModTypeForLoc>();
        private readonly ObservableCollection<ModTypeForTreeView> gptmdModTypeForTreeViewObservableCollection = new ObservableCollection<ModTypeForTreeView>();

        #endregion Private Fields

        #region Public Constructors

        public GptmdTaskWindow()
        {
            InitializeComponent();
            PopulateChoices();

            TheTask = new GptmdTask();
            UpdateFieldsFromTask(TheTask);

            this.saveButton.Content = "Add the GPTMD Task";
        }

        public GptmdTaskWindow(GptmdTask myGPTMDtask)
        {
            InitializeComponent();
            PopulateChoices();

            TheTask = myGPTMDtask;
            UpdateFieldsFromTask(TheTask);
        }

        #endregion Public Constructors

        #region Internal Properties

        internal GptmdTask TheTask { get; private set; }

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

        private void UpdateFieldsFromTask(GptmdTask task)
        {
            missedCleavagesTextBox.Text = task.CommonParameters.DigestionParams.MaxMissedCleavages.ToString(CultureInfo.InvariantCulture);
            txtMinPeptideLength.Text = task.CommonParameters.DigestionParams.MinPeptideLength.HasValue ? task.CommonParameters.DigestionParams.MinPeptideLength.Value.ToString(CultureInfo.InvariantCulture) : "";
            txtMaxPeptideLength.Text = task.CommonParameters.DigestionParams.MaxPeptideLength.HasValue ? task.CommonParameters.DigestionParams.MaxPeptideLength.Value.ToString(CultureInfo.InvariantCulture) : "";
            proteaseComboBox.SelectedItem = task.CommonParameters.DigestionParams.Protease;
            maxModificationIsoformsTextBox.Text = task.CommonParameters.DigestionParams.MaxModificationIsoforms.ToString(CultureInfo.InvariantCulture);
            initiatorMethionineBehaviorComboBox.SelectedIndex = (int)task.CommonParameters.DigestionParams.InitiatorMethionineBehavior;

            productMassToleranceTextBox.Text = task.CommonParameters.ProductMassTolerance.Value.ToString(CultureInfo.InvariantCulture);
            productMassToleranceComboBox.SelectedIndex = task.CommonParameters.ProductMassTolerance is AbsoluteTolerance ? 0 : 1;
            precursorMassToleranceTextBox.Text = task.CommonParameters.PrecursorMassTolerance.Value.ToString(CultureInfo.InvariantCulture);
            precursorMassToleranceComboBox.SelectedIndex = task.CommonParameters.PrecursorMassTolerance is AbsoluteTolerance ? 0 : 1;

            //maxDegreesOfParallelism.Text = task.CommonParameters.MaxParallelFilesToAnalyze.ToString();
            bCheckBox.IsChecked = task.CommonParameters.BIons;
            yCheckBox.IsChecked = task.CommonParameters.YIons;
            cCheckBox.IsChecked = task.CommonParameters.CIons;
            zdotCheckBox.IsChecked = task.CommonParameters.ZdotIons;
            //conserveMemoryCheckBox.IsChecked = task.CommonParameters.ConserveMemory;
            minScoreAllowed.Text = task.CommonParameters.ScoreCutoff.ToString(CultureInfo.InvariantCulture);

            OutputFileNameTextBox.Text = task.CommonParameters.TaskDescriptor;

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

            //localizeAllCheckBox.IsChecked = task.CommonParameters.LocalizeAll;
            //if (task.CommonParameters.LocalizeAll)
            {
                foreach (var heh in localizeModTypeForTreeViewObservableCollection)
                    heh.Use = false;
            }
            //else
            //{
            //    foreach (var heh in localizeModTypeForTreeViewObservableCollection)
            //    {
            //        if (task.CommonParameters.ListOfModTypesLocalize.Contains(heh.DisplayName))
            //            heh.Use = true;
            //        else
            //            heh.Use = false;
            //    }
            //}

            foreach (var mod in task.GptmdParameters.ListOfModsGptmd)
            {
                var theModType = gptmdModTypeForTreeViewObservableCollection.FirstOrDefault(b => b.DisplayName.Equals(mod.Item1));
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
                    gptmdModTypeForTreeViewObservableCollection.Add(theModType);
                    theModType.Children.Add(new ModForTreeView("UNKNOWN MODIFICATION!", true, mod.Item2, true, theModType));
                }
            }

            foreach (var ye in variableModTypeForTreeViewObservableCollection)
                ye.VerifyCheckState();
            foreach (var ye in fixedModTypeForTreeViewObservableCollection)
                ye.VerifyCheckState();
            foreach (var ye in gptmdModTypeForTreeViewObservableCollection)
                ye.VerifyCheckState();
        }

        private void PopulateChoices()
        {
            foreach (Protease protease in GlobalVariables.ProteaseDictionary.Values)
                proteaseComboBox.Items.Add(protease);
            proteaseComboBox.SelectedIndex = 12;

            foreach (string initiatior_methionine_behavior in Enum.GetNames(typeof(InitiatorMethionineBehavior)))
                initiatorMethionineBehaviorComboBox.Items.Add(initiatior_methionine_behavior);

            productMassToleranceComboBox.Items.Add("Absolute");
            productMassToleranceComboBox.Items.Add("ppm");
            precursorMassToleranceComboBox.Items.Add("Absolute");
            precursorMassToleranceComboBox.Items.Add("ppm");

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

            //foreach (var hm in GlobalVariables.AllModsKnown.GroupBy(b => b.modificationType))
            //    localizeModTypeForTreeViewObservableCollection.Add(new ModTypeForLoc(hm.Key));
            //localizeModsTreeView.DataContext = localizeModTypeForTreeViewObservableCollection;

            foreach (var hm in GlobalVariables.AllModsKnown.GroupBy(b => b.modificationType))
            {
                var theModType = new ModTypeForTreeView(hm.Key, false);
                gptmdModTypeForTreeViewObservableCollection.Add(theModType);
                foreach (var uah in hm)
                    theModType.Children.Add(new ModForTreeView(uah.ToString(), false, uah.id, false, theModType));
            }
            gptmdModsTreeView.DataContext = gptmdModTypeForTreeViewObservableCollection;
        }

        private void CancelButton_Click(object sender, RoutedEventArgs e)
        {
            DialogResult = false;
        }

        private void SaveButton_Click(object sender, RoutedEventArgs e)
        {
            CommonParameters CommonParamsToSave = new CommonParameters();

            DigestionParams digestionParamsToSave = new DigestionParams();
            digestionParamsToSave.MaxMissedCleavages = int.Parse(missedCleavagesTextBox.Text, CultureInfo.InvariantCulture);
            digestionParamsToSave.MinPeptideLength = int.TryParse(txtMinPeptideLength.Text, NumberStyles.Any, CultureInfo.InvariantCulture, out int temp) ? (int?)temp : null;
            digestionParamsToSave.MaxPeptideLength = int.TryParse(txtMaxPeptideLength.Text, NumberStyles.Any, CultureInfo.InvariantCulture, out temp) ? (int?)temp : null;
            digestionParamsToSave.Protease = (Protease)proteaseComboBox.SelectedItem;
            digestionParamsToSave.MaxModificationIsoforms = int.Parse(maxModificationIsoformsTextBox.Text, CultureInfo.InvariantCulture);
            digestionParamsToSave.InitiatorMethionineBehavior = (InitiatorMethionineBehavior)initiatorMethionineBehaviorComboBox.SelectedIndex;
            CommonParamsToSave.DigestionParams = digestionParamsToSave;

            if (OutputFileNameTextBox.Text != "")
                CommonParamsToSave.TaskDescriptor = OutputFileNameTextBox.Text;
            else
                CommonParamsToSave.TaskDescriptor = "GPTMDTask";

            if (productMassToleranceComboBox.SelectedIndex == 0)
                CommonParamsToSave.ProductMassTolerance = new AbsoluteTolerance(double.Parse(productMassToleranceTextBox.Text, CultureInfo.InvariantCulture));
            else
                CommonParamsToSave.ProductMassTolerance = new PpmTolerance(double.Parse(productMassToleranceTextBox.Text, CultureInfo.InvariantCulture));

            if (precursorMassToleranceComboBox.SelectedIndex == 0)
                CommonParamsToSave.PrecursorMassTolerance = new AbsoluteTolerance(double.Parse(precursorMassToleranceTextBox.Text, CultureInfo.InvariantCulture));
            else
                CommonParamsToSave.PrecursorMassTolerance = new PpmTolerance(double.Parse(precursorMassToleranceTextBox.Text, CultureInfo.InvariantCulture));

            CommonParamsToSave.BIons = bCheckBox.IsChecked.Value;
            CommonParamsToSave.YIons = yCheckBox.IsChecked.Value;
            CommonParamsToSave.CIons = cCheckBox.IsChecked.Value;
            CommonParamsToSave.ZdotIons = zdotCheckBox.IsChecked.Value;
            //CommonParamsToSave.ConserveMemory = conserveMemoryCheckBox.IsChecked.Value;
            CommonParamsToSave.ScoreCutoff = double.Parse(minScoreAllowed.Text, CultureInfo.InvariantCulture);

            var listOfModsVariable = new List<(string, string)>();
            foreach (var heh in variableModTypeForTreeViewObservableCollection)
                listOfModsVariable.AddRange(heh.Children.Where(b => b.Use).Select(b => (b.Parent.DisplayName, b.DisplayName)));
            CommonParamsToSave.ListOfModsVariable = listOfModsVariable;

            var listOfModsFixed = new List<(string, string)>();
            foreach (var heh in fixedModTypeForTreeViewObservableCollection)
                listOfModsFixed.AddRange(heh.Children.Where(b => b.Use).Select(b => (b.Parent.DisplayName, b.DisplayName)));
            CommonParamsToSave.ListOfModsFixed = listOfModsFixed;

            //if (localizeAllCheckBox.IsChecked.Value)
            {
                CommonParamsToSave.ListOfModTypesLocalize = null;
                CommonParamsToSave.LocalizeAll = true;
            }
            //else
            //{
            //    CommonParamsToSave.LocalizeAll = false;
            //    CommonParamsToSave.ListOfModTypesLocalize = localizeModTypeForTreeViewObservableCollection.Where(b => b.Use.HasValue && b.Use.Value).Select(b => b.DisplayName).ToList();
            //}
            TheTask.GptmdParameters.ListOfModsGptmd = new List<(string, string)>();
            foreach (var heh in gptmdModTypeForTreeViewObservableCollection)
                TheTask.GptmdParameters.ListOfModsGptmd.AddRange(heh.Children.Where(b => b.Use).Select(b => (b.Parent.DisplayName, b.DisplayName)));

            //if (int.TryParse(maxDegreesOfParallelism.Text, out int jsakdf))
            //    CommonParamsToSave.MaxParallelFilesToAnalyze = jsakdf;

            TheTask.CommonParameters = CommonParamsToSave;

            DialogResult = true;
        }

        #endregion Private Methods
    }
}