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
    /// Interaction logic for CalibrateTaskWindow.xaml
    /// </summary>
    public partial class CalibrateTaskWindow : Window
    {
        #region Private Fields

        private readonly ObservableCollection<ModTypeForTreeView> fixedModTypeForTreeViewObservableCollection = new ObservableCollection<ModTypeForTreeView>();
        private readonly ObservableCollection<ModTypeForTreeView> variableModTypeForTreeViewObservableCollection = new ObservableCollection<ModTypeForTreeView>();
        private readonly ObservableCollection<ModTypeForLoc> localizeModTypeForTreeViewObservableCollection = new ObservableCollection<ModTypeForLoc>();

        #endregion Private Fields

        #region Public Constructors

        public CalibrateTaskWindow()
        {
            InitializeComponent();
            PopulateChoices();

            TheTask = new CalibrationTask();
            UpdateFieldsFromTask(TheTask);

            this.saveButton.Content = "Add the Calibration Task";
        }

        public CalibrateTaskWindow(CalibrationTask myCalibrateTask)
        {
            InitializeComponent();
            PopulateChoices();

            TheTask = myCalibrateTask;
            UpdateFieldsFromTask(TheTask);
        }

        #endregion Public Constructors

        #region Internal Properties

        internal CalibrationTask TheTask { get; private set; }

        #endregion Internal Properties

        #region Private Methods

        private void UpdateFieldsFromTask(CalibrationTask task)
        {
            missedCleavagesTextBox.Text = task.CommonParameters.DigestionParams.MaxMissedCleavages == int.MaxValue ? "" : task.CommonParameters.DigestionParams.MaxMissedCleavages.ToString(CultureInfo.InvariantCulture);
            txtMinPeptideLength.Text = task.CommonParameters.DigestionParams.MinPeptideLength.ToString(CultureInfo.InvariantCulture);
            txtMaxPeptideLength.Text = task.CommonParameters.DigestionParams.MaxPeptideLength == int.MaxValue ? "" : task.CommonParameters.DigestionParams.MaxPeptideLength.ToString(CultureInfo.InvariantCulture);
            proteaseComboBox.SelectedItem = task.CommonParameters.DigestionParams.Protease;
            maxModificationIsoformsTextBox.Text = task.CommonParameters.DigestionParams.MaxModificationIsoforms.ToString(CultureInfo.InvariantCulture);
            initiatorMethionineBehaviorComboBox.SelectedIndex = (int)task.CommonParameters.DigestionParams.InitiatorMethionineBehavior;

            productMassToleranceTextBox.Text = task.CommonParameters.ProductMassTolerance.Value.ToString(CultureInfo.InvariantCulture);
            productMassToleranceComboBox.SelectedIndex = task.CommonParameters.ProductMassTolerance is AbsoluteTolerance ? 0 : 1;
            precursorMassToleranceTextBox.Text = task.CommonParameters.PrecursorMassTolerance.Value.ToString(CultureInfo.InvariantCulture);
            precursorMassToleranceComboBox.SelectedIndex = task.CommonParameters.PrecursorMassTolerance is AbsoluteTolerance ? 0 : 1;

            bCheckBox.IsChecked = task.CommonParameters.BIons;
            yCheckBox.IsChecked = task.CommonParameters.YIons;
            zdotCheckBox.IsChecked = task.CommonParameters.ZdotIons;
            cCheckBox.IsChecked = task.CommonParameters.CIons;

            //writeIntermediateFilesCheckBox.IsChecked = task.CalibrationParameters.WriteIntermediateFiles;

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
        }

        private void CancelButton_Click(object sender, RoutedEventArgs e)
        {
            DialogResult = false;
        }

        private void SaveButton_Click(object sender, RoutedEventArgs e)
        {
            #region Check Task Validity

            if (missedCleavagesTextBox.Text.Length == 0)
            {
                missedCleavagesTextBox.Text = int.MaxValue.ToString();
            }
            if (txtMinPeptideLength.Text.Length == 0)
            {
                MessageBox.Show("The minimum peptide length must be a positive integer");
                return;
            }
            if (txtMaxPeptideLength.Text.Length == 0)
            {
                txtMaxPeptideLength.Text = int.MaxValue.ToString();
            }
            if (!double.TryParse(productMassToleranceTextBox.Text, out double pmt) || pmt <= 0)
            {
                MessageBox.Show("The product mass tolerance contains unrecognized characters. \n You entered " + '"' + productMassToleranceTextBox.Text + '"' + "\n Please enter a positive number.");
                return;
            }
            if (!double.TryParse(precursorMassToleranceTextBox.Text, out double premt) || premt <= 0)
            {
                MessageBox.Show("The precursor mass tolerance contains unrecognized characters. \n You entered " + '"' + precursorMassToleranceTextBox.Text + '"' + "\n Please enter a positive number.");
                return;
            }
            if (!double.TryParse(minScoreAllowed.Text, out double msa) || msa < 1)
            {
                MessageBox.Show("The minimum score allowed contains unrecognized characters. \n You entered " + '"' + minScoreAllowed.Text + '"' + "\n Please enter a positive, non-zero number.");
                return;
            }
            if (!int.TryParse(maxModificationIsoformsTextBox.Text, out int mmi) || mmi < 1)
            {
                MessageBox.Show("The maximum number of modification isoforms contains unrecognized characters. \n You entered " + '"' + maxModificationIsoformsTextBox.Text + '"' + "\n Please enter a positive, non-zero number.");
                return;
            }

            #endregion Check Task Validity

            Protease protease = (Protease)proteaseComboBox.SelectedItem;
            int MaxMissedCleavages = int.Parse(missedCleavagesTextBox.Text, CultureInfo.InvariantCulture);
            int MinPeptideLength = int.Parse(txtMinPeptideLength.Text, NumberStyles.Any, CultureInfo.InvariantCulture);
            int MaxPeptideLength = int.Parse(txtMaxPeptideLength.Text, NumberStyles.Any, CultureInfo.InvariantCulture);
            int MaxModificationIsoforms = int.Parse(maxModificationIsoformsTextBox.Text, CultureInfo.InvariantCulture);
            DigestionParams digestionParamsToSave = new DigestionParams(protease: protease.Name, MaxMissedCleavages: MaxMissedCleavages, MinPeptideLength: MinPeptideLength, MaxPeptideLength: MaxPeptideLength, MaxModificationIsoforms: MaxModificationIsoforms);

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
            CommonParameters CommonParamsToSave = new CommonParameters(
                DigestionParams: digestionParamsToSave,
                BIons: bCheckBox.IsChecked.Value,
                YIons: yCheckBox.IsChecked.Value,
                CIons: cCheckBox.IsChecked.Value,
                ZdotIons: zdotCheckBox.IsChecked.Value,
                ScoreCutoff: double.Parse(minScoreAllowed.Text, CultureInfo.InvariantCulture),
                ListOfModsFixed: listOfModsFixed,
                ListOfModsVariable: listOfModsVariable,
                ProductMassTolerance: ProductMassTolerance,
                PrecursorMassTolerance: PrecursorMassTolerance);

            if (OutputFileNameTextBox.Text != "")
            {
                CommonParamsToSave.TaskDescriptor = OutputFileNameTextBox.Text;
            }
            else
            {
                CommonParamsToSave.TaskDescriptor = "CalibrateTask";
            }

            TheTask.CommonParameters = CommonParamsToSave;

            DialogResult = true;
        }

        private void Row_DoubleClick(object sender, MouseButtonEventArgs e)
        {
            var ye = sender as DataGridCell;
            if (ye.Content is TextBlock hm && !string.IsNullOrEmpty(hm.Text))
            {
                System.Diagnostics.Process.Start(hm.Text);
            }
        }

        private void PreviewIfInt(object sender, TextCompositionEventArgs e)
        {
            e.Handled = !TextBoxIntAllowed(e.Text);
        }

        private static bool TextBoxIntAllowed(String Text2)
        {
            return Array.TrueForAll(Text2.ToCharArray(),
                delegate (Char c) { return Char.IsDigit(c) || Char.IsControl(c); });
        }

        #endregion Private Methods
    }
}