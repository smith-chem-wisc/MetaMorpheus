using Chemistry;
using EngineLayer;
using GuiFunctions.Models;
using GuiFunctions.Util;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Windows;

namespace MetaMorpheusGUI
{
    /// <summary>
    /// Interaction logic for CustomMIonLossWindow.xaml
    /// </summary>
    public partial class CustomMIonLossWindow : Window
    {
        public List<CustomMIonLoss> CreatedLosses { get; private set; }

        public CustomMIonLossWindow()
        {
            InitializeComponent();
            
            // Set default selection based on current mode
            if (GuiFunctions.GuiGlobalParamsViewModel.Instance.IsRnaMode)
            {
                rnaModeCheckBox.IsChecked = true;
            }
            else
            {
                proteinModeCheckBox.IsChecked = true;
            }
        }

        private void SaveCustomLoss_Click(object sender, RoutedEventArgs e)
        {
            // Validate inputs
            string name = nameTextBox.Text?.Trim();
            string annotation = annotationTextBox.Text?.Trim();
            string formulaText = chemicalFormulaTextBox.Text?.Trim();

            if (string.IsNullOrEmpty(name))
            {
                MessageBox.Show("Please enter a name for the M-Ion loss.", "Validation Error", MessageBoxButton.OK, MessageBoxImage.Warning);
                return;
            }

            if (string.IsNullOrEmpty(annotation))
            {
                MessageBox.Show("Please enter an annotation for the M-Ion loss.", "Validation Error", MessageBoxButton.OK, MessageBoxImage.Warning);
                return;
            }

            if (string.IsNullOrEmpty(formulaText))
            {
                MessageBox.Show("Please enter a chemical formula.", "Validation Error", MessageBoxButton.OK, MessageBoxImage.Warning);
                return;
            }

            // Check that at least one mode is selected
            if (proteinModeCheckBox.IsChecked != true && rnaModeCheckBox.IsChecked != true)
            {
                MessageBox.Show("Please select at least one applicable mode.", "Validation Error", MessageBoxButton.OK, MessageBoxImage.Warning);
                return;
            }

            // Parse chemical formula
            ChemicalFormula formula;
            try
            {
                formula = ChemicalFormula.ParseFormula(formulaText);
            }
            catch (Exception ex)
            {
                MessageBox.Show($"Invalid chemical formula: {ex.Message}", "Validation Error", MessageBoxButton.OK, MessageBoxImage.Error);
                return;
            }

            // Create list of analyte types based on selections
            var analyteTypes = new List<AnalyteType>();
            if (proteinModeCheckBox.IsChecked == true)
            {
                analyteTypes.Add(AnalyteType.Peptide);
            }
            if (rnaModeCheckBox.IsChecked == true)
            {
                analyteTypes.Add(AnalyteType.Oligo);
            }

            // Create custom M-Ion losses (one per selected mode)
            CreatedLosses = new List<CustomMIonLoss>();
            try
            {
                foreach (var analyteType in analyteTypes)
                {
                    var customLoss = new CustomMIonLoss(name, annotation, formula, analyteType);
                    CustomMIonLossManager.AddCustomMIonLoss(customLoss);
                    CreatedLosses.Add(customLoss);
                }
                
                string modesMessage = analyteTypes.Count > 1 
                    ? $"Protein and RNA modes" 
                    : analyteTypes[0] == AnalyteType.Peptide ? "Protein mode" : "RNA mode";
                
                MessageBox.Show($"Successfully added custom M-Ion loss '{name}' for {modesMessage}.", 
                    "Success", MessageBoxButton.OK, MessageBoxImage.Information);
                DialogResult = true;
                Close();
            }
            catch (InvalidOperationException ex)
            {
                MessageBox.Show(ex.Message, "Duplicate Loss", MessageBoxButton.OK, MessageBoxImage.Warning);
            }
            catch (Exception ex)
            {
                MessageBox.Show($"Error saving custom M-Ion loss: {ex.Message}", "Error", MessageBoxButton.OK, MessageBoxImage.Error);
            }
        }

        private void CancelCustomLoss_Click(object sender, RoutedEventArgs e)
        {
            DialogResult = false;
            Close();
        }
    }
}
