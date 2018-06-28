using Chemistry;
using EngineLayer;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Windows;

namespace MetaMorpheusGUI
{
    /// <summary>
    /// Interaction logic for AddCustomModWindow.xaml
    /// </summary>
    public partial class AddCustomModWindow : Window
    {
        public static Dictionary<string, string> termTypeToParsableTermType;

        public AddCustomModWindow()
        {
            InitializeComponent();

            if (termTypeToParsableTermType == null)
            {
                termTypeToParsableTermType = new Dictionary<string, string>();
                termTypeToParsableTermType.Add("Any", "Anywhere.");
                termTypeToParsableTermType.Add("Peptide N-terminus", "Peptide N-terminal.");
                termTypeToParsableTermType.Add("Peptide C-terminus", "Peptide C-terminal.");
                termTypeToParsableTermType.Add("Protein N-terminus", "N-terminal.");
                termTypeToParsableTermType.Add("Protein C-terminus", "C-terminal.");
            }

            foreach (var kvp in termTypeToParsableTermType)
            {
                terminusTypeComboBox.Items.Add(kvp.Key);
            }

            terminusTypeComboBox.SelectedIndex = 0;
        }

        public void SaveCustomMod(object sender, RoutedEventArgs e)
        {
            string modsDirectory = Path.Combine(GlobalVariables.DataDir, @"Mods");
            var modsFiles = Directory.GetFiles(modsDirectory);
            string customModsPath = Path.Combine(modsDirectory, @"CustomMods.txt");
            List<string> customModsText = new List<string>();

            if (!File.Exists(customModsPath))
            {
                customModsText.Add("Custom Modifications");
            }
            else
            {
                customModsText = File.ReadAllLines(customModsPath).ToList();
            }

            string myModName = modNameTextBox.Text;
            string myAminoAcidMotifText = aminoAcidMotifTextBox.Text;
            string chemicalFormulaText = chemicalFormulaTextBox.Text;
            string modMassText = modMassTextBox.Text;
            string neutralLossText = neutralLossTextBox.Text;
            string categoryText = categoryTextBox.Text;
            string termType = termTypeToParsableTermType[terminusTypeComboBox.Text];

            if (ErrorsDetected(myModName, myAminoAcidMotifText, termType, modMassText, chemicalFormulaText, neutralLossText, categoryText))
            {
                return;
            }

            // write custom mod to mods file
            customModsText.Add("ID   " + myModName);
            customModsText.Add("TG   " + myAminoAcidMotifText);
            customModsText.Add("PP   " + termType);
            customModsText.Add("MT   " + categoryText);

            if (!string.IsNullOrEmpty(neutralLossText))
            {
                customModsText.Add("NL   " + neutralLossText);
            }
            if (!string.IsNullOrEmpty(chemicalFormulaText))
            {
                customModsText.Add("CF   " + chemicalFormulaText);
            }
            if (!string.IsNullOrEmpty(modMassText))
            {
                customModsText.Add("MM   " + modMassText);
            }

            customModsText.Add(@"//");

            // write/read temp file to make sure the mod is readable, then delete it
            string tempPath = Path.Combine(modsDirectory, @"temp.txt");
            try
            {
                File.WriteAllLines(tempPath, customModsText);
                GlobalVariables.AddMods(UsefulProteomicsDatabases.PtmListLoader.ReadModsFromFile(tempPath));
                File.Delete(tempPath);
            }
            catch (Exception ex)
            {
                MessageBox.Show("Problem parsing custom mod: " + ex.Message, "Error", MessageBoxButton.OK, MessageBoxImage.Hand);
                File.Delete(tempPath);
                return;
            }

            // delete old custom mods file, write new one
            try
            {
                File.Delete(customModsPath);
                File.WriteAllLines(customModsPath, customModsText);
            }
            catch (Exception ex)
            {
                MessageBox.Show("Problem saving custom mod to file: " + ex.Message, "Error", MessageBoxButton.OK, MessageBoxImage.Hand);
                return;
            }

            DialogResult = true;
        }

        public void CancelCustomMod(object sender, RoutedEventArgs e)
        {
            DialogResult = false;
        }

        private bool ErrorsDetected(string myModName, string aminoAcidOrMotif, string termType, string mass, string chemFormula, string neutralLoss, string category)
        {
            // parse input
            if (string.IsNullOrEmpty(myModName) || string.IsNullOrEmpty(category))
            {
                MessageBox.Show("Both the mod name and category need to be specified", "Error", MessageBoxButton.OK, MessageBoxImage.Hand);
                return true;
            }
            if (string.IsNullOrEmpty(mass) && string.IsNullOrEmpty(chemFormula))
            {
                MessageBox.Show("Either the mass or chemical formula needs to be specified", "Error", MessageBoxButton.OK, MessageBoxImage.Hand);
                return true;
            }
            if (!string.IsNullOrEmpty(mass) && !double.TryParse(mass, out double dmass))
            {
                MessageBox.Show("Could not parse modification mass", "Error", MessageBoxButton.OK, MessageBoxImage.Hand);
                return true;
            }
            if (!string.IsNullOrEmpty(neutralLoss) && !double.TryParse(neutralLoss, out double dNeutralLoss))
            {
                MessageBox.Show("Could not parse neutral loss mass", "Error", MessageBoxButton.OK, MessageBoxImage.Hand);
                return true;
            }
            try
            {
                ChemicalFormula cf = ChemicalFormula.ParseFormula(chemFormula);
            }
            catch
            {
                MessageBox.Show("Could not parse chemical formula", "Error", MessageBoxButton.OK, MessageBoxImage.Hand);
                return true;
            }

            return false;
        }
    }
}
