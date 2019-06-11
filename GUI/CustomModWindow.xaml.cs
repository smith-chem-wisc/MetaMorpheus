using Chemistry;
using EngineLayer;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Windows;
using MassSpectrometry;
using Proteomics;
using System.Globalization;

namespace MetaMorpheusGUI
{
    /// <summary>
    /// Interaction logic for AddCustomModWindow.xaml
    /// </summary>
    public partial class CustomModWindow : Window
    {
        public static Dictionary<string, string> locationRestrictions;

        public CustomModWindow()
        {
            InitializeComponent();

            if (locationRestrictions == null)
            {
                locationRestrictions = new Dictionary<string, string>();
                locationRestrictions.Add("Anywhere", "Anywhere.");
                locationRestrictions.Add("Peptide N-Terminus", "Peptide N-terminal.");
                locationRestrictions.Add("Peptide C-Terminus", "Peptide C-terminal.");
                locationRestrictions.Add("Protein N-Terminus", "N-terminal.");
                locationRestrictions.Add("Protein C-Terminus", "C-terminal.");
            }

            foreach (string locationRestriction in locationRestrictions.Keys)
            {
                locationRestrictionComboBox.Items.Add(locationRestriction);
            }

            foreach (DissociationType type in GlobalVariables.AllSupportedDissociationTypes.Values)
            {
                if (!type.Equals(DissociationType.Custom))
                {
                    dissociationTypeComboBox.Items.Add(type);
                }
            }

            locationRestrictionComboBox.SelectedItem = "Anywhere";
            dissociationTypeComboBox.SelectedItem = DissociationType.HCD;
        }

        public void SaveCustomMod_Click(object sender, RoutedEventArgs e)
        {
            string modsDirectory = Path.Combine(GlobalVariables.DataDir, @"Mods");
            string customModsPath = Path.Combine(modsDirectory, @"CustomModifications.txt");
            List<string> customModsText = new List<string>();

            if (!File.Exists(customModsPath))
            {
                customModsText.Add("Custom Modifications");
            }
            else
            {
                customModsText = File.ReadAllLines(customModsPath).ToList();
            }

            string idText = originalIdTextBox.Text;
            string motifText = motifTextBox.Text;
            string chemicalFormulaText = chemicalFormulaTextBox.Text;
            string modMassText = modMassTextBox.Text;
            string neutralLossText = neutralLossTextBox.Text;
            string diagnosticIonText = diagnosticIonTextBox.Text;
            string modificationTypeText = modificationTypeTextBox.Text;
            string locationRestriction = locationRestrictions[locationRestrictionComboBox.Text];
            DissociationType disType = GlobalVariables.AllSupportedDissociationTypes[dissociationTypeComboBox.Text];

            if (ErrorsDetected(idText, motifText, modMassText, chemicalFormulaText, neutralLossText, modificationTypeText, diagnosticIonText))
            {
                return;
            }

            // create custom mod
            Dictionary<DissociationType, List<double>> neutralLosses = null;
            if (!string.IsNullOrEmpty(neutralLossText))
            {
                neutralLosses = new Dictionary<DissociationType, List<double>>
                {
                    { disType, neutralLossText.Split(',').Select(p => double.Parse(p, CultureInfo.InvariantCulture)).ToList() }
                };
            }

            Dictionary<DissociationType, List<double>> diagnosticIons = null;
            if (!string.IsNullOrEmpty(diagnosticIonText))
            {
                diagnosticIons = new Dictionary<DissociationType, List<double>>()
                {
                    { disType, diagnosticIonText.Split(',').Select(p => double.Parse(p, CultureInfo.InvariantCulture).ToMass(1)).ToList() }
                };
            }

            ModificationMotif.TryGetMotif(motifText, out ModificationMotif finalMotif);

            ChemicalFormula chemicalFormula = null;
            if (!string.IsNullOrEmpty(chemicalFormulaText))
            {
                chemicalFormula = ChemicalFormula.ParseFormula(chemicalFormulaText);
            }

            double? modMass = null;
            if (!string.IsNullOrEmpty(modMassText))
            {
                modMass = double.Parse(modMassText, CultureInfo.InvariantCulture);
            }

            Modification modification = new Modification(
                _originalId: idText,
                _modificationType: modificationTypeText,
                _target: finalMotif,
                _locationRestriction: locationRestriction,
                _chemicalFormula: chemicalFormula,
                _monoisotopicMass: modMass,
                _neutralLosses: neutralLosses,
                _diagnosticIons: diagnosticIons);

            if (GlobalVariables.AllModsKnownDictionary.ContainsKey(modification.IdWithMotif))
            {
                MessageBox.Show("A modification already exists with the name: " + modification.IdWithMotif, "Error", MessageBoxButton.OK, MessageBoxImage.Hand);
                return;
            }

            // write custom mod to mods file

            // write/read temp file to make sure the mod is readable, then delete it
            string tempPath = Path.Combine(modsDirectory, @"temp.txt");
            try
            {
                List<string> temp = new List<string> { modification.ToString(), @"//" };
                File.WriteAllLines(tempPath, temp);
                var parsedMods = UsefulProteomicsDatabases.PtmListLoader.ReadModsFromFile(tempPath, out var errors);

                if (parsedMods.Count() != 1)
                {
                    MessageBox.Show("Problem parsing custom mod: One mod was expected, a different number was generated", "Error", MessageBoxButton.OK, MessageBoxImage.Hand);
                    return;
                }

                if (errors.Any())
                {
                    string concatErrors = string.Join(Environment.NewLine, errors.Select(p => p.Item2));
                    MessageBox.Show("Problem(s) parsing custom mod: " + Environment.NewLine + concatErrors, "Error", MessageBoxButton.OK, MessageBoxImage.Hand);
                    return;
                }

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
                customModsText.Add(modification.ToString());
                customModsText.Add(@"//");
                File.Delete(customModsPath);
                File.WriteAllLines(customModsPath, customModsText);
            }
            catch (Exception ex)
            {
                MessageBox.Show("Problem saving custom mod to file: " + ex.Message, "Error", MessageBoxButton.OK, MessageBoxImage.Hand);
                return;
            }

            GlobalVariables.AddMods(new List<Modification> { modification }, false);

            DialogResult = true;
        }

        private void CancelCustomMod_Click(object sender, RoutedEventArgs e)
        {
            DialogResult = false;
        }

        private bool ErrorsDetected(string myModName, string motif, string mass, string chemFormula, string neutralLoss, string modType, string diagnosticIon)
        {
            // parse input
            if (string.IsNullOrEmpty(myModName) || string.IsNullOrEmpty(modType))
            {
                MessageBox.Show("The mod name and type need to be specified", "Error", MessageBoxButton.OK, MessageBoxImage.Hand);
                return true;
            }

            if (modType.Contains(':'))
            {
                MessageBox.Show("Modification Type cannot contain ':'", "Error", MessageBoxButton.OK, MessageBoxImage.Hand);
                return true;
            }

            if (!string.IsNullOrEmpty(motif))
            {
                if (motif.Count(char.IsUpper) != 1)
                {
                    MessageBox.Show("Motif must contain exactly one uppercase letter", "Error", MessageBoxButton.OK, MessageBoxImage.Hand);
                    return true;
                }
            }
            else
            {
                MessageBox.Show("Motif must be defined", "Error", MessageBoxButton.OK, MessageBoxImage.Hand);
                return true;
            }

            if (string.IsNullOrEmpty(mass) && string.IsNullOrEmpty(chemFormula))
            {
                MessageBox.Show("Either the mass or chemical formula needs to be specified", "Error", MessageBoxButton.OK, MessageBoxImage.Hand);
                return true;
            }

            if (!string.IsNullOrEmpty(chemFormula))
            {
                try
                {
                    ChemicalFormula cf = ChemicalFormula.ParseFormula(chemFormula);
                }
                catch
                {
                    MessageBox.Show("Could not parse chemical formula", "Error", MessageBoxButton.OK, MessageBoxImage.Hand);
                    return true;
                }
            }

            if (!string.IsNullOrEmpty(mass) && !double.TryParse(mass, NumberStyles.Any, CultureInfo.InvariantCulture, out double dmass))
            {
                MessageBox.Show("Could not parse modification mass", "Error", MessageBoxButton.OK, MessageBoxImage.Hand);
                return true;
            }

            try
            {
                if (!string.IsNullOrEmpty(neutralLoss))
                {
                    neutralLoss.Split(',').Select(p => double.Parse(p, CultureInfo.InvariantCulture)).ToList();
                }
                if (!string.IsNullOrEmpty(diagnosticIon))
                {
                    diagnosticIon.Split(',').Select(p => double.Parse(p, CultureInfo.InvariantCulture)).ToList();
                }
            }
            catch
            {
                MessageBox.Show("Neutral losses and diagnostic ions must be entered as numbers separated by ','", "Error", MessageBoxButton.OK, MessageBoxImage.Hand);
                return true;
            }

            return false;
        }
    }
}