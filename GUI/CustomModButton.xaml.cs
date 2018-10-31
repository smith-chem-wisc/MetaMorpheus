using Chemistry;
using EngineLayer;
using Proteomics;
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
    public partial class CustomModButtonWindow : Window
    {
        public static Dictionary<string, string> termTypeToParsableTermType;
        public static Dictionary<string, MassSpectrometry.DissociationType> dissociationType;

        public CustomModButtonWindow()
        {
            InitializeComponent();

            if (termTypeToParsableTermType == null)
            {
                termTypeToParsableTermType = new Dictionary<string, string>();
                termTypeToParsableTermType.Add("Anywhere", "Anywhere.");
                termTypeToParsableTermType.Add("Peptide N-Terminus", "Peptide N-terminal.");
                termTypeToParsableTermType.Add("Peptide C-Terminus", "Peptide C-terminal.");
                termTypeToParsableTermType.Add("Protein N-Terminus", "N-terminal.");
                termTypeToParsableTermType.Add("Protein C-Terminus", "C-terminal.");
            }

            foreach (var kvp in termTypeToParsableTermType)
            {
                locationRestrictionComboBox.Items.Add(kvp.Key);
            }

            dissociationType = GlobalVariables.AllSupportedDissociationTypes;

            foreach (var type in dissociationType.Values)
            {
                dissociationTypeComboBox.Items.Add(type);
            }

            locationRestrictionComboBox.SelectedIndex = 0;
            dissociationTypeComboBox.SelectedValue = MassSpectrometry.DissociationType.HCD;
        }

        public void SaveCustomMod_Click(object sender, RoutedEventArgs e)
        {
            string modsDirectory = Path.Combine(GlobalVariables.DataDir, @"Mods");
            var modsFiles = Directory.GetFiles(modsDirectory);
            string customModsPath = Path.Combine(modsDirectory, @"UserCustomModifications.txt");
            List<string> customModsText = new List<string>();

            if (!File.Exists(customModsPath))
            {
                customModsText.Add("Custom Modifications");
            }
            else
            {
                customModsText = File.ReadAllLines(customModsPath).ToList();
            }

            string myModName = originalIdTextBox.Text;
            string motifText = motifTextBox.Text;
            string chemicalFormulaText = chemicalFormulaTextBox.Text;
            string modMassText = modMassTextBox.Text;
            string neutralLossText = neutralLossTextBox.Text;
            string diagnosticIonText = diagnosticIonTextBox.Text;
            string modificationTypeText = modificationTypeTextBox.Text;
            string locationRestriction = termTypeToParsableTermType[locationRestrictionComboBox.Text];
            MassSpectrometry.DissociationType disType = dissociationType[dissociationTypeComboBox.Text];

            if (ErrorsDetected(myModName, motifText, modMassText, chemicalFormulaText, neutralLossText, modificationTypeText, diagnosticIonText))
            {
                return;
            }

            // write custom mod to mods file
            customModsText.Add("ID   " + myModName);
            customModsText.Add("TG   " + motifText);
            customModsText.Add("PP   " + locationRestriction);
            customModsText.Add("MT   " + modificationTypeText);

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
                //GlobalVariables.AddMods(UsefulProteomicsDatabases.PtmListLoader.ReadModsFromFile(tempPath));
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

            Dictionary<MassSpectrometry.DissociationType, List<double>> neutralLosses;
            if (!string.IsNullOrEmpty(neutralLossText))
            {
                neutralLosses = new Dictionary<MassSpectrometry.DissociationType, List<double>>()
            {
                { disType, neutralLossText.Split(',').Select(double.Parse).ToList() }
            };
            }
            else
            {
                neutralLosses = null;
            };

            Dictionary<MassSpectrometry.DissociationType, List<double>> diagnosticIons;
            if (!string.IsNullOrEmpty(diagnosticIonText))
            {
                diagnosticIons = new Dictionary<MassSpectrometry.DissociationType, List<double>>()
            {
                { disType, diagnosticIonText.Split(',').Select(double.Parse).ToList() }
            };
            }
            else
            {
                diagnosticIons = null;
            }

            var motif = Proteomics.ModificationMotif.TryGetMotif(motifText, out Proteomics.ModificationMotif finalMotif);
            ChemicalFormula chemicalFormula = null;
            if (!string.IsNullOrEmpty(chemicalFormulaText))
            {
                chemicalFormula = ChemicalFormula.ParseFormula(chemicalFormulaText);
            }
            else
            {
                chemicalFormula = null;
            }
            double? modMass = null;
            if (!string.IsNullOrEmpty(modMassText))
            {
                modMass = double.Parse(modMassText);
            }

            Modification newMod = new Modification(myModName, null, modificationTypeText, null, finalMotif, locationRestriction,
                chemicalFormula, modMass, null, null, null, neutralLosses, diagnosticIons, null);

            GlobalVariables.AddMods(new List<Modification> { newMod });

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
                MessageBox.Show("Both the mod name and type need to be specified", "Error", MessageBoxButton.OK, MessageBoxImage.Hand);
                return true;
            }
            if (modType.Contains(':'))
            {
                MessageBox.Show("Modification Type cannot contain ':'", "Error", MessageBoxButton.OK, MessageBoxImage.Hand);
                return true;
            }
            if (!string.IsNullOrEmpty(motif))
            {
                var list = motif.Where(char.IsUpper);
                if (!(list.Count() == 1))
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
            if (!string.IsNullOrEmpty(mass) && !double.TryParse(mass, out double dmass))
            {
                MessageBox.Show("Could not parse modification mass", "Error", MessageBoxButton.OK, MessageBoxImage.Hand);
                return true;
            }
            try
            {
                if (!string.IsNullOrEmpty(neutralLoss))
                {
                    neutralLoss.Split(',').Select(double.Parse).ToList();
                }
                if (!string.IsNullOrEmpty(diagnosticIon))
                {
                    diagnosticIon.Split(',').Select(double.Parse).ToList();
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