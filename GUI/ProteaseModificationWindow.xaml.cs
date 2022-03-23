using Proteomics;
using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Navigation;
using System.Windows.Shapes;
using Proteomics.ProteolyticDigestion;
using System.IO;
using System.Globalization;
using MzLibUtil;
using System.ComponentModel;
using MassSpectrometry;
using Chemistry;
using EngineLayer;

namespace GUI
{
    /// <summary>
    /// Interaction logic for ProteaseModificationWindow.xaml
    /// Allows users to make custom proteases for digestion that have associated modifications
    /// </summary>
    public partial class ProteaseModificationWindow : Window
    {
        public bool proteaseModAdded = false;
        public static Dictionary<string, string> locationRestrictions;
        public string modName = "";
        public ProteaseModificationWindow()
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

            locationRestrictionComboBox.SelectedItem = "Anywhere";
        }
                
        //Save all the user provided information in the proteases file for future use
        private void SaveCustomProteaseMod_Click(object sender, RoutedEventArgs e)
        {
            string proteaseModDirectory = System.IO.Path.Combine(GlobalVariables.DataDir, @"Mods");
            string proteaseModFilePath = System.IO.Path.Combine(proteaseModDirectory, @"ProteaseMods.txt");
            List<string> customModsText = new List<string>();
            customModsText = File.ReadAllLines(proteaseModFilePath).ToList();
            //all of the protease properties that the user provided
            string idText = proteaseModIDTextBox.Text;
            string motifText = proteaseModificationMotifTextBox.Text;
            string chemicalFormulaText = chemicalFormulaTextBox.Text;
            string modMassText = monoisotopicMassTextBox.Text;
            string neutralLossText = null;
            string diagnosticIonText = null;
            string modificationTypeText = "Protease";
            string locationRestriction = locationRestrictions[locationRestrictionComboBox.Text];
            DissociationType disType = DissociationType.HCD;

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

            modName = modification.IdWithMotif;
            // write/read temp file to make sure the mod is readable, then delete it
            string tempPath = System.IO.Path.Combine(proteaseModDirectory, @"temp.txt");
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
                File.Delete(proteaseModFilePath);
                File.WriteAllLines(proteaseModFilePath, customModsText);
            }
            catch (Exception ex)
            {
                MessageBox.Show("Problem saving custom mod to file: " + ex.Message, "Error", MessageBoxButton.OK, MessageBoxImage.Hand);
                return;
            }

            GlobalVariables.AddMods(new List<Modification> { modification }, false);
            GlobalVariables.ProteaseMods = UsefulProteomicsDatabases.PtmListLoader.ReadModsFromFile(proteaseModFilePath, out var errorList).ToList();
            DialogResult = true;
            proteaseModAdded = true;
        }

        private void ClearCustomProteaseMod_Click(object sender, RoutedEventArgs e)
        {
            proteaseModIDTextBox.Clear();
            proteaseModificationMotifTextBox.Clear();
            chemicalFormulaTextBox.Clear();
            locationRestrictions.Clear();            
        }

    }
}
