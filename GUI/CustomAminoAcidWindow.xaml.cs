using Chemistry;
using EngineLayer;
using Proteomics.AminoAcidPolymer;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Windows;

namespace MetaMorpheusGUI
{
    /// <summary>
    /// Interaction logic for CustomAminoAcidWindow.xaml
    /// </summary>
    public partial class CustomAminoAcidWindow : Window
    {
        public CustomAminoAcidWindow()
        {
            InitializeComponent();
        }

        private void SaveCustomAminoAcid_Click(object sender, RoutedEventArgs e)
        {
            //VALIDATE INPUT
            if (AminoAcidTextBox.Text.Length == 0)
            {
                MessageBox.Show("Please specify the character that represents a synthetic amino acid in the database");
                return;
            }
            char aminoAcidLetter = AminoAcidTextBox.Text.First();

            ChemicalFormula formula;
            try
            {
                formula = ChemicalFormula.ParseFormula(ChemicalFormulaTextBox.Text);
            }
            catch
            { 
                MessageBox.Show("The checmical formula '"+ ChemicalFormulaTextBox.Text +"' could not be parsed. Please try again.");
                return;
            }

            if (GlobalVariables.InvalidAminoAcids.Contains(aminoAcidLetter))
            {
                MessageBox.Show("The amino acid '" + aminoAcidLetter + "' cannot be assigned. " +
                    "\nThis character is used for modification motifs or as result delimiters. " +
                    "\nThe following amino acids are not allowed:" +
                     string.Join(", ", GlobalVariables.InvalidAminoAcids.Select(x => x.ToString())) + ")");
                return;
            }

            //check if the specified amino acid already exists
            if (Residue.TryGetResidue(aminoAcidLetter, out Residue residue))
            {
                MessageBox.Show("The amino acid '" + aminoAcidLetter + "' already exists." +
                    "\nMonoisotopic Mass: " + residue.MonoisotopicMass +
                    "\nChemical Formula: " + residue.ThisChemicalFormula.Formula +

                    "\n\nYou may overwrite this amino acid by manually deleting/modifying the current entry. " +
                    "\nThis can be done in MetaMorpheus by navigating to 'Data' in the top-left corner, " +
                    "selecting 'Open folder with mods/data files' from the drop down menu, " +
                    "opening the folder 'CustomAminoAcids', and opening the file 'CustomAminoAcids.txt." +
                    "\nMetaMorpheus will need to be restarted for these changes to take effect." +
                    
                    "\n\nAmino acids can be reset to their default values by deleting the file 'CustomAminoAcids.txt' and restarting MetaMorpheus.");
                return;
            }

            //Alright, the entry is valid
            //Append the entry to the CustomAminoAcids.txt file
            string aminoAcidDirectory = Path.Combine(GlobalVariables.DataDir, @"CustomAminoAcids");
            string customAminoAcidPath = Path.Combine(aminoAcidDirectory, @"CustomAminoAcids.txt");

            //check that the file exists and create it if it doesn't
            if (!File.Exists(customAminoAcidPath))
            {
                GlobalVariables.WriteAminoAcidsFile();
            }

            //save it in the amino acid file
            List<string> customAminoAcidsText = File.ReadAllLines(customAminoAcidPath).ToList();
            customAminoAcidsText.Add(AminoAcidTextBox.Text + '\t' + aminoAcidLetter + '\t' + formula.MonoisotopicMass.ToString() + '\t' + formula.Formula); //tsv Name, one letter, monoisotopic, chemical formula
            File.WriteAllLines(customAminoAcidPath, customAminoAcidsText);

            //add the mod to the residue dictionary
            Residue.AddNewResiduesToDictionary(new List<Residue> { new Residue(AminoAcidTextBox.Text, aminoAcidLetter, AminoAcidTextBox.Text, formula, ModificationSites.Any) });
            MessageBox.Show("Success! Amino Acid '" + aminoAcidLetter + "' has been added to the dictionary." +
                "\nMonoisotopic Mass: " + formula.MonoisotopicMass.ToString() +
                "\nChemical Formula: " + formula.Formula);
            DialogResult = true;
        }

        private void CancelCustomAminoAcid_Click(object sender, RoutedEventArgs e)
        {
            DialogResult = false;
        }
    }
}
