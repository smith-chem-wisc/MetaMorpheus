using Proteomics.ProteolyticDigestion;
using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.IO;
using EngineLayer;
using Omics.Digestion;
using Proteomics.AminoAcidPolymer;

namespace MetaMorpheusGUI.Views
{
    /// <summary>
    /// Interaction logic for CustomProteaseWindow.xaml
    /// </summary>
    public partial class CustomProteaseWindow : Window
    {
        public CustomProteaseWindow()
        {
            InitializeComponent();

            foreach (CleavageSpecificity c in Enum.GetValues(typeof(CleavageSpecificity)))
            {
                CleavageSpecComboBox.Items.Add(c);
            }
        }

        private void SaveProteaseButton_Click(object sender, RoutedEventArgs e)
        {
            if (ProteaseNameTextBox.Text.Length == 0)
            {
                MessageBox.Show("Please specify the name of the protease");
                return;
            }
            if (CleavageSpecComboBox.Text.Length == 0)
            {
                MessageBox.Show("Please specify the cleavage specificity of the protease");
                return;
            }
            if (DigestionMotifTextBox.Text.Length == 0)
            {
                MessageBox.Show("Please specify the digestion motif(s) of the protease");
                return;
            }
            string proteaseName = ProteaseNameTextBox.Text;
            string psiMsAccessionNum = PsiMsAccessionNumTextBox.Text;
            string psiMsName = PsiNameTextBox.Text;
            string motifString = DigestionMotifTextBox.Text;

            CleavageSpecificity cleavageSpecificity;
            if (Enum.TryParse(CleavageSpecComboBox.Text, out CleavageSpecificity _cleaveageSpec))
            {
                cleavageSpecificity = _cleaveageSpec;
            }
            else
            {
                MessageBox.Show("Please specify the digestion motif(s) of the protease");
                return;
            }

            List<DigestionMotif> digestionMotifs;
            var aminoAcidsInMotifString = motifString.Where(char.IsLetter).ToList();
            foreach(var aa in aminoAcidsInMotifString)
            {
                if (!Residue.ResiduesDictionary.ContainsKey(aa.ToString()))
                {
                    MessageBox.Show("Motif(s) contains invalid amino acids.");
                    return;
                }
            }
            try
            {
                digestionMotifs = DigestionMotif.ParseDigestionMotifsFromString(motifString);
            }
            catch
            {
                MessageBox.Show("The motif(s) could not be parsed.");
                return;
            }

            //create a new Protease and add to ProteaseDictionary
            Protease proteaseToAdd = new Protease(proteaseName, cleavageSpecificity, psiMsAccessionNum, psiMsName, digestionMotifs);

            if (ProteaseDictionary.Dictionary.ContainsKey(proteaseName))
            {
                MessageBox.Show("This protease already exists");
            }
            else
            {
                if (ProteaseDictionary.Dictionary.TryAdd(proteaseName, proteaseToAdd))
                {
                    MessageBox.Show("Success! Protease '" + proteaseName + "' has been added to the dictionary.");
                } 
                else
                {
                    MessageBox.Show("Failed to add Protease '" + proteaseName + "' to the dictionary.");
                }
            }
            
            //customProtease file path
            string ProtDirectory = Path.Combine(GlobalVariables.DataDir, @"ProteolyticDigestion");
            string customProteasePath = Path.Combine(ProtDirectory, @"CustomProtease.tsv");

            //check if the customProtease file exists, create it if not
            if (!File.Exists(customProteasePath))
            {
                File.Create(customProteasePath);
            }

            string lineToAdd = "\n" + proteaseName + "\t" + motifString + "\t\t\t" + CleavageSpecComboBox.Text + "\t" + psiMsAccessionNum + "\t" + psiMsName + "\t\t\t\t\t";
            var lines = new List<string> { lineToAdd };
            File.AppendAllLines(customProteasePath, lines);
            DialogResult = true;
        }

        private void CancelProteaseButton_Click(object sender, RoutedEventArgs e)
        {
            DialogResult = false;
        }

        private void CleavageSpecComboBox_SelectionChanged(object sender, SelectionChangedEventArgs e)
        {

        }
    }
}
