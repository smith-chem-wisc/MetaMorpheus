using Chemistry;
using EngineLayer;
using EngineLayer.GlycoSearch;
using Google.Protobuf.Collections;
using Org.BouncyCastle.Asn1.BC;
using System;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Security.Policy;
using System.Windows;
using static System.Runtime.InteropServices.JavaScript.JSType;

namespace MetaMorpheusGUI
{
    /// <summary>
    /// Interaction logic for CustomMonosaccharideWindow.xaml
    /// </summary>
    public partial class CustomMonosaccharideWindow : Window
    {
        public CustomMonosaccharideWindow()
        {
            InitializeComponent();
        }

        private void SaveCustomMonosaccharide_Click(object sender, RoutedEventArgs e)
        {
            string name = nameTextBox.Text.Trim();
            string codeText = codeTextBox.Text.Trim();
            string formulaText = chemicalFormulaTextBox.Text.Trim();
            string massText = massTextBox.Text.Trim();
            string ionsText = diagnosticIonsTextBox.Text.Trim();
            string descriptionText = descriptionTextBox.Text.Trim();

            if (ErrorsDetected(name, codeText, formulaText, massText, ionsText))
            {
                return;
            }

            char code = codeText[0];
            double massDa;
            if (!string.IsNullOrEmpty(formulaText))
            {
                massDa = ChemicalFormula.ParseFormula(formulaText).MonoisotopicMass;
            }
            else 
            {
                massDa = double.Parse(massText, NumberStyles.Float, CultureInfo.InvariantCulture);
            }
            int massScaled = (int)Math.Round(massDa * 1E5);

            int[] diagnosticIons = null;
            if(!string.IsNullOrEmpty(ionsText))
            {
                diagnosticIons = ionsText.Split(',')
                    .Select(p => (int)Math.Round(double.Parse(p.Trim(), NumberStyles.Float, CultureInfo.InvariantCulture) * 1E5))
                    .ToArray();
            }

            try 
            {
                Glycan.RegisterCustomMonosaccharide(name, code, massScaled, diagnosticIons);
            }
            catch (ArgumentException ex)
            {
                MessageBox.Show(ex.Message, "Error", MessageBoxButton.OK, MessageBoxImage.Hand);
                return;
            }

            //persist to file so the monosaccharide is still recognized on the next launch
            string glycanModsDirectory = Path.Combine(GlobalVariables.DataDir, @"Glycan_Mods");
            string customMonosaccharidePath = Path.Combine(glycanModsDirectory, "MonosaccharidesCustom.tsv");
            string line = string.Join("\t", name, code.ToString(), massDa.ToString(CultureInfo.InvariantCulture), ionsText, descriptionText);
            try
            {
                Directory.CreateDirectory(glycanModsDirectory);
                if (!File.Exists(customMonosaccharidePath))
                {
                    File.WriteAllLines(customMonosaccharidePath, new[] { "Name\tSingleCharCode\tMonoisotopicMass\tDiagnosticIonMasses\tDescription", line });
                }
                else
                {
                    File.AppendAllLines(customMonosaccharidePath, new[] { line });
                }
            }
            catch (Exception ex) 
            {
                MessageBox.Show("The monosaccharide '" + name + "' is available for this session, but could not be saved to file for future sessions: " + ex.Message, "Warning", MessageBoxButton.OK, MessageBoxImage.Warning);
            }
            DialogResult = true;
        }

        private void CancelCustomMonosaccharide_Click(object sender, RoutedEventArgs e)
        {
            DialogResult = false;
        }

        private bool ErrorsDetected(string name, string code, string formula, string mass, string ions)
        {
            if (string.IsNullOrWhiteSpace(name))
            {
                MessageBox.Show("The monosaccharide name needs to be specified", "Error", MessageBoxButton.OK, MessageBoxImage.Hand);
                return true;
            }
            else if(Glycan.NameCharDic.ContainsKey(name))
            {
                MessageBox.Show("A monosaccharide already exists with the name: " + name, "Error", MessageBoxButton.OK, MessageBoxImage.Hand);
                return true;
            }
            if (string.IsNullOrWhiteSpace(code))
            {
                MessageBox.Show("Please enter a code for the monosaccharide.", "Error", MessageBoxButton.OK, MessageBoxImage.Hand);
                return true;
            }
            else if (code.Length != 1)
            {
                MessageBox.Show("Single - Char Code must be exactly one character", "Error", MessageBoxButton.OK, MessageBoxImage.Hand);
                return true;
            }
            else if (!char.IsAsciiLetter(code[0]))
            {
                MessageBox.Show("Single-Char Code must be an ASCII letter", "Error", MessageBoxButton.OK, MessageBoxImage.Hand);
                return true;
            }
            else if (Glycan.CharMassDic.ContainsKey(code[0])) 
            {
                MessageBox.Show("A monosaccharide already exists with the code: " + code, "Error", MessageBoxButton.OK, MessageBoxImage.Hand);
                return true;
            }

            if (string.IsNullOrWhiteSpace(formula) && string.IsNullOrWhiteSpace(mass)) 
            {
                MessageBox.Show("Either the monoisotopic mass or chemical formula needs to be specified", "Error", MessageBoxButton.OK, MessageBoxImage.Hand);
                return true;
            }
            else if (!string.IsNullOrWhiteSpace(formula))
            {
                try
                {
                    ChemicalFormula.ParseFormula(formula);
                }
                catch 
                {
                    MessageBox.Show("Could not parse chemical formula", "Error", MessageBoxButton.OK, MessageBoxImage.Hand);
                    return true;
                } 
            }

            if (!string.IsNullOrEmpty(mass) && !double.TryParse(mass, NumberStyles.Float, CultureInfo.InvariantCulture, out _))
            {
                MessageBox.Show("Could not parse monoisotopic mass", "Error", MessageBoxButton.OK, MessageBoxImage.Hand);
                return true;
            }

            if (!string.IsNullOrWhiteSpace(ions)) 
            {
                try 
                {
                    ions.Split(',').Select(p => double.Parse(p.Trim(), NumberStyles.Float, CultureInfo.InvariantCulture)).ToList();
                }
                catch 
                {
                    MessageBox.Show("Diagnostic ions must be entered as numbers separated by ','", "Error", MessageBoxButton.OK, MessageBoxImage.Hand);
                    return true;
                }
            }
            
            return false;
        }
    }
}
