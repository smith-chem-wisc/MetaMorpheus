using Chemistry;
using EngineLayer;
using System;
using System.Globalization;
using System.Linq;
using System.Windows;

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
            try 
            {
                string warning = GlycanDatabase.PersistCustomMonosaccharide(name, codeText, formulaText, massText, ionsText, descriptionText);
                if (warning != null)
                {
                    MessageBox.Show(warning, "Warning", MessageBoxButton.OK, MessageBoxImage.Warning);
                }
                else 
                {
                    MessageBox.Show(
                          $"The monosaccharide \"{name}\" was added to your dictionary. It can now be recognized when you build your glycan database.",
                          "Success", MessageBoxButton.OK, MessageBoxImage.Information);
                }
            }
            catch (Exception ex)
            {
                MessageBox.Show("Error saving custom monosaccharide: " + ex.Message, "Error", MessageBoxButton.OK, MessageBoxImage.Hand);
                return;
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
            else if (Glycan.NameCharDic.ContainsKey(name))
            {
                MessageBox.Show("A monosaccharide already exists with the name: " + name, "Error", MessageBoxButton.OK, MessageBoxImage.Hand);
                return true;
            }
            // Due to our tsv reader, any glycan name "Name" will be discarded on the next launch, so we don't want to allow it to be added in the glycan menu.
            else if (name.Equals("Name", StringComparison.OrdinalIgnoreCase)) 
            {
                MessageBox.Show("'Name' is reserved for the column-header row and cannot be used as a monosaccharide name", "Error", MessageBoxButton.OK, MessageBoxImage.Hand);
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


            double? formulaMassDa = null;
            if (!string.IsNullOrWhiteSpace(formula))
            {
                try { formulaMassDa = ChemicalFormula.ParseFormula(formula).MonoisotopicMass; }
                catch { MessageBox.Show("Could not parse chemical formula", "Error", MessageBoxButton.OK, MessageBoxImage.Hand); return true; }

                if (formulaMassDa <= 0 || formulaMassDa > 20000)
                {
                    MessageBox.Show("The chemical formula's mass must be a positive number below 20000 Da", "Error", MessageBoxButton.OK, MessageBoxImage.Hand);
                    return true;
                }
            }

            double? typedMassDa = null;
            if (!string.IsNullOrWhiteSpace(mass))
            {
                if (!double.TryParse(mass, NumberStyles.Float, CultureInfo.InvariantCulture, out double parsed) || parsed <= 0 || parsed > 20000)
                {
                    MessageBox.Show("Monoisotopic mass must be a positive number below 20000 Da", "Error", MessageBoxButton.OK, MessageBoxImage.Hand);
                    return true;
                }
                typedMassDa = parsed;
            }

            if (formulaMassDa is null && typedMassDa is null)
            {
                MessageBox.Show("Either the monoisotopic mass or chemical formula needs to be specified", "Error", MessageBoxButton.OK, MessageBoxImage.Hand);
                return true;
            }

            if (formulaMassDa.HasValue && typedMassDa.HasValue && Math.Abs(formulaMassDa.Value - typedMassDa.Value) > 0.01)
            {
                var result = MessageBox.Show(
                    $"The chemical formula's mass ({formulaMassDa:F5} Da) does not match the entered monoisotopic mass ({typedMassDa:F5} Da).\n\n" +
                    "The chemical formula's mass will be used and the typed mass will be discarded. Continue?",
                    "Mass mismatch", MessageBoxButton.YesNo, MessageBoxImage.Warning);
                if (result == MessageBoxResult.No) return true;
            }


            if (!string.IsNullOrWhiteSpace(ions))
            {
                double[] parsedIons;
                try
                {
                    parsedIons = ions.Split(',').Select(p => double.Parse(p.Trim(), NumberStyles.Float, CultureInfo.InvariantCulture)).ToArray();
                }
                catch
                {
                    MessageBox.Show("Diagnostic ions must be entered as numbers separated by ','", "Error", MessageBoxButton.OK, MessageBoxImage.Hand);
                    return true;
                }

                if (parsedIons.Any(ionDa => ionDa <= 0 || ionDa > 20000))
                {
                    MessageBox.Show("Diagnostic ions must be positive numbers below 20000 Da", "Error", MessageBoxButton.OK, MessageBoxImage.Hand);
                    return true;
                }
            }

            return false;
        }
    }
}
