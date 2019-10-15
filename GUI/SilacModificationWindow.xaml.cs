using System;
using System.Linq;
using System.Windows;
using System.Windows.Input;
using Proteomics;
using Proteomics.AminoAcidPolymer;
using Chemistry;
using System.Text.RegularExpressions;
using System.Globalization;

namespace MetaMorpheusGUI
{
    /// <summary>
    /// Interaction logic for ExperimentalDesignWindow.xaml
    /// </summary>
    public partial class SilacModificationWindow : Window
    {
        internal SilacInfoForDataGrid SilacLabel;
        private static readonly Regex FormulaRegex = new Regex(@"\s*([A-Z][a-z]*)(?:\{([0-9]+)\})?(-)?([0-9]+)?\s*", RegexOptions.Compiled);
        private static readonly Regex ValidateFormulaRegex = new Regex("^(" + FormulaRegex + ")+$", RegexOptions.Compiled);
        private int CarbonCount;
        private int NitrogenCount;
        private int OxygenCount;
        private int HydrogenCount;
        private int SulfurCount;
        private double AminoAcidMonoisotopicMass; //starts 0
        private char AminoAcid; //starts '\0'
        private bool ModifyingFormula;

        public SilacModificationWindow()
        {
            InitializeComponent();
        }

        private void KeyPressed(object sender, KeyEventArgs e)
        {
            if (e.Key == Key.Return)
            {
                SaveSilacLabelsButton_Click(sender, e);
            }
            else if (e.Key == Key.Escape)
            {
                CancelSilacButton_Click(sender, e);
            }
        }

        private void AminoAcidLookup(object sender, System.Windows.Controls.TextChangedEventArgs e)
        {
            System.Windows.Controls.TextBox textBox = ((System.Windows.Controls.TextBox)sender);
            char letter = textBox.Text.ToUpper().FirstOrDefault();
            if (letter == 0)
            {
                C12.IsEnabled = false;
                C13.IsEnabled = false;
                N14.IsEnabled = false;
                N15.IsEnabled = false;
                O16.IsEnabled = false;
                O18.IsEnabled = false;
                H1.IsEnabled = false;
                H2.IsEnabled = false;
                S32.IsEnabled = false;
                S34.IsEnabled = false;
            }
            else if (textBox.Text.Length != 1) //if too long, don't let them add stuff and don't change anything
            {
                textBox.Text = textBox.Text[0].ToString();
            }
            else if (Residue.TryGetResidue(letter, out Residue residue))
            {
                AminoAcid = letter;
                AminoAcidMonoisotopicMass = residue.MonoisotopicMass;
                ChemicalFormula.Text = residue.ThisChemicalFormula.Formula;
                ParseChemicalFormula();
                MassDifference.Text = "0.000";
                C12.IsEnabled = true;
                C13.IsEnabled = true;
                N14.IsEnabled = true;
                N15.IsEnabled = true;
                O16.IsEnabled = true;
                O18.IsEnabled = true;
                H1.IsEnabled = true;
                H2.IsEnabled = true;
                S32.IsEnabled = true;
                S34.IsEnabled = true;
            }
            else
            {
                MessageBox.Show('"' + letter.ToString() + '"' + " is not a valid amino acid. Please enter a single amino acid letter.");
                txtBoxAminoAcidLookup.Text = "";
            }
        }

        private void AminoAcidCountsModified(object sender, System.Windows.Controls.TextChangedEventArgs e)
        {
            if (!ModifyingFormula) //if we're changing the formula, don't let the textboxes mess with eachother
            {
                ModifyingFormula = true; //lock
                //figure out which element it is
                System.Windows.Controls.TextBox box = (System.Windows.Controls.TextBox)sender;
                switch (box.Name)
                {
                    case "C12":
                        ElementIsotopeTextBoxHelper(CarbonCount, C12, C13);
                        break;
                    case "C13":
                        ElementIsotopeTextBoxHelper(CarbonCount, C13, C12);
                        break;
                    case "N14":
                        ElementIsotopeTextBoxHelper(NitrogenCount, N14, N15);
                        break;
                    case "N15":
                        ElementIsotopeTextBoxHelper(NitrogenCount, N15, N14);
                        break;
                    case "O16":
                        ElementIsotopeTextBoxHelper(OxygenCount, O16, O18);
                        break;
                    case "O18":
                        ElementIsotopeTextBoxHelper(OxygenCount, O18, O16);
                        break;
                    case "H1":
                        ElementIsotopeTextBoxHelper(HydrogenCount, H1, H2);
                        break;
                    case "H2":
                        ElementIsotopeTextBoxHelper(HydrogenCount, H2, H1);
                        break;
                    case "S32":
                        ElementIsotopeTextBoxHelper(SulfurCount, S32, S34);
                        break;
                    case "S34":
                        ElementIsotopeTextBoxHelper(SulfurCount, S34, S32);
                        break;
                    default:
                        MessageBox.Show("Element " + box.Name + " has not been implemented.");
                        break;
                }

                //update the chemical formula
                ChemicalFormula f = new ChemicalFormula();
                if (!C12.Text.Equals(""))
                {
                    f.Add(PeriodicTable.GetElement("C"), Convert.ToInt16(C12.Text));
                }
                if (!C13.Text.Equals(""))
                {
                    f.Add(PeriodicTable.GetElement("C")[13], Convert.ToInt16(C13.Text));
                }
                if (!N14.Text.Equals(""))
                {
                    f.Add(PeriodicTable.GetElement("N"), Convert.ToInt16(N14.Text));
                }
                if (!N15.Text.Equals(""))
                {
                    f.Add(PeriodicTable.GetElement("N")[15], Convert.ToInt16(N15.Text));
                }
                if (!O16.Text.Equals(""))
                {
                    f.Add(PeriodicTable.GetElement("O"), Convert.ToInt16(O16.Text));
                }
                if (!O18.Text.Equals(""))
                {
                    f.Add(PeriodicTable.GetElement("O")[18], Convert.ToInt16(O18.Text));
                }
                if (!H1.Text.Equals(""))
                {
                    f.Add(PeriodicTable.GetElement("H"), Convert.ToInt16(H1.Text));
                }
                if (!H2.Text.Equals(""))
                {
                    f.Add(PeriodicTable.GetElement("H")[2], Convert.ToInt16(H2.Text));
                }
                if (!S32.Text.Equals(""))
                {
                    f.Add(PeriodicTable.GetElement("S"), Convert.ToInt16(S32.Text));
                }
                if (!S34.Text.Equals(""))
                {
                    f.Add(PeriodicTable.GetElement("S")[34], Convert.ToInt16(S34.Text));
                }
                MassDifference.Text = Math.Round(f.MonoisotopicMass - AminoAcidMonoisotopicMass, 3).ToString("F3");
                ChemicalFormula.Text = f.Formula;
                ModifyingFormula = false; //unlock
            }
        }

        private void ElementIsotopeTextBoxHelper(int count, System.Windows.Controls.TextBox modifiedTextBox, System.Windows.Controls.TextBox complementaryTextBox)
        {
            if (modifiedTextBox.Text.Equals("")) //if deleted, set the complement to be the max
            {
                complementaryTextBox.Text = count.ToString();
            }
            else
            {
                int complementaryValue = count - Convert.ToInt16(modifiedTextBox.Text);
                if (complementaryValue > count || complementaryValue < 0) //if an invalid value was provided (too high or too low)
                {
                    //modify the modified textbox
                    modifiedTextBox.Text = (count - Convert.ToInt16(complementaryTextBox.Text)).ToString();
                }
                else
                {
                    //modify the complementary textbox
                    complementaryTextBox.Text = (count - Convert.ToInt16(modifiedTextBox.Text)).ToString();
                }
            }
        }

        private void ParseChemicalFormula()
        {
            ModifyingFormula = true; //lock
            string formula = ChemicalFormula.Text;
            if (!ValidateFormulaRegex.IsMatch(formula))
            {
                MessageBox.Show("Input string for chemical formula was in an incorrect format: " + formula);
            }
            else if (formula.Contains("-"))
            {
                MessageBox.Show("Element numbers cannot be negative.");
            }
            else
            {
                //Clear all textboxes before populating them
                C12.Text = "0";
                C13.Text = "0";
                N14.Text = "0";
                N15.Text = "0";
                O16.Text = "0";
                O18.Text = "0";
                H1.Text = "0";
                H2.Text = "0";
                S32.Text = "0";
                S34.Text = "0";

                //populate the textboxes
                foreach (Match match in FormulaRegex.Matches(formula))
                {
                    string chemsym = match.Groups[1].Value; // Group 1: Chemical Symbol

                    Element element = PeriodicTable.GetElement(chemsym);

                    int quantityOfElement = match.Groups[4].Success ? // Group 4 (optional): Number of Elements
                        int.Parse(match.Groups[4].Value, CultureInfo.InvariantCulture) :
                        1;

                    if (match.Groups[2].Success) // Group 2 (optional): Isotope Mass Number
                    {
                        // Adding heavy isotope!
                        switch (element.AtomicSymbol)
                        {
                            case "C":
                                C13.Text = quantityOfElement.ToString();
                                break;
                            case "N":
                                N15.Text = quantityOfElement.ToString();
                                break;
                            case "O":
                                O18.Text = quantityOfElement.ToString();
                                break;
                            case "H":
                                H2.Text = quantityOfElement.ToString();
                                break;
                            case "S":
                                S34.Text = quantityOfElement.ToString();
                                break;
                            default:
                                MessageBox.Show("Element " + element.AtomicSymbol + " has not been implemented.");
                                break;
                        }
                    }
                    else
                    {
                        // Adding element!
                        switch (element.AtomicSymbol)
                        {
                            case "C":
                                C12.Text = quantityOfElement.ToString();
                                break;
                            case "N":
                                N14.Text = quantityOfElement.ToString();
                                break;
                            case "O":
                                O16.Text = quantityOfElement.ToString();
                                break;
                            case "H":
                                H1.Text = quantityOfElement.ToString();
                                break;
                            case "S":
                                S32.Text = quantityOfElement.ToString();
                                break;
                            default:
                                MessageBox.Show("Element " + element.AtomicSymbol + " has not been implemented.");
                                break;
                        }
                    }
                }
                CarbonCount = Convert.ToInt16(C12.Text) + Convert.ToInt16(C13.Text);
                NitrogenCount = Convert.ToInt16(N14.Text) + Convert.ToInt16(N15.Text);
                OxygenCount = Convert.ToInt16(O16.Text) + Convert.ToInt16(O18.Text);
                HydrogenCount = Convert.ToInt16(H1.Text) + Convert.ToInt16(H2.Text);
                SulfurCount = Convert.ToInt16(S32.Text) + Convert.ToInt16(S34.Text);
            }
            ModifyingFormula = false; //unlock
        }

        private void CheckIfNumber(object sender, TextCompositionEventArgs e)
        {
            e.Handled = !GlobalGuiSettings.CheckIsNumber(e.Text);
        }

        private void SaveSilacLabelsButton_Click(object sender, RoutedEventArgs e)
        {
            if (ValidateClick())
            {
                SilacLabel = SilacLabelStorageForMultiLabeledConditions ?? //use if labels have been added
                    new SilacInfoForDataGrid(new SilacLabel(AminoAcid, '\0', ChemicalFormula.Text, Convert.ToDouble(MassDifference.Text)), DetermineExperimentType());
                DialogResult = true;
            }
        }

        private SilacInfoForDataGrid SilacLabelStorageForMultiLabeledConditions;

        private void AddAdditionalButton_Click(object sender, RoutedEventArgs e)
        {
            if (ValidateClick())
            {
                //save current label
                SilacLabel savedLabel = new SilacLabel(AminoAcid, '\0', ChemicalFormula.Text, Convert.ToDouble(MassDifference.Text));
                SilacLabelStorageForMultiLabeledConditions = new SilacInfoForDataGrid(savedLabel, DetermineExperimentType());

                //open a new window to allow for more mods on this condition
                var dialog = new SilacModificationWindow();
                //update the radiobutton fields
                dialog.MultiplexRadioButton.IsChecked = MultiplexRadioButton.IsChecked;
                //dialog.StartConditionRadioButton.IsChecked = StartConditionRadioButton.IsChecked;
                //dialog.EndConditionRadioButton.IsChecked = EndConditionRadioButton.IsChecked;

                if (dialog.ShowDialog() == true)
                {
                    char originalAminoAcid = SilacLabelStorageForMultiLabeledConditions.SilacLabel[0].OriginalAminoAcid;
                    //As you loop throught this, the process looks something like addAdditional->addAdditional->Save (for 3 labels).
                    //But the saving actually happens as Save->add->add by storing the list of labels in each dialog.SilacLabel
                    //So you need to search against the dialog.silacLabel instead of the storage, which is only temporarily populated at the end.
                    if (dialog.SilacLabel.SilacLabel.Any(x => x.OriginalAminoAcid == originalAminoAcid))
                    {
                        MessageBox.Show("The amino acid '" + originalAminoAcid.ToString() + "' was labeled multiple times for this condition. Please try again.");
                        CancelSilacButton_Click(sender, e);
                    }
                    else
                    {
                        SilacLabelStorageForMultiLabeledConditions.AddAdditionalLabel(dialog.SilacLabel);
                        //update the ExperimentType to the most recent (allows the user to change their mind)
                        SilacLabelStorageForMultiLabeledConditions.LabelType = dialog.SilacLabel.LabelType;
                        SaveSilacLabelsButton_Click(sender, e);
                    }
                }
            }
        }

        private bool ValidateClick()
        {
            if (AminoAcid == '\0')
            {
                MessageBox.Show("Please select an amino acid to be labeled.");
                return false;
            }
            else if (MassDifference.Text.Equals("0.000"))
            {
                MessageBox.Show("No heavy isotopes were specified.\n" +
                    "Please type the number of each heavy isotope in the corresponding boxes.\n" +
                    "Unlabeled " + '"' + "light" + '"' + "amino acids are automatically searched and should not be specified.");
                return false;
            }
            else
            {
                return true;
            }
        }

        private void CancelSilacButton_Click(object sender, RoutedEventArgs e)
        {
            DialogResult = false;
        }

        private ExperimentType DetermineExperimentType()
        {
            return MultiplexRadioButton.IsChecked.Value ? ExperimentType.Multiplex :
             //StartConditionRadioButton.IsChecked.Value ? ExperimentType.Start :
             ExperimentType.End;
        }

        public enum ExperimentType { Multiplex, Start, End }
    }
}