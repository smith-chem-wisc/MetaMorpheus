using Proteomics.AminoAcidPolymer;
using System;
using System.Collections.Generic;
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
            if (AminoAcidTextBox.Text.Length == 0)
            {
                MessageBox.Show("Please specify the character that represents a synthetic amino acid in the database");
                return;
            }
            char aminoAcidLetter = AminoAcidTextBox.Text.First();

            //check if it already exists
            if(Residue.TryGetResidue(aminoAcidLetter, out Residue residue))
            {
                MessageBox.Show("The amino acid " + aminoAcidLetter + " already exists." +
                    "\nMass: " + residue.MonoisotopicMass +
                    "\nChemical Formula: " + residue.ThisChemicalFormula +
                    "\nYou may overwrite this amino acid by deleting the current entry. " +
                    "\nThis can be done in MetaMorpheus by navigating to 'Data' in the top-left corner, " +
                    "\nselecting 'Open folder with mods/data files' from the drop down menu," +
                    "\nopening the folder 'CustomAminoAcids', and opening the file 'CustomAminoAcids.txt");
            }
        }

        private void CancelCustomAminoAcid_Click(object sender, RoutedEventArgs e)
        {
            DialogResult = false;
        }
    }
}
