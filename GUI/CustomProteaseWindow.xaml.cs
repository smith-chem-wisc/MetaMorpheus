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
using EngineLayer;

namespace GUI
{
    /// <summary>
    /// Interaction logic for CustomProteaseWindow.xaml
    /// Allows users to make custom proteases for digestion
    /// </summary>
    public partial class CustomProteaseWindow : Window
    {
        public bool proteaseAdded = false;
        public string modName = "";
        public CustomProteaseWindow()
        {
            InitializeComponent();
            PopulateListBoxes();            
        }
        //Fill in list boxes with options
        private void PopulateListBoxes()
        {
            cleavageSpecificityListBox.Items.Add("full");
            cleavageSpecificityListBox.Items.Add("semi");
            cleavageTerminusListBox.Items.Add("C");
            cleavageTerminusListBox.Items.Add("N");
        }

        private void OpenProteaseModification_Click(object sender, RoutedEventArgs e)
        {
            var dialog = new ProteaseModificationWindow();
            dialog.ShowDialog();
            if (dialog.proteaseModAdded == true)
            {
                modName = dialog.modName;
                UsefulProteomicsDatabases.PtmListLoader.ReadModsFromFile(System.IO.Path.Combine(GlobalVariables.DataDir, @"Mods", @"ProteaseMods.txt"), out List<(Modification,string)> filteredModificationsWithWarnings);                
            }

        }

        //Save all the user provided information in the proteases file for future use
        private void SaveCustomProtease_Click(object sender, RoutedEventArgs e)
        {
            string proteaseDirectory = System.IO.Path.Combine(GlobalVariables.DataDir, @"ProteolyticDigestion");
            string proteaseFilePath = System.IO.Path.Combine(proteaseDirectory, @"proteases.tsv");
            List<string> proteaseFileText = new List<string>();
            proteaseFileText = File.ReadAllLines(proteaseFilePath).ToList();

            //all of the protease properties that the user provided
            string name = proteaseNameTextBox.Text;
            string allCleavageResidues = sequencesInducingCleavageTextBox.Text;
            string allResiduesStoppingCleavage = sequencesPreventingCleavageBox.Text;
            var cleavageTerminus = (string)cleavageTerminusListBox.SelectedItem;
            var cleavageSpecificity = (string)cleavageSpecificityListBox.SelectedItem;
            string psiAccession = psiAccessionNumber.Text;
            string psiNames = psiName.Text;            
            
            //formatting these properties for writing to the protease file, so they can be read in each time ProteaseGuru is used
            string proteaseInfo = name + "\t" ;

            var singleCleavageSites = new List<string>();
            var singlePreventionSites = new List<string>();

            if (allCleavageResidues != "")
            {
                singleCleavageSites = allCleavageResidues.Split(',').ToList();
            }

            if (allResiduesStoppingCleavage != "")
            {
                singlePreventionSites = allResiduesStoppingCleavage.Split(',').ToList();
            }

            if (cleavageTerminus == "C")
            {
                string cleavageMotif = "";
                var residues = singleCleavageSites.Count();
                var count = 1;
                foreach (var residue in singleCleavageSites)
                {
                    cleavageMotif += residue;
                    if (singlePreventionSites.Count() != 0)
                    {
                        foreach (var prevent in singlePreventionSites)
                        {
                            cleavageMotif += "[" + prevent + "]";
                        }
                    }
                    if (count < residues)
                    {
                        cleavageMotif += "|,";
                        count++;
                    }
                    else if (count == residues)
                    {
                        cleavageMotif += "|";
                    }
                                        
                }
                proteaseInfo += cleavageMotif;
            }
            else if (cleavageTerminus == "N")
            {
                string cleavageMotif = "";
                var residues = singleCleavageSites.Count();
                var count = 1;
                foreach (var residue in singleCleavageSites)
                {
                    cleavageMotif += "|"+ residue;
                    if (singlePreventionSites.Count() != 0)
                    {
                        foreach (var prevent in singlePreventionSites)
                        {
                            cleavageMotif += "[" + prevent + "]";
                        }
                    }
                    if (count < residues)
                    {
                        cleavageMotif += ",";
                        count++;
                    }
                }
                proteaseInfo += cleavageMotif;
            }
            if (modName != "")
            {
                proteaseInfo += "\t" + "\t" + "\t" + cleavageSpecificity + "\t" + psiAccession + "\t" + psiNames + "\t" + "\t" + modName;
            }
            else
            {
                proteaseInfo += "\t" + "\t" + "\t" + cleavageSpecificity + "\t" + psiAccession + "\t" + psiNames + "\t" + "\t";
            }
            proteaseFileText.Add(proteaseInfo);
            File.WriteAllLines(proteaseFilePath, proteaseFileText);
            ProteaseDictionary.Dictionary = ProteaseDictionary.LoadProteaseDictionary(proteaseFilePath, GlobalVariables.ProteaseMods);
            proteaseAdded = true;
        }

        private void ClearCustomProtease_Click(object sender, RoutedEventArgs e)
        {
            proteaseNameTextBox.Clear();
            sequencesInducingCleavageTextBox.Clear();
            sequencesPreventingCleavageBox.Clear();
            cleavageTerminusListBox.SelectedIndex = -1;
            cleavageSpecificityListBox.SelectedIndex = -1;
            psiAccessionNumber.Clear();
            psiName.Clear();
        }

    }
}
