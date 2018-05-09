using EngineLayer;
using Proteomics;
using System;
using System.Collections.Generic;
using System.IO;
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
using System.Windows.Shapes;

namespace MetaMorpheusGUI
{
    /// <summary>
    /// Interaction logic for AddCustomModWindow.xaml
    /// </summary>
    public partial class AddCustomModWindow : Window
    {
        public AddCustomModWindow()
        {
            InitializeComponent();
        }

        public void SaveCustomMod(object sender, RoutedEventArgs e)
        {
            string modsDirectory = System.IO.Path.Combine(GlobalVariables.DataDir, @"Mods");
            var modsFiles = Directory.GetFiles(modsDirectory);
            string customModsPath = System.IO.Path.Combine(modsDirectory, @"CustomMods.txt");
            List<string> customModsText = new List<string>();

            if(!File.Exists(customModsPath))
            {
                customModsText.Add("Custom Modifications");
            }
            else
            {
                customModsText = File.ReadAllLines(customModsPath).ToList();
            }

            string myModName = modNameTextBox.Text;
            string myAminoAcidMotifText = aminoAcidMotifTextBox.Text;
            string chemicalFormulaText = chemicalFormulaTextBox.Text;
            string modMassText = modMassTextBox.Text;
            string neutralLossText = neutralLossTextBox.Text;
            string categoryText = categoryTextBox.Text;

            // write custom mod to mods file
            customModsText.Add("ID   " + myModName);
            customModsText.Add("TG   " + myAminoAcidMotifText);
            customModsText.Add("PP   " + "Anywhere.");
            customModsText.Add("NL   " + neutralLossText);
            customModsText.Add("CF   " + chemicalFormulaText);
            customModsText.Add("MT   " + categoryText);
            customModsText.Add("MM   " + modMassText);
            customModsText.Add(@"//");

            // write/read temp file to make sure the mod is readable, then delete it
            string tempPath = System.IO.Path.Combine(modsDirectory, @"temp.txt");
            try
            {
                File.WriteAllLines(tempPath, customModsText);
                GlobalVariables.AddMods(UsefulProteomicsDatabases.PtmListLoader.ReadModsFromFile(tempPath));
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

            DialogResult = true;
        }

        public void CancelCustomMod(object sender, RoutedEventArgs e)
        {
            DialogResult = false;
        }
    }
}
