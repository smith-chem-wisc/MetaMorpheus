using Proteomics.ProteolyticDigestion;
using System;
using System.Collections;
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
using System.IO;
using EngineLayer;
using System.Xml;
using Easy.Common.Extensions;
using System.Linq;
using MassSpectrometry;
using Omics.Digestion;

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
            string proteaseName = ProteaseNameTextBox.Text;
            string psiMsAccessionNum = PsiMsAccessionNumTextBox.Text;
            string psiMsName = PsiNameTextBox.Text;
            string motifString = DigestionMotifTextBox.Text;
            Dictionary<string, CleavageSpecificity> CleavageSpecDic = new Dictionary<string, CleavageSpecificity>();
            foreach (CleavageSpecificity c in Enum.GetValues(typeof(CleavageSpecificity))) {
                CleavageSpecDic.Add(c.ToString(), c);
            }
            CleavageSpecificity cleavageSpecificity = CleavageSpecDic[CleavageSpecComboBox.Text];
            List<DigestionMotif> digestionMotifs;

            try
            {
                digestionMotifs = DigestionMotif.ParseDigestionMotifsFromString(motifString);
            }
            catch
            {
                MessageBox.Show("The motif(s) could not be parsed.");
                return;
            }

            //load the default ProteaseDictionary
            string folderPath = Environment.GetFolderPath(Environment.SpecialFolder.ProgramFiles);
            string path = ((!string.IsNullOrWhiteSpace(folderPath) && AppDomain.CurrentDomain.BaseDirectory.Contains(folderPath) && !AppDomain.CurrentDomain.BaseDirectory.Contains("Jenkins")) ? System.IO.Path.Combine(Environment.GetFolderPath(Environment.SpecialFolder.LocalApplicationData), "MetaMorpheus") : AppDomain.CurrentDomain.BaseDirectory);
            string path2 = System.IO.Path.Combine(path, "ProteolyticDigestion", "proteases.tsv");
            ProteaseDictionary.LoadProteaseDictionary(path2);

            //create a new Protease and add to ProteaseDictionary
            Protease proteaseToAdd = new Protease(proteaseName, cleavageSpecificity, psiMsAccessionNum, psiMsName, digestionMotifs);

            if (ProteaseDictionary.Dictionary.ContainsKey(proteaseName))
            {
                MessageBox.Show("This protease already exists");
            }
            else
            {
                ProteaseDictionary.Dictionary.Add(proteaseName, proteaseToAdd);
            }

            if (!ProteaseDictionary.Dictionary.TryAdd(proteaseName, proteaseToAdd))
            {
                MessageBox.Show("Success! Protease '" + proteaseName + "' has been added to the dictionary.");
            }
            
            //customProtease file path
            string ProtDirectory = Path.Combine(GlobalVariables.DataDir, @"ProteolyticDigestion");
            string customProteasePath = Path.Combine(ProtDirectory, @"CustomProtease.tsv");

            //check if the customProtease file exists, create it if not
            if (!File.Exists(customProteasePath))
            {
                File.Create(customProteasePath);
                //copy the original file to custom file
                File.Copy(path2, customProteasePath, true);
            }

            string lineToAdd = "\n" + proteaseName + "\t" + motifString + "\t\t\t" + CleavageSpecComboBox.Text + "\t" + psiMsAccessionNum + "\t" + psiMsName + "\t\t\t\t\t";
            File.AppendAllText(customProteasePath, lineToAdd);
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
