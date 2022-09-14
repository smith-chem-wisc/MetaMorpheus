using System.Windows;
using System.Windows.Input;
using System.IO;
using EngineLayer;
using System.Collections.Generic;
using System.Linq;
using System;
using System.Globalization;

namespace MetaMorpheusGUI
{
    /// <summary>
    /// Interaction logic for CustomCrosslinkerWindow.xaml
    /// </summary>
    public partial class CustomCrosslinkerWindow : Window
    {
        public CustomCrosslinkerWindow()
        {
            InitializeComponent();
            CommonChoice();
        }

        private void CommonChoice()
        {
            txtUdXLkerAminoAcids.Text = "K";
            txtUdXLkerAminoAcids2.Text = "K";
        }

        private void SaveButton_Click(object sender, RoutedEventArgs e)
        {
            string crosslinkerPath = Path.Combine(GlobalVariables.DataDir, @"Data");
            string customCrosslinkPath = Path.Combine(crosslinkerPath, @"CustomCrosslinkers.tsv");

            List<string> customCrosslinkerText = new List<string>();

            if (File.Exists(customCrosslinkPath))
            {
                customCrosslinkerText = File.ReadAllLines(customCrosslinkPath).ToList();
            }
            else
            {
                customCrosslinkerText.Add("Name\tCrosslinkAminoAcid\tCrosslinkerAminoAcid2\tCleavable\tDissociationType\tCrosslinkerTotalMass\tCrosslinkerShortMass\tCrosslinkerLongMass\tQuenchMassH2O\tQuenchMassNH2\tQuenchMassTris");
            }

            string name = txtUdXLKerName.Text;
            string aminoAcid1 = txtUdXLkerAminoAcids.Text;
            string aminoAcid2 = txtUdXLkerAminoAcids2.Text;
            bool isCleavable = ckbUdXLkerCleavable.IsChecked.Value;
            string dissociationTypes = txtDissociationType.Text.ToUpper();
            double mass =  double.Parse(txtUdXLkerTotalMs.Text == "" ? "0" : txtUdXLkerTotalMs.Text, CultureInfo.InvariantCulture);
            double shortMass = double.Parse(txtUdXLkerShortMass.Text == "" ? "0" : txtUdXLkerShortMass.Text, CultureInfo.InvariantCulture);
            double longMass = double.Parse(txtUdXLkerLongMass.Text == "" ? "0" : txtUdXLkerLongMass.Text, CultureInfo.InvariantCulture);
            double H2OQuenchMass = txtH2OQuenchMass.Text == "" ? (mass + 18.01052) : double.Parse(txtH2OQuenchMass.Text, CultureInfo.InvariantCulture);
            double NH2QuenchMass = txtNH2QuenchMass.Text == "" ? (mass + 17.02652) : double.Parse(txtNH2QuenchMass.Text, CultureInfo.InvariantCulture);
            double TrisQuenchMass = txtTrisQuenchMass.Text == "" ? (mass + 121.07392): double.Parse(txtTrisQuenchMass.Text, CultureInfo.InvariantCulture);

            if (GlobalVariables.Crosslinkers.Any(p => p.CrosslinkerName.Contains(name)))
            {
                MessageBox.Show("A crosslinker already exists with the name: " + name, "Error", MessageBoxButton.OK, MessageBoxImage.Hand);
                return;
            }

            var newCrosslinker = new Crosslinker(crosslinkerName: name, crosslinkerModSites: aminoAcid1, crosslinkerModSites2: aminoAcid2, totalMass: mass, cleavable: isCleavable, dissociationTypes: dissociationTypes,
                    cleaveMassShort: shortMass, cleaveMassLong: longMass, loopMass: mass, deadendMassH2O: H2OQuenchMass, deadendMassNH2: NH2QuenchMass, deadendMassTris: TrisQuenchMass);

            customCrosslinkerText.Add(newCrosslinker.ToString(true));

            try
            {                                 
                File.Delete(customCrosslinkPath);
                File.WriteAllLines(customCrosslinkPath, customCrosslinkerText);
            }
            catch (Exception ex)
            {
                MessageBox.Show("Problem saving custom crosslinker to file: " + ex.Message, "Error", MessageBoxButton.OK, MessageBoxImage.Hand);
                return;
            }

            GlobalVariables.AddCrosslinkers(new List<Crosslinker> { newCrosslinker });
            DialogResult = true;
        }

        private void CancelButton_Click(object sender, RoutedEventArgs e)
        {
            DialogResult = false;
        }

        private void CheckIfNumber(object sender, TextCompositionEventArgs e)
        {
            e.Handled = !GlobalGuiSettings.CheckIsNumber(e.Text);
        }
    }
}
