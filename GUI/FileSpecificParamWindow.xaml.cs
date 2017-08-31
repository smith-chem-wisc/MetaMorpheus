using EngineLayer;
using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.Globalization;
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
using TaskLayer;

namespace MetaMorpheusGUI
{
    /// <summary>
    /// Interaction logic for FileSpecificParamWindow.xaml
    /// </summary>
    /// 


    public partial class FileSpecificParamWindow : Window
    {
        internal FileSpecificSettings[] FileSpecificSettingsList { get; private set; }
        internal ObservableCollection<RawDataForDataGrid> SelectedRaw { get; private set; }

        public FileSpecificParamWindow(ObservableCollection<RawDataForDataGrid> selectedRaw)
        {
            InitializeComponent();
            SelectedRaw = selectedRaw;
            FileSpecificSettingsList = new FileSpecificSettings[selectedRaw.Count];

            foreach (Protease protease in GlobalEngineLevelSettings.ProteaseDictionary.Values)
                proteaseComboBox.Items.Add(protease);
            proteaseComboBox.SelectedIndex = 12;

            foreach (string initiatior_methionine_behavior in Enum.GetNames(typeof(InitiatorMethionineBehavior)))
                initiatorMethionineBehaviorComboBox.Items.Add(initiatior_methionine_behavior);
            string path = selectedRaw[0].FilePath;

            productMassToleranceComboBox.Items.Add("Absolute");
            productMassToleranceComboBox.Items.Add("Ppm");

        }

        public void WriteTomls(ObservableCollection<RawDataForDataGrid> a)
        {



        }

        private void PreviewIfInt(object sender, TextCompositionEventArgs e)
        {
            e.Handled = !TextBoxIntAllowed(e.Text);
        }

        private static Boolean TextBoxIntAllowed(String Text2)
        {
            return Array.TrueForAll<Char>(Text2.ToCharArray(),
                delegate (Char c) { return Char.IsDigit(c) || Char.IsControl(c); });
        }

        private void Cancel_Click(object sender, RoutedEventArgs e)
        {
            DialogResult = false;
        }

        private void Save_Click(object sender, RoutedEventArgs e)
        {
            for (int i = 0; i < FileSpecificSettingsList.Count(); i++)
            {
                FileSpecificSettingsList[i] = new FileSpecificSettings();
                FileSpecificSettingsList[i].DeconvolutionIntensityRatio = null;
                FileSpecificSettingsList[i].ConserveMemory = conserveMemoryCheckBox.IsChecked.Value;
                FileSpecificSettingsList[i].DeconvolutionIntensityRatio = TryParseNullable(DeconvolutionIntensityRatioTextBox.Text);
                if (TryParseNullableDouble((DeconvolutionMassToleranceInPpmTextBox.Text)) == null)
                {
                    FileSpecificSettingsList[i].DeconvolutionMassTolerance = null;
                }
                else
                    FileSpecificSettingsList[i].DeconvolutionMassTolerance = new PpmTolerance(double.Parse(DeconvolutionMassToleranceInPpmTextBox.Text, CultureInfo.InvariantCulture));
                FileSpecificSettingsList[i].DeconvolutionMaxAssumedChargeState = TryParseNullable(DeconvolutionMaxAssumedChargeStateTextBox.Text);
                FileSpecificSettingsList[i].DoPrecursorDeconvolution = DoPrecursorDeconvolution.IsChecked.Value;
                FileSpecificSettingsList[i].InitiatorMethionineBehavior = (InitiatorMethionineBehavior)initiatorMethionineBehaviorComboBox.SelectedIndex;
                FileSpecificSettingsList[i].MaxMissedCleavages = TryParseNullable(missedCleavagesTextBox.Text);
                FileSpecificSettingsList[i].MaxModificationIsoforms = TryParseNullable(maxModificationIsoformsTextBox.Text);
                FileSpecificSettingsList[i].MinPeptideLength = int.TryParse(txtMinPeptideLength.Text, NumberStyles.Any, CultureInfo.InvariantCulture, out int temp) ? (int?)temp : null;
                FileSpecificSettingsList[i].MaxPeptideLength = int.TryParse(txtMaxPeptideLength.Text, NumberStyles.Any, CultureInfo.InvariantCulture, out temp) ? (int?)temp : null;
                if (TryParseNullableDouble(productMassToleranceTextBox.Text) != null)
                {
                    if (productMassToleranceComboBox.SelectedIndex == 0)
                        FileSpecificSettingsList[i].ProductMassTolerance = new AbsoluteTolerance(double.Parse(productMassToleranceTextBox.Text, CultureInfo.InvariantCulture));
                    else
                        FileSpecificSettingsList[i].ProductMassTolerance = new PpmTolerance(double.Parse(productMassToleranceTextBox.Text, CultureInfo.InvariantCulture));

                }
                else
                    FileSpecificSettingsList[i].ProductMassTolerance = null;

                FileSpecificSettingsList[i].Protease = (Protease)proteaseComboBox.SelectedItem;
                FileSpecificSettingsList[i].ScoreCutoff = TryParseNullableDouble(minScoreAllowed.Text);
                FileSpecificSettingsList[i].TotalPartitions = TryParseNullable(numberOfDatabaseSearchesTextBox.Text);
                FileSpecificSettingsList[i].UseProvidedPrecursorInfo = UseProvidedPrecursorInfo.IsChecked.Value;

            }

            DialogResult = true;
        
        }

        //helper methods to return null if field is empty
        private int? TryParseNullable(string val)
        {
            return int.TryParse(val, out int outValue) ? (int?)outValue : null;
        }

        private double? TryParseNullableDouble(string val)
        {
            return double.TryParse(val, out double outValue) ? (double?)outValue : null;
        }


    }
}



