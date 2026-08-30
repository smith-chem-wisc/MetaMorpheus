using System;
using System.Windows.Controls;
using EngineLayer.DatabaseLoading;

namespace MetaMorpheusGUI
{
    /// <summary>
    /// Interaction logic for FastaHeaderParsingControl.xaml
    /// </summary>
    public partial class FastaHeaderParsingControl : UserControl
    {
        public FastaHeaderParsingControl()
        {
            InitializeComponent();

            foreach (FastaHeaderFormat format in Enum.GetValues<FastaHeaderFormat>())
                FastaHeaderFormatComboBox.Items.Add(format);

            SetFromParameters(new FastaHeaderParsingParameters());
        }

        public void SetFromParameters(FastaHeaderParsingParameters parameters)
        {
            parameters ??= new FastaHeaderParsingParameters();

            FastaHeaderFormatComboBox.SelectedItem = parameters.HeaderFormat;
            CustomAccessionRegexTextBox.Text = parameters.CustomAccessionRegex;
            CustomFullNameRegexTextBox.Text = parameters.CustomFullNameRegex;
            CustomNameRegexTextBox.Text = parameters.CustomNameRegex;
            CustomGeneNameRegexTextBox.Text = parameters.CustomGeneNameRegex;
            CustomOrganismRegexTextBox.Text = parameters.CustomOrganismRegex;
            CustomOrganismIdRegexTextBox.Text = parameters.CustomOrganismIdRegex;

            UpdateCustomEnabledState();
        }

        public FastaHeaderParsingParameters ToParameters()
        {
            var format = FastaHeaderFormatComboBox.SelectedItem is FastaHeaderFormat selected
                ? selected
                : FastaHeaderFormat.UniProt;

            return new FastaHeaderParsingParameters(format,
                CustomAccessionRegexTextBox.Text,
                CustomFullNameRegexTextBox.Text,
                CustomNameRegexTextBox.Text,
                CustomGeneNameRegexTextBox.Text,
                CustomOrganismRegexTextBox.Text,
                CustomOrganismIdRegexTextBox.Text);
        }

        private void FastaHeaderFormatComboBox_SelectionChanged(object sender, SelectionChangedEventArgs e)
        {
            UpdateCustomEnabledState();
        }

        private void UpdateCustomEnabledState()
        {
            if (CustomRegexGrid == null || PresetSummaryTextBlock == null)
                return;

            bool isCustom = FastaHeaderFormatComboBox.SelectedItem is FastaHeaderFormat.Custom;
            CustomRegexGrid.IsEnabled = isCustom;
            PresetSummaryTextBlock.Text = isCustom
                ? ""
                : DescribePreset(FastaHeaderFormatComboBox.SelectedItem as FastaHeaderFormat?);
        }

        private static string DescribePreset(FastaHeaderFormat? format) => format switch
        {
            FastaHeaderFormat.UniProt => "sp|P12345|AATM_RABIT ... OS= GN=",
            FastaHeaderFormat.Ensembl => "ENSP00000001 pep:... gene:ENSG00000001",
            FastaHeaderFormat.Gencode => "ENSP0000001.1|ENST...|...|gene symbol|...",
            FastaHeaderFormat.Ncbi => "gi|16128008|ref|NP_414555.1| description [organism]",
            _ => ""
        };
    }
}
