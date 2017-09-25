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
using System.Windows.Shapes;
using MetaMorpheusGUI;
using EngineLayer;
using TaskLayer;
using System.IO;
using Nett;
using MzLibUtil;

namespace MetaMorpheusGUI
{
    /// <summary>
    /// Interaction logic for ChangeParametersWindow.xaml
    /// </summary>
    public partial class ChangeParametersWindow : Window
    {
        internal ObservableCollection<RawDataForDataGrid> SelectedRaw { get; private set; }
        internal FileSpecificSettings[] FileSpecificSettingsList { get; private set; }
        Parameter[] paramList;


        public ChangeParametersWindow(ObservableCollection<RawDataForDataGrid> selectedRaw)
        {


            SelectedRaw = selectedRaw;
            FileSpecificSettingsList = new FileSpecificSettings[selectedRaw.Count];
            InitializeComponent();

            paramList = InitializeParameterList();
            ParameterGrid.Items.RemoveAt(0);
            for (int i = 0; i < paramList.Count(); i++)
            {
                ParameterGrid.Items.Add(paramList[i]);
            }
            //ParameterGrid.Items.get
            // ParameterGrid.Items.Add(new Parameter("asdf", "comboBox"));



            FileSpecificSettingsList = new FileSpecificSettings[selectedRaw.Count];
            if (selectedRaw.Count == 1)
            {
                string tomlFileName = System.IO.Path.ChangeExtension(SelectedRaw[0].FilePath, ".toml").ToString();
                Dictionary<string, KeyValuePair<string, Nett.TomlObject>> tomlSettingsList;
                if (File.Exists(tomlFileName))
                {
                    var paramFile = Toml.ReadFile(tomlFileName, MetaMorpheusTask.tomlConfig);
                    tomlSettingsList = paramFile.ToDictionary(p => p.Key);
                    FileSpecificSettings settings = new FileSpecificSettings(tomlSettingsList);
                    //UpdateAndPopulateFields(settings);
                }
                else
                {
                    tomlSettingsList = new Dictionary<string, KeyValuePair<string, TomlObject>>();
                    FileSpecificSettings settings = new FileSpecificSettings(tomlSettingsList);
                    settings.ConserveMemory = null;
                    settings.DoPrecursorDeconvolution = null;
                    settings.UseProvidedPrecursorInfo = null;
                    //UpdateAndPopulateFields(settings);
                }
            }

            //if more than one file selected have to do processing to check equality.
            else
            {
                string[] tomlFileNames = new string[SelectedRaw.Count];
                List<string> nonEqualValueNames = new List<string>();
                List<Dictionary<string, KeyValuePair<string, Nett.TomlObject>>> tomlSettingsListList = new List<Dictionary<string, KeyValuePair<string, Nett.TomlObject>>>();
                
                //Set tomlSettingListList (list of tomlSettings)
                for (int i = 0; i < selectedRaw.Count; i++)
                {
                    
                    tomlFileNames[i] = System.IO.Path.ChangeExtension(SelectedRaw[i].FilePath, ".toml").ToString();
                    if (File.Exists(tomlFileNames[i]))
                    {
                        var paramFile = Toml.ReadFile(tomlFileNames[i], MetaMorpheusTask.tomlConfig);
                        tomlSettingsListList.Add(paramFile.ToDictionary(p => p.Key));
                    }
                    else
                    {
                        tomlSettingsListList[i] = new Dictionary<string, KeyValuePair<string, TomlObject>>();
                    }


                }

                //check every value against every other value; if same, then add to list to set 'different in multiple files' later in Parameter Object
                Console.WriteLine(tomlSettingsListList[0]);
                for (int i = 0; i < tomlSettingsListList.Count; i++)
                {
                    for (int j = 0; j < tomlSettingsListList.Count; j++)
                    {
                        foreach (var key in tomlSettingsListList[i].Keys)
                        {
                            if (tomlSettingsListList[j].ContainsKey(tomlSettingsListList[i][key].Key))
                            {
                                if (!tomlSettingsListList[i][key].Value.Equals(tomlSettingsListList[j][key].Value)) 
                                {
                                    nonEqualValueNames.Add(tomlSettingsListList[i][key].Key);
                                }
                            }
                        }
                    }

                }
                //now set settings so that it reflects if values are different
                FileSpecificSettings settings = new FileSpecificSettings(tomlSettingsListList[0]);
                foreach(string key in nonEqualValueNames)
                {
                    if (key == "ConserveMemory")
                    {
                        paramList[1].Different = true;
                    }
                    if (key == "Protease")
                    {
                        paramList[0].Different = true;
                    }
                    if (key == "ConserveMemory")
                    {
                        paramList[1].Different = true;
                    }
                    if (key == "ConserveMemory")
                    {
                        paramList[1].Different = true;
                    }
                    if (key == "ConserveMemory")
                    {
                        paramList[1].Different = true;
                    }
                    if (key == "ConserveMemory")
                    {
                        paramList[1].Different = true;
                    }
                }

            }

            //DialogResult = true;



            /* private void UpdateAndPopulateFields(FileSpecificSettings settings)
              {
                  CommonParameters commonparams = new CommonParameters();
                  foreach (Protease protease in GlobalEngineLevelSettings.ProteaseDictionary.Values)
                      proteaseComboBox.Items.Add(protease);
                  proteaseComboBox.SelectedIndex = 12;

                  foreach (string initiatior_methionine_behavior in Enum.GetNames(typeof(InitiatorMethionineBehavior)))
                      initiatorMethionineBehaviorComboBox.Items.Add(initiatior_methionine_behavior);

                  productMassToleranceComboBox.Items.Add("Absolute");
                  productMassToleranceComboBox.Items.Add("Ppm");

                  DoPrecursorDeconvolution.IsChecked = settings.DoPrecursorDeconvolution ?? commonparams.DoPrecursorDeconvolution;
                  conserveMemoryCheckBox.IsChecked = settings.ConserveMemory ?? commonparams.ConserveMemory; ;
                  UseProvidedPrecursorInfo.IsChecked = settings.UseProvidedPrecursorInfo ?? commonparams.UseProvidedPrecursorInfo;
                  if (settings.DeconvolutionIntensityRatio != null)
                      DeconvolutionIntensityRatioTextBox.Text = settings.DeconvolutionIntensityRatio.ToString();
                  if (settings.DeconvolutionMaxAssumedChargeState != null)
                      DeconvolutionMaxAssumedChargeStateTextBox.Text = settings.DeconvolutionMaxAssumedChargeState.ToString();
                  if (settings.DeconvolutionMassTolerance != null)
                      DeconvolutionMassToleranceInPpmTextBox.Text = settings.DeconvolutionMassTolerance.ToString();
                  initiatorMethionineBehaviorComboBox.SelectedIndex = 1;
                  if (settings.MaxMissedCleavages != null)
                      missedCleavagesTextBox.Text = settings.MaxMissedCleavages.ToString();
                  txtMinPeptideLength.Text = settings.MinPeptideLength.HasValue ? settings.MinPeptideLength.Value.ToString(CultureInfo.InvariantCulture) : "";
                  txtMaxPeptideLength.Text = settings.MaxPeptideLength.HasValue ? settings.MaxPeptideLength.Value.ToString(CultureInfo.InvariantCulture) : "";
                  if (settings.MaxModificationIsoforms != null)
                      maxModificationIsoformsTextBox.Text = settings.MaxModificationIsoforms.ToString();
                  if (settings.TotalPartitions != null)
                      numberOfDatabaseSearchesTextBox.Text = settings.TotalPartitions.ToString();
                  if (settings.Protease != null)
                      proteaseComboBox.SelectedItem = settings.Protease;
                  if (settings.ScoreCutoff != null)
                      minScoreAllowed.Text = settings.ScoreCutoff.ToString();
                  if (settings.Max_mods_for_peptide != null)
                      MaxModsForPeptides.Text = settings.Max_mods_for_peptide.ToString();
              }*/
        }
        private void Save_Click(object sender, RoutedEventArgs e)
        {
            for (int i = 0; i < FileSpecificSettingsList.Count(); i++)
            {
                FileSpecificSettingsList[i] = new FileSpecificSettings();
                int? index = paramList[0].Value as int?;
                if (index.HasValue)
                {
                    FileSpecificSettingsList[i].Protease = GlobalEngineLevelSettings.ProteaseDictionary.ElementAt(index.Value).Value;    //paramList[0].Value;
                }
                FileSpecificSettingsList[i].ConserveMemory = paramList[1].Value as bool?;
                FileSpecificSettingsList[i].Max_mods_for_peptide = TryParseNullable(paramList[2].Value as string);
                FileSpecificSettingsList[i].DeconvolutionIntensityRatio = paramList[3].Value as double?;
                FileSpecificSettingsList[i].DoPrecursorDeconvolution = paramList[4].Value as bool?;
                FileSpecificSettingsList[i].UseProvidedPrecursorInfo = paramList[5].Value as bool?;
                FileSpecificSettingsList[i].ScoreCutoff = paramList[6].Value as double?;
                FileSpecificSettingsList[i].ProductMassTolerance = paramList[7].Value as Tolerance;
                FileSpecificSettingsList[i].DeconvolutionMaxAssumedChargeState = paramList[8].Value as int?;
                FileSpecificSettingsList[i].TotalPartitions = paramList[9].Value as int?;
                FileSpecificSettingsList[i].MaxModificationIsoforms = paramList[10].Value as int?;
                FileSpecificSettingsList[i].MaxPeptideLength = paramList[11].Value as int?;
                FileSpecificSettingsList[i].MinPeptideLength = paramList[12].Value as int?;
                FileSpecificSettingsList[i].MaxMissedCleavages = paramList[13].Value as int?;
                if (paramList[14].Value != null)
                    FileSpecificSettingsList[i].InitiatorMethionineBehavior = (InitiatorMethionineBehavior)paramList[14].Value;
                FileSpecificSettingsList[i].DeconvolutionMassTolerance = paramList[15].Value as Tolerance;
                FileSpecificSettingsList[i].TrimMsMsPeaks = paramList[16].Value as bool?;
                FileSpecificSettingsList[i].TrimMs1Peaks = paramList[17].Value as bool?;
                FileSpecificSettingsList[i].MinRatio = paramList[18].Value as double?;
                FileSpecificSettingsList[i].TopNpeaks = paramList[19].Value as int?;


            }

            DialogResult = true;

        }

        private void CheckIfEqual()
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

        private Parameter[] InitializeParameterList()
        {



            Parameter[] paramList = new Parameter[20];

            paramList[0] = new Parameter("Protease", "ComboBoxProtease");
            paramList[1] = new Parameter("Conserve Memory", "Bool");
            paramList[2] = new Parameter("Max Mods For Peptides", "TextBox");
            paramList[3] = new Parameter("Deconvolution Intensity Ratio", "TextBox");
            paramList[4] = new Parameter("Do Precursor Deconvolution", "Bool");
            paramList[5] = new Parameter("Use Provided Precursor Info", "Bool");
            paramList[6] = new Parameter("Minimum MM Score Allowed", "TextBox");
            paramList[7] = new Parameter("Product Mass Tolerance in PPM", "TextBox");
            paramList[8] = new Parameter("Deconvolution Max Assumed Charge State", "TextBox");
            paramList[9] = new Parameter("Number of Database Partitions", "TextBox");
            paramList[10] = new Parameter("Max Modification Isoforms", "TextBox");
            paramList[11] = new Parameter("Max Peptide Length", "TextBox");
            paramList[12] = new Parameter("Min Peptide Length", "TextBox");
            paramList[13] = new Parameter("Max Missed Cleavages", "TextBox");
            paramList[14] = new Parameter("Initiator Methonine", "ComboBoxInit");
            paramList[15] = new Parameter("Deconvolution Mass Tolerance in PPM", "TextBox");
            paramList[16] = new Parameter("TrimMsMsPeaks", "Bool");
            paramList[17] = new Parameter("TrimMs1Peaks", "Bool");
            paramList[18] = new Parameter("MinRatio", "TextBox");
            paramList[19] = new Parameter("TopNpeaks", "TextBox");

            return paramList;
        }

        void Unset(object sender, RoutedEventArgs e)
        {
            int index = ParameterGrid.SelectedIndex;
            MessageBox.Show(index.ToString());
            Parameter a = ParameterGrid.Items[index] as Parameter;
            a.Value = null;
            ParameterGrid.Items.Refresh();
        }

        private void Cancel_Click(object sender, RoutedEventArgs e)
        {
            DialogResult = false;
        }

        static private int? TryParseNullable(string val)
        {
            return int.TryParse(val, out int outValue) ? (int?)outValue : null;
        }
    }
}


