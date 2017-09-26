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
            ParameterGrid.Items.Refresh();


            FileSpecificSettingsList = new FileSpecificSettings[selectedRaw.Count];


            if (selectedRaw.Count == 1)
            {
                string tomlFileName = System.IO.Path.ChangeExtension(SelectedRaw[0].FilePath, ".toml");
                Dictionary<string, KeyValuePair<string, Nett.TomlObject>> tomlSettingsList;
                if (File.Exists(tomlFileName))
                {
                    var paramFile = Toml.ReadFile(tomlFileName, MetaMorpheusTask.tomlConfig);
                    tomlSettingsList = paramFile.ToDictionary(p => p.Key);
                    FileSpecificSettings settings = new FileSpecificSettings(tomlSettingsList);

                    UpdateAndPopulateFields(settings);
                }
                else
                {
                    tomlSettingsList = new Dictionary<string, KeyValuePair<string, TomlObject>>();
                    FileSpecificSettings settings = new FileSpecificSettings(tomlSettingsList);
                    settings.ConserveMemory = null;
                    settings.DoPrecursorDeconvolution = null;
                    settings.UseProvidedPrecursorInfo = null;
                    UpdateAndPopulateFields(settings);
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

                    tomlFileNames[i] = System.IO.Path.ChangeExtension(SelectedRaw[i].FilePath, ".toml");
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
                foreach (string key in nonEqualValueNames)
                {
                    if (key == "Protease")
                    {
                        paramList[0].Different = true;
                    }

                    if (key == "ConserveMemory")
                    {
                        paramList[1].Different = true;
                    }

                    if (key == "Max_mods_for_peptide")
                    {
                        paramList[2].Different = true;
                    }
                    if (key == "ConserveMemory")
                    {
                        paramList[3].Different = true;
                    }
                    if (key == "DeconvolutionIntensityRatio")
                    {
                        paramList[4].Different = true;
                    }
                    if (key == "UseProvidedPrecursorInfo")
                    {
                        paramList[5].Different = true;
                    }

                    if (key == "ScoreCutoff")
                    {
                        paramList[6].Different = true;
                    }

                    if (key == "ProductMassTolerance")
                    {
                        paramList[7].Different = true;
                    }

                    if (key == "DeconvolutionMaxAssumedChargeState")
                    {
                        paramList[8].Different = true;
                    }

                    if (key == "TotalPartitions")
                    {
                        paramList[9].Different = true;
                    }

                    if (key == "MaxModificationIsoforms")
                    {
                        paramList[10].Different = true;
                    }

                    if (key == "MaxPeptideLength")
                    {
                        paramList[11].Different = true;
                    }

                    if (key == "MinPeptideLength")
                    {
                        paramList[12].Different = true;
                    }

                    if (key == "MaxMissedCleavages")
                    {
                        paramList[13].Different = true;
                    }

                    if (key == "InitiatorMethionineBehavior")
                    {
                        paramList[14].Different = true;
                    }

                    if (key == "DeconvolutionMassTolerance")
                    {
                        paramList[15].Different = true;
                    }

                    if (key == "TrimMsMsPeaks")
                    {
                        paramList[16].Different = true;
                    }

                    if (key == "TrimMs1Peaks")
                    {
                        paramList[17].Different = true;
                    }

                    if (key == "MinRatio")
                    {
                        paramList[18].Different = true;
                    }

                    if (key == "TopNpeaks")
                    {
                        paramList[19].Different = true;
                    }
                }
                UpdateAndPopulateFields(settings);
            }




        }
        private void UpdateAndPopulateFields(FileSpecificSettings settings)
        {

            var a = ParameterGrid.Items.GetItemAt(0) as Parameter;
            int? index = paramList[0].ProtList.IndexOf(settings.Protease);
            if (index.HasValue)
                a.Value = index;
            //GlobalEngineLevelSettings.ProteaseDictionary
            var b = ParameterGrid.Items.GetItemAt(1) as Parameter;
            b.Value = settings.ConserveMemory;
            var c = ParameterGrid.Items.GetItemAt(2) as Parameter;
            c.Value = settings.Max_mods_for_peptide;
            paramList[2].Value = c.Value;
            var d = ParameterGrid.Items.GetItemAt(3) as Parameter;
            d.Value = settings.DeconvolutionIntensityRatio;
            paramList[3].Value = d.Value;
            var e = ParameterGrid.Items.GetItemAt(4) as Parameter;
            e.Value = settings.DoPrecursorDeconvolution;
            var f = ParameterGrid.Items.GetItemAt(5) as Parameter;
            f.Value = settings.UseProvidedPrecursorInfo;
            var g = ParameterGrid.Items.GetItemAt(6) as Parameter;
            g.Value = settings.ScoreCutoff;
            var h = ParameterGrid.Items.GetItemAt(7) as Parameter;
            h.Value = settings.ProductMassTolerance.Value;
            var i = ParameterGrid.Items.GetItemAt(8) as Parameter;
            i.Value = settings.DeconvolutionMaxAssumedChargeState;
            var j = ParameterGrid.Items.GetItemAt(9) as Parameter;
            j.Value = settings.TotalPartitions;
            var k = ParameterGrid.Items.GetItemAt(10) as Parameter;
            k.Value = settings.MaxModificationIsoforms;
            var l = ParameterGrid.Items.GetItemAt(11) as Parameter;
            l.Value = settings.MaxPeptideLength;
            var m = ParameterGrid.Items.GetItemAt(12) as Parameter;
            m.Value = settings.MinPeptideLength;
            var n = ParameterGrid.Items.GetItemAt(13) as Parameter;
            n.Value = settings.MaxMissedCleavages;
            int? index2 = paramList[0].InitList.IndexOf(settings.InitiatorMethionineBehavior.ToString());
            var o = ParameterGrid.Items.GetItemAt(14) as Parameter;
            if (index2.HasValue && index2.Value != 0)
                o.Value = index2;
            else
                o.Value = null;
            var p = ParameterGrid.Items.GetItemAt(15) as Parameter;
            p.Value = settings.DeconvolutionMassTolerance.Value;
            var q = ParameterGrid.Items.GetItemAt(16) as Parameter;
            q.Value = settings.TrimMsMsPeaks;
            var r = ParameterGrid.Items.GetItemAt(17) as Parameter;
            r.Value = settings.TrimMs1Peaks;
            var s = ParameterGrid.Items.GetItemAt(18) as Parameter;
            s.Value = settings.MinRatio;
            var t = ParameterGrid.Items.GetItemAt(19) as Parameter;
            t.Value = settings.TopNpeaks;
            var u = ParameterGrid.Items.GetItemAt(20) as Parameter;
            

            ParameterGrid.Items.Refresh();

        }

        private void Row_DoubleClick(object sender, MouseButtonEventArgs e)
        {
            var ye = sender as DataGridCell;
            if (ye.Content is TextBlock hm &&
            !string.IsNullOrEmpty(hm.Text)) { }
        }

        private void Save_Click(object sender, RoutedEventArgs e)
        {
            ParameterGrid.Items.Refresh();

            for (int i = 0; i < FileSpecificSettingsList.Count(); i++)
            {
                FileSpecificSettingsList[i] = new FileSpecificSettings();
                int? index = paramList[0].Value as int?;

                string toleranceType = "PPM";
                if (paramList[20].Value != null)
                {
                    if ((int)paramList[20].Value == 0)
                        toleranceType = "Absolute";
                    if ((int)paramList[20].Value == 1)
                        toleranceType = "PPM";
                }
                    if (index.HasValue && index >= 0)
                {
                    FileSpecificSettingsList[i].Protease = paramList[0].ProtList[index.Value];
                }


                FileSpecificSettingsList[i].ConserveMemory = paramList[1].Value as bool?;
                if (paramList[2].Value != null)
                {
                    int.TryParse(paramList[2].Value.ToString(), out var a);
                    FileSpecificSettingsList[i].Max_mods_for_peptide = a;
                }
                if (paramList[3].Value != null)
                {
                    int.TryParse(paramList[3].Value.ToString(), out var a);
                    FileSpecificSettingsList[i].DeconvolutionIntensityRatio = a;
                }
                FileSpecificSettingsList[i].DoPrecursorDeconvolution = paramList[4].Value as bool?;
                FileSpecificSettingsList[i].UseProvidedPrecursorInfo = paramList[5].Value as bool?;
                if (paramList[6].Value != null)
                {
                    double.TryParse(paramList[6].Value.ToString(), out var a);
                    FileSpecificSettingsList[i].ScoreCutoff = a;
                }
                if (paramList[7].Value != null)
                    FileSpecificSettingsList[i].ProductMassTolerance = Tolerance.ParseToleranceString(paramList[7].Value + " " + toleranceType); ;

                if (paramList[8].Value != null)
                {
                    int.TryParse(paramList[8].Value.ToString(), out var a);
                    FileSpecificSettingsList[i].DeconvolutionMaxAssumedChargeState = a;
                }
                if (paramList[9].Value != null)
                {
                    int.TryParse(paramList[9].Value.ToString(), out var a);
                    FileSpecificSettingsList[i].TotalPartitions = a;
                }
                if (paramList[10].Value != null)
                {
                    int.TryParse(paramList[10].Value.ToString(), out var a);
                    FileSpecificSettingsList[i].MaxModificationIsoforms = a;
                }
                if (paramList[11].Value != null)
                {
                    int.TryParse(paramList[11].Value.ToString(), out var a);
                    FileSpecificSettingsList[i].MaxPeptideLength = a;
                }
                if (paramList[12].Value != null)
                {
                    int.TryParse(paramList[12].Value.ToString(), out var a);
                    FileSpecificSettingsList[i].MinPeptideLength = a;
                }
                if (paramList[13].Value != null)
                {
                    int.TryParse(paramList[13].Value.ToString(), out var a);
                    FileSpecificSettingsList[i].MaxMissedCleavages = a;
                }

                if (paramList[14].Value != null)
                    if (!paramList[14].Value.Equals("Undefined"))
                        FileSpecificSettingsList[i].InitiatorMethionineBehavior = (InitiatorMethionineBehavior)paramList[14].Value;
                    else
                        paramList[14].Value = null;
                if (paramList[15].Value != null)
                    FileSpecificSettingsList[i].DeconvolutionMassTolerance = Tolerance.ParseToleranceString(paramList[15].Value + " " + toleranceType);
                FileSpecificSettingsList[i].TrimMsMsPeaks = paramList[16].Value as bool?;
                FileSpecificSettingsList[i].TrimMs1Peaks = paramList[17].Value as bool?;
                if (paramList[18].Value != null)
                {
                    int.TryParse(paramList[18].Value.ToString(), out var a);
                    FileSpecificSettingsList[i].MinRatio = a;
                }

                if (paramList[19].Value != null)
                {
                    int.TryParse(paramList[19].Value.ToString(), out var a);
                    FileSpecificSettingsList[i].TopNpeaks = a;
                }

            }

            DialogResult = true;

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
            paramList = new Parameter[21];

            paramList[0] = new Parameter("Protease", "ComboBoxProtease");
            paramList[1] = new Parameter("Conserve Memory", "Bool");
            paramList[2] = new Parameter("Max Mods For Peptides", "TextBox");
            paramList[3] = new Parameter("Deconvolution Intensity Ratio", "TextBox");
            paramList[4] = new Parameter("Do Precursor Deconvolution", "Bool");
            paramList[5] = new Parameter("Use Provided Precursor Info", "Bool");
            paramList[6] = new Parameter("Minimum MM Score Allowed", "TextBox");
            paramList[7] = new Parameter("Product Mass Tolerance", "TextBox");
            paramList[8] = new Parameter("Deconvolution Max Assumed Charge State", "TextBox");
            paramList[9] = new Parameter("Number of Database Partitions", "TextBox");
            paramList[10] = new Parameter("Max Modification Isoforms", "TextBox");
            paramList[11] = new Parameter("Max Peptide Length", "TextBox");
            paramList[12] = new Parameter("Min Peptide Length", "TextBox");
            paramList[13] = new Parameter("Max Missed Cleavages", "TextBox");
            paramList[14] = new Parameter("Initiator Methonine", "ComboBoxInit");
            paramList[15] = new Parameter("Deconvolution Mass Tolerance", "TextBox");
            paramList[16] = new Parameter("TrimMsMsPeaks", "Bool");
            paramList[17] = new Parameter("TrimMs1Peaks", "Bool");
            paramList[18] = new Parameter("MinRatio", "TextBox");
            paramList[19] = new Parameter("TopNpeaks", "TextBox");
            paramList[20] = new Parameter("Tolerance Type", "ProductMassToleranceList");

            return paramList;
        }

        void Unset(object sender, RoutedEventArgs e)
        {
            int index = ParameterGrid.SelectedIndex;
            Parameter a = ParameterGrid.Items[index] as Parameter;
            a.Value = null;
            ParameterGrid.Items.Refresh();
        }

        private void Cancel_Click(object sender, RoutedEventArgs e)
        {
            DialogResult = false;
        }

    }
}


