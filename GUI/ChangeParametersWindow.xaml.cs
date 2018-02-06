using EngineLayer;
using MzLibUtil;
using Nett;
using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.IO;
using System.Linq;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Input;
using TaskLayer;

namespace MetaMorpheusGUI
{
    /// <summary>
    /// Interaction logic for ChangeParametersWindow.xaml
    /// </summary>
    public partial class ChangeParametersWindow : Window
    {
        #region Private Fields

        private Parameter[] paramList;
        private Parameter[] tempParamList;

        #endregion Private Fields

        #region Public Constructors

        public ChangeParametersWindow(ObservableCollection<RawDataForDataGrid> selectedRaw)
        {
            SelectedRaw = selectedRaw;
            FileSpecificSettingsList = new FileSpecificSettings[selectedRaw.Count];
            InitializeComponent();
            tempParamList = InitializeParameterList();
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
                    Array.Copy(paramList, tempParamList, 20);
                    UpdateAndPopulateFields(settings);
                }
                else
                {
                    tomlSettingsList = new Dictionary<string, KeyValuePair<string, TomlObject>>();
                    FileSpecificSettings settings = new FileSpecificSettings(tomlSettingsList);
                    settings.ConserveMemory = null;
                    settings.DoPrecursorDeconvolution = null;
                    settings.UseProvidedPrecursorInfo = null;
                    Array.Copy(paramList, tempParamList, 20);
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
                        tomlSettingsListList.Add(new Dictionary<string, KeyValuePair<string, TomlObject>>());
                    }
                }

                //check every value against every other value; if same, then add to list to set 'different in multiple files' later in Parameter Object
                for (int i = 0; i < tomlSettingsListList.Count; i++)
                {
                    for (int j = 0; j < tomlSettingsListList.Count; j++)
                    {
                        foreach (var key in tomlSettingsListList[j].Keys)
                        {
                            if (tomlSettingsListList[i].ContainsKey(key))
                            {
                                if (tomlSettingsListList[j].ContainsKey(tomlSettingsListList[i][key].Key))
                                {
                                    //Need to convert first letter yo uppercase to match Type
                                    string h = tomlSettingsListList[j][key].Value.ReadableTypeName;
                                    //convert string to appropritate format
                                    if (!string.IsNullOrEmpty(h))
                                    {
                                        h = h.First().ToString().ToUpper() + h.Substring(1);
                                    }
                                    if (h == "Int")
                                    {
                                        h = "Int32";
                                    }
                                    if (h == "Bool")
                                    {
                                        h = "Boolean";
                                    }
                                    if (h == "Float")
                                    {
                                        h = "Single";
                                    }

                                    Type typeAsString = Type.GetType("System." + h);
                                    var a = tomlSettingsListList[j][key].Value.Get(typeAsString);

                                    var b = tomlSettingsListList[i][key].Value.Get(typeAsString);

                                    if (!a.Equals(b))
                                    {
                                        nonEqualValueNames.Add(tomlSettingsListList[j][key].Key);
                                    }
                                }
                            }
                            //add different checked if null and others are not
                            else nonEqualValueNames.Add(key);
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

                    if (key == "DeconvolutionIntensityRatio")
                    {
                        paramList[3].Different = true;
                    }

                    if (key == "DoPrecursorDeconvolution")
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
                TempSettings = new FileSpecificSettings[selectedRaw.Count];
                Array.Copy(paramList, tempParamList, 20);
                for (int i = 0; i < tomlSettingsListList.Count(); i++)
                {
                    TempSettings[i] = new FileSpecificSettings(tomlSettingsListList[i]);
                }

                UpdateAndPopulateFields(settings);
            }
        }

        #endregion Public Constructors

        #region Internal Properties

        internal ObservableCollection<RawDataForDataGrid> SelectedRaw { get; private set; }
        internal FileSpecificSettings[] FileSpecificSettingsList { get; private set; }

        #endregion Internal Properties

        #region Private Properties

        private FileSpecificSettings[] TempSettings { get; set; }

        #endregion Private Properties

        #region Private Methods

        private static Boolean TextBoxIntAllowed(String Text2)
        {
            return Array.TrueForAll<Char>(Text2.ToCharArray(),
                delegate (Char c) { return Char.IsDigit(c) || Char.IsControl(c); });
        }

        //Updates fields of display so that it reflects current settings
        private void UpdateAndPopulateFields(FileSpecificSettings settings)
        {
            var a = ParameterGrid.Items.GetItemAt(0) as Parameter;
            int? index = paramList[0].ProtList.IndexOf(settings.Protease);
            if (index.HasValue)
                a.Value = index;
            if (paramList[0].Different)
                a.Value = null;
            var b = ParameterGrid.Items.GetItemAt(1) as Parameter;
            b.Value = settings.ConserveMemory;
            if (paramList[1].Different)
                b.Value = null;
            var c = ParameterGrid.Items.GetItemAt(2) as Parameter;
            c.Value = settings.Max_mods_for_peptide;
            if (paramList[2].Different)
                c.Value = null;
            var d = ParameterGrid.Items.GetItemAt(3) as Parameter;
            d.Value = settings.DeconvolutionIntensityRatio;
            paramList[3].Value = d.Value;
            if (paramList[3].Different)
                d.Value = null;
            var e = ParameterGrid.Items.GetItemAt(4) as Parameter;
            e.Value = settings.DoPrecursorDeconvolution;
            if (paramList[4].Different)
                e.Value = null;
            var f = ParameterGrid.Items.GetItemAt(5) as Parameter;
            f.Value = settings.UseProvidedPrecursorInfo;
            if (paramList[5].Different)
                f.Value = null;
            var g = ParameterGrid.Items.GetItemAt(6) as Parameter;
            g.Value = settings.ScoreCutoff;
            if (paramList[6].Different)
                g.Value = null;
            var h = ParameterGrid.Items.GetItemAt(7) as Parameter;
            if (settings.ProductMassTolerance != null)
                h.Value = settings.ProductMassTolerance.Value;
            if (paramList[7].Different)
                h.Value = null;
            var i = ParameterGrid.Items.GetItemAt(8) as Parameter;
            i.Value = settings.DeconvolutionMaxAssumedChargeState;
            if (paramList[8].Different)
                i.Value = null;
            var j = ParameterGrid.Items.GetItemAt(9) as Parameter;
            j.Value = settings.TotalPartitions;
            if (paramList[9].Different)
                j.Value = null;
            var k = ParameterGrid.Items.GetItemAt(10) as Parameter;
            k.Value = settings.MaxModificationIsoforms;
            if (paramList[10].Different)
                k.Value = null;
            var l = ParameterGrid.Items.GetItemAt(11) as Parameter;
            l.Value = settings.MaxPeptideLength;
            if (paramList[11].Different)
                l.Value = null;
            var m = ParameterGrid.Items.GetItemAt(12) as Parameter;
            m.Value = settings.MinPeptideLength;
            if (paramList[12].Different)
                m.Value = null;
            var n = ParameterGrid.Items.GetItemAt(13) as Parameter;
            n.Value = settings.MaxMissedCleavages;
            if (paramList[13].Different)
                n.Value = null;
            int? index2 = paramList[0].InitList.IndexOf(settings.InitiatorMethionineBehavior.ToString());
            var o = ParameterGrid.Items.GetItemAt(14) as Parameter;
            if (index2.HasValue && index2.Value != 0)
                o.Value = index2;
            else
                o.Value = null;
            if (paramList[14].Different)
                o.Value = null;
            var p = ParameterGrid.Items.GetItemAt(15) as Parameter;
            if (settings.DeconvolutionMassTolerance != null)
                p.Value = settings.DeconvolutionMassTolerance.Value;
            if (paramList[15].Different)
                p.Value = null;
            var q = ParameterGrid.Items.GetItemAt(16) as Parameter;
            q.Value = settings.TrimMsMsPeaks;
            if (paramList[16].Different)
                q.Value = null;
            var r = ParameterGrid.Items.GetItemAt(17) as Parameter;
            r.Value = settings.TrimMs1Peaks;
            if (paramList[17].Different)
                r.Value = null;
            var s = ParameterGrid.Items.GetItemAt(18) as Parameter;
            s.Value = settings.MinRatio;
            if (paramList[18].Different)
                s.Value = null;
            var t = ParameterGrid.Items.GetItemAt(19) as Parameter;
            t.Value = settings.TopNpeaks;
            if (paramList[19].Different)
                t.Value = null;
            var v = ParameterGrid.Items.GetItemAt(0) as Parameter;
            //tolerance type
            int? index3 = paramList[20].ProductMassToleranceList.IndexOf(settings.ToleranceType);
            var w = ParameterGrid.Items.GetItemAt(0) as Parameter;
            if (index.HasValue)
                w.Value = index;
            if (paramList[20].Different)
                w.Value = null;

            ParameterGrid.Items.Refresh();
        }

        private void Row_DoubleClick(object sender, MouseButtonEventArgs e)
        {
            var ye = sender as DataGridCell;
            if (ye.Content is TextBlock hm &&
            !string.IsNullOrEmpty(hm.Text)) { }
        }

        //upon clicking save button, settings are saved to memory in FileSettingsList in order to be
        //written to a toml file later
        private void Save_Click(object sender, RoutedEventArgs e)
        {
            ParameterGrid.Items.Refresh();
            Array.Copy(tempParamList, paramList, 20);

            for (int i = 0; i < FileSpecificSettingsList.Count(); i++)
            {
                FileSpecificSettingsList[i] = new FileSpecificSettings();
                int? index = paramList[0].Value as int?;
               

                //Cases for each Parameter:
                //1. If one file with Changed value: write value to File settings
                //2. If one file with unchanged Value: Do nothing; will be null
                //3. If multiple files editing the same setting: write same setting to each file
                /* 4. If multiple files with different values for an attribute that is unchanged, then get each paramter setting
                from Temporary Settings. If this parameter is changed, then it will be the same so going to temp settings
                is Not necesarry
                */

                //Protease
                if (index.HasValue && index >= 0 && paramList[0].HasChanged)
                {
                    FileSpecificSettingsList[i].Protease = paramList[0].ProtList[index.Value];
                }
                else if (FileSpecificSettingsList.Count() > 1 && TempSettings[i].Protease != null && FileSpecificSettingsList[i].Protease == null)
                {
                    FileSpecificSettingsList[i].Protease = TempSettings[i].Protease;
                }

                //Conserve Memory
                if (paramList[1].Value != null)
                    FileSpecificSettingsList[i].ConserveMemory = paramList[1].Value as bool?;
                else if (FileSpecificSettingsList.Count() > 1 && TempSettings[i].ConserveMemory != null && FileSpecificSettingsList[i].ConserveMemory == null)
                    FileSpecificSettingsList[i].ConserveMemory = TempSettings[i].ConserveMemory;

                //Max Mods for Peptide
                if (paramList[2].Value != null)
                {
                    int.TryParse(paramList[2].Value.ToString(), out var a);
                    FileSpecificSettingsList[i].Max_mods_for_peptide = a;
                }
                else if (FileSpecificSettingsList.Count() > 1 && TempSettings[i].Max_mods_for_peptide != null && FileSpecificSettingsList[i].Max_mods_for_peptide == null)
                    FileSpecificSettingsList[i].Max_mods_for_peptide = TempSettings[i].Max_mods_for_peptide;

                //Deconvolution Intensity Ratio
                if (paramList[3].Value != null)
                {
                    int.TryParse(paramList[3].Value.ToString(), out var a);
                    FileSpecificSettingsList[i].DeconvolutionIntensityRatio = a;
                }
                else if (FileSpecificSettingsList.Count() > 1 && TempSettings[i].DeconvolutionIntensityRatio != null && FileSpecificSettingsList[i].DeconvolutionIntensityRatio == null)
                    FileSpecificSettingsList[i].DeconvolutionIntensityRatio = TempSettings[i].DeconvolutionIntensityRatio;

                //Precursor Deconvolution
                if (paramList[4].Value != null)
                    FileSpecificSettingsList[i].DoPrecursorDeconvolution = paramList[4].Value as bool?;
                else if (FileSpecificSettingsList.Count() > 1 && TempSettings[i].DoPrecursorDeconvolution != null && FileSpecificSettingsList[i].DoPrecursorDeconvolution == null)
                    FileSpecificSettingsList[i].DoPrecursorDeconvolution = TempSettings[i].DoPrecursorDeconvolution;

                //Use Provided Precursor Info
                if (paramList[5].Value != null)
                    FileSpecificSettingsList[i].UseProvidedPrecursorInfo = paramList[5].Value as bool?;
                else if (FileSpecificSettingsList.Count() > 1 && TempSettings[i].UseProvidedPrecursorInfo != null && FileSpecificSettingsList[i].UseProvidedPrecursorInfo == null)
                    FileSpecificSettingsList[i].UseProvidedPrecursorInfo = TempSettings[i].UseProvidedPrecursorInfo;

                //Score Cutoff
                if (paramList[6].Value != null)
                {
                    double.TryParse(paramList[6].Value.ToString(), out var a);
                    FileSpecificSettingsList[i].ScoreCutoff = a;
                }
                else if (FileSpecificSettingsList.Count() > 1 && TempSettings[i].ScoreCutoff != null && FileSpecificSettingsList[i].ScoreCutoff == null)
                    FileSpecificSettingsList[i].ScoreCutoff = TempSettings[i].ScoreCutoff;

                //Product Mass Tolerance
                if (paramList[7].Value != null)
                    FileSpecificSettingsList[i].ProductMassTolerance = Tolerance.ParseToleranceString(paramList[7].Value + " " + FileSpecificSettingsList[i].ToleranceType);
                else if (FileSpecificSettingsList.Count() > 1 && TempSettings[i].ProductMassTolerance != null && FileSpecificSettingsList[i].ProductMassTolerance == null)
                    FileSpecificSettingsList[i].ProductMassTolerance = TempSettings[i].ProductMassTolerance;

                //Deconvolution Max ASsumed Charge State
                if (paramList[8].Value != null)
                {
                    int.TryParse(paramList[8].Value.ToString(), out var a);
                    FileSpecificSettingsList[i].DeconvolutionMaxAssumedChargeState = a;
                }
                else if (FileSpecificSettingsList.Count() > 1 && TempSettings[i].DeconvolutionMaxAssumedChargeState != null && FileSpecificSettingsList[i].DeconvolutionMaxAssumedChargeState == null)
                    FileSpecificSettingsList[i].DeconvolutionMaxAssumedChargeState = TempSettings[i].DeconvolutionMaxAssumedChargeState;

                //Total Partitions
                if (paramList[9].Value != null)
                {
                    int.TryParse(paramList[9].Value.ToString(), out var a);
                    FileSpecificSettingsList[i].TotalPartitions = a;
                }
                else if (FileSpecificSettingsList.Count() > 1 && TempSettings[i].TotalPartitions != null && FileSpecificSettingsList[i].TotalPartitions == null)
                    FileSpecificSettingsList[i].TotalPartitions = TempSettings[i].TotalPartitions;

                //Max Modification Isoforms
                if (paramList[10].Value != null)
                {
                    int.TryParse(paramList[10].Value.ToString(), out var a);
                    FileSpecificSettingsList[i].MaxModificationIsoforms = a;
                }
                else if (FileSpecificSettingsList.Count() > 1 && TempSettings[i].MaxModificationIsoforms != null && FileSpecificSettingsList[i].MaxModificationIsoforms == null)
                    FileSpecificSettingsList[i].MaxModificationIsoforms = TempSettings[i].MaxModificationIsoforms;

                //Max Pep Length
                if (paramList[11].Value != null)
                {
                    int.TryParse(paramList[11].Value.ToString(), out var a);
                    FileSpecificSettingsList[i].MaxPeptideLength = a;
                }
                else if (FileSpecificSettingsList.Count() > 1 && TempSettings[i].MaxPeptideLength != null && FileSpecificSettingsList[i].MaxPeptideLength == null)
                    FileSpecificSettingsList[i].MaxPeptideLength = TempSettings[i].MaxPeptideLength;

                //Min Pep Length
                if (paramList[12].Value != null)
                {
                    int.TryParse(paramList[12].Value.ToString(), out var a);
                    FileSpecificSettingsList[i].MinPeptideLength = a;
                }
                else if (FileSpecificSettingsList.Count() > 1 && TempSettings[i].MinPeptideLength != null && FileSpecificSettingsList[i].MinPeptideLength == null)
                    FileSpecificSettingsList[i].MinPeptideLength = TempSettings[i].MinPeptideLength;

                //Max Missed Cleavages
                if (paramList[13].Value != null)
                {
                    int.TryParse(paramList[13].Value.ToString(), out var a);
                    FileSpecificSettingsList[i].MaxMissedCleavages = a;
                }
                else if (FileSpecificSettingsList.Count() > 1 && TempSettings[i].MaxMissedCleavages != null && FileSpecificSettingsList[i].MaxMissedCleavages == null)
                    FileSpecificSettingsList[i].MaxMissedCleavages = TempSettings[i].MaxMissedCleavages;

                //Init Methonine Behavior
                if (paramList[14].Value != null)
                    if (!paramList[14].Value.Equals("Undefined"))
                        FileSpecificSettingsList[i].InitiatorMethionineBehavior = (InitiatorMethionineBehavior)paramList[14].Value;
                    else
                        paramList[14].Value = null;
                else if (FileSpecificSettingsList.Count() > 1 && TempSettings[i].InitiatorMethionineBehavior != 0 && FileSpecificSettingsList[i].InitiatorMethionineBehavior == 0)
                    FileSpecificSettingsList[i].InitiatorMethionineBehavior = TempSettings[i].InitiatorMethionineBehavior;

                //Deconvolution Mass Tolerance
                if (paramList[15].Value != null)
                    FileSpecificSettingsList[i].DeconvolutionMassTolerance = Tolerance.ParseToleranceString(paramList[15].Value + " " + FileSpecificSettingsList[i].ToleranceType);
                else if (FileSpecificSettingsList.Count() > 1 && TempSettings[i].DeconvolutionMassTolerance != null && FileSpecificSettingsList[i].DeconvolutionMassTolerance == null)
                    FileSpecificSettingsList[i].DeconvolutionMassTolerance = TempSettings[i].DeconvolutionMassTolerance;

                //Trim Ms/Ms Peaks
                if (paramList[16].Value != null)
                    FileSpecificSettingsList[i].TrimMsMsPeaks = paramList[16].Value as bool?;
                else if (FileSpecificSettingsList.Count() > 1 && TempSettings[i].TrimMsMsPeaks != null && FileSpecificSettingsList[i].TrimMsMsPeaks == null)
                    FileSpecificSettingsList[i].TrimMsMsPeaks = TempSettings[i].TrimMsMsPeaks;

                //TrimMs1Peaks
                if (paramList[17].Value != null)
                    FileSpecificSettingsList[i].TrimMs1Peaks = paramList[17].Value as bool?;
                else if (FileSpecificSettingsList.Count() > 1 && TempSettings[i].TrimMs1Peaks != null && FileSpecificSettingsList[i].TrimMs1Peaks == null)
                    FileSpecificSettingsList[i].TrimMs1Peaks = TempSettings[i].TrimMs1Peaks;

                //min Ratio
                if (paramList[18].Value != null)
                {
                    int.TryParse(paramList[18].Value.ToString(), out var a);
                    FileSpecificSettingsList[i].MinRatio = a;
                }
                else if (FileSpecificSettingsList.Count() > 1 && TempSettings[i].MinRatio != null && FileSpecificSettingsList[i].MinRatio == null)
                    FileSpecificSettingsList[i].MinRatio = TempSettings[i].MinRatio;

                //Top N Peaks
                if (paramList[19].Value != null)
                {
                    int.TryParse(paramList[19].Value.ToString(), out var a);
                    FileSpecificSettingsList[i].TopNpeaks = a;
                }
                else if (FileSpecificSettingsList.Count() > 1 && TempSettings[i].TopNpeaks != null && FileSpecificSettingsList[i].TopNpeaks == null)
                    FileSpecificSettingsList[i].TopNpeaks = TempSettings[i].TopNpeaks;

                int? index2 = paramList[20].Value as int?;


                if (index2.HasValue && index >= 0 && paramList[20].HasChanged)
                {
                    FileSpecificSettingsList[i].ToleranceType = paramList[20].ProductMassToleranceList[index2.Value];
                }
                else if (FileSpecificSettingsList.Count() > 1 && TempSettings[i].ToleranceType != null && FileSpecificSettingsList[i].ToleranceType == null)
                {
                    FileSpecificSettingsList[i].ToleranceType = TempSettings[i].ToleranceType;
                }
            }

            DialogResult = true;
        }

        private void PreviewIfInt(object sender, TextCompositionEventArgs e)
        {
            e.Handled = !TextBoxIntAllowed(e.Text);
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

        //Undoes value(sets to null)
        private void Unset(object sender, RoutedEventArgs e)
        {
            int index = ParameterGrid.SelectedIndex;
            Parameter a = ParameterGrid.Items[index] as Parameter;
            a.Value = null;
            if (TempSettings != null)
            {
                for (int i = 0; i < TempSettings.Length; i++)
                {
                    FileSpecificSettings b = TempSettings.GetValue(i) as FileSpecificSettings;
                    switch (index)
                    {
                        case 0:
                            b.Protease = null;
                            break;

                        case 1:
                            b.ConserveMemory = null;
                            break;

                        case 2:
                            b.Max_mods_for_peptide = null;
                            break;

                        case 3:
                            b.DeconvolutionIntensityRatio = null;
                            break;

                        case 4:
                            b.DoPrecursorDeconvolution = null;
                            break;

                        case 5:
                            b.UseProvidedPrecursorInfo = null;
                            break;

                        case 6:
                            b.ScoreCutoff = null;
                            break;

                        case 7:
                            b.ProductMassTolerance = null;
                            break;

                        case 8:
                            b.DeconvolutionMaxAssumedChargeState = null;
                            break;

                        case 9:
                            b.TotalPartitions = null;
                            break;

                        case 10:
                            b.MaxModificationIsoforms = null;
                            break;

                        case 11:
                            b.MaxPeptideLength = null;
                            break;

                        case 12:
                            b.MinPeptideLength = null;
                            break;

                        case 13:
                            b.MaxMissedCleavages = null;
                            break;

                        case 14:
                            b.InitiatorMethionineBehavior = 0;
                            break;

                        case 15:
                            b.DeconvolutionMassTolerance = null;
                            break;

                        case 16:
                            b.TrimMsMsPeaks = null;
                            break;

                        case 17:
                            b.TrimMs1Peaks = null;
                            break;

                        case 18:
                            b.MinRatio = null;
                            break;

                        case 19:
                            b.TopNpeaks = null;
                            break;
                        case 20:
                            b.ToleranceType = null;
                            break;
                    }
                }
            }

            ParameterGrid.Items.Refresh();
        }

        private void Cancel_Click(object sender, RoutedEventArgs e)
        {
            DialogResult = false;
        }

        #endregion Private Methods
    }
}