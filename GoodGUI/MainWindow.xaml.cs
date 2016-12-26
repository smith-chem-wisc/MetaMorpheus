using MetaMorpheus;
using mzCal;
using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Threading;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;

namespace GoodGUI
{
    /// <summary>
    /// Interaction logic for MainWindow.xaml
    /// </summary>
    public partial class MainWindow : Window
    {
        private static ObservableCollection<RawDataAndResults> rawDataAndResultslist = new ObservableCollection<RawDataAndResults>();
        private static ObservableCollection<XMLdb> xMLdblist = new ObservableCollection<XMLdb>();
        private static ObservableCollection<GptmdPTMlist> GPTMDfileList = new ObservableCollection<GptmdPTMlist>();
        private static ObservableCollection<ModList> SearchfileList = new ObservableCollection<ModList>();
        private ParamsObject po;

        private enum Task
        {
            Search,
            GPTMD,
            Calibrate
        }

        public MainWindow()
        {
            InitializeComponent();

            mzCalIO.mzCalIO.Load();

            dataGridGPTMD.DataContext = GPTMDfileList;
            dataGridSearch.DataContext = SearchfileList;
            dataGridXMLs.DataContext = xMLdblist;
            dataGridDatafilesAndResults.DataContext = rawDataAndResultslist;

            GPTMDfileList.Add(new GptmdPTMlist("m.txt", true));
            GPTMDfileList.Add(new GptmdPTMlist("r.txt", false));
            GPTMDfileList.Add(new GptmdPTMlist("s.txt", false));

            SearchfileList.Add(new ModList("v.txt", false, true, false));
            SearchfileList.Add(new ModList("f.txt", false, false, true));
            SearchfileList.Add(new ModList("p.txt", true, false, false));
            SearchfileList.Add(new ModList("m.txt", true, false, false));
            SearchfileList.Add(new ModList("r.txt", false, false, false));
            SearchfileList.Add(new ModList("s.txt", false, false, false));

            // RAW FILES
            //addFile(@"C:\Users\stepa\Data\CalibrationPaperData\OrigData\Mouse\04-29-13_B6_Frac1_9uL.raw");
            //addFile(@"C:\Users\stepa\Data\CalibrationPaperData\OrigData\Mouse\04-29-13_B6_Frac2_9p5uL.raw");

            // MZID FILES
            //addFile(@"C:\Users\stepa\Data\CalibrationPaperData\OrigData\Mouse\2016-10-20-09-17\04-29-13_B6_Frac1_9uL.mzid");

            // XML
            addFile(@"C:\Users\stepa\Data\CalibrationPaperData\OrigData\uniprot-mouse-reviewed-10-3-2016.xml");
            //addFile(@"C:\Users\stepa\Data\CalibrationPaperData\OrigData\uniprot-human-reviewed-10-3-2016.xml");
            //addFile(@"C:\Users\stepa\Data\CalibrationPaperData\OrigData\cRAP-11-11-2016.xml");

            // Calib FILES
            //addFile(@"C:\Users\stepa\Data\CalibrationPaperData\OrigData\Mouse\04-29-13_B6_Frac1_9uL-Calibrated.mzML");
            //addFile(@"C:\Users\stepa\Data\CalibrationPaperData\Step2\Mouse\Calib-0.1.2\04-29-13_B6_Frac1_9uL-Calibrated.mzML");
            addFile(@"C:\Users\stepa\Data\CalibrationPaperData\Step2\Mouse\Calib-0.1.2\04-29-13_B6_Frac9_9p5uL-Calibrated.mzML");
            //addFile(@"C:\Users\stepa\Data\CalibrationPaperData\Step2\Jurkat\Calib-0.1.2\120426_Jurkat_highLC_Frac16-Calibrated.mzML");

            // TSV file
            //addFile(@"C:\Users\stepa\Data\CalibrationPaperData\OrigData\Mouse\2016-10-20-13-24\04-29-13_B6_Frac1_9uL-Calibrated.psmtsv");
            //addFile(@"C:\Users\stepa\Data\CalibrationPaperData\Step2\Mouse\Calib-0.1.2\2016-10-21-11-07\04-29-13_B6_Frac1_9uL-Calibrated.psmtsv");
            //addFile(@"C:\Users\stepa\Data\CalibrationPaperData\Step2\Mouse\Calib-0.1.2\2016-10-21-15-18\aggregate.psmtsv");

            foreach (Protease protease in ProteaseDictionary.Instance.Values)
                comboBox.Items.Add(protease);
            comboBox.SelectedIndex = 12;
            RegOutput(ProteaseDictionary.Instance.Count + " proteases loaded from proteases.tsv");
            foreach (string initiatior_methionine_behavior in Enum.GetNames(typeof(InitiatorMethionineBehavior)))
                cboInitiatorMethionineBehavior.Items.Add(initiatior_methionine_behavior.ToLower());
            cboInitiatorMethionineBehavior.SelectedIndex = 2;
            comboBox1.Items.Add("Da");
            comboBox1.Items.Add("ppm");
            comboBox1.SelectedIndex = 1;
            comboBox1_Copy.Items.Add("Da");
            comboBox1_Copy.Items.Add("ppm");
            comboBox1_Copy.SelectedIndex = 0;
            AminoAcidMasses.LoadAminoAcidMasses();
            RegOutput("Amino acid masses loaded from amino_acids.tsv");

            try
            {
                using (StreamReader amino_acids = new StreamReader(Path.Combine(Path.GetDirectoryName(Environment.GetCommandLineArgs()[0]), "LimitScans.txt")))
                {
                    while (amino_acids.Peek() != -1)
                    {
                        SpectraLimits.Add(Convert.ToInt32(amino_acids.ReadLine()));
                    }
                    RegOutput("Will work with " + SpectraLimits.Count + " spectra read from LimitScans.txt");
                }
            }
            catch (FileNotFoundException)
            {
                RegOutput("LimitScans.txt not found, looking at all scans");
            }

            try
            {
                using (StreamReader amino_acids = new StreamReader(Path.Combine(Path.GetDirectoryName(Environment.GetCommandLineArgs()[0]), "LimitProteins.txt")))
                {
                    while (amino_acids.Peek() != -1)
                    {
                        SpecificProteinSelection.Add(amino_acids.ReadLine());
                    }
                    RegOutput("Will work with " + SpecificProteinSelection.Count + " proteins read from LimitProteins.txt");
                }
            }
            catch (FileNotFoundException)
            {
                RegOutput("LimitProteins.txt not found, looking at all proteins");
            }

            //myListOfEntries.Add(new DataPath(@"C:\Users\stepa\Data\CalibrationPaperData\Step2\Mouse\Calib-0.1.2\04-29-13_B6_Frac1_9uL-Calibrated.mzML"));
            //myListOfEntries.Add(new DataPath(@"C:\Users\stepa\Data\CalibrationPaperData\Step2\Mouse\Calib-0.1.2\04-29-13_B6_Frac6_5uL-Calibrated.mzML"));
            //myListOfEntries.Add(new DataPath(@"C:\Users\stepa\Data\CalibrationPaperData\Step2\Mouse\Calib-0.1.2\04-29-13_B6_Frac9_9p5uL-Calibrated.mzML"));

            //XMLtextBox.Text = @"C:\Users\stepa\Data\CalibrationPaperData\Step4\Mouse.XML";
            //textBox3.Text = @" - 121.95; -120.95; -119.95; -118.95; -117.95; -116.95; -115.95; -114.95; -113.95; -112.95; -111.95; -110.95; -109.95; -108.95; -107.95; -106.95; -105.95; -104.95; -103.95; -102.95; -101.95; -100.95; -99.95; -98.95; -97.95; -96.95; -95.95; -94.95; -93.95; -92.95; -91.95; -90.95; -89.95; -88.95; -87.95; -86.95; -85.95; -84.95; -83.95; -82.95; -81.95; -80.95; -79.95; -78.95; -77.95; -76.95; -75.95; -74.95; -73.95; -72.95; -71.95; -70.95; -69.95; -68.95; -67.95; -66.95; -65.95; -64.95; -63.95; -62.95; -61.95; -60.95; -59.95; -58.95; -57.95; -56.95; -55.95; -54.95; -53.95; -52.95; -51.95; -50.95; -49.95; -48.95; -47.95; -46.95; -45.95; -44.95; -43.95; -42.95; -41.95; -40.95; -39.95; -38.95; -37.95; -36.95; -35.95; -34.95; -33.95; -32.95; -31.95; -30.95; -29.95; -28.95; -27.95; -26.95; -25.95; -24.95; -23.95; -22.95; -21.95; -20.95; -19.95; -18.95; -17.95; -16.95; -15.95; -14.95; -13.95; -12.95; -11.95; -10.95; -9.95; -8.95; -7.95; -6.95; -5.95; -4.95; -3.95; -2.95; -1.95; -0.95; 0.05; 1.05; 2.05; 3.05; 4.05; 5.05; 6.05; 7.05; 8.05; 9.05; 10.05; 11.05; 12.05; 13.05; 14.05; 15.05; 16.05; 17.05; 18.05; 19.05; 20.05; 21.05; 22.05; 23.05; 24.05; 25.05; 26.05; 27.05; 28.05; 29.05; 30.05; 31.05; 32.05; 33.05; 34.05; 35.05; 36.05; 37.05; 38.05; 39.05; 40.05; 41.05; 42.05; 43.05; 44.05; 45.05; 46.05; 47.05; 48.05; 49.05; 50.05; 51.05; 52.05; 53.05; 54.05; 55.05; 56.05; 57.05; 58.05; 59.05; 60.05; 61.05; 62.05; 63.05; 64.05; 65.05; 66.05; 67.05; 68.05; 69.05; 70.05; 71.05; 72.05; 73.05; 74.05; 75.05; 76.05; 77.05; 78.05; 79.05; 80.05; 81.05; 82.05; 83.05; 84.05; 85.05; 86.05; 87.05; 88.05; 89.05; 90.05; 91.05; 92.05; 93.05; 94.05; 95.05; 96.05; 97.05; 98.05; 99.05; 100.05; 101.05; 102.05; 103.05; 104.05; 105.05; 106.05; 107.05; 108.05; 109.05; 110.05; 111.05; 112.05; 113.05; 114.05; 115.05; 116.05; 117.05; 118.05; 119.05; 120.05; 121.05; 122.05; 123.05; 124.05; 125.05; 126.05; 127.05; 128.05; 129.05; 130.05; 131.05; 132.05; 133.05; 134.05; 135.05; 136.05; 137.05; 138.05; 139.05; 140.05; 141.05; 142.05; 143.05; 144.05; 145.05; 146.05; 147.05; 148.05; 149.05; 150.05; 151.05; 152.05; 153.05; 154.05; 155.05; 156.05; 157.05; 158.05; 159.05; 160.05; 161.05; 162.05; 163.05; 164.05; 165.05; 166.05; 167.05; 168.05; 169.05; 170.05; 171.05; 172.05; 173.05; 174.05; 175.05; 176.05; 177.05; 178.05; 179.05; 180.05; 181.05; 182.05; 183.05; 184.05; 185.05; 186.05; 187.05; 188.05; 189.05; 190.05; 191.05; 192.05; 193.05; 194.05; 195.05; 196.05; 197.05; 198.05; 199.05; 200.05; 201.05; 202.05; 203.05; 204.05; 205.05; 206.05; 207.05; 208.05; 209.05; 210.05; 211.05; 212.05; 213.05; 214.05; 215.05; 216.05; 217.05; 218.05; 219.05; 220.05";
            //textBox3.Text = @"-15.95; 0.05; 16.05";
            //textBox3.Text = @"-17.026549";

            textBox3.Text = @"0";
            //textBox3.Text = @"-18; -17; 0; 1; 2; 12; 16; 22; 32; 42; 80; 113";
            //textBox3.Text = @"0.05; 1.05; 12.05; 16.05; 80.05";

            //textBox2.Text = @"0.15";
            //textBox2.Text = @"5";
            //textBox2.Text = @"10";
            textBox2.Text = @"1000";
            comboBox1.SelectedIndex = 0; // Da
            //comboBox1.SelectedIndex = 1; // ppm

            //checkBoxDecoy.IsChecked = false;
            checkBoxDecoy.IsChecked = true;

            //maxModificationIsoformsTextBox.Text = 10000.ToString();

            po = new ParamsObject();
            po.outLabelStatusHandler += NewoutLabelStatus;
            po.outProgressBarHandler += NewoutProgressBar;
            po.outRichTextBoxHandler += NewoutRichTextBox;
            po.outSuccessfullyFinishedTaskHandler += NewSuccessfullyFinishedTask;
            po.outSuccessfullyStartingTaskHandler += NewSuccessfullyStartingTask;
            po.SuccessfullyFinishedFileHandler += NewSuccessfullyFinishedFile;
            po.outTextBoxHandler += NewoutTextBox;
        }

        private void addFile(string filepath)
        {
            var theExtension = Path.GetExtension(filepath);
            filepath = Path.GetFullPath(filepath);
            RawDataAndResults res;
            switch (theExtension)
            {
                case ".raw":
                case ".mzML":
                    res = GetCorrespondingRawDataAndResultsEntry(filepath);
                    if (res != null)
                        res.AddFilePath(filepath);
                    else
                        rawDataAndResultslist.Add(new RawDataAndResults(filepath, null, null));
                    break;

                case ".psmtsv":
                    res = GetCorrespondingRawDataAndResultsEntry(filepath);
                    if (res != null)
                        res.AddTSV(filepath);
                    else
                        rawDataAndResultslist.Add(new RawDataAndResults(null, null, filepath));
                    break;

                case ".mzid":
                    res = GetCorrespondingRawDataAndResultsEntry(filepath);
                    if (res != null)
                        res.AddMZID(filepath);
                    else
                        rawDataAndResultslist.Add(new RawDataAndResults(null, filepath, null));
                    break;

                case ".xml":
                    xMLdblist.Add(new XMLdb(filepath));
                    break;
            }
            dataGridDatafilesAndResults.Items.Refresh();
        }

        private RawDataAndResults GetCorrespondingRawDataAndResultsEntry(string filepath)
        {
            var fileNameNoExtension = Path.GetFileNameWithoutExtension(filepath);
            foreach (var a in rawDataAndResultslist)
            {
                if (a.FileName != null)
                {
                    var aNoExtension = Path.GetFileNameWithoutExtension(a.FileName);
                    if (aNoExtension.Equals(fileNameNoExtension))
                        return a;
                }
                if (a.mzidName != null)
                {
                    var aNoExtension = Path.GetFileNameWithoutExtension(a.mzidName);
                    if (aNoExtension.Equals(fileNameNoExtension))
                        return a;
                }
                if (a.psmsTSVName != null)
                {
                    var aNoExtension = Path.GetFileNameWithoutExtension(a.psmsTSVName);
                    if (aNoExtension.Equals(fileNameNoExtension))
                        return a;
                }
            }
            return null;
        }

        private void ClearRaw_Click(object sender, RoutedEventArgs e)
        {
            rawDataAndResultslist.Clear();
        }

        private void ButtonCalibrate_Click(object sender, RoutedEventArgs e)
        {
            string errorMessage = VerifyState(Task.Calibrate);
            if (errorMessage != null)
            {
                ErrorOutput(errorMessage);
                return;
            }

            var t = new Thread(() => DoTask(Task.Calibrate, po, null));
            t.IsBackground = true;
            t.Start();
        }

        private void ButtonSearch_Click(object sender, RoutedEventArgs e)
        {
            string errorMessage = VerifyState(Task.Search);
            if (errorMessage != null)
            {
                ErrorOutput(errorMessage);
                return;
            }

            List<double> accepted_precursor_mass_errors = new List<double>();
            foreach (string accepted_precursor_mass_error_text in textBox3.Text.Split(';'))
            {
                double accepted_precursor_mass_error;
                double.TryParse(accepted_precursor_mass_error_text, NumberStyles.Float | NumberStyles.AllowThousands, CultureInfo.InvariantCulture, out accepted_precursor_mass_error);
                accepted_precursor_mass_errors.Add(accepted_precursor_mass_error);
            }

            SearchParamsObject spo = new SearchParamsObject(
                checkBoxDecoy.IsChecked == true ? true : false,
                (Protease)comboBox.SelectedItem,
                Convert.ToInt32(missedCleavagesTextBox.Text), (InitiatorMethionineBehavior)Enum.Parse(typeof(InitiatorMethionineBehavior), cboInitiatorMethionineBehavior.Text, true),
                new MassTolerance(Convert.ToDouble(textBox2.Text), (MassToleranceUnits)comboBox1.SelectedIndex),
                accepted_precursor_mass_errors,
                new MassTolerance(Convert.ToDouble(textBox2_Copy.Text), (MassToleranceUnits)comboBox1_Copy.SelectedIndex),
                bCheckBox.IsChecked == true ? true : false,
                yCheckBox.IsChecked == true ? true : false,
                Convert.ToInt32(maxModificationIsoformsTextBox.Text),
                Convert.ToInt32(maxModificationsTextBox.Text)
                );

            var t = new Thread(() => DoTask(Task.Search, po, spo));
            t.IsBackground = true;
            t.Start();
        }

        private static void DoTask(Task task, ParamsObject po, object customObject)
        {
            po.startingTask();
            po.statusLabel("Working...");
            po.ReportProgress(0);
            switch (task)
            {
                case Task.Search:
                    po.RTBoutput("Starting database search");
                    SearchParamsObject paramsHere = (SearchParamsObject)customObject;
                    bool assign_charge_states = true;
                    bool deisotope = false;
                    bool on_the_fly_decoys = paramsHere.lookAtDecoys;
                    Protease protease = paramsHere.protease;
                    int max_missed_cleavages = paramsHere.missedCleavages;
                    InitiatorMethionineBehavior initiator_methionine_behavior = paramsHere.initiatorMethionineBehavior;
                    int max_variable_mod_isoforms = paramsHere.max_variable_mod_isoforms;
                    int max_mods_for_peptide = paramsHere.max_mods_for_peptide;
                    int min_assumed_precursor_charge_state = 2;
                    int max_assumed_precursor_charge_state = 4;
                    int max_peaks = 400;
                    MassTolerance precursor_mass_tolerance = paramsHere.precursorMassTolerance;

                    MassTolerance product_mass_tolerance = paramsHere.productMassTolerance;
                    double max_false_discovery_rate = 0.01;
                    List<string> filesToSearch = rawDataAndResultslist.Where(b => b.Use).Select(b => b.FileName).ToList();
                    string output_folder = Path.Combine(Path.GetDirectoryName(filesToSearch[0]), DateTime.Now.ToString("yyyy-MM-dd-HH-mm-ss", CultureInfo.InvariantCulture));

                    po.RTBoutput("output folder " + output_folder);
                    if (!Directory.Exists(output_folder))
                        Directory.CreateDirectory(output_folder);

                    Dictionary<string, List<MorpheusModification>> modsToLocalize = new Dictionary<string, List<MorpheusModification>>();
                    var modsInXML = ProteomeDatabaseReader.ReadXMLmodifications(xMLdblist.Where(b => b.Use).Select(b => b.FileName));
                    var modsInXMLtoTrim = ProteomeDatabaseReader.ReadXMLmodifications(xMLdblist.Where(b => b.Use).Select(b => b.FileName));
                    var modsKnown = SearchfileList.Where(b => b.Localize).SelectMany(b => b.getMods()).ToList();
                    foreach (var knownMod in modsKnown)
                        if (modsInXML.Contains(knownMod.NameInXML))
                        {
                            if (modsToLocalize.ContainsKey(knownMod.NameInXML))
                                modsToLocalize[knownMod.NameInXML].Add(knownMod);
                            else
                                modsToLocalize.Add(knownMod.NameInXML, new List<MorpheusModification>() { knownMod });
                            modsInXMLtoTrim.Remove(knownMod.NameInXML);
                        }
                    foreach (var ok in modsInXMLtoTrim)
                        modsToLocalize.Add(ok, new List<MorpheusModification>() { new MorpheusModification(ok) });

                    po.RTBoutput("The following modifications are in the XML databases, but not in the selected Localize PTM lists:");
                    po.RTBoutput(string.Join("\n", modsInXMLtoTrim));

                    po.RTBoutput("The following modifications are VARIABLE, i.e. we attempt to place them at every appropriate residue:");
                    po.RTBoutput(string.Join(", ", SearchfileList.Where(b => b.Variable).SelectMany(b => b.getMods()).Select(b => b.NameInXML)));

                    po.RTBoutput("The following modifications are FIXED, i.e. we assign them at every appropriate residue:");
                    po.RTBoutput(string.Join(", ", SearchfileList.Where(b => b.Fixed).SelectMany(b => b.getMods()).Select(b => b.NameInXML)));

                    string extraLogStuff = "";
                    extraLogStuff += "PTM files: " + string.Join(", ", SearchfileList.Where(b => b.Localize).Select(b => b.FileName)) + "\n";
                    extraLogStuff += "XML databases: " + string.Join(", ", xMLdblist.Where(b => b.Use).Select(b => b.FileName)) + "\n";
                    extraLogStuff += "Localized mods: " + string.Join(", ", modsToLocalize.Select(b => b.Key)) + "\n";
                    extraLogStuff += "Variable mods: " + string.Join(", ", SearchfileList.Where(b => b.Variable).SelectMany(b => b.getMods()).Select(b => b.NameInXML)) + "\n";
                    extraLogStuff += "Fixed mods: " + string.Join(", ", SearchfileList.Where(b => b.Fixed).SelectMany(b => b.getMods()).Select(b => b.NameInXML)) + "\n";

                    DatabaseSearcher database_searcher = new DatabaseSearcher(GetDatas(filesToSearch),
                        min_assumed_precursor_charge_state, max_assumed_precursor_charge_state,
                        max_peaks,
                        assign_charge_states, deisotope, on_the_fly_decoys,
                        protease, max_missed_cleavages, initiator_methionine_behavior,
                        max_variable_mod_isoforms,
                        max_mods_for_peptide,
                        precursor_mass_tolerance,
                        paramsHere.accepted_precursor_mass_errors,
                        product_mass_tolerance,
                        max_false_discovery_rate,
                        output_folder,
                        paramsHere.bions,
                        paramsHere.yions,
                        xMLdblist.Where(b => b.Use).SelectMany(b => b.getProteins(on_the_fly_decoys, modsToLocalize)).ToList(),
                        SearchfileList.Where(b => b.Variable).SelectMany(b => b.getMods()).ToList(),
                        SearchfileList.Where(b => b.Fixed).SelectMany(b => b.getMods()).ToList(),
                        extraLogStuff,
                        filesToSearch.Count);
                    po.AttachoutRichTextBoxHandlerr(handler => database_searcher.outputHandler += handler);
                    po.AttachoSuccessfullyFinishedFileHandler(handler => database_searcher.finishedFileHandler += handler);
                    po.AttachoutProgressBarHandler(handler => database_searcher.progressHandler += handler);
                    po.AttachoutLabelStatusHandler(handler => database_searcher.labelStatusHandler += handler);
                    database_searcher.DoSearch();
                    break;

                case Task.Calibrate:
                    po.RTBoutput("Starting calibrate task");
                    foreach (var anEntry in rawDataAndResultslist.Where(b => b.FileName != null && b.Use == true && b.mzidName != null && b.UseMzid == true && b.UsePsmsTsv == false).ToList())
                    {
                        SoftwareLockMassParams a = mzCalIO.mzCalIO.GetReady(anEntry.FileName, anEntry.mzidName);
                        po.AttachoutRichTextBoxHandlerr(handler => a.outputHandler += handler);
                        po.AttachoSuccessfullyFinishedFileHandler(handler => a.finishedFileHandler += handler);
                        po.AttachoutProgressBarHandler(handler => a.progressHandler += handler);
                        po.AttachoutRichTextBoxHandlerr(handler => a.watchHandler += handler);
                        SoftwareLockMassRunner.Run(a);
                    }
                    break;

                case Task.GPTMD:
                    po.RTBoutput("Starting GPTMD");
                    GPTMDParamsObject paramsHeree = (GPTMDParamsObject)customObject;

                    IEnumerable<Tuple<double, double>> combos = ReadCombos(@"combos.txt");

                    modsToLocalize = new Dictionary<string, List<MorpheusModification>>();
                    modsInXML = ProteomeDatabaseReader.ReadXMLmodifications(xMLdblist.Where(b => b.Use).Select(b => b.FileName));
                    modsInXMLtoTrim = ProteomeDatabaseReader.ReadXMLmodifications(xMLdblist.Where(b => b.Use).Select(b => b.FileName));
                    var gptmdmodsKnown = GPTMDfileList.Where(b => b.Use).SelectMany(b => b.getMods()).ToList();
                    foreach (var knownMod in gptmdmodsKnown)
                        if (modsInXML.Contains(knownMod.NameInXML))
                        {
                            if (modsToLocalize.ContainsKey(knownMod.NameInXML))
                                modsToLocalize[knownMod.NameInXML].Add(knownMod);
                            else
                                modsToLocalize.Add(knownMod.NameInXML, new List<MorpheusModification>() { knownMod });
                            modsInXMLtoTrim.Remove(knownMod.NameInXML);
                        }

                    foreach (var ok in modsInXMLtoTrim)
                        modsToLocalize.Add(ok, new List<MorpheusModification>() { new MorpheusModification(ok) });

                    var okk = xMLdblist.First(b => b.Use);
                    var outputFileName = Path.Combine(Path.GetDirectoryName(okk.FileName), Path.GetFileNameWithoutExtension(okk.FileName) + "-GPTMD-" + DateTime.Now.ToString("yyyy-MM-dd-HH-mm", CultureInfo.InvariantCulture) + ".xml");
                    GPTMD gptmd = new GPTMD();
                    po.AttachoutRichTextBoxHandlerr(handler => gptmd.outputHandler += handler);
                    po.AttachoutProgressBarHandler(handler => gptmd.progressHandler += handler);
                    po.AttachoutLabelStatusHandler(handler => gptmd.labelStatusHandler += handler);

                    gptmd.GPTMDD(
                        rawDataAndResultslist.First(b => b.UsePsmsTsv),
                        xMLdblist.Where(b => b.Use).SelectMany(b => b.getProteins(false, modsToLocalize)).ToList(),
                        combos.ToList(),
                        GPTMDfileList.Where(b => b.Use).SelectMany(b => b.getMods()).ToList(),
                        paramsHeree.monoisotopic,
                        outputFileName);
                    po.FinishedFile(outputFileName);
                    break;
            }
            po.ReportProgress(100);
            po.RTBoutput("SUCCESS!");
            po.statusLabel("Ready");
            po.FinishedTask();
        }

        private static IEnumerable<Tuple<double, double>> ReadCombos(string v)
        {
            var lines = File.ReadLines(v);
            foreach (var line in lines)
            {
                var ye = line.Split(' ');
                yield return new Tuple<double, double>(Convert.ToDouble(ye[0]), Convert.ToDouble(ye[1]));
            }
        }

        private static IEnumerable<TandemMassSpectra> GetDatas(List<string> filesToSearch)
        {
            // convert all paths to absolute for outputs
            for (int i = 0; i < filesToSearch.Count; i++)
            {
                //Console.WriteLine(Path.GetExtension(myListOfEntries[i].filepath));
                if (Path.GetExtension(filesToSearch[i]).Equals(".mzML"))
                {
                    //var ok = new MzMLTandemMassSpectra();
                    //ok.filename = filesToSearch[i];
                    //yield return ok;
                }
                else if (Path.GetExtension(filesToSearch[i]).Equals(".raw"))
                {
                    //var ok = new ThermoTandemMassSpectra();
                    //ok.filename = filesToSearch[i];
                    //yield return ok;
                }
            }
            return null;
        }

        private string VerifyState(Task task)
        {
            RegOutput("Verifying state...");
            switch (task)
            {
                case Task.Search:
                    if (xMLdblist.Where(b => b.Use).Count() == 0)
                        return "No XML database selected";
                    foreach (var ye in rawDataAndResultslist)
                    {
                        if (ye.UseMzid || ye.UsePsmsTsv)
                            return "Not searching " + ye.Use + "because corresponding result is checked";
                        if (ye.Use && (ye.mzidName != null || ye.psmsTSVName != null))
                            return "Not searching " + ye.Use + "because there are existing results";
                    }
                    if (rawDataAndResultslist.Where(b => b.Use).Count() == 0)
                        return "No files selected for search";
                    break;

                case Task.Calibrate:
                    foreach (var ye in rawDataAndResultslist)
                    {
                        if (ye.Use == true)
                        {
                            // So, seems that want to calibrate this one...
                            if (ye.UseMzid == false)
                                return "Please select mzid file for " + ye.FileName;
                            if (ye.UsePsmsTsv == true)
                                return "Please unselect tsv file for " + ye.FileName;
                            if (ye.FileName == null)
                                return "Filename is null";
                            if (ye.mzidName == null)
                                return "mzid is null";
                        }
                        if (ye.UseMzid == true)
                        {
                            // So, seems that want to calibrate this one...
                            if (ye.Use == false)
                                return "Please select spectra file for " + ye.FileName;
                            if (ye.UsePsmsTsv == true)
                                return "Please unselect tsv file for " + ye.FileName;
                            if (ye.FileName == null)
                                return "Filename is null";
                            if (ye.mzidName == null)
                                return "mzid is null";
                        }
                    }
                    break;

                case Task.GPTMD:
                    if (rawDataAndResultslist.Where(b => b.Use).Count() > 0
                        || rawDataAndResultslist.Where(b => b.UseMzid).Count() > 0)
                        return "Please uncheck spectra or mzid files";
                    if (rawDataAndResultslist.Where(b => b.UsePsmsTsv).Count() != 1
                        || rawDataAndResultslist.Where(b => b.UseMzid).Count() > 0)
                        return "GPTMD works on a single tsv file";
                    if (xMLdblist.Where(b => b.Use).Count() == 0
                        || rawDataAndResultslist.Where(b => b.UseMzid).Count() > 0)
                        return "Please select XML database file(s)";
                    break;
            }
            RegOutput("Verified successfully");
            return null;
        }

        private void ErrorOutput(string v)
        {
            if (outRichTextBox.Document.Blocks.LastBlock.Foreground != Brushes.Red)
            {
                Paragraph p = new Paragraph(new Run(v + "\n"));
                p.Margin = new Thickness(0);
                p.Foreground = Brushes.Red;
                p.FontSize = 20;
                outRichTextBox.Document.Blocks.Add(p);
            }
            else
                outRichTextBox.AppendText(v + "\n");
            outRichTextBox.ScrollToEnd();
        }

        private void RegOutput(string v)
        {
            if (outRichTextBox.Document.Blocks.LastBlock.Foreground != Brushes.Black)
            {
                Paragraph p = new Paragraph(new Run(v + "\n"));
                p.Margin = new Thickness(0);
                p.Foreground = Brushes.Black;
                p.FontSize = 12;
                outRichTextBox.Document.Blocks.Add(p);
            }
            else
                outRichTextBox.AppendText(v + "\n");
            outRichTextBox.ScrollToEnd();
        }

        private void ButtonGPTMD_Click(object sender, RoutedEventArgs e)
        {
            string errorMessage = VerifyState(Task.GPTMD);
            if (errorMessage != null)
            {
                ErrorOutput(errorMessage);
                return;
            }
            GPTMDParamsObject gpo = new GPTMDParamsObject(checkBoxMonoisotopic.IsChecked.Value);
            var t = new Thread(() => DoTask(Task.GPTMD, po, gpo));
            t.IsBackground = true;
            t.Start();
        }

        private void AddPTMlist_Click(object sender, RoutedEventArgs e)
        {
            // Create the OpenFIleDialog object
            Microsoft.Win32.OpenFileDialog openPicker = new Microsoft.Win32.OpenFileDialog();
            openPicker.Filter = "TXT Files|*.txt";
            openPicker.FilterIndex = 1;
            openPicker.RestoreDirectory = true;
            if (openPicker.ShowDialog() == true)
            {
                addFile(openPicker.FileName);
            }
        }

        private void AddXML_Click(object sender, RoutedEventArgs e)
        {
            // Create the OpenFIleDialog object
            Microsoft.Win32.OpenFileDialog openPicker = new Microsoft.Win32.OpenFileDialog();
            openPicker.Filter = "XML Files|*.xml";
            openPicker.FilterIndex = 1;
            openPicker.RestoreDirectory = true;
            if (openPicker.ShowDialog() == true)
            {
                addFile(openPicker.FileName);
            }
        }

        private void AddRaw_Click(object sender, RoutedEventArgs e)
        {
            Microsoft.Win32.OpenFileDialog openFileDialog1 = new Microsoft.Win32.OpenFileDialog();

            openFileDialog1.Filter = "Spectra Files(*.raw;*.mzML;*.mzid;*.psmtsv)|*.raw;*.mzML;*.mzid;*.psmtsv";
            openFileDialog1.FilterIndex = 1;
            openFileDialog1.RestoreDirectory = true;
            openFileDialog1.Multiselect = true;

            if (openFileDialog1.ShowDialog() == true)
                foreach (var filepath in openFileDialog1.FileNames)
                    addFile(filepath);
        }

        private void Window_Drop(object sender, DragEventArgs e)
        {
            string[] files = (string[])e.Data.GetData(DataFormats.FileDrop);
            foreach (var file in files)
            {
                addFile(file);
            }
        }

        private void RemoveRaw_Click(object sender, RoutedEventArgs e)
        {
        }

        private void button_Click(object sender, RoutedEventArgs e)
        {
            System.Diagnostics.Process.Start(@"help.txt");
        }

        private void Row_DoubleClick(object sender, MouseButtonEventArgs e)
        {
            var ye = sender as DataGridCell;
            var hm = ye.Content as TextBlock;
            if (hm != null && !hm.Text.Equals(""))
            {
                System.Diagnostics.Process.Start(hm.Text);
            }
        }
    }
}