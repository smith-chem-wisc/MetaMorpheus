using InternalLogicEngineLayer;
using InternalLogicTaskLayer;
using OldInternalLogic;
using Spectra;
using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.IO;
using System.Linq;
using System.Threading;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;

namespace MetaMorpheusGUI
{
    /// <summary>
    /// Interaction logic for MainWindow.xaml
    /// </summary>
    public partial class MainWindow : Window
    {
        private static ObservableCollection<RawData> rawDataObservableCollection = new ObservableCollection<RawData>();
        private static ObservableCollection<XMLdb> xmlDBobservableCollection = new ObservableCollection<XMLdb>();
        private static ObservableCollection<ModList> modListObservableCollection = new ObservableCollection<ModList>();
        private static ObservableCollection<SearchMode> searchModeObservableCollection = new ObservableCollection<SearchMode>();
        private static ObservableCollection<FinishedFile> finishedFileObservableCollection = new ObservableCollection<FinishedFile>();
        private static ObservableCollection<MyTaskEngine> taskEngineObservableCollection = new ObservableCollection<MyTaskEngine>();

        public static string elementsLocation = @"elements.dat";

        public MainWindow()
        {
            InitializeComponent();

            UsefulProteomicsDatabases.Loaders.LoadElements(elementsLocation);

            // modificationsDataGrid.DataContext = ModFileList;
            dataGridXMLs.DataContext = xmlDBobservableCollection;
            dataGridDatafiles.DataContext = rawDataObservableCollection;
            tasksDataGrid.DataContext = taskEngineObservableCollection;
            outputFilesDataGrid.DataContext = finishedFileObservableCollection;

            modListObservableCollection.Add(new ModList("f.txt"));
            modListObservableCollection.Add(new ModList("v.txt"));
            modListObservableCollection.Add(new ModList("p.txt"));
            modListObservableCollection.Add(new ModList("m.txt"));
            modListObservableCollection.Add(new ModList("r.txt"));
            modListObservableCollection.Add(new ModList("s.txt"));

            LoadSearchModesFromFile();

            xmlDBobservableCollection.Add(new XMLdb(@"C:\Users\stepa\Data\CalibrationPaperData\OrigData\uniprot-mouse-reviewed-12-23-2016.xml"));

            rawDataObservableCollection.Add(new RawData(@"C:\Users\stepa\Data\CalibrationPaperData\Step2\Mouse\Calib-0.1.2\04-29-13_B6_Frac9_9p5uL-Calibrated.mzML"));

            RegOutput(ProteaseDictionary.Instance.Count + " proteases loaded from proteases.tsv");
            AminoAcidMasses.LoadAminoAcidMasses();
            RegOutput("Amino acid masses loaded from amino_acids.tsv");

            //po.newDbsHandler += AddNewDB;
            //po.newSpectrasHandler += AddNewSpectra;

            AllTasksEngine.startingAllTasksEngineHandler += NewSuccessfullyStartingAllTasks;
            AllTasksEngine.finishedAllTasksEngineHandler += NewSuccessfullyFinishedAllTasks;

            MyTaskEngine.startingSingleTaskHander += Po_startingSingleTaskHander;
            MyTaskEngine.finishedSingleTaskHandler += Po_finishedSingleTaskHandler;
            MyTaskEngine.finishedWritingFileHandler += NewSuccessfullyFinishedWritingFile;

            MyEngine.outProgressHandler += NewoutProgressBar;

            MyEngine.outLabelStatusHandler += NewoutLabelStatus;
            MyEngine.outRichTextBoxHandler += NewoutRichTextBox;

            UpdateTaskGuiStuff();
        }

        private void LoadSearchModesFromFile()
        {
            searchModeObservableCollection.Add(new DotSearchMode("5ppm", new double[] { 0 }, new Tolerance(ToleranceUnit.PPM, 5)));
            searchModeObservableCollection.Add(new IntervalSearchMode("twoPointOneDalton", new List<DoubleRange>() { new DoubleRange(-2.1, 2.1) }));
        }

        private void AddNewDB(object sender, List<string> e)
        {
            if (!Dispatcher.CheckAccess())
            {
                Dispatcher.BeginInvoke(new Action(() => AddNewDB(sender, e)));
            }
            else
            {
                foreach (var uu in xmlDBobservableCollection)
                    uu.Use = false;
                foreach (var uu in e)
                    xmlDBobservableCollection.Add(new XMLdb(uu));
            }
        }

        private void AddNewSpectra(object sender, List<string> e)
        {
            if (!Dispatcher.CheckAccess())
            {
                Dispatcher.BeginInvoke(new Action(() => AddNewSpectra(sender, e)));
            }
            else
            {
                foreach (var uu in rawDataObservableCollection)
                    uu.Use = false;
                foreach (var uu in e)
                    rawDataObservableCollection.Add(new RawData(uu));
            }
        }

        private void Po_startingSingleTaskHander(object sender, SingleTaskEventArgs s)
        {
            if (!Dispatcher.CheckAccess())
            {
                Dispatcher.BeginInvoke(new Action(() => Po_startingSingleTaskHander(sender, s)));
            }
            else
            {
                s.theTask.IsMySelected = true;
                statusLabel.Content = "Running " + s.theTask.taskType + " task";
                outProgressBar.IsIndeterminate = true;

                // Update highlighting

                tasksDataGrid.Items.Refresh();
                dataGridDatafiles.Items.Refresh();
                dataGridXMLs.Items.Refresh();
            }
        }

        private void Po_finishedSingleTaskHandler(object sender, SingleTaskEventArgs s)
        {
            if (!Dispatcher.CheckAccess())
            {
                Dispatcher.BeginInvoke(new Action(() => Po_finishedSingleTaskHandler(sender, s)));
            }
            else
            {
                s.theTask.IsMySelected = false;
                statusLabel.Content = "Finished " + s.theTask.taskType + " task";
                outProgressBar.Value = 100;

                // Update highlighting and possibly add to the datafiles/xml list
                tasksDataGrid.Items.Refresh();
                dataGridDatafiles.Items.Refresh();
                dataGridXMLs.Items.Refresh();
            }
        }

        //private void TasksDataGrid_SelectionChanged(object sender, SelectionChangedEventArgs e)
        //{
        //    // Need to decide what happened:
        //    // If user clicked on one existing one, do something

        //    DataGrid datagrid = sender as DataGrid; // e.Source could have been used instead of sender as well
        //    if (datagrid.SelectedIndex >= 0 && datagrid.SelectedIndex != highlightedTaskIndex)
        //    {
        //        ((MyTask)datagrid.Items[highlightedTaskIndex]).IsMySelected = false;
        //        highlightedTaskIndex = datagrid.SelectedIndex;
        //        ((MyTask)datagrid.Items[highlightedTaskIndex]).IsMySelected = true;
        //        datagrid.Items.Refresh();
        //        //UpdateTaskStuff();
        //    }
        //}

        //private void UpdateTaskStuff()
        //{
        //    if (tasksDataGrid.Items.Count > 0)
        //    {
        //        var selectedTask = tasksDataGrid.Items[highlightedTaskIndex] as MyTask;
        //        if (selectedTask != null)
        //        {
        //            switch (selectedTask.taskType)
        //            {
        //                case MyTaskEnum.Calibrate:
        //                    taskTabControl.SelectedIndex = 0;
        //                    ((Control)taskTabControl.Items[0]).IsEnabled = true;
        //                    ((Control)taskTabControl.Items[1]).IsEnabled = false;
        //                    ((Control)taskTabControl.Items[2]).IsEnabled = false;
        //                    addTaskButton.IsEnabled = false;
        //                    RemoveTask.IsEnabled = true;
        //                    break;
        //                case MyTaskEnum.Search:
        //                    taskTabControl.SelectedIndex = 1;
        //                    ((Control)taskTabControl.Items[0]).IsEnabled = false;
        //                    ((Control)taskTabControl.Items[1]).IsEnabled = true;
        //                    ((Control)taskTabControl.Items[2]).IsEnabled = false;
        //                    addTaskButton.IsEnabled = false;
        //                    RemoveTask.IsEnabled = true;
        //                    break;
        //                case MyTaskEnum.GPTMD:
        //                    taskTabControl.SelectedIndex = 2;
        //                    ((Control)taskTabControl.Items[0]).IsEnabled = false;
        //                    ((Control)taskTabControl.Items[1]).IsEnabled = false;
        //                    ((Control)taskTabControl.Items[2]).IsEnabled = true;
        //                    addTaskButton.IsEnabled = false;
        //                    RemoveTask.IsEnabled = true;
        //                    break;
        //                case MyTaskEnum.NewTask:
        //                    ((Control)taskTabControl.Items[0]).IsEnabled = true;
        //                    ((Control)taskTabControl.Items[1]).IsEnabled = true;
        //                    ((Control)taskTabControl.Items[2]).IsEnabled = true;
        //                    addTaskButton.IsEnabled = true;
        //                    RemoveTask.IsEnabled = false;
        //                    break;
        //            }
        //            SetTextBoxValues(selectedTask);
        //        }
        //    }
        //}

        private void addFinishedFile(string filepath)
        {
            finishedFileObservableCollection.Add(new FinishedFile(filepath));
            outputFilesDataGrid.Items.Refresh();
        }

        private RawData GetCorrespondingRawDataAndResultsEntry(string filepath)
        {
            var fileNameNoExtension = Path.GetFileNameWithoutExtension(filepath);
            foreach (var a in rawDataObservableCollection)
            {
                if (a.FileName != null)
                {
                    var aNoExtension = Path.GetFileNameWithoutExtension(a.FileName);
                    if (aNoExtension.Equals(fileNameNoExtension))
                        return a;
                }
                //if (a.mzidName != null)
                //{
                //    var aNoExtension = Path.GetFileNameWithoutExtension(a.mzidName);
                //    if (aNoExtension.Equals(fileNameNoExtension))
                //        return a;
                //}
                //if (a.psmsTSVName != null)
                //{
                //    var aNoExtension = Path.GetFileNameWithoutExtension(a.psmsTSVName);
                //    if (aNoExtension.Equals(fileNameNoExtension))
                //        return a;
                //}
            }
            return null;
        }

        private void ClearRaw_Click(object sender, RoutedEventArgs e)
        {
            rawDataObservableCollection.Clear();
        }

        //private void ButtonCalibrate_Click(object sender, RoutedEventArgs e)
        //{
        //    if (errorMessage != null)
        //    {
        //        ErrorOutput(errorMessage);
        //        return;
        //    }

        //    var t = new Thread(() => DoTask(MyTaskEnum.Calibrate, po, null));
        //    t.IsBackground = true;
        //    t.Start();
        //}

        //private void ButtonSearch_Click(object sender, RoutedEventArgs e)
        //{
        //    if (errorMessage != null)
        //    {
        //        ErrorOutput(errorMessage);
        //        return;
        //    }

        //    List<double> accepted_precursor_mass_errors = new List<double>();
        //    foreach (string accepted_precursor_mass_error_text in textBox3.Text.Split(';'))
        //    {
        //        double accepted_precursor_mass_error;
        //        double.TryParse(accepted_precursor_mass_error_text, NumberStyles.Float | NumberStyles.AllowThousands, CultureInfo.InvariantCulture, out accepted_precursor_mass_error);
        //        accepted_precursor_mass_errors.Add(accepted_precursor_mass_error);
        //    }

        //    //SearchParamsObject spo = new SearchParamsObject(
        //    //    checkBoxDecoy.IsChecked == true ? true : false,
        //    //    (Protease)comboBox.SelectedItem,
        //    //    Convert.ToInt32(missedCleavagesTextBox.Text), (InitiatorMethionineBehavior)Enum.Parse(typeof(InitiatorMethionineBehavior), cboInitiatorMethionineBehavior.Text, true),
        //    //    new MassTolerance(Convert.ToDouble(textBox2.Text), (MassToleranceUnits)comboBox1.SelectedIndex),
        //    //    accepted_precursor_mass_errors,
        //    //    new MassTolerance(Convert.ToDouble(textBox2_Copy.Text), (MassToleranceUnits)comboBox1_Copy.SelectedIndex),
        //    //    bCheckBox.IsChecked == true ? true : false,
        //    //    yCheckBox.IsChecked == true ? true : false,
        //    //    Convert.ToInt32(maxModificationIsoformsTextBox.Text)
        //    //    );

        //    //var t = new Thread(() => DoTask(Task.Search, po, spo));
        //    //t.IsBackground = true;
        //    //t.Start();
        //}

        //private static void DoTask(MyTaskEnum task, ParamsObject po, object customObject)
        //{
        //    po.startingAllTasks();
        //    po.statusLabel("Working...");
        //    po.ReportProgress(0);
        //    switch (task)
        //    {
        //        case MyTaskEnum.Search:
        //            po.RTBoutput("Starting database search");
        //            SearchParamsObject paramsHere = (SearchParamsObject)customObject;
        //            bool assign_charge_states = true;
        //            bool deisotope = false;
        //            bool on_the_fly_decoys = paramsHere.lookAtDecoys;
        //            Protease protease = paramsHere.protease;
        //            int max_missed_cleavages = paramsHere.missedCleavages;
        //            InitiatorMethionineBehavior initiator_methionine_behavior = paramsHere.initiatorMethionineBehavior;
        //            int max_variable_mod_isoforms = paramsHere.max_variable_mod_isoforms;
        //            int max_mods_for_peptide = paramsHere.max_mods_for_peptide;
        //            int min_assumed_precursor_charge_state = 2;
        //            int max_assumed_precursor_charge_state = 4;
        //            int max_peaks = 400;
        //            MassTolerance precursor_mass_tolerance = paramsHere.precursorMassTolerance;

        //            MassTolerance product_mass_tolerance = paramsHere.productMassTolerance;
        //            double max_false_discovery_rate = 0.01;
        //            List<string> filesToSearch = rawDataAndResultslist.Where(b => b.Use).Select(b => b.FileName).ToList();
        //            string output_folder = Path.Combine(Path.GetDirectoryName(filesToSearch[0]), DateTime.Now.ToString("yyyy-MM-dd-HH-mm-ss", CultureInfo.InvariantCulture));

        //            po.RTBoutput("output folder " + output_folder);
        //            if (!Directory.Exists(output_folder))
        //                Directory.CreateDirectory(output_folder);

        //            Dictionary<string, List<MorpheusModification>> modsToLocalize = new Dictionary<string, List<MorpheusModification>>();
        //            var modsInXML = ProteomeDatabaseReader.ReadXMLmodifications(xMLdblist.Where(b => b.Use).Select(b => b.FileName));
        //            var modsInXMLtoTrim = ProteomeDatabaseReader.ReadXMLmodifications(xMLdblist.Where(b => b.Use).Select(b => b.FileName));
        //            var modsKnown = ModFileList.Where(b => b.Localize).SelectMany(b => b.getMods()).ToList();
        //            foreach (var knownMod in modsKnown)
        //                if (modsInXML.Contains(knownMod.NameInXML))
        //                {
        //                    if (modsToLocalize.ContainsKey(knownMod.NameInXML))
        //                        modsToLocalize[knownMod.NameInXML].Add(knownMod);
        //                    else
        //                        modsToLocalize.Add(knownMod.NameInXML, new List<MorpheusModification>() { knownMod });
        //                    modsInXMLtoTrim.Remove(knownMod.NameInXML);
        //                }
        //            foreach (var ok in modsInXMLtoTrim)
        //                modsToLocalize.Add(ok, new List<MorpheusModification>() { new MorpheusModification(ok) });

        //            po.RTBoutput("The following modifications are in the XML databases, but not in the selected Localize PTM lists:");
        //            po.RTBoutput(string.Join("\n", modsInXMLtoTrim));

        //            po.RTBoutput("The following modifications are VARIABLE, i.e. we attempt to place them at every appropriate residue:");
        //            po.RTBoutput(string.Join(", ", ModFileList.Where(b => b.Variable).SelectMany(b => b.getMods()).Select(b => b.NameInXML)));

        //            po.RTBoutput("The following modifications are FIXED, i.e. we assign them at every appropriate residue:");
        //            po.RTBoutput(string.Join(", ", ModFileList.Where(b => b.Fixed).SelectMany(b => b.getMods()).Select(b => b.NameInXML)));

        //            string extraLogStuff = "";
        //            extraLogStuff += "PTM files: " + string.Join(", ", ModFileList.Where(b => b.Localize).Select(b => b.FileName)) + "\n";
        //            extraLogStuff += "XML databases: " + string.Join(", ", xMLdblist.Where(b => b.Use).Select(b => b.FileName)) + "\n";
        //            extraLogStuff += "Localized mods: " + string.Join(", ", modsToLocalize.Select(b => b.Key)) + "\n";
        //            extraLogStuff += "Variable mods: " + string.Join(", ", ModFileList.Where(b => b.Variable).SelectMany(b => b.getMods()).Select(b => b.NameInXML)) + "\n";
        //            extraLogStuff += "Fixed mods: " + string.Join(", ", ModFileList.Where(b => b.Fixed).SelectMany(b => b.getMods()).Select(b => b.NameInXML)) + "\n";

        //            DatabaseSearcher database_searcher = new DatabaseSearcher(GetDatas(filesToSearch),
        //                min_assumed_precursor_charge_state, max_assumed_precursor_charge_state,
        //                max_peaks,
        //                assign_charge_states, deisotope, on_the_fly_decoys,
        //                protease, max_missed_cleavages, initiator_methionine_behavior,
        //                max_variable_mod_isoforms,
        //                max_mods_for_peptide,
        //                precursor_mass_tolerance,
        //                paramsHere.accepted_precursor_mass_errors,
        //                product_mass_tolerance,
        //                max_false_discovery_rate,
        //                output_folder,
        //                paramsHere.bions,
        //                paramsHere.yions,
        //                xMLdblist.Where(b => b.Use).SelectMany(b => b.getProteins(on_the_fly_decoys, modsToLocalize)).ToList(),
        //                ModFileList.Where(b => b.Variable).SelectMany(b => b.getMods()).ToList(),
        //                ModFileList.Where(b => b.Fixed).SelectMany(b => b.getMods()).ToList(),
        //                extraLogStuff,
        //                filesToSearch.Count);
        //            po.AttachoutRichTextBoxHandlerr(handler => database_searcher.outputHandler += handler);
        //            po.AttachoSuccessfullyFinishedFileHandler(handler => database_searcher.finishedFileHandler += handler);
        //            po.AttachoutProgressBarHandler(handler => database_searcher.progressHandler += handler);
        //            po.AttachoutLabelStatusHandler(handler => database_searcher.labelStatusHandler += handler);
        //            database_searcher.DoSearch();
        //            break;

        //        case MyTaskEnum.Calibrate:

        //            po.RTBoutput("Starting calibrate task");
        //            foreach (var anEntry in rawDataAndResultslist.Where(b => b.FileName != null && b.Use == true).ToList())
        //            {
        //                //ClassicSearchParams searchParams = new ClassicSearchParams(myMsDataFile, spectraFileIndex, variableModifications, fixedModifications, localizeableModifications, proteinList, fragmentTolerance, protease, searchModes[0]);
        //                //ClassicSearchEngine searchEngine = new ClassicSearchEngine(searchParams);
        //                //ClassicSearchResults searchResults = (ClassicSearchResults)searchEngine.Run();

        //                //SoftwareLockMassParams a = mzCalIO.mzCalIO.GetReady(anEntry.FileName, anEntry.mzidName);
        //                //po.AttachoutRichTextBoxHandlerr(handler => a.outputHandler += handler);
        //                //po.AttachoSuccessfullyFinishedFileHandler(handler => a.finishedFileHandler += handler);
        //                //po.AttachoutProgressBarHandler(handler => a.progressHandler += handler);
        //                //po.AttachoutRichTextBoxHandlerr(handler => a.watchHandler += handler);
        //                //SoftwareLockMassRunner.Run(a);
        //            }
        //            break;

        //        case MyTaskEnum.GPTMD:
        //            po.RTBoutput("Starting GPTMD");
        //            GPTMDParamsObject paramsHeree = (GPTMDParamsObject)customObject;

        //            IEnumerable<Tuple<double, double>> combos = ReadCombos(@"combos.txt");

        //            modsToLocalize = new Dictionary<string, List<MorpheusModification>>();
        //            modsInXML = ProteomeDatabaseReader.ReadXMLmodifications(xMLdblist.Where(b => b.Use).Select(b => b.FileName));
        //            modsInXMLtoTrim = ProteomeDatabaseReader.ReadXMLmodifications(xMLdblist.Where(b => b.Use).Select(b => b.FileName));
        //            var gptmdmodsKnown = GPTMDfileList.Where(b => b.Use).SelectMany(b => b.getMods()).ToList();
        //            foreach (var knownMod in gptmdmodsKnown)
        //                if (modsInXML.Contains(knownMod.NameInXML))
        //                {
        //                    if (modsToLocalize.ContainsKey(knownMod.NameInXML))
        //                        modsToLocalize[knownMod.NameInXML].Add(knownMod);
        //                    else
        //                        modsToLocalize.Add(knownMod.NameInXML, new List<MorpheusModification>() { knownMod });
        //                    modsInXMLtoTrim.Remove(knownMod.NameInXML);
        //                }

        //            foreach (var ok in modsInXMLtoTrim)
        //                modsToLocalize.Add(ok, new List<MorpheusModification>() { new MorpheusModification(ok) });

        //            var okk = xMLdblist.First(b => b.Use);
        //            var outputFileName = Path.Combine(Path.GetDirectoryName(okk.FileName), Path.GetFileNameWithoutExtension(okk.FileName) + "-GPTMD-" + DateTime.Now.ToString("yyyy-MM-dd-HH-mm", CultureInfo.InvariantCulture) + ".xml");
        //            GPTMD gptmd = new GPTMD();
        //            po.AttachoutRichTextBoxHandlerr(handler => gptmd.outputHandler += handler);
        //            po.AttachoutProgressBarHandler(handler => gptmd.progressHandler += handler);
        //            po.AttachoutLabelStatusHandler(handler => gptmd.labelStatusHandler += handler);

        //            gptmd.GPTMDD(
        //                rawDataAndResultslist.First(),
        //                xMLdblist.Where(b => b.Use).SelectMany(b => b.getProteins(false, modsToLocalize)).ToList(),
        //                combos.ToList(),
        //                GPTMDfileList.Where(b => b.Use).SelectMany(b => b.getMods()).ToList(),
        //                paramsHeree.monoisotopic,
        //                outputFileName);
        //            po.FinishedFile(outputFileName);
        //            break;
        //    }
        //    po.ReportProgress(100);
        //    po.RTBoutput("SUCCESS!");
        //    po.statusLabel("Ready");
        //    po.FinishedAllTasks();
        //}

        private static IEnumerable<Tuple<double, double>> ReadCombos(string v)
        {
            var lines = File.ReadLines(v);
            foreach (var line in lines)
            {
                var ye = line.Split(' ');
                yield return new Tuple<double, double>(Convert.ToDouble(ye[0]), Convert.ToDouble(ye[1]));
            }
        }

        //private static IEnumerable<TandemMassSpectra> GetDatas(List<string> filesToSearch)
        //{
        //    // convert all paths to absolute for outputs
        //    for (int i = 0; i < filesToSearch.Count; i++)
        //    {
        //        //p.po.RTBoutput(Path.GetExtension(myListOfEntries[i].filepath));
        //        if (Path.GetExtension(filesToSearch[i]).Equals(".mzML"))
        //        {
        //            //var ok = new MzMLTandemMassSpectra();
        //            //ok.filename = filesToSearch[i];
        //            //yield return ok;
        //        }
        //        else if (Path.GetExtension(filesToSearch[i]).Equals(".raw"))
        //        {
        //            //var ok = new ThermoTandemMassSpectra();
        //            //ok.filename = filesToSearch[i];
        //            //yield return ok;
        //        }
        //    }
        //    return null;
        //}

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
            if (outRichTextBox != null)
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
        }

        //private void ButtonGPTMD_Click(object sender, RoutedEventArgs e)
        //{
        //    GPTMDParamsObject gpo = new GPTMDParamsObject(checkBoxMonoisotopic.IsChecked.Value);
        //    var t = new Thread(() => DoTask(MyTaskEnum.GPTMD, po, gpo));
        //    t.IsBackground = true;
        //    t.Start();
        //}

        private void AddPTMlist_Click(object sender, RoutedEventArgs e)
        {
            // Create the OpenFIleDialog object
            Microsoft.Win32.OpenFileDialog openPicker = new Microsoft.Win32.OpenFileDialog();
            openPicker.Filter = "TXT Files|*.txt";
            openPicker.FilterIndex = 1;
            openPicker.RestoreDirectory = true;
            if (openPicker.ShowDialog() == true)
            {
                addFinishedFile(openPicker.FileName);
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
                xmlDBobservableCollection.Add(new XMLdb(openPicker.FileName));
            }
            dataGridXMLs.Items.Refresh();
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
                    rawDataObservableCollection.Add(new RawData(filepath));
            dataGridDatafiles.Items.Refresh();
        }

        private void Window_Drop(object sender, DragEventArgs e)
        {
            string[] files = (string[])e.Data.GetData(DataFormats.FileDrop);
            foreach (var file in files)
            {
                var theExtension = Path.GetExtension(file);
                switch (theExtension)
                {
                    case ".raw":
                    case ".mzML":
                        rawDataObservableCollection.Add(new RawData(file));
                        break;

                    case ".xml":
                        xmlDBobservableCollection.Add(new XMLdb(file));
                        break;
                }
                dataGridDatafiles.Items.Refresh();
            }
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

        //private void TabSelectionChanges(object sender, SelectionChangedEventArgs e)
        //{
        //    TabControl tabControl = sender as TabControl; // e.Source could have been used instead of sender as well
        //    if (dataGridDatafiles != null)
        //    {
        //        if (tabControl != null && dataGridDatafiles.Columns.Count > 0)
        //        {
        //            if (tabControl.SelectedIndex == 2)
        //                //GPTMD - show extra column
        //                modificationsDataGrid.Columns[0].Visibility = Visibility.Visible;

        //            else
        //                // Hide extra column
        //                modificationsDataGrid.Columns[0].Visibility = Visibility.Hidden;
        //        }
        //        dataGridDatafiles.Items.Refresh();
        //    }
        //}

        private void RunAllTasks_Click(object sender, RoutedEventArgs e)
        {
            AllTasksEngine a = new AllTasksEngine(taskEngineObservableCollection.ToList(), rawDataObservableCollection.Where(b => b.Use).Select(b => b.FileName).ToList(), xmlDBobservableCollection.Where(b => b.Use).Select(b => b.FileName).ToList());
            var t = new Thread(() => a.Run());
            t.IsBackground = true;
            t.Start();
        }

        private void ClearTasks_Click(object sender, RoutedEventArgs e)
        {
            taskEngineObservableCollection.Clear();
            UpdateTaskGuiStuff();
        }

        private void UpdateTaskGuiStuff()
        {
            if (taskEngineObservableCollection.Count == 0)
            {
                RunTasksButton.IsEnabled = false;
                RemoveLastTaskButton.IsEnabled = false;
                ClearTasksButton.IsEnabled = false;
            }
            else
            {
                RunTasksButton.IsEnabled = true;
                RemoveLastTaskButton.IsEnabled = true;
                ClearTasksButton.IsEnabled = true;
            }
        }

        //private void addTaskButton_Click(object sender, RoutedEventArgs e)
        //{
        //    var he = new MyTask(taskTabControl.SelectedIndex);
        //    SetTaskValues(he);
        //    taskList.Insert(tasksDataGrid.Items.Count - 1, he);
        //    highlightedTaskIndex = tasksDataGrid.Items.Count - 1;
        //    RunTasksButton.IsEnabled = true;
        //}

        private void addSearchTaskButton_Click(object sender, RoutedEventArgs e)
        {
            var dialog = new SearchTaskWindow(modListObservableCollection, searchModeObservableCollection);
            if (dialog.ShowDialog() == true)
            {
                taskEngineObservableCollection.Add(dialog.TheTask);
                UpdateTaskGuiStuff();
            }
            else
            {
            }
        }

        private void addCalibrateTaskButton_Click(object sender, RoutedEventArgs e)
        {
            var dialog = new CalibrateTaskWindow(modListObservableCollection);
            if (dialog.ShowDialog() == true)
            {
                taskEngineObservableCollection.Add(dialog.TheTask);
                UpdateTaskGuiStuff();
            }
            else
            {
            }
        }

        private void addGPTMDTaskButton_Click(object sender, RoutedEventArgs e)
        {
            var dialog = new GPTMDTaskWindow(modListObservableCollection);
            if (dialog.ShowDialog() == true)
            {
                taskEngineObservableCollection.Add(dialog.TheTask);
                UpdateTaskGuiStuff();
            }
            else
            {
            }
        }

        private void RemoveLastTask_Click(object sender, RoutedEventArgs e)
        {
            taskEngineObservableCollection.RemoveAt(taskEngineObservableCollection.Count - 1);
            UpdateTaskGuiStuff();
        }

        private void tasksDataGrid_MouseDoubleClick(object sender, MouseButtonEventArgs e)
        {
            var a = sender as DataGrid;
            var ok = (MyTaskEngine)a.SelectedItem;
            switch (ok.taskType)
            {
                case MyTaskEnum.Search:
                    var searchDialog = new SearchTaskWindow(ok as SearchTask, modListObservableCollection, searchModeObservableCollection);
                    searchDialog.ShowDialog();
                    break;

                case MyTaskEnum.GPTMD:
                    var gptmddialog = new GPTMDTaskWindow(ok as GPTMDTask, modListObservableCollection);
                    gptmddialog.ShowDialog();
                    break;

                case MyTaskEnum.Calibrate:
                    var calibratedialog = new CalibrateTaskWindow(ok as CalibrationTask, modListObservableCollection);
                    calibratedialog.ShowDialog();
                    break;
            }
        }
        private void NewoutLabelStatus(object sender, string s)
        {
            if (!Dispatcher.CheckAccess())
            {
                Dispatcher.BeginInvoke(new Action(() => NewoutLabelStatus(sender, s)));
            }
            else
            {
                outProgressBar.IsIndeterminate = true;
                statusLabel.Content = s;
            }
        }

        private void NewoutProgressBar(object sender, ProgressEventArgs s)
        {
            if (!Dispatcher.CheckAccess())
            {
                Dispatcher.BeginInvoke(new Action(() => NewoutProgressBar(sender, s)));
            }
            else
            {
                outProgressBar.IsIndeterminate = false;
                outProgressBar.Value = s.new_progress;
                statusLabel.Content = s.v;
            }
        }

        private void NewoutRichTextBox(object sender, string tup)
        {
            if (!Dispatcher.CheckAccess())
            {
                Dispatcher.BeginInvoke(new Action(() => NewoutRichTextBox(sender, tup)));
            }
            else
            {
                RegOutput(tup);
            }
        }

        private void NewRefreshBetweenTasks(object sender, EventArgs e)
        {
            if (!Dispatcher.CheckAccess())
            {
                Dispatcher.BeginInvoke(new Action(() => NewRefreshBetweenTasks(sender, e)));
            }
            else
            {
                dataGridDatafiles.Items.Refresh();
                dataGridXMLs.Items.Refresh();
            }
        }

        private void NewSuccessfullyStartingAllTasks(object sender, EventArgs e)
        {
            if (!Dispatcher.CheckAccess())
            {
                Dispatcher.BeginInvoke(new Action(() => NewSuccessfullyStartingAllTasks(sender, e)));
            }
            else
            {
                //TODO: Check those
                XMLdbPanel.IsEnabled = false;
                DatafilesStackPanel.IsEnabled = false;
                addSearchTaskButton.IsEnabled = false;
                addCalibrateTaskButton.IsEnabled = false;
                addGPTMDTaskButton.IsEnabled = false;
                tasksPanel.IsEnabled = false;

                statusLabel.Content = "Starting all tasks...";
                outProgressBar.IsIndeterminate = true;

                dataGridDatafiles.Items.Refresh();
            }
        }

        private void NewSuccessfullyFinishedAllTasks(object sender, EventArgs e)
        {
            if (!Dispatcher.CheckAccess())
            {
                Dispatcher.BeginInvoke(new Action(() => NewSuccessfullyFinishedAllTasks(sender, e)));
            }
            else
            {
                //TODO: Check those
                XMLdbPanel.IsEnabled = true;
                DatafilesStackPanel.IsEnabled = true;
                addSearchTaskButton.IsEnabled = true;
                addCalibrateTaskButton.IsEnabled = true;
                addGPTMDTaskButton.IsEnabled = true;
                tasksPanel.IsEnabled = true;

                statusLabel.Content = "Finished all tasks!";
                outProgressBar.IsIndeterminate = false;
                outProgressBar.Value = 100;

                dataGridDatafiles.Items.Refresh();
            }
        }

        private void NewSuccessfullyFinishedWritingFile(object sender, SingleFileEventArgs v)
        {
            if (!Dispatcher.CheckAccess())
            {
                Dispatcher.BeginInvoke(new Action(() => NewSuccessfullyFinishedWritingFile(sender, v)));
            }
            else
            {
                addFinishedFile(v.writtenFile);
            }
        }

        private void NewUnSuccessfullyFinishedTask(object sender, string v)
        {
            if (!Dispatcher.CheckAccess())
            {
                Dispatcher.BeginInvoke(new Action(() => NewUnSuccessfullyFinishedTask(sender, v)));
            }
            else
            {
                ErrorOutput(v);
            }
        }
    }
}
