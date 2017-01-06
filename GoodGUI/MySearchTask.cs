using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using MetaMorpheus;
using static GoodGUI.MainWindow;
using IndexSearchAndAnalyze;
using Spectra;
using System.IO;
using MassSpectrometry;
using IO.MzML;
using IO.Thermo;
using System.Linq;
using mzCal;
using Proteomics;
using FragmentGeneration;
using System.Globalization;

namespace GoodGUI
{
    public class MySearchTask : MyTask
    {

        public MySearchTask(ObservableCollection<ModListForSearch> modFileList) : base(1)
        {
            this.modFileList = modFileList;

            modFileList[0].Fixed = true;
            modFileList[1].Variable = true;
            modFileList[2].Localize = true;

            // Set default values here:
            searchDecoy = true;
            maxMissedCleavages = 2;
            protease = ProteaseDictionary.Instance["trypsin (no proline rule)"];
            maxModificationIsoforms = 4096;
            initiatorMethionineBehavior = InitiatorMethionineBehavior.Variable;
            productMassTolerance = new Tolerance(ToleranceUnit.Absolute, 0.01);
            bIons = true;
            yIons = true;
        }

        public string precursorMassToleranceTextBox { get; internal set; }
        public int precursorMassToleranceComboBox { get; internal set; }
        public int maxMissedCleavages { get; internal set; }
        public Protease protease { get; internal set; }
        public int maxModificationIsoforms { get; internal set; }
        public InitiatorMethionineBehavior initiatorMethionineBehavior { get; internal set; }
        public Tolerance productMassTolerance { get; internal set; }
        public bool bIons { get; internal set; }
        public bool yIons { get; internal set; }
        public bool searchDecoy { get; internal set; }
        public ObservableCollection<ModListForSearch> modFileList { get; private set; }

        internal override void DoTask(ObservableCollection<RawData> completeRawFileListCollection, ObservableCollection<XMLdb> completeXmlDbList, ParamsObject po)
        {
            var currentRawFileList = completeRawFileListCollection.Where(b => b.Use).Select(b => b.FileName).ToList();
            for (int spectraFileIndex = 0; spectraFileIndex < currentRawFileList.Count; spectraFileIndex++)
            {
                var origDataFile = currentRawFileList[spectraFileIndex];
                po.RTBoutput("Loading spectra file...");
                IMsDataFile<IMzSpectrum<MzPeak>> myMsDataFile;
                if (Path.GetExtension(origDataFile).Equals(".mzML"))
                    myMsDataFile = new Mzml(origDataFile, 400);
                else
                    myMsDataFile = new ThermoRawFile(origDataFile, 400);
                po.RTBoutput("Opening spectra file...");
                myMsDataFile.Open();
                po.RTBoutput("Finished opening spectra file " + Path.GetFileName(origDataFile));

                List<MorpheusModification> variableModifications = modFileList.Where(b => b.Variable).SelectMany(b => b.getMods()).ToList();
                List<MorpheusModification> fixedModifications = modFileList.Where(b => b.Fixed).SelectMany(b => b.getMods()).ToList();
                List<MorpheusModification> localizeableModifications = modFileList.Where(b => b.Localize).SelectMany(b => b.getMods()).ToList();
                Dictionary<string, List<MorpheusModification>> identifiedModsInXML;
                HashSet<string> unidentifiedModStrings;
                GenerateModsFromStrings(completeXmlDbList.Select(b => b.FileName).ToList(), localizeableModifications, out identifiedModsInXML, out unidentifiedModStrings);

                var proteinList = completeXmlDbList.SelectMany(b => b.getProteins(searchDecoy, identifiedModsInXML)).ToList();

                List<SearchMode> searchModes = new List<SearchMode>();
                searchModes.Add(new IntervalSearchMode("withinhalfAdaltonOfZero", new List<double>() { 0 }, new Tolerance(ToleranceUnit.PPM, 5)));

                ClassicSearchParams searchParams = new ClassicSearchParams(myMsDataFile, spectraFileIndex, variableModifications, fixedModifications, localizeableModifications, proteinList, productMassTolerance, protease, searchModes[0], po.RTBoutput, po.ReportProgress);
                ClassicSearchEngine searchEngine = new ClassicSearchEngine(searchParams);
                ClassicSearchResults searchResults = (ClassicSearchResults)searchEngine.Run();
                po.RTBoutput(searchResults.ToString());

            }
        }
    }
}