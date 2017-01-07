using IndexSearchAndAnalyze;
using IO.MzML;
using IO.Thermo;
using MassSpectrometry;
using MetaMorpheus;
using Spectra;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.IO;
using System.Linq;

namespace GoodGUI
{
    public class MySearchTask : MyTask
    {
        public MySearchTask(IEnumerable<ModList> modList) : base(1)
        {
            // Set default values here:
            searchDecoy = true;
            maxMissedCleavages = 2;
            protease = ProteaseDictionary.Instance["trypsin (no proline rule)"];
            maxModificationIsoforms = 4096;
            initiatorMethionineBehavior = InitiatorMethionineBehavior.Variable;
            productMassTolerance = new Tolerance(ToleranceUnit.Absolute, 0.01);
            bIons = true;
            yIons = true;
            listOfModListsForSearch = new List<ModListForSearch>();
            foreach (var uu in modList)
                listOfModListsForSearch.Add(new ModListForSearch(uu));
            listOfModListsForSearch[0].Fixed = true;
            listOfModListsForSearch[1].Variable = true;
            listOfModListsForSearch[2].Localize = true;
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
        public List<ModListForSearch> listOfModListsForSearch { get; internal set; }

        internal override void DoTask(ObservableCollection<RawData> completeRawFileListCollection, ObservableCollection<XMLdb> completeXmlDbList, AllTasksParams po)
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

                List<MorpheusModification> variableModifications = listOfModListsForSearch.Where(b => b.Variable).SelectMany(b => b.getMods()).ToList();
                List<MorpheusModification> fixedModifications = listOfModListsForSearch.Where(b => b.Fixed).SelectMany(b => b.getMods()).ToList();
                List<MorpheusModification> localizeableModifications = listOfModListsForSearch.Where(b => b.Localize).SelectMany(b => b.getMods()).ToList();
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