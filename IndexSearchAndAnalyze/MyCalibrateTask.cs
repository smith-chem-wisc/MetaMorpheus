using IO.MzML;
using IO.Thermo;
using MassSpectrometry;
using MetaMorpheus;
using Spectra;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.IO;
using System.Linq;

namespace IndexSearchAndAnalyze
{
    public class MyCalibrateTask : MyTask
    {
        public MyCalibrateTask() : base(0)
        {
        }

        public Tolerance precursorMassTolerance { get; set; }

        public MyCalibrateTask(ObservableCollection<ModList> modList) : base(0)
        {
            // Set default values here:
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
            precursorMassTolerance = new Tolerance(ToleranceUnit.PPM, 10);
        }

        public override MyTaskResults DoTask(AllTasksParams po)
        {
            MyTaskResults myTaskResults = new MyTaskResults();
            myTaskResults.newSpectra = new List<string>();
            var currentRawFileList = po.rawDataAndResultslist;

            Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>> compactPeptideToProteinPeptideMatching = new Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>>();

            Dictionary<CompactPeptide, PeptideWithSetModifications> fullSequenceToProteinSingleMatch = new Dictionary<CompactPeptide, PeptideWithSetModifications>();

            SearchMode searchMode = new DotSearchMode("withinhalfAdaltonOfZero", new double[] { 0 }, new Tolerance(ToleranceUnit.PPM, 5));
            List<SearchMode> searchModes = new List<SearchMode>() { searchMode };

            List<ParentSpectrumMatch>[] allPsms = new List<ParentSpectrumMatch>[searchModes.Count];
            for (int j = 0; j < searchModes.Count; j++)
                allPsms[j] = new List<ParentSpectrumMatch>();

            po.status("Loading modifications...");
            List<MorpheusModification> variableModifications = listOfModListsForSearch.Where(b => b.Variable).SelectMany(b => b.getMods()).ToList();
            List<MorpheusModification> fixedModifications = listOfModListsForSearch.Where(b => b.Fixed).SelectMany(b => b.getMods()).ToList();
            List<MorpheusModification> localizeableModifications = listOfModListsForSearch.Where(b => b.Localize).SelectMany(b => b.getMods()).ToList();
            Dictionary<string, List<MorpheusModification>> identifiedModsInXML;
            HashSet<string> unidentifiedModStrings;
            GenerateModsFromStrings(po.xMLdblist, localizeableModifications, out identifiedModsInXML, out unidentifiedModStrings);

            po.status("Loading proteins...");
            var proteinList = po.xMLdblist.SelectMany(b => getProteins(true, identifiedModsInXML, b)).ToList();

            for (int spectraFileIndex = 0; spectraFileIndex < currentRawFileList.Count; spectraFileIndex++)
            {
                var origDataFile = currentRawFileList[spectraFileIndex];
                po.status("Loading spectra file...");
                IMsDataFile<IMzSpectrum<MzPeak>> myMsDataFile;
                if (Path.GetExtension(origDataFile).Equals(".mzML"))
                    myMsDataFile = new Mzml(origDataFile, 400);
                else
                    myMsDataFile = new ThermoRawFile(origDataFile, 400);
                po.status("Opening spectra file...");
                myMsDataFile.Open();
                po.output("Finished opening spectra file " + Path.GetFileName(origDataFile));

                ClassicSearchParams searchParams = new ClassicSearchParams(myMsDataFile, spectraFileIndex, variableModifications, fixedModifications, localizeableModifications, proteinList, productMassTolerance, protease, searchModes, po);
                ClassicSearchEngine searchEngine = new ClassicSearchEngine(searchParams);
                ClassicSearchResults searchResults = (ClassicSearchResults)searchEngine.Run();
                po.output(searchResults.ToString());

                for (int i = 0; i < searchModes.Count; i++)
                    allPsms[i].AddRange(searchResults.outerPsms[i]);

                // Run analysis on single file results
                AnalysisParams analysisParams = new AnalysisParams(searchResults.outerPsms, compactPeptideToProteinPeptideMatching, proteinList, variableModifications, fixedModifications, localizeableModifications, protease, searchModes, myMsDataFile, productMassTolerance, (BinTreeStructure myTreeStructure, string s) => Writing.WriteTree(myTreeStructure, output_folder, Path.GetFileNameWithoutExtension(origDataFile) + s, po), (List<NewPsmWithFDR> h, string s) => Writing.WriteToTabDelimitedTextFileWithDecoys(h, output_folder, Path.GetFileNameWithoutExtension(origDataFile) + s, po), po);
                AnalysisEngine analysisEngine = new AnalysisEngine(analysisParams);
                AnalysisResults analysisResults = (AnalysisResults)analysisEngine.Run();

                po.output(analysisResults.ToString());

                var identifications = analysisResults.allResultingIdentifications[0];

                myMsDataFile.Close();
                myMsDataFile = null;

                //Now can calibrate!!!

                SoftwareLockMassParams a = mzCalIO.GetReady(origDataFile, identifications, productMassTolerance, po, myTaskResults);

                SoftwareLockMassRunner.Run(a);
            }
            return myTaskResults;
        }
    }
}