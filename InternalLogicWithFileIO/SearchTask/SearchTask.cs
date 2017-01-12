using InternalLogicEngineLayer;
using IO.MzML;
using IO.Thermo;
using MassSpectrometry;
using OldInternalLogic;
using Spectra;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace InternalLogicTaskLayer
{
    public class SearchTask : MyTaskEngine
    {

        public bool searchDecoy { get; set; }

        public bool classicSearch { get; set; }

        public bool doParsimony { get; set; }
        public List<SearchModeFoSearch> searchModes { get; set; }

        public List<ModListForSearchTask> listOfModListsForSearch { get; set; }

        public SearchTask(IEnumerable<ModList> modList, IEnumerable<SearchMode> inputSearchModes)
        {
            // Set default values here:
            classicSearch = true;
            doParsimony = true;
            searchDecoy = true;
            maxMissedCleavages = 2;
            protease = ProteaseDictionary.Instance["trypsin (no proline rule)"];
            maxModificationIsoforms = 4096;
            initiatorMethionineBehavior = InitiatorMethionineBehavior.Variable;
            productMassTolerance = new Tolerance(ToleranceUnit.Absolute, 0.01);
            bIons = true;
            yIons = true;
            listOfModListsForSearch = new List<ModListForSearchTask>();
            foreach (var uu in modList)
                listOfModListsForSearch.Add(new ModListForSearchTask(uu));
            listOfModListsForSearch[0].Fixed = true;
            listOfModListsForSearch[1].Variable = true;
            listOfModListsForSearch[2].Localize = true;

            searchModes = new List<SearchModeFoSearch>();
            foreach (var uu in inputSearchModes)
                searchModes.Add(new SearchModeFoSearch(uu));
            searchModes[0].Use = true;
            this.taskType = MyTaskEnum.Search;
        }


        protected override MyResults RunSpecific()
        {
            Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>> compactPeptideToProteinPeptideMatching = new Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>>();

            Dictionary<CompactPeptide, PeptideWithSetModifications> fullSequenceToProteinSingleMatch = new Dictionary<CompactPeptide, PeptideWithSetModifications>();

            status("Loading modifications...");
            List<MorpheusModification> variableModifications = listOfModListsForSearch.Where(b => b.Variable).SelectMany(b => b.getMods()).ToList();
            List<MorpheusModification> fixedModifications = listOfModListsForSearch.Where(b => b.Fixed).SelectMany(b => b.getMods()).ToList();
            List<MorpheusModification> localizeableModifications = listOfModListsForSearch.Where(b => b.Localize).SelectMany(b => b.getMods()).ToList();
            Dictionary<string, List<MorpheusModification>> identifiedModsInXML;
            HashSet<string> unidentifiedModStrings;
            GenerateModsFromStrings(xmlDbFilenameList, localizeableModifications, out identifiedModsInXML, out unidentifiedModStrings);

            List<SearchMode> searchModesS = searchModes.Where(b => b.Use).Select(b => b.sm).ToList();

            List<ParentSpectrumMatch>[] allPsms = new List<ParentSpectrumMatch>[searchModesS.Count];
            for (int j = 0; j < searchModesS.Count; j++)
                allPsms[j] = new List<ParentSpectrumMatch>();

            status("Loading proteins...");
            var proteinList = xmlDbFilenameList.SelectMany(b => getProteins(searchDecoy, identifiedModsInXML, b)).ToList();

            List<CompactPeptide> peptideIndex = null;
            Dictionary<float, List<int>> fragmentIndexDict = null;
            float[] keys = null;
            List<int>[] fragmentIndex = null;

            if (!classicSearch)
            {
                status("Making fragment dictionary...");

                GetPeptideAndFragmentIndices(out peptideIndex, out fragmentIndexDict, listOfModListsForSearch, searchDecoy, variableModifications, fixedModifications, localizeableModifications, proteinList, protease, output_folder);

                keys = fragmentIndexDict.OrderBy(b => b.Key).Select(b => b.Key).ToArray();
                fragmentIndex = fragmentIndexDict.OrderBy(b => b.Key).Select(b => b.Value).ToArray();
            }

            var currentRawFileList = rawDataFilenameList;
            for (int spectraFileIndex = 0; spectraFileIndex < currentRawFileList.Count; spectraFileIndex++)
            {
                var origDataFile = currentRawFileList[spectraFileIndex];
                status("Loading spectra file...");
                IMsDataFile<IMzSpectrum<MzPeak>> myMsDataFile;
                if (Path.GetExtension(origDataFile).Equals(".mzML"))
                    myMsDataFile = new Mzml(origDataFile, 400);
                else
                    myMsDataFile = new ThermoRawFile(origDataFile, 400);
                status("Opening spectra file...");
                myMsDataFile.Open();
                output("Finished opening spectra file " + Path.GetFileName(origDataFile));

                ClassicSearchEngine classicSearchEngine = null;
                ClassicSearchResults classicSearchResults = null;

                ModernSearchEngine modernSearchEngine = null;
                ModernSearchResults modernSearchResults = null;

                // run classic search
                if (classicSearch)
                {
                    classicSearchEngine = new ClassicSearchEngine(myMsDataFile, spectraFileIndex, variableModifications, fixedModifications, localizeableModifications, proteinList, productMassTolerance, protease, searchModesS);

                    classicSearchResults = (ClassicSearchResults)classicSearchEngine.Run();
                    output(classicSearchResults.ToString());
                    for (int i = 0; i < searchModesS.Count; i++)
                        allPsms[i].AddRange(classicSearchResults.outerPsms[i]);

                    AnalysisEngine analysisEngine = new AnalysisEngine(classicSearchResults.outerPsms, compactPeptideToProteinPeptideMatching, proteinList, variableModifications, fixedModifications, localizeableModifications, protease, searchModesS, myMsDataFile, productMassTolerance, (BinTreeStructure myTreeStructure, string s) => WriteTree(myTreeStructure, output_folder, Path.GetFileNameWithoutExtension(origDataFile) + s), (List<NewPsmWithFDR> h, string s) => WriteToTabDelimitedTextFileWithDecoys(h, output_folder, Path.GetFileNameWithoutExtension(origDataFile) + s), doParsimony);

                    AnalysisResults analysisResults = (AnalysisResults)analysisEngine.Run();

                    output(analysisResults.ToString());
                }

                // run modern search
                else
                {
                    modernSearchEngine = new ModernSearchEngine(myMsDataFile, spectraFileIndex, peptideIndex, keys, fragmentIndex, variableModifications, fixedModifications, localizeableModifications, proteinList, productMassTolerance.Value, protease, searchModesS);

                    modernSearchResults = (ModernSearchResults)modernSearchEngine.Run();
                    output(modernSearchResults.ToString());
                    for (int i = 0; i < searchModesS.Count; i++)
                        allPsms[i].AddRange(modernSearchResults.newPsms[i]);

                    AnalysisEngine analysisEngine = new AnalysisEngine(modernSearchResults.newPsms.Select(b => b.ToArray()).ToArray(), compactPeptideToProteinPeptideMatching, proteinList, variableModifications, fixedModifications, localizeableModifications, protease, searchModesS, myMsDataFile, productMassTolerance, (BinTreeStructure myTreeStructure, string s) => WriteTree(myTreeStructure, output_folder, Path.GetFileNameWithoutExtension(origDataFile) + s), (List<NewPsmWithFDR> h, string s) => WriteToTabDelimitedTextFileWithDecoys(h, output_folder, Path.GetFileNameWithoutExtension(origDataFile) + s), doParsimony);

                    AnalysisResults analysisResults = (AnalysisResults)analysisEngine.Run();

                    output(analysisResults.ToString());
                }
            }

            if (currentRawFileList.Count > 1)
            {
                AnalysisEngine analysisEngine = new AnalysisEngine(allPsms.Select(b => b.ToArray()).ToArray(), compactPeptideToProteinPeptideMatching, proteinList, variableModifications, fixedModifications, localizeableModifications, protease, searchModesS, null, productMassTolerance, (BinTreeStructure myTreeStructure, string s) => WriteTree(myTreeStructure, output_folder, "aggregate"), (List<NewPsmWithFDR> h, string s) => WriteToTabDelimitedTextFileWithDecoys(h, output_folder, "aggregate" + s), doParsimony);

                AnalysisResults analysisResults = (AnalysisResults)analysisEngine.Run();
                output(analysisResults.ToString());
            }
            return new MySearchTaskResults(this);
        }

        public override void ValidateParams()
        {
            foreach (var huh in listOfModListsForSearch)
            {
                if (huh.Fixed && huh.Localize)
                    throw new EngineValidationException("Not allowed to set same modifications to both fixed and localize");
                if (huh.Fixed && huh.Variable)
                    throw new EngineValidationException("Not allowed to set same modifications to both fixed and variable");
                if (huh.Localize && huh.Variable)
                    throw new EngineValidationException("Not allowed to set same modifications to both localize and variable");
            }
        }
    }
}