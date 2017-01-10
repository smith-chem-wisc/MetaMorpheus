using IO.MzML;
using IO.Thermo;
using MassSpectrometry;
using MetaMorpheus;
using Spectra;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.IO;
using System.Linq;

namespace IndexSearchAndAnalyze
{
    public class MySearchTask : MyTask
    {
        public bool searchDecoy { get; set; }

        public bool classicSearch { get; set; }

        public MySearchTask(IEnumerable<ModList> modList) : base(1)
        {
            // Set default values here:
            classicSearch = true;
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

        public override MyTaskResults DoTask(AllTasksParams po)
        {
            Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>> compactPeptideToProteinPeptideMatching = new Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>>();

            Dictionary<CompactPeptide, PeptideWithSetModifications> fullSequenceToProteinSingleMatch = new Dictionary<CompactPeptide, PeptideWithSetModifications>();

            SearchMode searchMode = new IntervalSearchMode("withinhalfAdaltonOfZero", new List<double>() { 0 }, new Tolerance(ToleranceUnit.PPM, 5));
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
            var proteinList = po.xMLdblist.SelectMany(b => getProteins(searchDecoy, identifiedModsInXML, b)).ToList();

            po.status("Making fragment dictionary...");
            List<CompactPeptide> peptideIndex;
            Dictionary<float, List<int>> fragmentIndexDict;

            Indices.GetPeptideAndFragmentIndices(out peptideIndex, out fragmentIndexDict, null, listOfModListsForSearch, searchDecoy, variableModifications, fixedModifications, localizeableModifications, proteinList, protease, po);

            var keys = fragmentIndexDict.OrderBy(b => b.Key).Select(b => b.Key).ToArray();
            var fragmentIndex = fragmentIndexDict.OrderBy(b => b.Key).Select(b => b.Value).ToArray();

            var currentRawFileList = po.rawDataAndResultslist;
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
                po.RTBoutput("Finished opening spectra file " + Path.GetFileName(origDataFile));

                ClassicSearchParams classicSearchParams = null;
                ClassicSearchEngine classicSearchEngine = null;
                ClassicSearchResults classicSearchResults = null;

                ModernSearchParams modernSearchParams = null;
                ModernSearchEngine modernSearchEngine = null;
                ModernSearchResults modernSearchResults = null;

                // run classic search
                if (classicSearch)
                {
                    // classic
                    // IMsDataFile<IMzSpectrum<MzPeak>> myMsDataFile, int spectraFileIndex, List<MorpheusModification> variableModifications, List<MorpheusModification> fixedModifications, List<MorpheusModification> localizeableModifications, List<Protein> proteinList, Tolerance fragmentTolerance, Protease protease, SearchMode searchMode, AllTasksParams a2) : base(a2)
                    classicSearchParams = new ClassicSearchParams(myMsDataFile, spectraFileIndex, variableModifications, fixedModifications, localizeableModifications, proteinList, productMassTolerance, protease, searchModes, po);
                    classicSearchEngine = new ClassicSearchEngine(classicSearchParams);
                    classicSearchResults = (ClassicSearchResults)classicSearchEngine.Run();
                    po.RTBoutput(classicSearchResults.ToString());
                    for (int i = 0; i < searchModes.Count; i++)
                        allPsms[i].AddRange(classicSearchResults.outerPsms[i]);

                    AnalysisParams analysisParams = new AnalysisParams(classicSearchResults.outerPsms, compactPeptideToProteinPeptideMatching, proteinList, variableModifications, fixedModifications, localizeableModifications, protease, searchModes, myMsDataFile, productMassTolerance, (BinTreeStructure myTreeStructure, string s) => Writing.WriteTree(myTreeStructure, output_folder, Path.GetFileNameWithoutExtension(origDataFile) + s, po), (List<NewPsmWithFDR> h, string s) => Writing.WriteToTabDelimitedTextFileWithDecoys(h, output_folder, Path.GetFileNameWithoutExtension(origDataFile) + s, po), po);
                    AnalysisEngine analysisEngine = new AnalysisEngine(analysisParams);
                    AnalysisResults analysisResults = (AnalysisResults)analysisEngine.Run();

                    po.RTBoutput(analysisResults.ToString());
                }

                // run modern search
                else
                {
                    // modern
                    //                                          myMsDataFile, spectraFileIndex, peptideIndex, keys, fragmentIndex, variableModifications, fixedModifications, localizeableModifications, proteinList, fragmentTolerance,    protease, searchModes,    a2)
                    modernSearchParams = new ModernSearchParams(myMsDataFile, spectraFileIndex, peptideIndex, keys, fragmentIndex, variableModifications, fixedModifications, localizeableModifications, proteinList, productMassTolerance.Value, protease, null, po);
                    modernSearchEngine = new ModernSearchEngine(modernSearchParams);
                    modernSearchResults = (ModernSearchResults)modernSearchEngine.Run();
                    po.RTBoutput(modernSearchResults.ToString());
                    for (int i = 0; i < searchModes.Count; i++)
                        allPsms[i].AddRange(modernSearchResults.newPsms[i]);

                    AnalysisParams analysisParams = new AnalysisParams(modernSearchResults.newPsms.Select(b => b.ToArray()).ToArray(), compactPeptideToProteinPeptideMatching, proteinList, variableModifications, fixedModifications, localizeableModifications, protease, searchModes, myMsDataFile, productMassTolerance, (BinTreeStructure myTreeStructure, string s) => Writing.WriteTree(myTreeStructure, output_folder, Path.GetFileNameWithoutExtension(origDataFile) + s, po), (List<NewPsmWithFDR> h, string s) => Writing.WriteToTabDelimitedTextFileWithDecoys(h, output_folder, Path.GetFileNameWithoutExtension(origDataFile) + s, po), po);
                    AnalysisEngine analysisEngine = new AnalysisEngine(analysisParams);
                    AnalysisResults analysisResults = (AnalysisResults)analysisEngine.Run();

                    po.RTBoutput(analysisResults.ToString());
                }
                

            }

            if (currentRawFileList.Count > 1)
            {
                AnalysisParams analysisParams = new AnalysisParams(allPsms.Select(b => b.ToArray()).ToArray(), compactPeptideToProteinPeptideMatching, proteinList, variableModifications, fixedModifications, localizeableModifications, protease, searchModes, null, productMassTolerance, (BinTreeStructure myTreeStructure, string s) => Writing.WriteTree(myTreeStructure, output_folder, "aggregate", po), (List<NewPsmWithFDR> h, string s) => Writing.WriteToTabDelimitedTextFileWithDecoys(h, output_folder, "aggregate" + s, po), po);
                AnalysisEngine analysisEngine = new AnalysisEngine(analysisParams);
                AnalysisResults analysisResults = (AnalysisResults)analysisEngine.Run();
                po.RTBoutput(analysisResults.ToString());
            }
            return new MyTaskResults();
        }
    }
}
