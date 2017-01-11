using InternalLogic;
using InternalLogicCalibration;
using IO.MzML;
using IO.Thermo;
using MassSpectrometry;
using OldInternalLogic;
using Spectra;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.IO;
using System.Linq;
using System;

namespace InternalLogicWithFileIO
{
    public class MyCalibrateTask : MyTaskEngine
    {
        public Tolerance precursorMassTolerance { get; set; }

        public MyCalibrateTask(ObservableCollection<ModList> modList)
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

        public override void ValidateParams()
        {
            throw new NotImplementedException();
        }

        protected override MyResults RunSpecific()
        {
            MyTaskResults myTaskResults = new MyTaskResults(this);
            myTaskResults.newSpectra = new List<string>();
            var currentRawFileList = rawDataAndResultslist;

            Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>> compactPeptideToProteinPeptideMatching = new Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>>();

            Dictionary<CompactPeptide, PeptideWithSetModifications> fullSequenceToProteinSingleMatch = new Dictionary<CompactPeptide, PeptideWithSetModifications>();

            SearchMode searchMode = new DotSearchMode("withinhalfAdaltonOfZero", new double[] { 0 }, new Tolerance(ToleranceUnit.PPM, 5));
            List<SearchMode> searchModes = new List<SearchMode>() { searchMode };

            List<ParentSpectrumMatch>[] allPsms = new List<ParentSpectrumMatch>[searchModes.Count];
            for (int j = 0; j < searchModes.Count; j++)
                allPsms[j] = new List<ParentSpectrumMatch>();

            status("Loading modifications...");
            List<MorpheusModification> variableModifications = listOfModListsForSearch.Where(b => b.Variable).SelectMany(b => b.getMods()).ToList();
            List<MorpheusModification> fixedModifications = listOfModListsForSearch.Where(b => b.Fixed).SelectMany(b => b.getMods()).ToList();
            List<MorpheusModification> localizeableModifications = listOfModListsForSearch.Where(b => b.Localize).SelectMany(b => b.getMods()).ToList();
            Dictionary<string, List<MorpheusModification>> identifiedModsInXML;
            HashSet<string> unidentifiedModStrings;
            GenerateModsFromStrings(xMLdblist, localizeableModifications, out identifiedModsInXML, out unidentifiedModStrings);

            status("Loading proteins...");
            var proteinList = xMLdblist.SelectMany(b => getProteins(true, identifiedModsInXML, b)).ToList();

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

                ClassicSearchEngine searchEngine = new ClassicSearchEngine(myMsDataFile, spectraFileIndex, variableModifications, fixedModifications, localizeableModifications, proteinList, productMassTolerance, protease, searchModes);

                ClassicSearchResults searchResults = (ClassicSearchResults)searchEngine.Run();
                output(searchResults.ToString());

                for (int i = 0; i < searchModes.Count; i++)
                    allPsms[i].AddRange(searchResults.outerPsms[i]);

                // Run analysis on single file results
                AnalysisEngine analysisEngine = new AnalysisEngine(searchResults.outerPsms, compactPeptideToProteinPeptideMatching, proteinList, variableModifications, fixedModifications, localizeableModifications, protease, searchModes, myMsDataFile, productMassTolerance, (BinTreeStructure myTreeStructure, string s) => Writing.WriteTree(myTreeStructure, output_folder, Path.GetFileNameWithoutExtension(origDataFile) + s), (List<NewPsmWithFDR> h, string s) => Writing.WriteToTabDelimitedTextFileWithDecoys(h, output_folder, Path.GetFileNameWithoutExtension(origDataFile) + s), false);

                AnalysisResults analysisResults = (AnalysisResults)analysisEngine.Run();

                output(analysisResults.ToString());

                var identifications = analysisResults.allResultingIdentifications[0];

                myMsDataFile.Close();
                myMsDataFile = null;

                //Now can calibrate!!!
                IMsDataFile<IMzSpectrum<MzPeak>> myMsDataFileForCalibration;
                if (Path.GetExtension(origDataFile).Equals(".mzML"))
                {
                    myMsDataFileForCalibration = new Mzml(origDataFile);
                    myMsDataFileForCalibration.Open();
                }
                else
                {
                    myMsDataFileForCalibration = new ThermoRawFile(origDataFile);
                    myMsDataFileForCalibration.Open();
                }
                int randomSeed = 1;

                // TODO: fix the tolerance calculation below
                var a = new CalibrationEngine(myMsDataFileForCalibration, randomSeed, productMassTolerance.Value * 2);

                a.identifications = identifications;
                a.mzRange = new DoubleRange(0, 0);

                //a.MS1spectraToWatch.Add(22557);

                //a.MS2spectraToWatch.Add(22564);

                a.matchesToExclude = new HashSet<int>();

                a.Run();
            }
            return myTaskResults;
        }
    }
}