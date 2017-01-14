using InternalLogicCalibration;
using InternalLogicEngineLayer;
using IO.MzML;
using IO.Thermo;
using MassSpectrometry;
using OldInternalLogic;
using Spectra;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.IO;
using System.Linq;
using System.Text;

namespace InternalLogicTaskLayer
{
    public class CalibrationTask : MyTaskEngine
    {

        #region Public Constructors

        public CalibrationTask(ObservableCollection<ModList> modList)
        {
            // Set default values here:
            maxMissedCleavages = 2;
            protease = ProteaseDictionary.Instance["trypsin (no proline rule)"];
            maxModificationIsoforms = 4096;
            initiatorMethionineBehavior = InitiatorMethionineBehavior.Variable;
            productMassTolerance = new Tolerance(ToleranceUnit.Absolute, 0.01);
            bIons = true;
            yIons = true;
            listOfModListsForCalibration = new List<ModListForCalibrationTask>();
            foreach (var uu in modList)
                listOfModListsForCalibration.Add(new ModListForCalibrationTask(uu));
            listOfModListsForCalibration[0].Fixed = true;
            listOfModListsForCalibration[1].Variable = true;
            listOfModListsForCalibration[2].Localize = true;
            precursorMassTolerance = new Tolerance(ToleranceUnit.PPM, 10);
            this.taskType = MyTaskEnum.Calibrate;
			maxNumPeaksPerScan = 400;
        }

        #endregion Public Constructors

        #region Public Properties

        public List<ModListForCalibrationTask> listOfModListsForCalibration { get; set; }
        public Tolerance precursorMassTolerance { get; set; }

        #endregion Public Properties

        #region Protected Methods

        protected override void ValidateParams()
        {
            if (xmlDbFilenameList == null)
                throw new EngineValidationException("xMLdblist cannot be null");
            if (xmlDbFilenameList.Count == 0)
                throw new EngineValidationException("xMLdblist cannot be empty");
            if (rawDataFilenameList == null)
                throw new EngineValidationException("rawDataAndResultslist cannot be null");
            if (rawDataFilenameList.Count == 0)
                throw new EngineValidationException("rawDataAndResultslist cannot be empty");
        }

        protected override string GetSpecificTaskInfo()
        {
            StringBuilder sb = new StringBuilder();
            sb.AppendLine("Fixed mod lists: " + string.Join(",", listOfModListsForCalibration.Where(b => b.Fixed).Select(b => b.FileName)));
            sb.AppendLine("Variable mod lists: " + string.Join(",", listOfModListsForCalibration.Where(b => b.Variable).Select(b => b.FileName)));
            sb.AppendLine("Localized mod lists: " + string.Join(",", listOfModListsForCalibration.Where(b => b.Localize).Select(b => b.FileName)));
            sb.AppendLine("precursorMassTolerance: " + precursorMassTolerance);
            return sb.ToString();
        }

        protected override MyResults RunSpecific()
        {
            MyTaskResults myTaskResults = new MyCalibrationTaskResults(this);
            myTaskResults.newSpectra = new List<string>();
            var currentRawFileList = rawDataFilenameList;

            Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>> compactPeptideToProteinPeptideMatching = new Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>>();

            Dictionary<CompactPeptide, PeptideWithSetModifications> fullSequenceToProteinSingleMatch = new Dictionary<CompactPeptide, PeptideWithSetModifications>();
            SearchMode searchMode;
            if (precursorMassTolerance.Unit == ToleranceUnit.PPM)
                searchMode = new SinglePpmAroundZeroSearchMode("", precursorMassTolerance.Value);
            else
                searchMode = new SingleAbsoluteAroundZeroSearchMode("", precursorMassTolerance.Value);
            List<SearchMode> searchModes = new List<SearchMode>() { searchMode };

            List<ParentSpectrumMatch>[] allPsms = new List<ParentSpectrumMatch>[1];
            allPsms[0] = new List<ParentSpectrumMatch>();

            status("Loading modifications...");
            List<MorpheusModification> variableModifications = listOfModListsForCalibration.Where(b => b.Variable).SelectMany(b => b.getMods()).ToList();
            List<MorpheusModification> fixedModifications = listOfModListsForCalibration.Where(b => b.Fixed).SelectMany(b => b.getMods()).ToList();
            List<MorpheusModification> localizeableModifications = listOfModListsForCalibration.Where(b => b.Localize).SelectMany(b => b.getMods()).ToList();
            Dictionary<string, List<MorpheusModification>> identifiedModsInXML;
            HashSet<string> unidentifiedModStrings;
            MatchXMLmodsToKnownMods(xmlDbFilenameList, localizeableModifications, out identifiedModsInXML, out unidentifiedModStrings);

            status("Loading proteins...");
            var proteinList = xmlDbFilenameList.SelectMany(b => getProteins(true, identifiedModsInXML, b)).ToList();

            for (int spectraFileIndex = 0; spectraFileIndex < currentRawFileList.Count; spectraFileIndex++)
            {
                var origDataFile = currentRawFileList[spectraFileIndex];
                status("Loading spectra file...");
                IMsDataFile<IMzSpectrum<MzPeak>> myMsDataFile;
                if (Path.GetExtension(origDataFile).Equals(".mzML"))
                    myMsDataFile = new Mzml(origDataFile, maxNumPeaksPerScan);
                else
                    myMsDataFile = new ThermoRawFile(origDataFile, maxNumPeaksPerScan);
                status("Opening spectra file...");
                myMsDataFile.Open();
                //output("Finished opening spectra file " + Path.GetFileName(origDataFile));

                ClassicSearchEngine searchEngine = new ClassicSearchEngine(myMsDataFile, spectraFileIndex, variableModifications, fixedModifications, localizeableModifications, proteinList, productMassTolerance, protease, searchModes);

                ClassicSearchResults searchResults = (ClassicSearchResults)searchEngine.Run();

                for (int i = 0; i < searchModes.Count; i++)
                    allPsms[i].AddRange(searchResults.outerPsms[i]);

                // Run analysis on single file results
                AnalysisEngine analysisEngine = new AnalysisEngine(searchResults.outerPsms, compactPeptideToProteinPeptideMatching, proteinList, variableModifications, fixedModifications, localizeableModifications, protease, searchModes, myMsDataFile, productMassTolerance, (BinTreeStructure myTreeStructure, string s) => WriteTree(myTreeStructure, output_folder, Path.GetFileNameWithoutExtension(origDataFile) + s), (List<NewPsmWithFDR> h, string s) => WritePSMsToTSV(h, output_folder, Path.GetFileNameWithoutExtension(origDataFile) + s), false);

                AnalysisResults analysisResults = (AnalysisResults)analysisEngine.Run();

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
                var a = new CalibrationEngine(myMsDataFileForCalibration, randomSeed, productMassTolerance.Value * 2, identifications);

                var result = (CalibrationResults)a.Run();

                status("Creating _indexedmzMLConnection, and putting data in it");
                var path = Path.Combine(Path.GetDirectoryName(output_folder), Path.GetFileNameWithoutExtension(origDataFile) + "-Calibrated.mzML");
                MzmlMethods.CreateAndWriteMyIndexedMZmlwithCalibratedSpectra(result.myMsDataFile, path);

                SucessfullyFinishedWritingFile(path);

                myTaskResults.newSpectra.Add(path);
            }
            return myTaskResults;
        }

        #endregion Protected Methods

    }
}