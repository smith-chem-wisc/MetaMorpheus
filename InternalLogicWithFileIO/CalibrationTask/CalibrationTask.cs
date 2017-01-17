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
using System.Threading.Tasks;

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
            productMassToleranceInDaltons = 0.01;
            bIons = true;
            yIons = true;
            listOfModListsForCalibration = new List<ModListForCalibrationTask>();
            foreach (var uu in modList)
                listOfModListsForCalibration.Add(new ModListForCalibrationTask(uu));
            listOfModListsForCalibration[0].Fixed = true;
            listOfModListsForCalibration[1].Variable = true;
            listOfModListsForCalibration[2].Localize = true;
            precursorMassTolerance = new Tolerance(ToleranceUnit.PPM, 10);
            taskType = MyTaskEnum.Calibrate;
            maxNumPeaksPerScan = 400;
        }

        #endregion Public Constructors

        #region Public Properties

        public List<ModListForCalibrationTask> listOfModListsForCalibration { get; set; }
        public Tolerance precursorMassTolerance { get; set; }

        public double productMassToleranceInDaltons { get; set; }

        #endregion Public Properties

        #region Protected Methods

        protected override string GetSpecificTaskInfo()
        {
            var sb = new StringBuilder();
            sb.AppendLine("Fixed mod lists: " + string.Join(",", listOfModListsForCalibration.Where(b => b.Fixed).Select(b => b.FileName)));
            sb.AppendLine("Variable mod lists: " + string.Join(",", listOfModListsForCalibration.Where(b => b.Variable).Select(b => b.FileName)));
            sb.AppendLine("Localized mod lists: " + string.Join(",", listOfModListsForCalibration.Where(b => b.Localize).Select(b => b.FileName)));
            sb.AppendLine("productMassToleranceInDaltons: " + productMassToleranceInDaltons);
            sb.Append("precursorMassTolerance: " + precursorMassTolerance);
            return sb.ToString();
        }

        protected override MyResults RunSpecific()
        {
            MyTaskResults myTaskResults = new MyCalibrationTaskResults(this);
            myTaskResults.newSpectra = new List<string>();
            var currentRawFileList = rawDataFilenameList;

            var compactPeptideToProteinPeptideMatching = new Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>>();

            SearchMode searchMode;
            if (precursorMassTolerance.Unit == ToleranceUnit.PPM)
                searchMode = new SinglePpmAroundZeroSearchMode("", precursorMassTolerance.Value);
            else
                searchMode = new SingleAbsoluteAroundZeroSearchMode("", precursorMassTolerance.Value);
            var searchModes = new List<SearchMode> { searchMode };

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

            Parallel.For(0, currentRawFileList.Count, spectraFileIndex =>
            {
                var origDataFileName = currentRawFileList[spectraFileIndex];
                LocalMs2Scan[] listOfSortedms2Scans;
                IMsDataFile<IMzSpectrum<MzPeak>> myMsDataFile;
                lock (myTaskResults)
                {
                    status("Loading spectra file " + origDataFileName + "...");
                    if (Path.GetExtension(origDataFileName).Equals(".mzML"))
                        myMsDataFile = new Mzml(origDataFileName, maxNumPeaksPerScan);
                    else
                        myMsDataFile = new ThermoRawFile(origDataFileName, maxNumPeaksPerScan);
                    status("Opening spectra file " + origDataFileName + "...");
                    myMsDataFile.Open();
                    listOfSortedms2Scans = myMsDataFile.Where(b => b.MsnOrder == 2).Select(b => new LocalMs2Scan(b)).OrderBy(b => b.precursorMass).ToArray();
                }

                var searchEngine = new ClassicSearchEngine(listOfSortedms2Scans, myMsDataFile.NumSpectra, spectraFileIndex, variableModifications, fixedModifications, proteinList, new Tolerance(ToleranceUnit.Absolute, productMassToleranceInDaltons), protease, searchModes);

                var searchResults = (ClassicSearchResults)searchEngine.Run();

                for (int i = 0; i < searchModes.Count; i++)
                    allPsms[i].AddRange(searchResults.outerPsms[i]);

                // Run analysis on single file results
                var analysisEngine = new AnalysisEngine(searchResults.outerPsms, compactPeptideToProteinPeptideMatching, proteinList, variableModifications, fixedModifications, localizeableModifications, protease, searchModes, myMsDataFile, new Tolerance(ToleranceUnit.Absolute, productMassToleranceInDaltons), (BinTreeStructure myTreeStructure, string s) => WriteTree(myTreeStructure, output_folder, Path.GetFileNameWithoutExtension(origDataFileName) + s), (List<NewPsmWithFDR> h, string s) => WritePSMsToTSV(h, output_folder, Path.GetFileNameWithoutExtension(origDataFileName) + s), null, false);

                var analysisResults = (AnalysisResults)analysisEngine.Run();

                var identifications = analysisResults.allResultingIdentifications[0];

                myMsDataFile.Close();
                myMsDataFile = null;

                //Now can calibrate!!!
                IMsDataFile<IMzSpectrum<MzPeak>> myMsDataFileForCalibration;
                if (Path.GetExtension(origDataFileName).Equals(".mzML"))
                {
                    myMsDataFileForCalibration = new Mzml(origDataFileName);
                    myMsDataFileForCalibration.Open();
                }
                else
                {
                    myMsDataFileForCalibration = new ThermoRawFile(origDataFileName);
                    myMsDataFileForCalibration.Open();
                }
                int randomSeed = 1;

                // TODO: fix the tolerance calculation below
                var a = new CalibrationEngine(myMsDataFileForCalibration, randomSeed, productMassToleranceInDaltons * 2, identifications);

                var result = (CalibrationResults)a.Run();

                status("Creating _indexedmzMLConnection, putting data in it, and writing!");
                var path = Path.Combine(output_folder, Path.GetFileNameWithoutExtension(origDataFileName) + "-Calibrated.mzML");
                MzmlMethods.CreateAndWriteMyIndexedMZmlwithCalibratedSpectra(result.myMsDataFile, path);

                SucessfullyFinishedWritingFile(path);

                lock (myTaskResults)
                {
                    myTaskResults.newSpectra.Add(path);
                }
            }
            );
            return myTaskResults;
        }

        #endregion Protected Methods

    }
}