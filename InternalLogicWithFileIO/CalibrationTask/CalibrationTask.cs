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
            MaxMissedCleavages = 2;
            Protease = ProteaseDictionary.Instance["trypsin (no proline rule)"];
            MaxModificationIsoforms = 4096;
            InitiatorMethionineBehavior = InitiatorMethionineBehavior.Variable;
            ProductMassToleranceInDaltons = 0.01;
            BIons = true;
            YIons = true;
            ListOfModListsForCalibration = new List<ModListForCalibrationTask>();
            foreach (var uu in modList)
                ListOfModListsForCalibration.Add(new ModListForCalibrationTask(uu));
            ListOfModListsForCalibration[0].Fixed = true;
            ListOfModListsForCalibration[1].Variable = true;
            ListOfModListsForCalibration[2].Localize = true;
            PrecursorMassTolerance = new Tolerance(ToleranceUnit.PPM, 10);
            TaskType = MyTask.Calibrate;
            MaxNumPeaksPerScan = 400;
        }

        #endregion Public Constructors

        #region Public Properties

        public List<ModListForCalibrationTask> ListOfModListsForCalibration { get; set; }
        public Tolerance PrecursorMassTolerance { get; set; }

        public double ProductMassToleranceInDaltons { get; set; }

        #endregion Public Properties

        #region Protected Properties

        protected override string SpecificTaskInfo
        {
            get
            {
                var sb = new StringBuilder();
                sb.AppendLine("Fixed mod lists: " + string.Join(",", ListOfModListsForCalibration.Where(b => b.Fixed).Select(b => b.FileName)));
                sb.AppendLine("Variable mod lists: " + string.Join(",", ListOfModListsForCalibration.Where(b => b.Variable).Select(b => b.FileName)));
                sb.AppendLine("Localized mod lists: " + string.Join(",", ListOfModListsForCalibration.Where(b => b.Localize).Select(b => b.FileName)));
                sb.AppendLine("productMassToleranceInDaltons: " + ProductMassToleranceInDaltons);
                sb.Append("precursorMassTolerance: " + PrecursorMassTolerance);
                return sb.ToString();
            }
        }

        #endregion Protected Properties

        #region Protected Methods

        protected override MyResults RunSpecific()
        {
            MyTaskResults myTaskResults = new MyCalibrationTaskResults(this);
            myTaskResults.newSpectra = new List<string>();
            var currentRawFileList = rawDataFilenameList;

            var compactPeptideToProteinPeptideMatching = new Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>>();

            SearchMode searchMode;
            if (PrecursorMassTolerance.Unit == ToleranceUnit.PPM)
                searchMode = new SinglePpmAroundZeroSearchMode("", PrecursorMassTolerance.Value);
            else
                searchMode = new SingleAbsoluteAroundZeroSearchMode("", PrecursorMassTolerance.Value);
            var searchModes = new List<SearchMode> { searchMode };

            List<ParentSpectrumMatch>[] allPsms = new List<ParentSpectrumMatch>[1];
            allPsms[0] = new List<ParentSpectrumMatch>();

            Status("Loading modifications...");
            List<MorpheusModification> variableModifications = ListOfModListsForCalibration.Where(b => b.Variable).SelectMany(b => b.Mods).ToList();
            List<MorpheusModification> fixedModifications = ListOfModListsForCalibration.Where(b => b.Fixed).SelectMany(b => b.Mods).ToList();
            List<MorpheusModification> localizeableModifications = ListOfModListsForCalibration.Where(b => b.Localize).SelectMany(b => b.Mods).ToList();
            Dictionary<string, List<MorpheusModification>> identifiedModsInXML;
            HashSet<string> unidentifiedModStrings;
            MatchXMLmodsToKnownMods(xmlDbFilenameList, localizeableModifications, out identifiedModsInXML, out unidentifiedModStrings);

            Status("Loading proteins...");
            var proteinList = xmlDbFilenameList.SelectMany(b => GetProteins(true, identifiedModsInXML, b)).ToList();

            Parallel.For(0, currentRawFileList.Count, spectraFileIndex =>
            {
                var origDataFileName = currentRawFileList[spectraFileIndex];
                LocalMS2Scan[] listOfSortedms2Scans;
                IMsDataFile<IMzSpectrum<MzPeak>> myMsDataFile;
                lock (myTaskResults)
                {
                    Status("Loading spectra file " + origDataFileName + "...");
                    if (Path.GetExtension(origDataFileName).Equals(".mzML"))
                        myMsDataFile = new Mzml(origDataFileName, MaxNumPeaksPerScan);
                    else
                        myMsDataFile = new ThermoRawFile(origDataFileName, MaxNumPeaksPerScan);
                    Status("Opening spectra file " + origDataFileName + "...");
                    myMsDataFile.Open();
                    listOfSortedms2Scans = myMsDataFile.Where(b => b.MsnOrder == 2).Select(b => new LocalMS2Scan(b)).OrderBy(b => b.PrecursorMass).ToArray();
                }

                var searchEngine = new ClassicSearchEngine(listOfSortedms2Scans, myMsDataFile.NumSpectra, spectraFileIndex, variableModifications, fixedModifications, proteinList, new Tolerance(ToleranceUnit.Absolute, ProductMassToleranceInDaltons), Protease, searchModes, MaxMissedCleavages, MaxModificationIsoforms, myMsDataFile.Name);

                var searchResults = (ClassicSearchResults)searchEngine.Run();

                for (int i = 0; i < searchModes.Count; i++)
                    allPsms[i].AddRange(searchResults.OuterPsms[i]);

                // Run analysis on single file results
                var analysisEngine = new AnalysisEngine(searchResults.OuterPsms, compactPeptideToProteinPeptideMatching, proteinList, variableModifications, fixedModifications, localizeableModifications, Protease, searchModes, myMsDataFile, new Tolerance(ToleranceUnit.Absolute, ProductMassToleranceInDaltons), (BinTreeStructure myTreeStructure, string s) => WriteTree(myTreeStructure, OutputFolder, Path.GetFileNameWithoutExtension(origDataFileName) + s), (List<NewPsmWithFdr> h, string s) => WritePsmsToTsv(h, OutputFolder, Path.GetFileNameWithoutExtension(origDataFileName) + s), null, false, MaxMissedCleavages, MaxModificationIsoforms);

                var analysisResults = (AnalysisResults)analysisEngine.Run();

                var identifications = analysisResults.AllResultingIdentifications[0];

                myMsDataFile.Close();

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

                int minMS1isotopicPeaksNeededForConfirmedIdentification = 3;
                int minMS2isotopicPeaksNeededForConfirmedIdentification = 2;
                int numFragmentsNeededForEveryIdentification = 10;
                double toleranceInMZforMS1Search = 0.01;

                // TODO: fix the tolerance calculation below
                var a = new CalibrationEngine(myMsDataFileForCalibration, randomSeed, ProductMassToleranceInDaltons * 2, identifications, minMS1isotopicPeaksNeededForConfirmedIdentification, minMS2isotopicPeaksNeededForConfirmedIdentification, numFragmentsNeededForEveryIdentification, toleranceInMZforMS1Search);

                var result = (CalibrationResults)a.Run();

                Status("Creating _indexedmzMLConnection, putting data in it, and writing!");
                var path = Path.Combine(OutputFolder, Path.GetFileNameWithoutExtension(origDataFileName) + "-Calibrated.mzML");
                MzmlMethods.CreateAndWriteMyIndexedMZmlwithCalibratedSpectra(result.MyMSDataFile, path);

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