using EngineLayer;
using EngineLayer.Analysis;
using EngineLayer.Calibration;
using EngineLayer.ClassicSearch;
using IO.MzML;
using IO.Thermo;
using MassSpectrometry;
using Proteomics;
using Spectra;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace TaskLayer
{
    public class CalibrationTask : MyTaskEngine
    {

        #region Public Constructors

        public CalibrationTask()
        {
            // Set default values here:
            MaxMissedCleavages = 2;
            Protease = ProteaseDictionary.Instance["trypsin"];
            MaxModificationIsoforms = 4096;
            InitiatorMethionineBehavior = InitiatorMethionineBehavior.Variable;
            ProductMassTolerance = new Tolerance(ToleranceUnit.Absolute, 0.01);
            PrecursorMassTolerance = new Tolerance(ToleranceUnit.PPM, 10);
            BIons = true;
            YIons = true;

            ListOfModListsFixed = new List<ModList> { AllModLists.First(b => b.FileName.EndsWith("f.txt")) };
            ListOfModListsVariable = new List<ModList> { AllModLists.First(b => b.FileName.EndsWith("v.txt")) };
            ListOfModListsLocalize = new List<ModList> { AllModLists.First(b => b.FileName.EndsWith("ptmlist.txt")) };

            TaskType = MyTask.Calibrate;
            MaxNumPeaksPerScan = 400;
        }

        #endregion Public Constructors

        #region Public Properties

        public List<ModList> ListOfModListsFixed { get; set; }
        public List<ModList> ListOfModListsVariable { get; set; }
        public List<ModList> ListOfModListsLocalize { get; set; }
        public Tolerance ProductMassTolerance { get; set; }
        public Tolerance PrecursorMassTolerance { get; set; }

        #endregion Public Properties

        #region Protected Properties

        protected override string SpecificTaskInfo
        {
            get
            {
                var sb = new StringBuilder();
                sb.AppendLine("Fixed mod lists: " + string.Join(",", ListOfModListsFixed.Select(b => b.FileName)));
                sb.AppendLine("Variable mod lists: " + string.Join(",", ListOfModListsVariable.Select(b => b.FileName)));
                sb.AppendLine("Localized mod lists: " + string.Join(",", ListOfModListsLocalize.Select(b => b.FileName)));
                sb.AppendLine("PrecursorMassTolerance: " + PrecursorMassTolerance);
                sb.Append("ProductMassTolerance: " + ProductMassTolerance);
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

            SearchMode searchMode;
            if (PrecursorMassTolerance.Unit == ToleranceUnit.PPM)
                searchMode = new SinglePpmAroundZeroSearchMode(PrecursorMassTolerance.Value);
            else
                searchMode = new SingleAbsoluteAroundZeroSearchMode(PrecursorMassTolerance.Value);
            var searchModes = new List<SearchMode> { searchMode };

            List<PsmParent>[] allPsms = new List<PsmParent>[1];
            allPsms[0] = new List<PsmParent>();

            Status("Loading modifications...");
            List<MetaMorpheusModification> variableModifications = ListOfModListsVariable.SelectMany(b => b.Mods).ToList();
            List<MetaMorpheusModification> fixedModifications = ListOfModListsFixed.SelectMany(b => b.Mods).ToList();
            List<MetaMorpheusModification> localizeableModifications = ListOfModListsLocalize.SelectMany(b => b.Mods).ToList();

            Dictionary<string, List<MetaMorpheusModification>> identifiedModsInXML;
            HashSet<string> unidentifiedModStrings;
            MatchXMLmodsToKnownMods(dbFilenameList, localizeableModifications, out identifiedModsInXML, out unidentifiedModStrings);

            Status("Loading proteins...");
            var proteinList = dbFilenameList.SelectMany(b => GetProteins(true, identifiedModsInXML, b)).ToList();

            List<ProductType> lp = new List<ProductType>();
            FragmentTypes fragmentTypesForCalibration = FragmentTypes.None;
            if (BIons)
            {
                fragmentTypesForCalibration = fragmentTypesForCalibration | FragmentTypes.b;
                lp.Add(ProductType.B);
            }
            if (YIons)
            {
                fragmentTypesForCalibration = fragmentTypesForCalibration | FragmentTypes.y;
                lp.Add(ProductType.Y);
            }

            Parallel.For(0, currentRawFileList.Count, spectraFileIndex =>
            {
                var compactPeptideToProteinPeptideMatching = new Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>>();
                var origDataFileName = currentRawFileList[spectraFileIndex];
                LocalMS2Scan[] listOfSortedms2Scans;
                IMsDataFile<IMzSpectrum<MzPeak>> myMsDataFile;
                lock (myTaskResults)
                {
                    StartingDataFile(origDataFileName);
                    Status("Loading spectra file " + origDataFileName + "...");
                    if (Path.GetExtension(origDataFileName).Equals(".mzML"))
                        myMsDataFile = new Mzml(origDataFileName, MaxNumPeaksPerScan);
                    else
                        myMsDataFile = new ThermoRawFile(origDataFileName, MaxNumPeaksPerScan);
                    Status("Opening spectra file " + origDataFileName + "...");
                    myMsDataFile.Open();
                    listOfSortedms2Scans = GetMs2Scans(myMsDataFile).OrderBy(b => b.PrecursorMass).ToArray();
                }

                var searchEngine = new ClassicSearchEngine(listOfSortedms2Scans, myMsDataFile.NumSpectra, variableModifications, fixedModifications, proteinList, ProductMassTolerance, Protease, searchModes, MaxMissedCleavages, MaxModificationIsoforms, myMsDataFile.Name, lp);

                var searchResults = (ClassicSearchResults)searchEngine.Run();

                for (int i = 0; i < searchModes.Count; i++)
                    allPsms[i].AddRange(searchResults.OuterPsms[i]);

                // Run analysis on single file results
                var analysisEngine = new AnalysisEngine(searchResults.OuterPsms, compactPeptideToProteinPeptideMatching, proteinList, variableModifications, fixedModifications, localizeableModifications, Protease, searchModes, myMsDataFile, ProductMassTolerance, (BinTreeStructure myTreeStructure, string s) => WriteTree(myTreeStructure, OutputFolder, Path.GetFileNameWithoutExtension(origDataFileName) + s), (List<NewPsmWithFdr> h, string s) => WritePsmsToTsv(h, OutputFolder, Path.GetFileNameWithoutExtension(origDataFileName) + s), null, false, MaxMissedCleavages, MaxModificationIsoforms, false, lp, double.NaN);

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

                // TODO: fix the tolerance calculation below

                var a = new CalibrationEngine(myMsDataFileForCalibration, randomSeed, ProductMassTolerance, identifications, minMS1isotopicPeaksNeededForConfirmedIdentification, minMS2isotopicPeaksNeededForConfirmedIdentification, numFragmentsNeededForEveryIdentification, PrecursorMassTolerance, fragmentTypesForCalibration);

                var result = a.Run();

                if (result is MyErroredResults)
                {
                    Warn(a.ToString());
                }
                else
                {
                    Status("Creating _indexedmzMLConnection, putting data in it, and writing!");
                    var path = Path.Combine(OutputFolder, Path.GetFileNameWithoutExtension(origDataFileName) + "-Calibrated.mzML");
                    lock (myTaskResults)
                    {
                        MzmlMethods.CreateAndWriteMyIndexedMZmlwithCalibratedSpectra(((CalibrationResults)result).MyMSDataFile, path);

                        SucessfullyFinishedWritingFile(path);

                        myTaskResults.newSpectra.Add(path);
                    }
                }
                FinishedDataFile(origDataFileName);
            }
            );
            return myTaskResults;
        }

        #endregion Protected Methods

    }
}