using EngineLayer;
using EngineLayer.Analysis;
using EngineLayer.Calibration;
using EngineLayer.ClassicSearch;
using IO.MzML;
using IO.Thermo;
using MassSpectrometry;
using MzLibUtil;
using Proteomics;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using UsefulProteomicsDatabases;

namespace TaskLayer
{
    public class CalibrationTask : MetaMorpheusTask

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

            MaxDegreeOfParallelism = -1;
        }

        #endregion Public Constructors

        #region Public Properties

        public List<ModList> ListOfModListsFixed { get; set; }
        public List<ModList> ListOfModListsVariable { get; set; }
        public List<ModList> ListOfModListsLocalize { get; set; }
        public Tolerance ProductMassTolerance { get; set; }
        public Tolerance PrecursorMassTolerance { get; set; }
        public int MaxDegreeOfParallelism { get; set; }

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

        #region Protected Internal Methods

        protected internal void WriteMs1DataPoints(List<LabeledMs1DataPoint> items, string outputFolder, string fileName)
        {
            var writtenFile = Path.Combine(outputFolder, fileName + ".ms1dptsv");
            using (StreamWriter output = new StreamWriter(writtenFile))
            {
                output.WriteLine(LabeledMs1DataPoint.TabSeparatedHeaderForMs1);
                foreach (var dp in items)
                {
                    output.Write(string.Join("\t", dp.inputs));
                    output.WriteLine("\t" + dp.label);
                }
            }
            SucessfullyFinishedWritingFile(writtenFile);
        }

        protected internal void WriteMs2DataPoints(List<LabeledMs2DataPoint> items, string outputFolder, string fileName)
        {
            var writtenFile = Path.Combine(outputFolder, fileName + ".ms2dptsv");
            using (StreamWriter output = new StreamWriter(writtenFile))
            {
                output.WriteLine(LabeledMs2DataPoint.TabSeparatedHeaderForMs1);
                foreach (var dp in items)
                {
                    output.Write(string.Join("\t", dp.inputs));
                    output.WriteLine("\t" + dp.label);
                }
            }
            SucessfullyFinishedWritingFile(writtenFile);
        }

        #endregion Protected Internal Methods

        #region Protected Methods

        protected override MyResults RunSpecific()
        {
            var myTaskResults = new MyCalibrationTaskResults(this);
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
            List<ModificationWithMass> variableModifications = ListOfModListsVariable.SelectMany(b => b.Mods).OfType<ModificationWithMass>().ToList();
            List<ModificationWithMass> fixedModifications = ListOfModListsFixed.SelectMany(b => b.Mods).OfType<ModificationWithMass>().ToList();
            List<ModificationWithMass> localizeableModifications = ListOfModListsLocalize.SelectMany(b => b.Mods).OfType<ModificationWithMass>().ToList();

            Status("Loading proteins...");
            Dictionary<string, Modification> um;
            var proteinList = dbFilenameList.SelectMany(b => ProteinDbLoader.LoadProteinDb(b.FileName, true, localizeableModifications, b.IsContaminant, out um)).ToList();

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

            object lock1 = new object();
            object lock2 = new object();
            ParallelOptions parallelOptions = new ParallelOptions();
            parallelOptions.MaxDegreeOfParallelism = MaxDegreeOfParallelism;
            Parallel.For(0, currentRawFileList.Count, parallelOptions, spectraFileIndex =>
            {
                var compactPeptideToProteinPeptideMatching = new Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>>();
                var origDataFileName = currentRawFileList[spectraFileIndex];
                LocalMS2Scan[] listOfSortedms2Scans;
                IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> myMsDataFile;
                lock (lock1) // Lock because reading is sequential
                {
                    StartingDataFile(origDataFileName);
                    Status("Loading spectra file " + origDataFileName + "...");
                    if (Path.GetExtension(origDataFileName).Equals(".mzML"))
                        myMsDataFile = Mzml.LoadAllStaticData(origDataFileName);
                    else
                        myMsDataFile = ThermoStaticData.LoadAllStaticData(origDataFileName);
                    Status("Opening spectra file " + origDataFileName + "...");
                    listOfSortedms2Scans = GetMs2Scans(myMsDataFile).OrderBy(b => b.MonoisotopicPrecursorMass).ToArray();
                }

                var searchEngine = new ClassicSearchEngine(listOfSortedms2Scans, myMsDataFile.NumSpectra, variableModifications, fixedModifications, proteinList, ProductMassTolerance, Protease, searchModes, MaxMissedCleavages, MaxModificationIsoforms, origDataFileName, lp);

                var searchResults = (ClassicSearchResults)searchEngine.Run();
                myTaskResults.AddResultText(searchResults);

                for (int i = 0; i < searchModes.Count; i++)
                    allPsms[i].AddRange(searchResults.OuterPsms[i]);

                InitiatorMethionineBehavior initiatorMethionineBehavior = InitiatorMethionineBehavior.Variable;
                // Run analysis on single file results
                var analysisEngine = new AnalysisEngine(searchResults.OuterPsms, compactPeptideToProteinPeptideMatching, proteinList, variableModifications, fixedModifications, localizeableModifications, Protease, searchModes, myMsDataFile, ProductMassTolerance, (BinTreeStructure myTreeStructure, string s) => WriteTree(myTreeStructure, OutputFolder, Path.GetFileNameWithoutExtension(origDataFileName) + s), (List<NewPsmWithFdr> h, string s) => WritePsmsToTsv(h, OutputFolder, Path.GetFileNameWithoutExtension(origDataFileName) + s), null, false, MaxMissedCleavages, MaxModificationIsoforms, false, lp, double.NaN, initiatorMethionineBehavior);

                var analysisResults = (AnalysisResults)analysisEngine.Run();
                myTaskResults.AddResultText(analysisResults);

                var identifications = analysisResults.AllResultingIdentifications[0];

                //Now can calibrate!!!
                IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> myMsDataFileForCalibration;
                if (Path.GetExtension(origDataFileName).Equals(".mzML"))
                {
                    myMsDataFileForCalibration = Mzml.LoadAllStaticData(origDataFileName);
                }
                else
                {
                    myMsDataFileForCalibration = ThermoStaticData.LoadAllStaticData(origDataFileName);
                }
                int randomSeed = 1;

                int minMS1isotopicPeaksNeededForConfirmedIdentification = 3;
                int minMS2isotopicPeaksNeededForConfirmedIdentification = 2;
                int numFragmentsNeededForEveryIdentification = 10;

                // TODO: fix the tolerance calculation below

                var a = new CalibrationEngine(myMsDataFileForCalibration, randomSeed, ProductMassTolerance, identifications, minMS1isotopicPeaksNeededForConfirmedIdentification, minMS2isotopicPeaksNeededForConfirmedIdentification, numFragmentsNeededForEveryIdentification, PrecursorMassTolerance, fragmentTypesForCalibration, (List<LabeledMs1DataPoint> theList, string s) => WriteMs1DataPoints(theList, OutputFolder, Path.GetFileNameWithoutExtension(origDataFileName) + s), (List<LabeledMs2DataPoint> theList, string s) => WriteMs2DataPoints(theList, OutputFolder, Path.GetFileNameWithoutExtension(origDataFileName) + s));

                var result = a.Run();
                myTaskResults.AddResultText(result);

                if (result is MyErroredResults)
                {
                    Warn(a.ToString());
                }
                else
                {
                    Status("Creating _indexedmzMLConnection, putting data in it, and writing!");
                    var path = Path.Combine(OutputFolder, Path.GetFileNameWithoutExtension(origDataFileName) + "-Calibrated.mzML");
                    lock (lock2) // Lock because writing is sequential
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