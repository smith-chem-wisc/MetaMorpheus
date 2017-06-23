using EngineLayer;
using EngineLayer.Analysis;
using EngineLayer.Calibration;
using EngineLayer.ClassicSearch;
using IO.MzML;
using IO.Thermo;
using MassSpectrometry;
using MzLibUtil;
using Proteomics;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace TaskLayer
{
    public class CalibrationTask : MetaMorpheusTask
    {

        #region Public Constructors

        public CalibrationTask() : base(MyTask.Calibrate)
        {
            // Set default values here:
            MaxMissedCleavages = 2;
            MinPeptideLength = null;
            MaxPeptideLength = null;
            Protease = GlobalTaskLevelSettings.ProteaseDictionary["trypsin"];
            MaxModificationIsoforms = 4096;
            InitiatorMethionineBehavior = InitiatorMethionineBehavior.Variable;
            ProductMassTolerance = new Tolerance(ToleranceUnit.Absolute, 0.01);
            PrecursorMassTolerance = new Tolerance(ToleranceUnit.PPM, 10);
            BIons = true;
            YIons = true;
            CIons = false;
            ZdotIons = false;

            LocalizeAll = true;

            NonLinearCalibration = true;

            ListOfModsVariable = new List<Tuple<string, string>> { new Tuple<string, string>("Common Variable", "Oxidation of M") };
            ListOfModsFixed = new List<Tuple<string, string>> { new Tuple<string, string>("Common Fixed", "Carbamidomethyl of C") };
            ListOfModsLocalize = GlobalTaskLevelSettings.AllModsKnown.Select(b => new Tuple<string, string>(b.modificationType, b.id)).ToList();

            MaxDegreeOfParallelism = -1;
        }

        #endregion Public Constructors

        #region Public Properties

        public static List<string> AllModLists { get; private set; }

        public InitiatorMethionineBehavior InitiatorMethionineBehavior { get; set; }

        public int MaxMissedCleavages { get; set; }

        public int? MinPeptideLength { get; set; }

        public int? MaxPeptideLength { get; set; }

        public int MaxModificationIsoforms { get; set; }

        public Protease Protease { get; set; }

        public bool BIons { get; set; }

        public bool YIons { get; set; }

        public bool ZdotIons { get; set; }

        public bool CIons { get; set; }
        public List<Tuple<string, string>> ListOfModsFixed { get; set; }
        public List<Tuple<string, string>> ListOfModsVariable { get; set; }
        public List<Tuple<string, string>> ListOfModsLocalize { get; set; }
        public Tolerance ProductMassTolerance { get; set; }
        public Tolerance PrecursorMassTolerance { get; set; }
        public int MaxDegreeOfParallelism { get; set; }
        public bool LocalizeAll { get; set; }
        public bool NonLinearCalibration { get; set; }

        #endregion Public Properties

        #region Public Methods

        public override string ToString()
        {
            var sb = new StringBuilder();
            sb.AppendLine(TaskType.ToString());
            sb.AppendLine("The initiator methionine behavior is set to "
                + InitiatorMethionineBehavior
                + " and the maximum number of allowed missed cleavages is "
                + MaxMissedCleavages);
            sb.AppendLine("MinPeptideLength: " + MinPeptideLength);
            sb.AppendLine("MaxPeptideLength: " + MaxPeptideLength);
            sb.AppendLine("maxModificationIsoforms: " + MaxModificationIsoforms);
            sb.AppendLine("protease: " + Protease);
            sb.AppendLine("bIons: " + BIons);
            sb.AppendLine("yIons: " + YIons);
            sb.AppendLine("cIons: " + CIons);
            sb.AppendLine("zdotIons: " + ZdotIons);
            //sb.AppendLine("Fixed mod lists: " + string.Join(",", ListOfModListsFixed));
            //sb.AppendLine("Variable mod lists: " + string.Join(",", ListOfModListsVariable));
            //sb.AppendLine("Localized mod lists: " + string.Join(",", ListOfModListsLocalize));
            sb.AppendLine("PrecursorMassTolerance: " + PrecursorMassTolerance);
            sb.Append("ProductMassTolerance: " + ProductMassTolerance);
            return sb.ToString();
        }

        #endregion Public Methods

        #region Protected Internal Methods

        protected internal void WriteMs1DataPoints(List<LabeledMs1DataPoint> items, string outputFolder, string fileName, List<string> nestedIDs)
        {
            var writtenFile = Path.Combine(outputFolder, fileName + ".ms1dptsv");
            using (StreamWriter output = new StreamWriter(writtenFile))
            {
                output.WriteLine(LabeledMs1DataPoint.TabSeparatedHeaderForMs1 + "\t" + NewPsmWithFdr.TabSeparatedHeader);
                foreach (var dp in items)
                {
                    output.Write(string.Join("\t", dp.Inputs));
                    output.Write("\t" + dp.Label);
                    output.WriteLine("\t" + dp.identification);
                }
            }
            SucessfullyFinishedWritingFile(writtenFile, nestedIDs);
        }

        protected internal void WriteMs2DataPoints(List<LabeledMs2DataPoint> items, string outputFolder, string fileName, List<string> nestedIDs)
        {
            var writtenFile = Path.Combine(outputFolder, fileName + ".ms2dptsv");
            using (StreamWriter output = new StreamWriter(writtenFile))
            {
                output.WriteLine(LabeledMs2DataPoint.TabSeparatedHeaderForMs1 + "\t" + NewPsmWithFdr.TabSeparatedHeader);
                foreach (var dp in items)
                {
                    output.Write(string.Join("\t", dp.Inputs));
                    output.Write("\t" + dp.Label);
                    output.WriteLine("\t" + dp.identification);
                }
            }
            SucessfullyFinishedWritingFile(writtenFile, nestedIDs);
        }

        #endregion Protected Internal Methods

        #region Protected Methods

        protected override MyTaskResults RunSpecific(string OutputFolder, List<DbForTask> dbFilenameList, List<string> currentRawFileList, string taskId)
        {
            myTaskResults = new MyTaskResults(this)
            {
                newSpectra = new List<string>()
            };
            SearchMode searchMode;
            if (PrecursorMassTolerance.Unit == ToleranceUnit.PPM)
                searchMode = new SinglePpmAroundZeroSearchMode(PrecursorMassTolerance.Value);
            else
                searchMode = new SingleAbsoluteAroundZeroSearchMode(PrecursorMassTolerance.Value);
            var searchModes = new List<SearchMode> { searchMode };

            Status("Loading modifications...", new List<string> { taskId });
            List<ModificationWithMass> variableModifications = GlobalTaskLevelSettings.AllModsKnown.OfType<ModificationWithMass>().Where(b => ListOfModsVariable.Contains(new Tuple<string, string>(b.modificationType, b.id))).ToList();
            List<ModificationWithMass> fixedModifications = GlobalTaskLevelSettings.AllModsKnown.OfType<ModificationWithMass>().Where(b => ListOfModsFixed.Contains(new Tuple<string, string>(b.modificationType, b.id))).ToList();

            List<ModificationWithMass> localizeableModifications;
            if (LocalizeAll)
                localizeableModifications = GlobalTaskLevelSettings.AllModsKnown.OfType<ModificationWithMass>().ToList();
            else
                localizeableModifications = GlobalTaskLevelSettings.AllModsKnown.OfType<ModificationWithMass>().Where(b => ListOfModsLocalize.Contains(new Tuple<string, string>(b.modificationType, b.id))).ToList();

            Dictionary<ModificationWithMass, ushort> modsDictionary = new Dictionary<ModificationWithMass, ushort>();
            foreach (var mod in fixedModifications)
                modsDictionary.Add(mod, 0);
            int i = 1;
            foreach (var mod in variableModifications)
            {
                modsDictionary.Add(mod, (ushort)i);
                i++;
            }
            foreach (var mod in localizeableModifications)
            {
                if (!modsDictionary.ContainsKey(mod))
                    modsDictionary.Add(mod, (ushort)i);
                i++;
            }

            Status("Loading proteins...", new List<string> { taskId });
            Dictionary<string, Modification> um;
            var proteinList = dbFilenameList.SelectMany(b => LoadProteinDb(b.FileName, true, localizeableModifications, b.IsContaminant, out um)).ToList();

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
            if (CIons)
            {
                fragmentTypesForCalibration = fragmentTypesForCalibration | FragmentTypes.c;
                lp.Add(ProductType.C);
            }
            if (ZdotIons)
            {
                fragmentTypesForCalibration = fragmentTypesForCalibration | FragmentTypes.zdot;
                lp.Add(ProductType.Zdot);
            }

            object lock1 = new object();
            object lock2 = new object();
            ParallelOptions parallelOptions = new ParallelOptions()
            {
                MaxDegreeOfParallelism = MaxDegreeOfParallelism
            };
            Status("Calibrating...", new List<string> { taskId });

            Parallel.For(0, currentRawFileList.Count, parallelOptions, spectraFileIndex =>
            {
                var compactPeptideToProteinPeptideMatching = new Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>>();
                var origDataFile = currentRawFileList[spectraFileIndex];
                Ms2ScanWithSpecificMass[] listOfSortedms2Scans;
                IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> myMsDataFile;
                lock (lock1) // Lock because reading is sequential
                {
                    NewCollection(Path.GetFileName(origDataFile), new List<string> { taskId, "Individual Searches", origDataFile });
                    StartingDataFile(origDataFile, new List<string> { taskId, "Individual Searches", origDataFile });
                    Status("Loading spectra file " + origDataFile + "...", new List<string> { taskId, "Individual Searches", origDataFile });
                    if (Path.GetExtension(origDataFile).Equals(".mzML"))
                        myMsDataFile = Mzml.LoadAllStaticData(origDataFile);
                    else
                        myMsDataFile = ThermoStaticData.LoadAllStaticData(origDataFile);
                    Status("Getting ms2 scans...", new List<string> { taskId, "Individual Searches", origDataFile });
                    bool findAllPrecursors = true;
                    bool useProvidedPrecursorInfo = true;
                    var intensityRatio = 4;
                    listOfSortedms2Scans = MetaMorpheusEngine.GetMs2Scans(myMsDataFile, findAllPrecursors, useProvidedPrecursorInfo, intensityRatio, origDataFile).OrderBy(b => b.PrecursorMass).ToArray();
                }

                StringBuilder sbForThisFile = new StringBuilder();
                sbForThisFile.AppendLine(origDataFile);

                var searchEngine = new ClassicSearchEngine(listOfSortedms2Scans, variableModifications, fixedModifications, proteinList, ProductMassTolerance, Protease, searchModes, MaxMissedCleavages, MinPeptideLength, MaxPeptideLength, MaxModificationIsoforms, lp, new List<string> { taskId, "Individual Searches", origDataFile }, false);

                var searchResults = (ClassicSearchResults)searchEngine.Run();

                InitiatorMethionineBehavior initiatorMethionineBehavior = InitiatorMethionineBehavior.Variable;
                // Run analysis on single file results
                var analysisEngine = new AnalysisEngine(searchResults.OuterPsms, compactPeptideToProteinPeptideMatching, proteinList, variableModifications, fixedModifications, Protease, searchModes, listOfSortedms2Scans, ProductMassTolerance, (BinTreeStructure myTreeStructure, string s) => WriteTree(myTreeStructure, OutputFolder, Path.GetFileNameWithoutExtension(origDataFile) + "_" + s, new List<string> { taskId, "Individual Searches", origDataFile }), (List<NewPsmWithFdr> h, string s, List<string> ss) => WritePsmsToTsv(h, OutputFolder, Path.GetFileNameWithoutExtension(origDataFile) + s, ss), null, (List<NewPsmWithFdr> h, List<ProteinGroup> g, SearchMode m, string s, List<string> ss) => WriteMzidentml(h, g, variableModifications, fixedModifications, new List<Protease> { Protease }, 0.01, m, ProductMassTolerance, MaxMissedCleavages, OutputFolder, Path.GetFileNameWithoutExtension(origDataFile) + s, ss), false, false, false, MaxMissedCleavages, MinPeptideLength, MaxPeptideLength, MaxModificationIsoforms, false, lp, double.NaN, initiatorMethionineBehavior, new List<string> { taskId, "Individual Searches", origDataFile }, false, false, 0, modsDictionary, myMsDataFile, null);

                var analysisResults = (AnalysisResults)analysisEngine.Run();

                var goodIdentifications = analysisResults.AllResultingIdentifications[0].Where(b => b.QValue < 0.01 && !b.IsDecoy).ToList();

                //Now can calibrate!!!

                int minMS1isotopicPeaksNeededForConfirmedIdentification = 3;
                int minMS2isotopicPeaksNeededForConfirmedIdentification = 2;
                int numFragmentsNeededForEveryIdentification = 10;

                // TODO: fix the tolerance calculation below

                var a = new CalibrationEngine(myMsDataFile, ProductMassTolerance, goodIdentifications, minMS1isotopicPeaksNeededForConfirmedIdentification, minMS2isotopicPeaksNeededForConfirmedIdentification, numFragmentsNeededForEveryIdentification, PrecursorMassTolerance, fragmentTypesForCalibration, (List<LabeledMs1DataPoint> theList, string s) => WriteMs1DataPoints(theList, OutputFolder, Path.GetFileNameWithoutExtension(origDataFile) + s, new List<string> { taskId, "Individual Searches", origDataFile }), (List<LabeledMs2DataPoint> theList, string s) => WriteMs2DataPoints(theList, OutputFolder, Path.GetFileNameWithoutExtension(origDataFile) + s, new List<string> { taskId, "Individual Searches", origDataFile }), false, new List<string> { taskId, "Individual Searches", origDataFile });

                MetaMorpheusEngineResults theResult = a.Run();

                CalibrationResults calibrationResult = theResult as CalibrationResults;

                if (calibrationResult == null)
                {
                    sbForThisFile.AppendLine(theResult.ToString());
                    Warn(theResult.ToString(), new List<string> { taskId, "Individual Searches", origDataFile });
                    Status("Errored", new List<string> { taskId, "Individual Searches", origDataFile });
                    return;
                }

                sbForThisFile.AppendLine("\t" + "Before Calib: MS1 (Th): " + calibrationResult.ms1meanSds.First());
                sbForThisFile.AppendLine("\t" + "Before Calib: MS2 (Th): " + calibrationResult.ms2meanSds.First());
                sbForThisFile.AppendLine("\t" + "After Linear: MS1 (Th): " + calibrationResult.ms1meanSds.Last());
                sbForThisFile.AppendLine("\t" + "After Linear: MS2 (Th): " + calibrationResult.ms2meanSds.Last());

                if (NonLinearCalibration)
                {
                    // Second search round

                    Status("Getting ms2 scans for second round...", new List<string> { taskId, "Individual Searches", origDataFile });
                    bool useProvidedPrecursorInfo = true;
                    bool findAllPrecursors = true;
                    var intensityRatio = 4;
                    var listOfSortedms2ScansTest = MetaMorpheusEngine.GetMs2Scans(myMsDataFile, findAllPrecursors, useProvidedPrecursorInfo, intensityRatio, origDataFile).OrderBy(b => b.PrecursorMass).ToArray();
                    var searchEngineTest = new ClassicSearchEngine(listOfSortedms2ScansTest, variableModifications, fixedModifications, proteinList, ProductMassTolerance, Protease, searchModes, MaxMissedCleavages, MinPeptideLength, MaxPeptideLength, MaxModificationIsoforms, lp, new List<string> { taskId, "Individual Searches", origDataFile }, false);
                    var searchResultsTest = (ClassicSearchResults)searchEngineTest.Run();
                    var analysisEngineTest = new AnalysisEngine(searchResultsTest.OuterPsms, compactPeptideToProteinPeptideMatching, proteinList, variableModifications, fixedModifications, Protease, searchModes, listOfSortedms2ScansTest, ProductMassTolerance, (BinTreeStructure myTreeStructure, string s) => WriteTree(myTreeStructure, OutputFolder, Path.GetFileNameWithoutExtension(origDataFile) + "_" + s + "test", new List<string> { taskId, "Individual Searches", origDataFile }), (List<NewPsmWithFdr> h, string s, List<string> ss) => WritePsmsToTsv(h, OutputFolder, Path.GetFileNameWithoutExtension(origDataFile) + "_" + s + "test", ss), null, (List<NewPsmWithFdr> h, List<ProteinGroup> g, SearchMode m, string s, List<string> ss) => WriteMzidentml(h, g, variableModifications, fixedModifications, new List<Protease> { Protease }, 0.01, m, ProductMassTolerance, MaxMissedCleavages, OutputFolder, Path.GetFileNameWithoutExtension(origDataFile) + "_" + s + "test", ss), false, false, false, MaxMissedCleavages, MinPeptideLength, MaxPeptideLength, MaxModificationIsoforms, false, lp, double.NaN, initiatorMethionineBehavior, new List<string> { taskId, "Individual Searches", origDataFile }, false, false, 0, modsDictionary, myMsDataFile, null);
                    var analysisResultsTest = (AnalysisResults)analysisEngineTest.Run();

                    theResult = new CalibrationEngine(myMsDataFile, ProductMassTolerance, analysisResultsTest.AllResultingIdentifications[0].Where(b => b.QValue < 0.01 && !b.IsDecoy).ToList(), minMS1isotopicPeaksNeededForConfirmedIdentification, minMS2isotopicPeaksNeededForConfirmedIdentification, numFragmentsNeededForEveryIdentification, PrecursorMassTolerance, fragmentTypesForCalibration, (List<LabeledMs1DataPoint> theList, string s) => WriteMs1DataPoints(theList, OutputFolder, Path.GetFileNameWithoutExtension(origDataFile) + s + "after", new List<string> { taskId, "Individual Searches", origDataFile }), (List<LabeledMs2DataPoint> theList, string s) => WriteMs2DataPoints(theList, OutputFolder, Path.GetFileNameWithoutExtension(origDataFile) + s + "after", new List<string> { taskId, "Individual Searches", origDataFile }), true, new List<string> { taskId, "Individual Searches", origDataFile }).Run() as CalibrationResults;

                    calibrationResult = theResult as CalibrationResults;

                    if (calibrationResult == null)
                    {
                        sbForThisFile.AppendLine(theResult.ToString());
                        Warn(theResult.ToString(), new List<string> { taskId, "Individual Searches", origDataFile });
                        Status("Errored", new List<string> { taskId, "Individual Searches", origDataFile });
                        return;
                    }

                    sbForThisFile.AppendLine("\t" + "After NonLin: MS1 (Th): " + calibrationResult.ms1meanSds.Last());
                    sbForThisFile.AppendLine("\t" + "After NonLin: MS2 (Th): " + calibrationResult.ms2meanSds.Last());
                }
                myTaskResults.AddNiceText(sbForThisFile.ToString());

                Status("Writing mzML!", new List<string> { taskId, "Individual Searches", origDataFile });
                var path = Path.Combine(OutputFolder, Path.GetFileNameWithoutExtension(origDataFile) + "-Calibrated.mzML");
                lock (lock2) // Lock because writing is sequential
                {
                    MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(myMsDataFile, path, false);

                    SucessfullyFinishedWritingFile(path, new List<string> { taskId, "Individual Searches", origDataFile });

                    myTaskResults.newSpectra.Add(path);
                }
                FinishedDataFile(origDataFile, new List<string> { taskId, "Individual Searches", origDataFile });
                ReportProgress(new ProgressEventArgs(100, "Done!", new List<string> { taskId, "Individual Searches", origDataFile }));
            }
            );
            ReportProgress(new ProgressEventArgs(100, "Done!", new List<string> { taskId, "Individual Searches" }));

            return myTaskResults;
        }

        #endregion Protected Methods

    }
}