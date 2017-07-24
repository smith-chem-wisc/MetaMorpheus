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
            MinPeptideLength = 5;
            MaxPeptideLength = null;
            Protease = GlobalTaskLevelSettings.ProteaseDictionary["trypsin"];
            MaxModificationIsoforms = 4096;
            InitiatorMethionineBehavior = InitiatorMethionineBehavior.Variable;
            ProductMassTolerance = new AbsoluteTolerance(0.01);
            PrecursorMassTolerance = new PpmTolerance(10);
            BIons = true;
            YIons = true;
            CIons = false;
            ZdotIons = false;

            LocalizeAll = true;

            NonLinearCalibration = true;

            ListOfModsVariable = new List<Tuple<string, string>> { new Tuple<string, string>("Common Variable", "Oxidation of M") };
            ListOfModsFixed = new List<Tuple<string, string>> { new Tuple<string, string>("Common Fixed", "Carbamidomethyl of C") };
            ListOfModsLocalize = GlobalTaskLevelSettings.AllModsKnown.Select(b => new Tuple<string, string>(b.modificationType, b.id)).ToList();

            MaxDegreeOfParallelism = 1;
            ConserveMemory = false;

            WriteIntermediateFiles = false;

            DoPrecursorDeconvolution = true;
            UseProvidedPrecursorInfo = true;
            DeconvolutionIntensityRatio = 4;
            DeconvolutionMaxAssumedChargeState = 10;
            DeconvolutionMassTolerance = new PpmTolerance(5);
        }

        #endregion Public Constructors

        #region Public Properties

        public static List<string> AllModLists { get; private set; }

        public InitiatorMethionineBehavior InitiatorMethionineBehavior { get; set; }

        public int MaxMissedCleavages { get; set; }
        public bool ConserveMemory { get; set; }

        public int? MinPeptideLength { get; set; }

        public int? MaxPeptideLength { get; set; }

        public int MaxModificationIsoforms { get; set; }

        public Protease Protease { get; set; }

        public bool BIons { get; set; }

        public bool YIons { get; set; }

        public bool ZdotIons { get; set; }

        public bool CIons { get; set; }
        public Tolerance ProductMassTolerance { get; set; }
        public Tolerance PrecursorMassTolerance { get; set; }
        public bool NonLinearCalibration { get; set; }
        public bool WriteIntermediateFiles { get; set; }

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
                output.WriteLine(LabeledMs1DataPoint.TabSeparatedHeaderForMs1 + "\t" + Psm.GetTabSeparatedHeader());
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
                output.WriteLine(LabeledMs2DataPoint.TabSeparatedHeaderForMs1 + "\t" + Psm.GetTabSeparatedHeader());
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
            MassDiffAcceptor searchMode;
            if (PrecursorMassTolerance is PpmTolerance)
                searchMode = new SinglePpmAroundZeroSearchMode(PrecursorMassTolerance.Value);
            else
                searchMode = new SingleAbsoluteAroundZeroSearchMode(PrecursorMassTolerance.Value);
            var searchModes = new List<MassDiffAcceptor> { searchMode };

            Status("Loading modifications...", new List<string> { taskId });
            List<ModificationWithMass> variableModifications = GlobalTaskLevelSettings.AllModsKnown.OfType<ModificationWithMass>().Where(b => ListOfModsVariable.Contains(new Tuple<string, string>(b.modificationType, b.id))).ToList();
            List<ModificationWithMass> fixedModifications = GlobalTaskLevelSettings.AllModsKnown.OfType<ModificationWithMass>().Where(b => ListOfModsFixed.Contains(new Tuple<string, string>(b.modificationType, b.id))).ToList();
            List<ModificationWithMass> localizeableModifications;
            if (LocalizeAll)
                localizeableModifications = GlobalTaskLevelSettings.AllModsKnown.OfType<ModificationWithMass>().ToList();
            else
                localizeableModifications = GlobalTaskLevelSettings.AllModsKnown.OfType<ModificationWithMass>().Where(b => ListOfModsLocalize.Contains(new Tuple<string, string>(b.modificationType, b.id))).ToList();

            Status("Loading proteins...", new List<string> { taskId });
            Dictionary<string, Modification> um;
            var proteinList = dbFilenameList.SelectMany(b => LoadProteinDb(b.FilePath, true, localizeableModifications, b.IsContaminant, out um)).ToList();

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
            ParallelOptions parallelOptions = new ParallelOptions();
            if (MaxDegreeOfParallelism.HasValue)
                parallelOptions.MaxDegreeOfParallelism = MaxDegreeOfParallelism.Value;
            Status("Calibrating...", new List<string> { taskId });

            Parallel.For(0, currentRawFileList.Count, parallelOptions, spectraFileIndex =>
            {
                var origDataFile = currentRawFileList[spectraFileIndex];
                Ms2ScanWithSpecificMass[] listOfSortedms2Scans;
                IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> myMsDataFile;
                lock (lock1) // Lock because reading is sequential
                {
                    NewCollection(Path.GetFileName(origDataFile), new List<string> { taskId, "Individual Spectra Files", origDataFile });
                    StartingDataFile(origDataFile, new List<string> { taskId, "Individual Spectra Files", origDataFile });
                    Status("Loading spectra file " + origDataFile + "...", new List<string> { taskId, "Individual Spectra Files", origDataFile });
                    if (Path.GetExtension(origDataFile).Equals(".mzML"))
                        myMsDataFile = Mzml.LoadAllStaticData(origDataFile);
                    else
                        myMsDataFile = ThermoStaticData.LoadAllStaticData(origDataFile);
                }

                Status("Getting ms2 scans...", new List<string> { taskId, "Individual Spectra Files", origDataFile });
                listOfSortedms2Scans = GetMs2Scans(myMsDataFile, origDataFile, DoPrecursorDeconvolution, UseProvidedPrecursorInfo, DeconvolutionIntensityRatio, DeconvolutionMaxAssumedChargeState, DeconvolutionMassTolerance).OrderBy(b => b.PrecursorMass).ToArray();

                StringBuilder sbForThisFile = new StringBuilder();
                sbForThisFile.AppendLine(origDataFile);

                int minMS1isotopicPeaksNeededForConfirmedIdentification = 3;
                int minMS2isotopicPeaksNeededForConfirmedIdentification = 2;
                int numFragmentsNeededForEveryIdentification = 10;

                // Linear calibration
                {
                    var searchEngine = new ClassicSearchEngine(listOfSortedms2Scans, variableModifications, fixedModifications, proteinList, ProductMassTolerance, Protease, searchModes, MaxMissedCleavages, MinPeptideLength, MaxPeptideLength, MaxModificationIsoforms, lp, new List<string> { taskId, "Individual Spectra Files", origDataFile }, ConserveMemory, InitiatorMethionineBehavior, false);

                    var searchResults = (SearchResults)searchEngine.Run();
                    List<Psm>[] allPsms = new List<Psm>[1];
                    allPsms[0] = searchResults.Psms[0].ToList();
                    // Group and order psms

                    SequencesToActualProteinPeptidesEngine sequencesToActualProteinPeptidesEngine = new SequencesToActualProteinPeptidesEngine(allPsms, proteinList, searchModes, Protease, MaxMissedCleavages, MinPeptideLength, MaxPeptideLength, InitiatorMethionineBehavior, fixedModifications, variableModifications, MaxModificationIsoforms, new List<string> { taskId, "Individual Spectra Files", origDataFile });
                    var res = (SequencesToActualProteinPeptidesEngineResults)sequencesToActualProteinPeptidesEngine.Run();
                    Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>> compactPeptideToProteinPeptideMatching = res.CompactPeptideToProteinPeptideMatching;

                    allPsms[0] = allPsms[0].Where(b => b != null).OrderByDescending(b => b.Score).ThenBy(b => Math.Abs(b.ScanPrecursorMass - b.MostProbableProteinInfo.PeptideMonoisotopicMass)).GroupBy(b => new Tuple<string, int, double>(b.FullFilePath, b.ScanNumber, b.MostProbableProteinInfo.PeptideMonoisotopicMass)).Select(b => b.First()).ToList();

                    // Run FDR analysis on single file results
                    new FdrAnalysisEngine(allPsms, searchModes, new List<string> { taskId, "Individual Spectra Files", origDataFile }).Run();

                    if (WriteIntermediateFiles)
                        WritePsmsToTsv(allPsms[0], OutputFolder, Path.GetFileNameWithoutExtension(origDataFile) + "PSMsBeforeLinearCalib", new List<string> { taskId, "Individual Spectra Files", origDataFile });

                    var goodIdentifications = allPsms[0].Where(b => b.FdrInfo.QValue < 0.01 && !b.MostProbableProteinInfo.IsDecoy).ToList();

                    Action<List<LabeledMs1DataPoint>, string> ms1Action = (List<LabeledMs1DataPoint> theList, string s) => {; };
                    Action<List<LabeledMs2DataPoint>, string> ms2Action = (List<LabeledMs2DataPoint> theList, string s) => {; };
                    if (WriteIntermediateFiles)
                    {
                        ms1Action = (List<LabeledMs1DataPoint> theList, string s) => WriteMs1DataPoints(theList, OutputFolder, Path.GetFileNameWithoutExtension(origDataFile) + "inLinear" + s, new List<string> { taskId, "Individual Spectra Files", origDataFile });
                        ms2Action = (List<LabeledMs2DataPoint> theList, string s) => WriteMs2DataPoints(theList, OutputFolder, Path.GetFileNameWithoutExtension(origDataFile) + "inLinear" + s, new List<string> { taskId, "Individual Spectra Files", origDataFile });
                    }

                    MetaMorpheusEngineResults theResult = new CalibrationEngine(myMsDataFile, ProductMassTolerance, goodIdentifications, minMS1isotopicPeaksNeededForConfirmedIdentification, minMS2isotopicPeaksNeededForConfirmedIdentification, numFragmentsNeededForEveryIdentification, PrecursorMassTolerance, fragmentTypesForCalibration, ms1Action, ms2Action, false, new List<string> { taskId, "Individual Spectra Files", origDataFile }).Run();

                    CalibrationResults calibrationResult = theResult as CalibrationResults;

                    if (calibrationResult == null)
                    {
                        sbForThisFile.AppendLine(theResult.ToString());
                        Warn(theResult.ToString(), new List<string> { taskId, "Individual Spectra Files", origDataFile });
                        Status("Errored", new List<string> { taskId, "Individual Spectra Files", origDataFile });
                        return;
                    }

                    sbForThisFile.AppendLine("\t" + "Before Calib: MS1 (Th): " + calibrationResult.ms1meanSds.First());
                    sbForThisFile.AppendLine("\t" + "Before Calib: MS2 (Th): " + calibrationResult.ms2meanSds.First());
                    sbForThisFile.AppendLine("\t" + "After Linear: MS1 (Th): " + calibrationResult.ms1meanSds.Last());
                    sbForThisFile.AppendLine("\t" + "After Linear: MS2 (Th): " + calibrationResult.ms2meanSds.Last());
                }
                if (NonLinearCalibration)
                {
                    // Second search round

                    Status("Getting ms2 scans for second round...", new List<string> { taskId, "Individual Spectra Files", origDataFile });
                    var listOfSortedms2ScansTest = GetMs2Scans(myMsDataFile, origDataFile, DoPrecursorDeconvolution, UseProvidedPrecursorInfo, DeconvolutionIntensityRatio, DeconvolutionMaxAssumedChargeState, DeconvolutionMassTolerance).OrderBy(b => b.PrecursorMass).ToArray();
                    var searchEngineTest = new ClassicSearchEngine(listOfSortedms2ScansTest, variableModifications, fixedModifications, proteinList, ProductMassTolerance, Protease, searchModes, MaxMissedCleavages, MinPeptideLength, MaxPeptideLength, MaxModificationIsoforms, lp, new List<string> { taskId, "Individual Spectra Files", origDataFile }, ConserveMemory, InitiatorMethionineBehavior, false);
                    var searchResultsTest = (SearchResults)searchEngineTest.Run();

                    // Group and order psms
                    List<Psm>[] allPsms = new List<Psm>[1];
                    allPsms[0] = searchResultsTest.Psms[0].ToList();

                    SequencesToActualProteinPeptidesEngine sequencesToActualProteinPeptidesEngineTest = new SequencesToActualProteinPeptidesEngine(allPsms, proteinList, searchModes, Protease, MaxMissedCleavages, MinPeptideLength, MaxPeptideLength, InitiatorMethionineBehavior, fixedModifications, variableModifications, MaxModificationIsoforms, new List<string> { taskId, "Individual Spectra Files", origDataFile });
                    var resTest = (SequencesToActualProteinPeptidesEngineResults)sequencesToActualProteinPeptidesEngineTest.Run();
                    Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>> compactPeptideToProteinPeptideMatchingTest = resTest.CompactPeptideToProteinPeptideMatching;

                    allPsms[0] = allPsms[0].Where(b => b != null).OrderByDescending(b => b.Score).ThenBy(b => Math.Abs(b.ScanPrecursorMass - b.MostProbableProteinInfo.PeptideMonoisotopicMass)).GroupBy(b => new Tuple<string, int, double>(b.FullFilePath, b.ScanNumber, b.MostProbableProteinInfo.PeptideMonoisotopicMass)).Select(b => b.First()).ToList();

                    new FdrAnalysisEngine(allPsms, searchModes, new List<string> { taskId, "Individual Spectra Files", origDataFile }).Run();

                    if (WriteIntermediateFiles)
                        WritePsmsToTsv(allPsms[0], OutputFolder, Path.GetFileNameWithoutExtension(origDataFile) + "PSMsBeforeNonLinearCalib", new List<string> { taskId, "Individual Spectra Files", origDataFile });

                    var goodIdentifications = allPsms[0].Where(b => b.FdrInfo.QValue < 0.01 && !b.MostProbableProteinInfo.IsDecoy).ToList();

                    Action<List<LabeledMs1DataPoint>, string> ms1Action = (List<LabeledMs1DataPoint> theList, string s) => {; };
                    Action<List<LabeledMs2DataPoint>, string> ms2Action = (List<LabeledMs2DataPoint> theList, string s) => {; };
                    if (WriteIntermediateFiles)
                    {
                        ms1Action = (List<LabeledMs1DataPoint> theList, string s) => WriteMs1DataPoints(theList, OutputFolder, Path.GetFileNameWithoutExtension(origDataFile) + "inNonLinear" + s, new List<string> { taskId, "Individual Spectra Files", origDataFile });
                        ms2Action = (List<LabeledMs2DataPoint> theList, string s) => WriteMs2DataPoints(theList, OutputFolder, Path.GetFileNameWithoutExtension(origDataFile) + "inNonLinear" + s, new List<string> { taskId, "Individual Spectra Files", origDataFile });
                    }

                    var theResult = new CalibrationEngine(myMsDataFile, ProductMassTolerance, goodIdentifications, minMS1isotopicPeaksNeededForConfirmedIdentification, minMS2isotopicPeaksNeededForConfirmedIdentification, numFragmentsNeededForEveryIdentification, PrecursorMassTolerance, fragmentTypesForCalibration, ms1Action, ms2Action, true, new List<string> { taskId, "Individual Spectra Files", origDataFile }).Run();
                    CalibrationResults calibrationResult = theResult as CalibrationResults;

                    if (calibrationResult == null)
                    {
                        sbForThisFile.AppendLine(theResult.ToString());
                        Warn(theResult.ToString(), new List<string> { taskId, "Individual Spectra Files", origDataFile });
                        Status("Errored", new List<string> { taskId, "Individual Spectra Files", origDataFile });
                        return;
                    }

                    sbForThisFile.AppendLine("\t" + "After NonLin: MS1 (Th): " + calibrationResult.ms1meanSds.Last());
                    sbForThisFile.AppendLine("\t" + "After NonLin: MS2 (Th): " + calibrationResult.ms2meanSds.Last());
                }

                if (WriteIntermediateFiles)
                {
                    var ms2ScansAfterCalib = GetMs2Scans(myMsDataFile, origDataFile, DoPrecursorDeconvolution, UseProvidedPrecursorInfo, DeconvolutionIntensityRatio, DeconvolutionMaxAssumedChargeState, DeconvolutionMassTolerance).OrderBy(b => b.PrecursorMass).ToArray();
                    var searchResultsAfterCalib = (SearchResults)new ClassicSearchEngine(ms2ScansAfterCalib, variableModifications, fixedModifications, proteinList, ProductMassTolerance, Protease, searchModes, MaxMissedCleavages, MinPeptideLength, MaxPeptideLength, MaxModificationIsoforms, lp, new List<string> { taskId, "Individual Spectra Files", origDataFile }, ConserveMemory, InitiatorMethionineBehavior, false).Run();

                    List<Psm>[] allPsms = new List<Psm>[1];
                    allPsms[0] = searchResultsAfterCalib.Psms[0].ToList();

                    SequencesToActualProteinPeptidesEngine sequencesToActualProteinPeptidesEngineTest = new SequencesToActualProteinPeptidesEngine(allPsms, proteinList, searchModes, Protease, MaxMissedCleavages, MinPeptideLength, MaxPeptideLength, InitiatorMethionineBehavior, fixedModifications, variableModifications, MaxModificationIsoforms, new List<string> { taskId, "Individual Spectra Files", origDataFile });
                    var resTest = (SequencesToActualProteinPeptidesEngineResults)sequencesToActualProteinPeptidesEngineTest.Run();
                    Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>> compactPeptideToProteinPeptideMatchingTest = resTest.CompactPeptideToProteinPeptideMatching;

                    // Group and order psms
                    allPsms[0] = searchResultsAfterCalib.Psms[0].Where(b => b != null).OrderByDescending(b => b.Score).ThenBy(b => Math.Abs(b.ScanPrecursorMass - b.MostProbableProteinInfo.PeptideMonoisotopicMass)).GroupBy(b => new Tuple<string, int, double>(b.FullFilePath, b.ScanNumber, b.MostProbableProteinInfo.PeptideMonoisotopicMass)).Select(b => b.First()).ToList();

                    new FdrAnalysisEngine(allPsms, searchModes, new List<string> { taskId, "Individual Spectra Files", origDataFile }).Run();

                    if (WriteIntermediateFiles)
                        WritePsmsToTsv(allPsms[0], OutputFolder, Path.GetFileNameWithoutExtension(origDataFile) + "PSMsAfterCalib", new List<string> { taskId, "Individual Spectra Files", origDataFile });
                }

                myTaskResults.AddNiceText(sbForThisFile.ToString());

                Status("Writing mzML!", new List<string> { taskId, "Individual Spectra Files", origDataFile });
                var path = Path.Combine(OutputFolder, Path.GetFileNameWithoutExtension(origDataFile) + "-Calibrated.mzML");

                MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(myMsDataFile, path, false);

                SucessfullyFinishedWritingFile(path, new List<string> { taskId, "Individual Spectra Files", origDataFile });

                myTaskResults.newSpectra.Add(path);
                FinishedDataFile(origDataFile, new List<string> { taskId, "Individual Spectra Files", origDataFile });
                ReportProgress(new ProgressEventArgs(100, "Done!", new List<string> { taskId, "Individual Spectra Files", origDataFile }));
            }
            );
            ReportProgress(new ProgressEventArgs(100, "Done!", new List<string> { taskId, "Individual Spectra Files" }));

            return myTaskResults;
        }

        #endregion Protected Methods

    }
}