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
            CommonParameters = new CommonParameters
            {
                ProductMassTolerance = new PpmTolerance(30),
                TrimMs1Peaks = false,
                TrimMsMsPeaks = false,
            };

            CalibrationParameters = new CalibrationParameters();
        }

        #endregion Public Constructors

        #region Public Properties

        public CalibrationParameters CalibrationParameters { get; set; }

        #endregion Public Properties

        #region Protected Internal Methods

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

        #endregion Protected Internal Methods

        #region Protected Methods

        protected override MyTaskResults RunSpecific(string OutputFolder, List<DbForTask> dbFilenameList, List<string> currentRawFileList, string taskId, FileSpecificSettings[] fileSettings)
        {
            myTaskResults = new MyTaskResults(this)
            {
                newSpectra = new List<string>()
            };
            MassDiffAcceptor searchMode;
            if (CalibrationParameters.PrecursorMassTolerance is PpmTolerance)
                searchMode = new SinglePpmAroundZeroSearchMode(CalibrationParameters.PrecursorMassTolerance.Value);
            else
                searchMode = new SingleAbsoluteAroundZeroSearchMode(CalibrationParameters.PrecursorMassTolerance.Value);
            var searchModes = searchMode;

            Status("Loading modifications...", new List<string> { taskId });
            List<ModificationWithMass> variableModifications = GlobalEngineLevelSettings.AllModsKnown.OfType<ModificationWithMass>().Where(b => CommonParameters.ListOfModsVariable.Contains(new Tuple<string, string>(b.modificationType, b.id))).ToList();
            List<ModificationWithMass> fixedModifications = GlobalEngineLevelSettings.AllModsKnown.OfType<ModificationWithMass>().Where(b => CommonParameters.ListOfModsFixed.Contains(new Tuple<string, string>(b.modificationType, b.id))).ToList();
            List<ModificationWithMass> localizeableModifications;
            if (CommonParameters.LocalizeAll)
                localizeableModifications = GlobalEngineLevelSettings.AllModsKnown.OfType<ModificationWithMass>().ToList();
            else
                localizeableModifications = GlobalEngineLevelSettings.AllModsKnown.OfType<ModificationWithMass>().Where(b => CommonParameters.ListOfModsLocalize.Contains(new Tuple<string, string>(b.modificationType, b.id))).ToList();

            Status("Loading proteins...", new List<string> { taskId });
            var proteinList = dbFilenameList.SelectMany(b => LoadProteinDb(b.FilePath, true, localizeableModifications, b.IsContaminant, out Dictionary<string, Modification> um)).ToList();

            List<ProductType> lp = new List<ProductType>();
            FragmentTypes fragmentTypesForCalibration = FragmentTypes.None;
            if (CommonParameters.BIons)
            {
                fragmentTypesForCalibration = fragmentTypesForCalibration | FragmentTypes.b;
                lp.Add(ProductType.B);
            }
            if (CommonParameters.YIons)
            {
                fragmentTypesForCalibration = fragmentTypesForCalibration | FragmentTypes.y;
                lp.Add(ProductType.Y);
            }
            if (CommonParameters.CIons)
            {
                fragmentTypesForCalibration = fragmentTypesForCalibration | FragmentTypes.c;
                lp.Add(ProductType.C);
            }
            if (CommonParameters.ZdotIons)
            {
                fragmentTypesForCalibration = fragmentTypesForCalibration | FragmentTypes.zdot;
                lp.Add(ProductType.Zdot);
            }

            proseCreatedWhileRunning.Append("The following calibration settings were used: ");
            proseCreatedWhileRunning.Append("protease = " + CommonParameters.DigestionParams.Protease + "; ");
            proseCreatedWhileRunning.Append("maximum missed cleavages = " + CommonParameters.DigestionParams.MaxMissedCleavages + "; ");
            proseCreatedWhileRunning.Append("minimum peptide length = " + CommonParameters.DigestionParams.MinPeptideLength + "; ");
            if (CommonParameters.DigestionParams.MaxPeptideLength == null)
            {
                proseCreatedWhileRunning.Append("maximum peptide length = unspecified; ");
            }
            else
            {
                proseCreatedWhileRunning.Append("maximum peptide length = " + CommonParameters.DigestionParams.MaxPeptideLength + "; ");
            }
            proseCreatedWhileRunning.Append("initiator methionine behavior = " + CommonParameters.DigestionParams.InitiatorMethionineBehavior + "; ");
            proseCreatedWhileRunning.Append("max modification isoforms = " + CommonParameters.DigestionParams.MaxModificationIsoforms + "; ");

            proseCreatedWhileRunning.Append("fixed modifications = " + string.Join(", ", fixedModifications.Select(m => m.id)) + "; ");
            proseCreatedWhileRunning.Append("variable modifications = " + string.Join(", ", variableModifications.Select(m => m.id)) + "; ");
            proseCreatedWhileRunning.Append("parent mass tolerance(s) = {" + searchModes.ToProseString() + "}; ");
            proseCreatedWhileRunning.Append("product mass tolerance = " + CommonParameters.ProductMassTolerance + " Da. ");
            proseCreatedWhileRunning.Append("The combined search database contained " + proteinList.Count + " total entries including " + proteinList.Where(p => p.IsContaminant).Count() + " contaminant sequences. ");

            object lock1 = new object();
            ParallelOptions parallelOptions = new ParallelOptions();
            if (CommonParameters.MaxDegreeOfParallelism.HasValue)
                parallelOptions.MaxDegreeOfParallelism = CommonParameters.MaxDegreeOfParallelism.Value;
            Status("Calibrating...", new List<string> { taskId });
            TerminusType terminusType = ProductTypeToTerminusType.IdentifyTerminusType(lp);
            Parallel.For(0, currentRawFileList.Count, parallelOptions, spectraFileIndex =>
            {
                var origDataFile = currentRawFileList[spectraFileIndex];
                IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> myMsDataFile;
                lock (lock1) // Lock because reading is sequential
                {
                    NewCollection(Path.GetFileName(origDataFile), new List<string> { taskId, "Individual Spectra Files", origDataFile });
                    StartingDataFile(origDataFile, new List<string> { taskId, "Individual Spectra Files", origDataFile });
                    Status("Loading spectra file " + origDataFile + "...", new List<string> { taskId, "Individual Spectra Files", origDataFile });
                    if (Path.GetExtension(origDataFile).Equals(".mzML", StringComparison.InvariantCultureIgnoreCase))
                        myMsDataFile = Mzml.LoadAllStaticData(origDataFile);
                    else
                        myMsDataFile = ThermoStaticData.LoadAllStaticData(origDataFile);
                }

                StringBuilder sbForThisFile = new StringBuilder();
                sbForThisFile.AppendLine(origDataFile);

                int minMS1isotopicPeaksNeededForConfirmedIdentification = 3;
                int minMS2isotopicPeaksNeededForConfirmedIdentification = 2;
                int numFragmentsNeededForEveryIdentification = 10;

                Random rnd = new Random(spectraFileIndex);
                // Linear calibration
                {
                    Status("Getting ms2 scans...", new List<string> { taskId, "Individual Spectra Files", origDataFile });
                    var listOfSortedms2Scans = GetMs2Scans(myMsDataFile, origDataFile, CommonParameters.DoPrecursorDeconvolution, CommonParameters.UseProvidedPrecursorInfo, CommonParameters.DeconvolutionIntensityRatio, CommonParameters.DeconvolutionMaxAssumedChargeState, CommonParameters.DeconvolutionMassTolerance).OrderBy(b => b.PrecursorMass).ToArray();

                    Psm[] allPsmsArray = new Psm[listOfSortedms2Scans.Length];

                    var searchEngine = new ClassicSearchEngine(allPsmsArray, listOfSortedms2Scans, variableModifications, fixedModifications, proteinList, lp, searchModes, false, CommonParameters, new List<string> { taskId, "Individual Spectra Files", origDataFile });
                    searchEngine.Run();
                    List<Psm> allPsms = allPsmsArray.ToList();
                    // Group and order psms

                    SequencesToActualProteinPeptidesEngine sequencesToActualProteinPeptidesEngine = new SequencesToActualProteinPeptidesEngine(allPsms, proteinList, fixedModifications, variableModifications, lp, new List<DigestionParams> { CommonParameters.DigestionParams }, new List<string> { taskId, "Individual Spectra Files", origDataFile });

                    var res = (SequencesToActualProteinPeptidesEngineResults)sequencesToActualProteinPeptidesEngine.Run();
                    Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>> compactPeptideToProteinPeptideMatching = res.CompactPeptideToProteinPeptideMatching;

                    foreach (var huh in allPsms)
                        if (huh != null && huh.MostProbableProteinInfo == null)
                            huh.MatchToProteinLinkedPeptides(compactPeptideToProteinPeptideMatching);

                    allPsms = allPsms.Where(b => b != null).OrderByDescending(b => b.Score).ThenBy(b => b.PeptideMonisotopicMass.HasValue ? Math.Abs(b.ScanPrecursorMass - b.PeptideMonisotopicMass.Value) : double.MaxValue).GroupBy(b => new Tuple<string, int, double?>(b.FullFilePath, b.ScanNumber, b.PeptideMonisotopicMass)).Select(b => b.First()).ToList();

                    // Run FDR analysis on single file results
                    new FdrAnalysisEngine(allPsms, searchModes, new List<string> { taskId, "Individual Spectra Files", origDataFile }).Run();

                    if (CalibrationParameters.WriteIntermediateFiles)
                    {
                        var writtenFile = Path.Combine(OutputFolder, Path.GetFileNameWithoutExtension(origDataFile) + "-PSMsBeforeLinearCalib.psmtsv");
                        WritePsmsToTsv(allPsms, writtenFile);
                        SucessfullyFinishedWritingFile(writtenFile, new List<string> { taskId, "Individual Spectra Files", origDataFile });
                    }

                    var goodIdentifications = allPsms.Where(b => b.FdrInfo.QValue < 0.01 && !b.IsDecoy).ToList();

                    Action<List<LabeledMs1DataPoint>, string> ms1Action = (List<LabeledMs1DataPoint> theList, string s) => {; };
                    Action<List<LabeledMs2DataPoint>, string> ms2Action = (List<LabeledMs2DataPoint> theList, string s) => {; };
                    if (CalibrationParameters.WriteIntermediateFiles)
                    {
                        ms1Action = (List<LabeledMs1DataPoint> theList, string s) => WriteMs1DataPoints(theList, OutputFolder, Path.GetFileNameWithoutExtension(origDataFile) + "inLinear" + s, new List<string> { taskId, "Individual Spectra Files", origDataFile });
                        ms2Action = (List<LabeledMs2DataPoint> theList, string s) => WriteMs2DataPoints(theList, OutputFolder, Path.GetFileNameWithoutExtension(origDataFile) + "inLinear" + s, new List<string> { taskId, "Individual Spectra Files", origDataFile });
                    }

                    MetaMorpheusEngineResults theResult = new CalibrationEngine(myMsDataFile, CommonParameters.ProductMassTolerance, goodIdentifications, minMS1isotopicPeaksNeededForConfirmedIdentification, minMS2isotopicPeaksNeededForConfirmedIdentification, numFragmentsNeededForEveryIdentification, CalibrationParameters.PrecursorMassTolerance, fragmentTypesForCalibration, ms1Action, ms2Action, false, rnd, new List<string> { taskId, "Individual Spectra Files", origDataFile }).Run();

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
                if (CalibrationParameters.NonLinearCalibration)
                {
                    // Second search round

                    Status("Getting ms2 scans for second round...", new List<string> { taskId, "Individual Spectra Files", origDataFile });
                    var listOfSortedms2ScansTest = GetMs2Scans(myMsDataFile, origDataFile, CommonParameters.DoPrecursorDeconvolution, CommonParameters.UseProvidedPrecursorInfo, CommonParameters.DeconvolutionIntensityRatio, CommonParameters.DeconvolutionMaxAssumedChargeState, CommonParameters.DeconvolutionMassTolerance).OrderBy(b => b.PrecursorMass).ToArray();
                    Psm[] allPsmsArray = new Psm[listOfSortedms2ScansTest.Length];
                    var searchEngineTest = new ClassicSearchEngine(allPsmsArray, listOfSortedms2ScansTest, variableModifications, fixedModifications, proteinList, lp, searchModes, false, CommonParameters, new List<string> { taskId, "Individual Spectra Files", origDataFile });

                    searchEngineTest.Run();

                    // Group and order psms
                    List<Psm> allPsms = allPsmsArray.ToList();

                    SequencesToActualProteinPeptidesEngine sequencesToActualProteinPeptidesEngineTest = new SequencesToActualProteinPeptidesEngine(allPsms, proteinList, fixedModifications, variableModifications, lp, new List<DigestionParams> { CommonParameters.DigestionParams }, new List<string> { taskId, "Individual Spectra Files", origDataFile });

                    var resTest = (SequencesToActualProteinPeptidesEngineResults)sequencesToActualProteinPeptidesEngineTest.Run();
                    Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>> compactPeptideToProteinPeptideMatchingTest = resTest.CompactPeptideToProteinPeptideMatching;

                    foreach (var huh in allPsms)
                        if (huh != null && huh.MostProbableProteinInfo == null)
                            huh.MatchToProteinLinkedPeptides(compactPeptideToProteinPeptideMatchingTest);

                    allPsms = allPsms.Where(b => b != null).OrderByDescending(b => b.Score).ThenBy(b => b.PeptideMonisotopicMass.HasValue ? Math.Abs(b.ScanPrecursorMass - b.PeptideMonisotopicMass.Value) : double.MaxValue).GroupBy(b => new Tuple<string, int, double?>(b.FullFilePath, b.ScanNumber, b.PeptideMonisotopicMass)).Select(b => b.First()).ToList();

                    new FdrAnalysisEngine(allPsms, searchModes, new List<string> { taskId, "Individual Spectra Files", origDataFile }).Run();

                    if (CalibrationParameters.WriteIntermediateFiles)
                    {
                        var writtenFile = Path.Combine(OutputFolder, Path.GetFileNameWithoutExtension(origDataFile) + "-PSMsBeforeNonLinearCalib.psmtsv");
                        WritePsmsToTsv(allPsms, writtenFile);
                        SucessfullyFinishedWritingFile(writtenFile, new List<string> { taskId, "Individual Spectra Files", origDataFile });
                    }
                    var goodIdentifications = allPsms.Where(b => b.FdrInfo.QValue < 0.01 && !b.IsDecoy).ToList();

                    Action<List<LabeledMs1DataPoint>, string> ms1Action = (List<LabeledMs1DataPoint> theList, string s) => {; };
                    Action<List<LabeledMs2DataPoint>, string> ms2Action = (List<LabeledMs2DataPoint> theList, string s) => {; };
                    if (CalibrationParameters.WriteIntermediateFiles)
                    {
                        ms1Action = (List<LabeledMs1DataPoint> theList, string s) => WriteMs1DataPoints(theList, OutputFolder, Path.GetFileNameWithoutExtension(origDataFile) + "inNonLinear" + s, new List<string> { taskId, "Individual Spectra Files", origDataFile });
                        ms2Action = (List<LabeledMs2DataPoint> theList, string s) => WriteMs2DataPoints(theList, OutputFolder, Path.GetFileNameWithoutExtension(origDataFile) + "inNonLinear" + s, new List<string> { taskId, "Individual Spectra Files", origDataFile });
                    }

                    var theResult = new CalibrationEngine(myMsDataFile, CommonParameters.ProductMassTolerance, goodIdentifications, minMS1isotopicPeaksNeededForConfirmedIdentification, minMS2isotopicPeaksNeededForConfirmedIdentification, numFragmentsNeededForEveryIdentification, CalibrationParameters.PrecursorMassTolerance, fragmentTypesForCalibration, ms1Action, ms2Action, true, rnd, new List<string> { taskId, "Individual Spectra Files", origDataFile }).Run();
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

                if (CalibrationParameters.WriteIntermediateFiles)
                {
                    var ms2ScansAfterCalib = GetMs2Scans(myMsDataFile, origDataFile, CommonParameters.DoPrecursorDeconvolution, CommonParameters.UseProvidedPrecursorInfo, CommonParameters.DeconvolutionIntensityRatio, CommonParameters.DeconvolutionMaxAssumedChargeState, CommonParameters.DeconvolutionMassTolerance).OrderBy(b => b.PrecursorMass).ToArray();
                    Psm[] allPsmsArray = new Psm[ms2ScansAfterCalib.Length];
                    new ClassicSearchEngine(allPsmsArray, ms2ScansAfterCalib, variableModifications, fixedModifications, proteinList, lp, searchModes, false, CommonParameters, new List<string> { taskId, "Individual Spectra Files", origDataFile }).Run();

                    List<Psm> allPsms = allPsmsArray.ToList();

                    SequencesToActualProteinPeptidesEngine sequencesToActualProteinPeptidesEngineTest = new SequencesToActualProteinPeptidesEngine(allPsms, proteinList, fixedModifications, variableModifications, lp, new List<DigestionParams> { CommonParameters.DigestionParams }, new List<string> { taskId, "Individual Spectra Files", origDataFile });
                    var resTest = (SequencesToActualProteinPeptidesEngineResults)sequencesToActualProteinPeptidesEngineTest.Run();
                    Dictionary<CompactPeptideBase, HashSet<PeptideWithSetModifications>> compactPeptideToProteinPeptideMatchingTest = resTest.CompactPeptideToProteinPeptideMatching;

                    foreach (var huh in allPsms)
                        if (huh != null && huh.MostProbableProteinInfo == null)
                            huh.MatchToProteinLinkedPeptides(compactPeptideToProteinPeptideMatchingTest);

                    // Group and order psms
                    allPsms = allPsmsArray.Where(b => b != null).OrderByDescending(b => b.Score).ThenBy(b => b.PeptideMonisotopicMass.HasValue ? Math.Abs(b.ScanPrecursorMass - b.PeptideMonisotopicMass.Value) : double.MaxValue).GroupBy(b => new Tuple<string, int, double?>(b.FullFilePath, b.ScanNumber, b.PeptideMonisotopicMass)).Select(b => b.First()).ToList();

                    new FdrAnalysisEngine(allPsms, searchModes, new List<string> { taskId, "Individual Spectra Files", origDataFile }).Run();

                    if (CalibrationParameters.WriteIntermediateFiles)
                    {
                        var writtenFile = Path.Combine(OutputFolder, Path.GetFileNameWithoutExtension(origDataFile) + "-PSMsAfterCalib.psmtsv");
                        WritePsmsToTsv(allPsms, writtenFile);
                        SucessfullyFinishedWritingFile(writtenFile, new List<string> { taskId, "Individual Spectra Files", origDataFile });
                    }
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