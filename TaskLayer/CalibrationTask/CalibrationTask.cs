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

            ListOfModsVariable = new List<Tuple<string, string>> { new Tuple<string, string>("Common Variable", "Oxidation of M") };
            ListOfModsFixed = new List<Tuple<string, string>> { new Tuple<string, string>("Common Fixed", "Carbamidomethyl of C") };
            ListOfModsLocalize = GlobalTaskLevelSettings.AllModsKnown.Select(b => new Tuple<string, string>(b.modificationType, b.id)).ToList();

            MaxDegreeOfParallelism = -1;
        }

        #endregion Public Constructors

        #region Public Properties

        public static List<string> AllModLists { get; private set; }

        public MyTask TaskType { get; internal set; }

        public InitiatorMethionineBehavior InitiatorMethionineBehavior { get; set; }

        public int MaxMissedCleavages { get; set; }

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
                output.WriteLine(LabeledMs1DataPoint.TabSeparatedHeaderForMs1);
                foreach (var dp in items)
                {
                    output.Write(string.Join("\t", dp.Inputs));
                    output.WriteLine("\t" + dp.Label);
                }
            }
            SucessfullyFinishedWritingFile(writtenFile, nestedIDs);
        }

        protected internal void WriteMs2DataPoints(List<LabeledMs2DataPoint> items, string outputFolder, string fileName, List<string> nestedIDs)
        {
            var writtenFile = Path.Combine(outputFolder, fileName + ".ms2dptsv");
            using (StreamWriter output = new StreamWriter(writtenFile))
            {
                output.WriteLine(LabeledMs2DataPoint.TabSeparatedHeaderForMs1);
                foreach (var dp in items)
                {
                    output.Write(string.Join("\t", dp.Inputs));
                    output.WriteLine("\t" + dp.Label);
                }
            }
            SucessfullyFinishedWritingFile(writtenFile, nestedIDs);
        }

        #endregion Protected Internal Methods

        #region Protected Methods

        protected override MyTaskResults RunSpecific(string OutputFolder, List<DbForTask> dbFilenameList, List<string> currentRawFileList, string taskId)
        {
            var myTaskResults = new MyTaskResults(this)
            {
                newSpectra = new List<string>()
            };
            SearchMode searchMode;
            if (PrecursorMassTolerance.Unit == ToleranceUnit.PPM)
                searchMode = new SinglePpmAroundZeroSearchMode(PrecursorMassTolerance.Value);
            else
                searchMode = new SingleAbsoluteAroundZeroSearchMode(PrecursorMassTolerance.Value);
            var searchModes = new List<SearchMode> { searchMode };

            List<PsmParent>[] allPsms = new List<PsmParent>[1];
            allPsms[0] = new List<PsmParent>();

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
            var proteinList = dbFilenameList.SelectMany(b => LoadProteinDb(b.FileName, true, localizeableModifications, b.IsContaminant, null, out um)).ToList();

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
                var origDataFileName = currentRawFileList[spectraFileIndex];
                LocalMS2Scan[] listOfSortedms2Scans;
                IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> myMsDataFile;
                lock (lock1) // Lock because reading is sequential
                {
                    NewCollection(Path.GetFileName(origDataFileName), new List<string> { taskId, "Individual Searches", origDataFileName });
                    StartingDataFile(origDataFileName, new List<string> { taskId, "Individual Searches", origDataFileName });
                    Status("Loading spectra file " + origDataFileName + "...", new List<string> { taskId, "Individual Searches", origDataFileName });
                    if (Path.GetExtension(origDataFileName).Equals(".mzML"))
                        myMsDataFile = Mzml.LoadAllStaticData(origDataFileName);
                    else
                        myMsDataFile = ThermoStaticData.LoadAllStaticData(origDataFileName);
                    Status("Opening spectra file " + origDataFileName + "...", new List<string> { taskId, "Individual Searches", origDataFileName });
                    listOfSortedms2Scans = MetaMorpheusEngine.GetMs2Scans(myMsDataFile).OrderBy(b => b.MonoisotopicPrecursorMass).ToArray();
                }

                var searchEngine = new ClassicSearchEngine(listOfSortedms2Scans, myMsDataFile.NumSpectra, variableModifications, fixedModifications, proteinList, ProductMassTolerance, Protease, searchModes, MaxMissedCleavages, MaxModificationIsoforms, origDataFileName, lp, new List<string> { taskId, "Individual Searches", origDataFileName }, false);

                var searchResults = (ClassicSearchResults)searchEngine.Run();
                myTaskResults.AddResultText(searchResults);

                for (int i = 0; i < searchModes.Count; i++)
                    allPsms[i].AddRange(searchResults.OuterPsms[i]);

                InitiatorMethionineBehavior initiatorMethionineBehavior = InitiatorMethionineBehavior.Variable;
                // Run analysis on single file results
                var analysisEngine = new AnalysisEngine(searchResults.OuterPsms, compactPeptideToProteinPeptideMatching, proteinList, variableModifications, fixedModifications, localizeableModifications, Protease, searchModes, myMsDataFile, ProductMassTolerance, (BinTreeStructure myTreeStructure, string s) => WriteTree(myTreeStructure, OutputFolder, Path.GetFileNameWithoutExtension(origDataFileName) + s, new List<string> { taskId, "Individual Searches", origDataFileName }), (List<NewPsmWithFdr> h, string s) => WritePsmsToTsv(h, OutputFolder, Path.GetFileNameWithoutExtension(origDataFileName) + s, new List<string> { taskId, "Individual Searches", origDataFileName }), null, false, false, false, MaxMissedCleavages, MaxModificationIsoforms, false, lp, double.NaN, initiatorMethionineBehavior, new List<string> { taskId, "Individual Searches", origDataFileName }, false, 0, 0);

                var analysisResults = (AnalysisResults)analysisEngine.Run();
                myTaskResults.AddResultText(analysisResults);

                var identifications = analysisResults.AllResultingIdentifications[0];

                //Now can calibrate!!!

                int minMS1isotopicPeaksNeededForConfirmedIdentification = 3;
                int minMS2isotopicPeaksNeededForConfirmedIdentification = 2;
                int numFragmentsNeededForEveryIdentification = 10;

                // TODO: fix the tolerance calculation below

                var a = new CalibrationEngine(myMsDataFile, ProductMassTolerance, identifications, minMS1isotopicPeaksNeededForConfirmedIdentification, minMS2isotopicPeaksNeededForConfirmedIdentification, numFragmentsNeededForEveryIdentification, PrecursorMassTolerance, fragmentTypesForCalibration, (List<LabeledMs1DataPoint> theList, string s) => WriteMs1DataPoints(theList, OutputFolder, Path.GetFileNameWithoutExtension(origDataFileName) + s, new List<string> { taskId, "Individual Searches", origDataFileName }), (List<LabeledMs2DataPoint> theList, string s) => WriteMs2DataPoints(theList, OutputFolder, Path.GetFileNameWithoutExtension(origDataFileName) + s, new List<string> { taskId, "Individual Searches", origDataFileName }), false, new List<string> { taskId, "Individual Searches", origDataFileName });

                var resultBeforeFC = a.Run();
                myTaskResults.AddResultText(resultBeforeFC);

                // Second search round

                var listOfSortedms2ScansTest = MetaMorpheusEngine.GetMs2Scans(myMsDataFile).OrderBy(b => b.MonoisotopicPrecursorMass).ToArray();
                var searchEngineTest = new ClassicSearchEngine(listOfSortedms2ScansTest, myMsDataFile.NumSpectra, variableModifications, fixedModifications, proteinList, ProductMassTolerance, Protease, searchModes, MaxMissedCleavages, MaxModificationIsoforms, origDataFileName, lp, new List<string> { taskId, "Individual Searches", origDataFileName }, false);
                var searchResultsTest = (ClassicSearchResults)searchEngineTest.Run();
                myTaskResults.AddResultText(searchResultsTest);
                var analysisEngineTest = new AnalysisEngine(searchResultsTest.OuterPsms, compactPeptideToProteinPeptideMatching, proteinList, variableModifications, fixedModifications, localizeableModifications, Protease, searchModes, myMsDataFile, ProductMassTolerance, (BinTreeStructure myTreeStructure, string s) => WriteTree(myTreeStructure, OutputFolder, Path.GetFileNameWithoutExtension(origDataFileName) + s + "test", new List<string> { taskId, "Individual Searches", origDataFileName }), (List<NewPsmWithFdr> h, string s) => WritePsmsToTsv(h, OutputFolder, Path.GetFileNameWithoutExtension(origDataFileName) + s + "test", new List<string> { taskId, "Individual Searches", origDataFileName }), null, false, false, false, MaxMissedCleavages, MaxModificationIsoforms, false, lp, double.NaN, initiatorMethionineBehavior, new List<string> { taskId, "Individual Searches", origDataFileName }, false, 0, 0);
                var analysisResultsTest = (AnalysisResults)analysisEngineTest.Run();
                myTaskResults.AddResultText(analysisResultsTest);

                //

                var resultAfterFC = new CalibrationEngine(myMsDataFile, ProductMassTolerance, analysisResultsTest.AllResultingIdentifications[0], minMS1isotopicPeaksNeededForConfirmedIdentification, minMS2isotopicPeaksNeededForConfirmedIdentification, numFragmentsNeededForEveryIdentification, PrecursorMassTolerance, fragmentTypesForCalibration, (List<LabeledMs1DataPoint> theList, string s) => WriteMs1DataPoints(theList, OutputFolder, Path.GetFileNameWithoutExtension(origDataFileName) + s + "after", new List<string> { taskId, "Individual Searches", origDataFileName }), (List<LabeledMs2DataPoint> theList, string s) => WriteMs2DataPoints(theList, OutputFolder, Path.GetFileNameWithoutExtension(origDataFileName) + s + "after", new List<string> { taskId, "Individual Searches", origDataFileName }), true, new List<string> { taskId, "Individual Searches", origDataFileName }).Run();
                myTaskResults.AddResultText(resultAfterFC);

                // Final search round - not required

                var listOfSortedms2ScansTest2 = MetaMorpheusEngine.GetMs2Scans(myMsDataFile).OrderBy(b => b.MonoisotopicPrecursorMass).ToArray();
                var searchEngineTest2 = new ClassicSearchEngine(listOfSortedms2ScansTest2, myMsDataFile.NumSpectra, variableModifications, fixedModifications, proteinList, ProductMassTolerance, Protease, searchModes, MaxMissedCleavages, MaxModificationIsoforms, origDataFileName, lp, new List<string> { taskId, "Individual Searches", origDataFileName }, false);
                var searchResultsTest2 = (ClassicSearchResults)searchEngineTest2.Run();
                myTaskResults.AddResultText(searchResultsTest2);
                var analysisEngineTest2 = new AnalysisEngine(searchResultsTest2.OuterPsms, compactPeptideToProteinPeptideMatching, proteinList, variableModifications, fixedModifications, localizeableModifications, Protease, searchModes, myMsDataFile, ProductMassTolerance, (BinTreeStructure myTreeStructure, string s) => WriteTree(myTreeStructure, OutputFolder, Path.GetFileNameWithoutExtension(origDataFileName) + s + "test2", new List<string> { taskId, "Individual Searches", origDataFileName }), (List<NewPsmWithFdr> h, string s) => WritePsmsToTsv(h, OutputFolder, Path.GetFileNameWithoutExtension(origDataFileName) + s + "test2", new List<string> { taskId, "Individual Searches", origDataFileName }), null, false, false, false, MaxMissedCleavages, MaxModificationIsoforms, false, lp, double.NaN, initiatorMethionineBehavior, new List<string> { taskId, "Individual Searches", origDataFileName }, false, 0, 0);

                var analysisResultsTest2 = (AnalysisResults)analysisEngineTest2.Run();
                myTaskResults.AddResultText(analysisResultsTest2);

                //

                if (resultBeforeFC is MyErroredResults)
                {
                    Warn(a.ToString(), new List<string> { taskId, "Individual Searches", origDataFileName });
                }
                else
                {
                    Status("Creating _indexedmzMLConnection, putting data in it, and writing!", new List<string> { taskId, "Individual Searches", origDataFileName });
                    var path = Path.Combine(OutputFolder, Path.GetFileNameWithoutExtension(origDataFileName) + "-Calibrated.mzML");
                    lock (lock2) // Lock because writing is sequential
                    {
                        MzmlMethods.CreateAndWriteMyMzmlWithCalibratedSpectra(myMsDataFile, path, false);

                        SucessfullyFinishedWritingFile(path, new List<string> { taskId, "Individual Searches", origDataFileName });

                        myTaskResults.newSpectra.Add(path);
                    }
                }
                FinishedDataFile(origDataFileName, new List<string> { taskId, "Individual Searches", origDataFileName });
                Status("Done!", new List<string> { taskId, "Individual Searches", origDataFileName });
            }
            );
            return myTaskResults;
        }

        #endregion Protected Methods

    }
}