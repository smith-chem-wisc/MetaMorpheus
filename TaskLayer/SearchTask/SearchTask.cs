using EngineLayer;
using EngineLayer.ClassicSearch;
using EngineLayer.Indexing;
using EngineLayer.ModernSearch;
using EngineLayer.NonSpecificEnzymeSearch;
using FlashLFQ;
using MassSpectrometry;
using MzLibUtil;
using Proteomics;
using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;

namespace TaskLayer
{
    public class SearchTask : MetaMorpheusTask
    {
        public SearchTask() : base(MyTask.Search)
        {
            CommonParameters = new CommonParameters();

            SearchParameters = new SearchParameters();
        }

        public SearchParameters SearchParameters { get; set; }

        public static MassDiffAcceptor GetMassDiffAcceptor(Tolerance precursorMassTolerance, MassDiffAcceptorType massDiffAcceptorType, string customMdac)
        {
            switch (massDiffAcceptorType)
            {
                case MassDiffAcceptorType.Exact:
                    if (precursorMassTolerance is PpmTolerance)
                        return new SinglePpmAroundZeroSearchMode(precursorMassTolerance.Value);
                    else
                        return new SingleAbsoluteAroundZeroSearchMode(precursorMassTolerance.Value);

                case MassDiffAcceptorType.OneMM:
                    return new DotMassDiffAcceptor("1mm", new List<double> { 0, 1.0029 }, precursorMassTolerance);

                case MassDiffAcceptorType.TwoMM:
                    return new DotMassDiffAcceptor("2mm", new List<double> { 0, 1.0029, 2.0052 }, precursorMassTolerance);

                case MassDiffAcceptorType.ThreeMM:
                    return new DotMassDiffAcceptor("3mm", new List<double> { 0, 1.0029, 2.0052, 3.0077 }, precursorMassTolerance);

                case MassDiffAcceptorType.ModOpen:
                    return new IntervalMassDiffAcceptor("-187andUp", new List<DoubleRange> { new DoubleRange(-187, double.PositiveInfinity) });

                case MassDiffAcceptorType.Open:
                    return new OpenSearchMode();

                case MassDiffAcceptorType.Custom:
                    return ParseSearchMode(customMdac);

                default:
                    throw new MetaMorpheusException("Unknown MassDiffAcceptorType");
            }
        }

        protected override MyTaskResults RunSpecific(string OutputFolder, List<DbForTask> dbFilenameList, List<string> currentRawFileList, string taskId, FileSpecificParameters[] fileSettingsList)
        {
            // disable quantification if a .mgf is being used
            if (SearchParameters.DoQuantification && currentRawFileList.Any(x => Path.GetExtension(x).Equals(".mgf", StringComparison.OrdinalIgnoreCase)))
            {
                SearchParameters.DoQuantification = false;
            }

            LoadModifications(taskId, out var variableModifications, out var fixedModifications, out var localizeableModificationTypes);

            // load proteins
            List<Protein> proteinList = LoadProteins(taskId, dbFilenameList, SearchParameters.SearchTarget, SearchParameters.DecoyType, localizeableModificationTypes, CommonParameters);

            // write prose settings
            ProseCreatedWhileRunning.Append("The following search settings were used: ");
            ProseCreatedWhileRunning.Append("protease = " + CommonParameters.DigestionParams.Protease + "; ");
            ProseCreatedWhileRunning.Append("maximum missed cleavages = " + CommonParameters.DigestionParams.MaxMissedCleavages + "; ");
            ProseCreatedWhileRunning.Append("minimum peptide length = " + CommonParameters.DigestionParams.MinPeptideLength + "; ");
            ProseCreatedWhileRunning.Append(CommonParameters.DigestionParams.MaxPeptideLength == int.MaxValue ?
                "maximum peptide length = unspecified; " :
                "maximum peptide length = " + CommonParameters.DigestionParams.MaxPeptideLength + "; ");
            ProseCreatedWhileRunning.Append("initiator methionine behavior = " + CommonParameters.DigestionParams.InitiatorMethionineBehavior + "; ");
            ProseCreatedWhileRunning.Append("fixed modifications = " + string.Join(", ", fixedModifications.Select(m => m.IdWithMotif)) + "; ");
            ProseCreatedWhileRunning.Append("variable modifications = " + string.Join(", ", variableModifications.Select(m => m.IdWithMotif)) + "; ");
            ProseCreatedWhileRunning.Append("max mods per peptide = " + CommonParameters.DigestionParams.MaxModsForPeptide + "; ");
            ProseCreatedWhileRunning.Append("max modification isoforms = " + CommonParameters.DigestionParams.MaxModificationIsoforms + "; ");
            ProseCreatedWhileRunning.Append("precursor mass tolerance = " + CommonParameters.PrecursorMassTolerance + "; ");
            ProseCreatedWhileRunning.Append("product mass tolerance = " + CommonParameters.ProductMassTolerance + "; ");
            ProseCreatedWhileRunning.Append("report PSM ambiguity = " + CommonParameters.ReportAllAmbiguity + ". ");
            ProseCreatedWhileRunning.Append("The combined search database contained " + proteinList.Count(p => !p.IsDecoy) + " non-decoy protein entries including " + proteinList.Count(p => p.IsContaminant) + " contaminant sequences. ");

            // start the search task
            MyTaskResults = new MyTaskResults(this);
            int numFdrCategories = (int)(Enum.GetValues(typeof(FdrCategory)).Cast<FdrCategory>().Last() + 1); //+1 because it starts at zero
            List<PeptideSpectralMatch>[] allPsms = new List<PeptideSpectralMatch>[numFdrCategories];
            foreach (int localFDRCategory in SearchParameters.LocalFdrCategories) //only add if we're using for FDR, else ignore it as null.
            {
                allPsms[localFDRCategory] = new List<PeptideSpectralMatch>();
            }

            FlashLfqResults flashLfqResults = null;

            MyFileManager myFileManager = new MyFileManager(SearchParameters.DisposeOfFileWhenDone);

            var fileSpecificCommonParams = fileSettingsList.Select(b => SetAllFileSpecificCommonParams(CommonParameters, b));

            int completedFiles = 0;
            object indexLock = new object();
            object psmLock = new object();

            Status("Searching files...", taskId);
            Status("Searching files...", new List<string> { taskId, "Individual Spectra Files" });

            Dictionary<string, int[]> numMs2SpectraPerFile = new Dictionary<string, int[]>();
            for (int spectraFileIndex = 0; spectraFileIndex < currentRawFileList.Count; spectraFileIndex++)
            {
                if (GlobalVariables.StopLoops) { break; }

                var origDataFile = currentRawFileList[spectraFileIndex];

                // mark the file as in-progress
                StartingDataFile(origDataFile, new List<string> { taskId, "Individual Spectra Files", origDataFile });

                CommonParameters combinedParams = SetAllFileSpecificCommonParams(CommonParameters, fileSettingsList[spectraFileIndex]);

                MassDiffAcceptor massDiffAcceptor = GetMassDiffAcceptor(combinedParams.PrecursorMassTolerance, SearchParameters.MassDiffAcceptorType, SearchParameters.CustomMdac);

                var thisId = new List<string> { taskId, "Individual Spectra Files", origDataFile };
                NewCollection(Path.GetFileName(origDataFile), thisId);
                Status("Loading spectra file...", thisId);
                MsDataFile myMsDataFile = myFileManager.LoadFile(origDataFile, combinedParams.TopNpeaks, combinedParams.MinRatio, combinedParams.TrimMs1Peaks, combinedParams.TrimMsMsPeaks, combinedParams);
                Status("Getting ms2 scans...", thisId);
                Ms2ScanWithSpecificMass[] arrayOfMs2ScansSortedByMass = GetMs2Scans(myMsDataFile, origDataFile, combinedParams.DoPrecursorDeconvolution, combinedParams.UseProvidedPrecursorInfo, combinedParams.DeconvolutionIntensityRatio, combinedParams.DeconvolutionMaxAssumedChargeState, combinedParams.DeconvolutionMassTolerance).OrderBy(b => b.PrecursorMass).ToArray();
                numMs2SpectraPerFile.Add(Path.GetFileNameWithoutExtension(origDataFile), new int[] { myMsDataFile.GetAllScansList().Count(p => p.MsnOrder == 2), arrayOfMs2ScansSortedByMass.Length });
                myFileManager.DoneWithFile(origDataFile);

                PeptideSpectralMatch[][] fileSpecificPsms = new PeptideSpectralMatch[numFdrCategories][]; //generate an array of all possible locals
                foreach (int localFDRCategory in SearchParameters.LocalFdrCategories) //only add if we're using for FDR, else ignore it as null.
                {
                    fileSpecificPsms[localFDRCategory] = new PeptideSpectralMatch[arrayOfMs2ScansSortedByMass.Length];
                }

                // modern search
                if (SearchParameters.SearchType == SearchType.Modern)
                {
                    for (int currentPartition = 0; currentPartition < combinedParams.TotalPartitions; currentPartition++)
                    {
                        List<PeptideWithSetModifications> peptideIndex = null;
                        List<Protein> proteinListSubset = proteinList.GetRange(currentPartition * proteinList.Count() / combinedParams.TotalPartitions, ((currentPartition + 1) * proteinList.Count() / combinedParams.TotalPartitions) - (currentPartition * proteinList.Count() / combinedParams.TotalPartitions));

                        Status("Getting fragment dictionary...", new List<string> { taskId });
                        var indexEngine = new IndexingEngine(proteinListSubset, variableModifications, fixedModifications, currentPartition, SearchParameters.DecoyType, combinedParams, SearchParameters.MaxFragmentSize, false, new List<string> { taskId });
                        List<int>[] fragmentIndex = null;
                        List<int>[] precursorIndex = null;
                        lock (indexLock)
                        {
                            GenerateIndexes(indexEngine, dbFilenameList, ref peptideIndex, ref fragmentIndex, ref precursorIndex, proteinList, GlobalVariables.AllModsKnown.ToList(), taskId);
                        }

                        Status("Searching files...", taskId);

                        new ModernSearchEngine(fileSpecificPsms, arrayOfMs2ScansSortedByMass, peptideIndex, fragmentIndex, currentPartition, combinedParams, massDiffAcceptor, SearchParameters.MaximumMassThatFragmentIonScoreIsDoubled, thisId).Run();

                        ReportProgress(new ProgressEventArgs(100, "Done with search " + (currentPartition + 1) + "/" + combinedParams.TotalPartitions + "!", thisId));
                    }
                }
                // nonspecific search
                else if (SearchParameters.SearchType == SearchType.NonSpecific)
                {
                    List<CommonParameters> paramsToUse = new List<CommonParameters> { combinedParams };
                    if (combinedParams.DigestionParams.SearchModeType == CleavageSpecificity.Semi) //if semi, we need to do both N and C to hit everything
                    {
                        paramsToUse.Clear();
                        List<FragmentationTerminus> terminiToUse = new List<FragmentationTerminus> { FragmentationTerminus.N, FragmentationTerminus.C };
                        foreach (FragmentationTerminus terminus in terminiToUse) //set both termini
                        {
                            paramsToUse.Add(combinedParams.CloneWithNewTerminus(terminus));
                        }
                    }
                    foreach (CommonParameters paramToUse in paramsToUse)
                    {
                        for (int currentPartition = 0; currentPartition < paramToUse.TotalPartitions; currentPartition++)
                        {
                            List<PeptideWithSetModifications> peptideIndex = null;

                            List<Protein> proteinListSubset = proteinList.GetRange(currentPartition * proteinList.Count() / paramToUse.TotalPartitions, ((currentPartition + 1) * proteinList.Count() / paramToUse.TotalPartitions) - (currentPartition * proteinList.Count() / paramToUse.TotalPartitions));

                            List<int>[] fragmentIndex = new List<int>[1];
                            List<int>[] precursorIndex = new List<int>[1];

                            Status("Getting fragment dictionary...", new List<string> { taskId });
                            var indexEngine = new IndexingEngine(proteinListSubset, variableModifications, fixedModifications, currentPartition, SearchParameters.DecoyType, paramToUse, SearchParameters.MaxFragmentSize, true, new List<string> { taskId });
                            lock (indexLock)
                                GenerateIndexes(indexEngine, dbFilenameList, ref peptideIndex, ref fragmentIndex, ref precursorIndex, proteinList, GlobalVariables.AllModsKnown.ToList(), taskId);
                            
                            Status("Searching files...", taskId);

                            new NonSpecificEnzymeSearchEngine(fileSpecificPsms, arrayOfMs2ScansSortedByMass, peptideIndex, fragmentIndex, precursorIndex, currentPartition, paramToUse, massDiffAcceptor, SearchParameters.MaximumMassThatFragmentIonScoreIsDoubled, thisId).Run();

                            ReportProgress(new ProgressEventArgs(100, "Done with search " + (currentPartition + 1) + "/" + paramToUse.TotalPartitions + "!", thisId));
                        }
                    }
                }
                // classic search
                else
                {
                    Status("Starting search...", thisId);
                    new ClassicSearchEngine(fileSpecificPsms, arrayOfMs2ScansSortedByMass, variableModifications, fixedModifications, proteinList, massDiffAcceptor, combinedParams, thisId).Run();

                    ReportProgress(new ProgressEventArgs(100, "Done with search!", thisId));
                }

                lock (psmLock)
                {
                    for (int i = 0; i < allPsms.Length; i++)
                    {
                        if (allPsms[i] != null)
                        {
                            allPsms[i].AddRange(fileSpecificPsms[i]);
                        }
                    }
                }

                completedFiles++;
                FinishedDataFile(origDataFile, new List<string> { taskId, "Individual Spectra Files", origDataFile });
                ReportProgress(new ProgressEventArgs(completedFiles / currentRawFileList.Count, "Searching...", new List<string> { taskId, "Individual Spectra Files" }));
            }

            ReportProgress(new ProgressEventArgs(100, "Done with all searches!", new List<string> { taskId, "Individual Spectra Files" }));

            PostSearchAnalysisParameters parameters = new PostSearchAnalysisParameters();
            parameters.SearchTaskResults = MyTaskResults;
            parameters.SearchTaskId = taskId;
            parameters.SearchParameters = SearchParameters;
            parameters.ProteinList = proteinList;
            parameters.AllPsms = allPsms;
            parameters.FixedModifications = fixedModifications;
            parameters.VariableModifications = variableModifications;
            parameters.ListOfDigestionParams = new HashSet<DigestionParams>(fileSpecificCommonParams.Select(p => p.DigestionParams));
            parameters.CurrentRawFileList = currentRawFileList;
            parameters.MyFileManager = myFileManager;
            parameters.NumNotches = GetNumNotches(SearchParameters.MassDiffAcceptorType, SearchParameters.CustomMdac);
            parameters.OutputFolder = OutputFolder;
            parameters.IndividualResultsOutputFolder = Path.Combine(OutputFolder, "Individual File Results");
            parameters.FlashLfqResults = flashLfqResults;
            parameters.FileSettingsList = fileSettingsList;
            parameters.NumMs2SpectraPerFile = numMs2SpectraPerFile;
            parameters.DatabaseFilenameList = dbFilenameList;
            PostSearchAnalysisTask postProcessing = new PostSearchAnalysisTask
            {
                Parameters = parameters,
                CommonParameters = CommonParameters
            };
            return postProcessing.Run();
        }

        private int GetNumNotches(MassDiffAcceptorType massDiffAcceptorType, string customMdac)
        {
            switch (massDiffAcceptorType)
            {
                case MassDiffAcceptorType.Exact: return 1;
                case MassDiffAcceptorType.OneMM: return 2;
                case MassDiffAcceptorType.TwoMM: return 3;
                case MassDiffAcceptorType.ThreeMM: return 4;
                case MassDiffAcceptorType.ModOpen: return 1;
                case MassDiffAcceptorType.Open: return 1;
                case MassDiffAcceptorType.Custom: return ParseSearchMode(customMdac).NumNotches;

                default: throw new MetaMorpheusException("Unknown MassDiffAcceptorType");
            }
        }

        private static MassDiffAcceptor ParseSearchMode(string text)
        {
            MassDiffAcceptor ye = null;

            var split = text.Split(' ');

            switch (split[1])
            {
                case "dot":

                    var massShifts = Array.ConvertAll(split[4].Split(','), Double.Parse);
                    var newString = split[2].Replace("�", "");
                    var toleranceValue = double.Parse(newString, CultureInfo.InvariantCulture);
                    if (split[3].ToUpperInvariant().Equals("PPM"))
                        ye = new DotMassDiffAcceptor(split[0], massShifts, new PpmTolerance(toleranceValue));
                    else if (split[3].ToUpperInvariant().Equals("DA"))
                        ye = new DotMassDiffAcceptor(split[0], massShifts, new AbsoluteTolerance(toleranceValue));
                    break;

                case "interval":
                    IEnumerable<DoubleRange> doubleRanges = Array.ConvertAll(split[2].Split(','), b => new DoubleRange(double.Parse(b.Trim(new char[] { '[', ']' }).Split(';')[0], CultureInfo.InvariantCulture), double.Parse(b.Trim(new char[] { '[', ']' }).Split(';')[1], CultureInfo.InvariantCulture)));
                    ye = new IntervalMassDiffAcceptor(split[0], doubleRanges);
                    break;

                case "OpenSearch":
                    ye = new OpenSearchMode();
                    break;

                case "daltonsAroundZero":
                    ye = new SingleAbsoluteAroundZeroSearchMode(double.Parse(split[2], CultureInfo.InvariantCulture));
                    break;

                case "ppmAroundZero":
                    ye = new SinglePpmAroundZeroSearchMode(double.Parse(split[2], CultureInfo.InvariantCulture));
                    break;

                default:
                    throw new MetaMorpheusException("Could not parse search mode string");
            }
            return ye;
        }
    }
}