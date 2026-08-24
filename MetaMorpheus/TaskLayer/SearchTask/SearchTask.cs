using Easy.Common.Extensions;
using EngineLayer;
using EngineLayer.ClassicSearch;
using EngineLayer.DatabaseLoading;
using EngineLayer.Indexing;
using EngineLayer.ModernSearch;
using EngineLayer.NonSpecificEnzymeSearch;
using EngineLayer.Util;
using FlashLFQ;
using MassSpectrometry;
using MzLibUtil;
using Omics;
using Omics.Digestion;
using Omics.Fragmentation;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using Readers;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Threading.Tasks;
using Chemistry;

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

        /// <summary>
        /// Builds the precursor search mode for a <see cref="MassDiffAcceptorType"/>.
        /// </summary>
        /// <remarks>
        /// Deliberately takes no <see cref="PrecursorMassMatchMode"/>. Whether candidates are matched on the
        /// monoisotopic mass or the most-abundant isotopic mass is now carried entirely by
        /// <paramref name="massDiffAcceptorType"/> (the <c>MostAbundant_*</c> members), so accepting the match
        /// mode as well would let a caller ask for something this method cannot honour.
        /// <see cref="RunSpecific"/> derives <see cref="CommonParameters.PrecursorMassMatchMode"/> from the
        /// acceptor type instead, which is what keeps the theory and observed sides in step.
        /// </remarks>
        public static MassDiffAcceptor GetMassDiffAcceptor(Tolerance precursorMassTolerance, MassDiffAcceptorType massDiffAcceptorType, string customMdac,
            AverageResidue averagineModel = null, double expectedIsotopeSpacing = Constants.C13MinusC12)
        {
            averagineModel ??= GlobalVariables.AnalyteType == AnalyteType.Oligo
                ? new OxyriboAveragine()
                : new Averagine();

            switch (massDiffAcceptorType)
            {
                case MassDiffAcceptorType.Exact:
                    if (precursorMassTolerance is PpmTolerance)
                    {
                        return new SinglePpmAroundZeroSearchMode(precursorMassTolerance.Value);
                    }
                    else
                    {
                        return new SingleAbsoluteAroundZeroSearchMode(precursorMassTolerance.Value);
                    }

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

                case MassDiffAcceptorType.PlusOrMinusThreeMM:
                    return new DotMassDiffAcceptor(
                        "PlusOrMinus3Da",
                        new List<double>
                        {
                            -3 * Chemistry.Constants.C13MinusC12,
                            -2 * Chemistry.Constants.C13MinusC12,
                            -1 * Chemistry.Constants.C13MinusC12,
                            0,
                            1 * Chemistry.Constants.C13MinusC12,
                            2 * Chemistry.Constants.C13MinusC12,
                            3 * Chemistry.Constants.C13MinusC12
                        },
                        precursorMassTolerance);

                case MassDiffAcceptorType.MostAbundant_Exact:
                    return new MostAbundantMassDiffAcceptor("mostAbundant", precursorMassTolerance, averagineModel, 0, expectedIsotopeSpacing);

                case MassDiffAcceptorType.MostAbundant_PlusMinusOne:
                    return new MostAbundantMassDiffAcceptor("mostAbundant_1", precursorMassTolerance, averagineModel, 1, expectedIsotopeSpacing);

                case MassDiffAcceptorType.MostAbundant_PlusMinusTwo:
                    return new MostAbundantMassDiffAcceptor("mostAbundant_2", precursorMassTolerance, averagineModel, 2, expectedIsotopeSpacing);

                default:
                    throw new MetaMorpheusException("Unknown MassDiffAcceptorType");
            }
        }

        /// <summary>
        /// The acceptor type equivalent to the pre-<c>MostAbundant_*</c> way of asking for a most-abundant
        /// search, i.e. <c>CommonParameters.PrecursorMassMatchMode = MostAbundant</c>. That setting used to
        /// override <see cref="MassDiffAcceptorType"/> outright and build a
        /// <see cref="MostAbundantMassDiffAcceptor"/> at <see cref="MostAbundantMassDiffAcceptor.DefaultMaxApexOffsetNeutrons"/>
        /// (2) neutrons of apex tolerance, which is exactly what <see cref="MassDiffAcceptorType.MostAbundant_PlusMinusTwo"/>
        /// builds today.
        /// </summary>
        private const MassDiffAcceptorType LegacyMostAbundantEquivalent = MassDiffAcceptorType.MostAbundant_PlusMinusTwo;

        /// <summary>
        /// Migrates a run that asks for most-abundant matching the old way — through
        /// <see cref="CommonParameters.PrecursorMassMatchMode"/> rather than a <c>MostAbundant_*</c>
        /// <see cref="MassDiffAcceptorType"/> — onto the acceptor type that now carries that request.
        /// </summary>
        /// <remarks>
        /// The GUI migrates on open, so this only reaches TOML and command-line runs — which is where
        /// most existing most-abundant configurations live, because the match mode was the only way to ask
        /// for one before the <c>MostAbundant_*</c> members existed. Without this, the acceptor built below
        /// would be monoisotopic and the run would come back silently changed. Migrating rather than merely
        /// warning is what keeps those runs searching the way they used to; the warning is what stops the
        /// change from being silent, and is the only record of it — <see cref="MetaMorpheusTask.RunTask"/>
        /// writes the settings toml before it calls <see cref="RunSpecific"/>, so that file still shows the
        /// original pairing.
        /// </remarks>
        private void MigrateLegacyMostAbundantRequest()
        {
            if (CommonParameters.PrecursorMassMatchMode != PrecursorMassMatchMode.MostAbundant
                || SearchParameters.MassDiffAcceptorType.IsMostAbundant())
            {
                return;
            }

            Warn($"PrecursorMassMatchMode is {PrecursorMassMatchMode.MostAbundant} but MassDiffAcceptorType is " +
                 $"{SearchParameters.MassDiffAcceptorType}, which is monoisotopic. Most-abundant matching is now selected " +
                 $"by the mass difference acceptor, so this search will use {LegacyMostAbundantEquivalent} — the equivalent " +
                 $"of the old setting. Set MassDiffAcceptorType explicitly to silence this.");

            SearchParameters.MassDiffAcceptorType = LegacyMostAbundantEquivalent;
        }

        protected override MyTaskResults RunSpecific(string OutputFolder, List<DbForTask> dbFilenameList, List<string> currentRawFileList, string taskId,
            FileSpecificParameters[] fileSettingsList)
        {
            MigrateLegacyMostAbundantRequest();

            MyTaskResults = new(this);
            MyFileManager myFileManager = new MyFileManager(SearchParameters.DisposeOfFileWhenDone);
            var fileSpecificCommonParams = fileSettingsList.Select(b => SetAllFileSpecificCommonParams(CommonParameters, b));

            // start loading first spectra file in the background
            string fileToLoad = currentRawFileList[0];
            Task<MsDataFile> nextFileLoadingTask = new(() => myFileManager.LoadFile(fileToLoad, SetAllFileSpecificCommonParams(CommonParameters, fileSettingsList[0])));
            nextFileLoadingTask.Start();

            if (SearchParameters.DoLabelFreeQuantification)
            {
                // disable quantification if a .mgf or Bruker file is being used. Bruker paths are normally the ".d"
                // folder, but an inner file (analysis.baf, analysis.tdf, ...) can reach here too, so check for both.
                if (currentRawFileList.Select(filepath => Path.GetExtension(filepath))
                    .Any(ext => ext.Equals(".mgf", StringComparison.OrdinalIgnoreCase) || ext.Equals(".d", StringComparison.OrdinalIgnoreCase) || ext.Equals(".msalign", StringComparison.OrdinalIgnoreCase) || BrukerDataDirectory.IsInnerFileExtension(ext)))
                {
                    SearchParameters.DoLabelFreeQuantification = false;
                }
                //if we're doing SILAC, assign and add the silac labels to the residue dictionary
                else if (SearchParameters.SilacLabels != null || SearchParameters.StartTurnoverLabel != null || SearchParameters.EndTurnoverLabel != null)
                {
                    char heavyLabel = 'a'; //char to assign
                    //add the Turnoverlabels to the silacLabels list. They weren't there before just to prevent duplication in the tomls
                    if (SearchParameters.StartTurnoverLabel != null || SearchParameters.EndTurnoverLabel != null)
                    {
                        //original silacLabels object is null, so we need to initialize it
                        SearchParameters.SilacLabels = new List<SilacLabel>();
                        if (SearchParameters.StartTurnoverLabel != null)
                        {
                            var updatedLabel = SilacConversions.UpdateAminoAcidLabel(SearchParameters.StartTurnoverLabel, heavyLabel);
                            heavyLabel = updatedLabel.nextHeavyLabel;
                            SearchParameters.StartTurnoverLabel = updatedLabel.updatedLabel;
                            SearchParameters.SilacLabels.Add(SearchParameters.StartTurnoverLabel);
                        }
                        if (SearchParameters.EndTurnoverLabel != null)
                        {
                            var updatedLabel = SilacConversions.UpdateAminoAcidLabel(SearchParameters.EndTurnoverLabel, heavyLabel);
                            heavyLabel = updatedLabel.nextHeavyLabel;
                            SearchParameters.EndTurnoverLabel = updatedLabel.updatedLabel;
                            SearchParameters.SilacLabels.Add(SearchParameters.EndTurnoverLabel);
                        }
                    }
                    else
                    {
                        //change the silac residues to lower case amino acids (currently null)
                        List<SilacLabel> updatedLabels = new List<SilacLabel>();
                        for (int i = 0; i < SearchParameters.SilacLabels.Count; i++)
                        {
                            var updatedLabel = SilacConversions.UpdateAminoAcidLabel(SearchParameters.SilacLabels[i], heavyLabel);
                            heavyLabel = updatedLabel.nextHeavyLabel;
                            updatedLabels.Add(updatedLabel.updatedLabel);
                        }
                        SearchParameters.SilacLabels = updatedLabels;
                    }
                }
            }
            //if no quant, remove any silac labels that may have been added, because they screw up downstream analysis
            if (!SearchParameters.DoLabelFreeQuantification) //using "if" instead of "else", because DoLabelFreeQuantification can change if it's an mgf
            {
                SearchParameters.SilacLabels = null;
            }

            // start loading all data in the background while task is being set up
            LoadModifications(taskId, out var variableModifications, out var fixedModifications, out var localizeableModificationTypes);

            // start loading proteins in the background
            var dbLoader = new DatabaseLoadingEngine(CommonParameters, this.FileSpecificParameters, [taskId], dbFilenameList, taskId, SearchParameters.DecoyType, SearchParameters.SearchTarget, localizeableModificationTypes, SearchParameters.TCAmbiguity, SearchParameters.WriteTargetDecoyFasta, OutputFolder);
            var proteinLoadingTask = dbLoader.RunAsync();
            List<IBioPolymer> bioPolymerList = null!;

            // load spectral libraries
            var spectralLibrary = LoadSpectralLibraries(taskId, dbFilenameList);

            // write prose settings
            ProseCreatedWhileRunning.Append("The following search settings were used: ");
            ProseCreatedWhileRunning.Append($"{GlobalVariables.AnalyteType.GetDigestionAgentLabel()} = " + CommonParameters.DigestionParams.DigestionAgent + "; ");
            ProseCreatedWhileRunning.Append($"search for truncated {GlobalVariables.AnalyteType.GetBioPolymerLabel().ToLower()} and proteolysis products = " + CommonParameters.AddTruncations + "; ");
            ProseCreatedWhileRunning.Append("maximum missed cleavages = " + CommonParameters.DigestionParams.MaxMissedCleavages + "; ");
            ProseCreatedWhileRunning.Append($"minimum {GlobalVariables.AnalyteType.GetUniqueFormLabel().ToLower()} length = " + CommonParameters.DigestionParams.MinLength + "; ");
            ProseCreatedWhileRunning.Append(CommonParameters.DigestionParams.MaxLength == int.MaxValue ?
                $"maximum {GlobalVariables.AnalyteType.GetUniqueFormLabel().ToLower()} length = unspecified; " :
                $"maximum {GlobalVariables.AnalyteType.GetUniqueFormLabel().ToLower()} length = " + CommonParameters.DigestionParams.MaxLength + "; ");
            if (CommonParameters.DigestionParams is DigestionParams digestionParams)
                ProseCreatedWhileRunning.Append("initiator methionine behavior = " + digestionParams.InitiatorMethionineBehavior + "; \n");
            ProseCreatedWhileRunning.Append("fixed modifications = " + string.Join(", ", fixedModifications.Select(m => m.IdWithMotif)) + "; ");
            ProseCreatedWhileRunning.Append("variable modifications = " + string.Join(", ", variableModifications.Select(m => m.IdWithMotif)) + "; ");
            ProseCreatedWhileRunning.Append($"max mods per {GlobalVariables.AnalyteType.GetUniqueFormLabel().ToLower()} = " + CommonParameters.DigestionParams.MaxMods + "; ");
            ProseCreatedWhileRunning.Append("max modification isoforms = " + CommonParameters.DigestionParams.MaxModificationIsoforms + "; ");
            ProseCreatedWhileRunning.Append("precursor mass tolerance = " + CommonParameters.PrecursorMassTolerance + "; ");
            ProseCreatedWhileRunning.Append("product mass tolerance = " + CommonParameters.ProductMassTolerance + "; ");
            ProseCreatedWhileRunning.Append($"report {GlobalVariables.AnalyteType.GetSpectralMatchLabel()} ambiguity = " + CommonParameters.ReportAllAmbiguity + ". ");

            // start the search task
            List<SpectralMatch> allPsms = new List<SpectralMatch>();

            //generate an array to store category specific fdr values (for speedy semi/nonspecific searches)
            int numFdrCategories = (int)(Enum.GetValues(typeof(FdrCategory)).Cast<FdrCategory>().Last() + 1); //+1 because it starts at zero
            List<SpectralMatch>[] allCategorySpecificPsms = new List<SpectralMatch>[numFdrCategories];
            for (int i = 0; i < numFdrCategories; i++)
            {
                allCategorySpecificPsms[i] = new List<SpectralMatch>();
            }

            FlashLfqResults flashLfqResults = null;

            int completedFiles = 0;
            object indexLock = new object();
            object psmLock = new object();
            int? decidedPartitions = null;

            Status("Searching files...", new List<string> { taskId });
            Status("Searching files...", new List<string> { taskId, "Individual Spectra Files" });

            Dictionary<string, int[]> numMs2SpectraPerFile = new Dictionary<string, int[]>(); // key is filename, value is an int array of length 2, where the first element is the number of MS2 spectra in the file, and the second element is the number of different deconvoluted precursors assigned to those scans
            bool collectedDigestionInformation = false;
            IDictionary<(string Accession, string BaseSequence), int> digestionCountDictionary = null;
            int numNotches = 0;
            for (int spectraFileIndex = 0; spectraFileIndex < currentRawFileList.Count; spectraFileIndex++)
            {
                if (GlobalVariables.StopLoops) { break; }

                var origDataFile = currentRawFileList[spectraFileIndex];

                // mark the file as in-progress
                StartingDataFile(origDataFile, new List<string> { taskId, "Individual Spectra Files", origDataFile });

                CommonParameters combinedParams = SetAllFileSpecificCommonParams(CommonParameters, fileSettingsList[spectraFileIndex]);

                // The theory side is governed by SearchParameters.MassDiffAcceptorType, while the observed precursor
                // mass comes from PrecursorMassMatchMode. Keep them aligned so a most-abundant acceptor actually
                // searches against most-abundant masses (and vice versa).
                combinedParams.PrecursorMassMatchMode = SearchParameters.MassDiffAcceptorType.IsMostAbundant()
                    ? PrecursorMassMatchMode.MostAbundant
                    : PrecursorMassMatchMode.Monoisotopic;

                MassDiffAcceptor massDiffAcceptor = GetMassDiffAcceptor(combinedParams.PrecursorMassTolerance, SearchParameters.MassDiffAcceptorType, SearchParameters.CustomMdac,
                    combinedParams.GetAverageResidue(), combinedParams.IsotopeSpacing());

                numNotches = massDiffAcceptor.NumNotches;

                var thisId = new List<string> { taskId, "Individual Spectra Files", origDataFile };
                NewCollection(Path.GetFileName(origDataFile), thisId);
                Status("Loading spectra file...", thisId);

                // ensure that the next file has finished loading from the async method
                nextFileLoadingTask.Wait();
                var myMsDataFile = nextFileLoadingTask.Result;

                // If the file is one which does not have precursor scans, but only precursor information, then we need to set the parameters accordingly
                // We do this by adjusting the transient combined params so that this can be done on a file by file basis. 
                if (myMsDataFile is Mgf or Ms2Align)
                {
                    combinedParams.DoPrecursorDeconvolution = false;
                    combinedParams.UseProvidedPrecursorInfo = true;
                }

                // If we're doing multiplex quantification, and there are MS3 scans, we assume that
                // MS3 was used for reporter ion detection, and adjust the parameters accordingly
                if (SearchParameters.DoMultiplexQuantification && myMsDataFile.Scans.Any(s => s.MsnOrder == 3))
                {
                    // In most experiments with MS3 scans for reporter ion detection, MS2ChildScanDissociationType is LowCID.
                    // However, we do not set it here to allow for flexibility in dissociation type selection.
                    combinedParams.MS3ChildScanDissociationType = DissociationType.HCD;
                }

                // if another file exists, then begin loading it in while the previous is being searched
                if (origDataFile != currentRawFileList.Last())
                {
                    int nextFileIndex = spectraFileIndex + 1;
                    nextFileLoadingTask = new Task<MsDataFile>(() => myFileManager.LoadFile(currentRawFileList[nextFileIndex], SetAllFileSpecificCommonParams(CommonParameters, fileSettingsList[nextFileIndex])));
                    nextFileLoadingTask.Start();
                }

                Status("Getting ms2 scans...", thisId);
                // Sort by the mass this search selects candidates on. ClassicSearchEngine binary-searches
                // that same mass, so the array and the search key must use the same quantity.
                Ms2ScanWithSpecificMass[] arrayOfMs2ScansSortedByMass = GetMs2Scans(myMsDataFile, origDataFile, combinedParams).OrderBy(b => b.GetPrecursorMassForSearch(combinedParams)).ToArray();
                // A most-abundant search falls back to the monoisotopic mass for any scan with no observed
                // apex. That is intended, but silent — report how often it happens so it is noticed.
                string fallbackWarning = arrayOfMs2ScansSortedByMass.GetMonoisotopicFallbackWarning(combinedParams, Path.GetFileName(origDataFile));
                if (fallbackWarning != null) { Warn(fallbackWarning); }
                numMs2SpectraPerFile.Add(Path.GetFileNameWithoutExtension(origDataFile), new int[] { myMsDataFile.GetAllScansList().Count(p => p.MsnOrder == 2), arrayOfMs2ScansSortedByMass.Length });
                myFileManager.DoneWithFile(origDataFile);

                SpectralMatch[] fileSpecificPsms = new SpectralMatch[arrayOfMs2ScansSortedByMass.Length];

                // ensure proteins are loaded in before proceeding with search
                switch (proteinLoadingTask.IsCompleted)
                {
                    case true when bioPolymerList is null: // has finished loading but not been set
                        bioPolymerList = (proteinLoadingTask.Result as DatabaseLoadingEngineResults).BioPolymers;
                        Status("Searching files...", new List<string> { taskId });
                        break;
                    case true when bioPolymerList.Any(): // has finished loading and already been set
                        break;
                    case false: // has not finished loading
                        proteinLoadingTask.Wait();
                        bioPolymerList = (proteinLoadingTask.Result as DatabaseLoadingEngineResults).BioPolymers;
                        Status("Searching files...", new List<string> { taskId });
                        break;
                }

                if (SearchParameters.DoMultiplexQuantification)
                {
                    IsobaricMassTag massTag = IsobaricMassTag.GetIsobaricMassTag(SearchParameters.MultiplexModId);
                    if (massTag == null) // Should probably warn/update results if null
                    {
                        throw new MetaMorpheusException("Could not find isobaric mass tag with the name " + SearchParameters.MultiplexModId);
                    }

                    foreach (var scan in arrayOfMs2ScansSortedByMass)
                    {
                        scan.SetIsobaricMassTagReporterIonIntensities(massTag);
                    }
                }

                // modern search
                if (SearchParameters.SearchType == SearchType.Modern)
                {
                    // Assume modern search is for proteins.
                    var proteinList = bioPolymerList.Cast<Protein>().ToList();
                    // scoped to indexing/searching only, so the settings the task reports stay as configured
                    CommonParameters indexParams = RaisePartitionsToFitMemory(proteinList, combinedParams, fixedModifications,
                        variableModifications, SearchParameters.SilacLabels, SearchParameters.StartTurnoverLabel,
                        SearchParameters.EndTurnoverLabel, SearchParameters.MaxFragmentSize, ref decidedPartitions);
                    for (int currentPartition = 0; currentPartition < indexParams.TotalPartitions; currentPartition++)
                    {
                        List<PeptideWithSetModifications> peptideIndex = null;
                        List<Protein> proteinListSubset = proteinList.GetRange(currentPartition * proteinList.Count / indexParams.TotalPartitions,
                            ((currentPartition + 1) * proteinList.Count / indexParams.TotalPartitions) - (currentPartition * proteinList.Count / indexParams.TotalPartitions));

                        Status("Getting fragment dictionary...", new List<string> { taskId });
                        var indexEngine = new IndexingEngine(proteinListSubset, variableModifications, fixedModifications, SearchParameters.SilacLabels,
                            SearchParameters.StartTurnoverLabel, SearchParameters.EndTurnoverLabel, currentPartition, SearchParameters.DecoyType, indexParams, FileSpecificParameters,
                            SearchParameters.MaxFragmentSize, false, dbFilenameList.Select(p => new FileInfo(p.FilePath)).ToList(), SearchParameters.TCAmbiguity, new List<string> { taskId });
                        FragmentIndex fragmentIndex = null;
                        List<int>[] precursorIndex = null;

                        lock (indexLock)
                        {
                            GenerateIndexes(indexEngine, dbFilenameList, ref peptideIndex, ref fragmentIndex, ref precursorIndex, proteinList, taskId);
                        }

                        Status("Searching files...", taskId);

                        new ModernSearchEngine(fileSpecificPsms, arrayOfMs2ScansSortedByMass, peptideIndex, fragmentIndex, currentPartition,
                            indexParams, this.FileSpecificParameters, massDiffAcceptor, SearchParameters.MaximumMassThatFragmentIonScoreIsDoubled, thisId).Run();

                        ReportProgress(new ProgressEventArgs(100, "Done with search " + (currentPartition + 1) + "/" + indexParams.TotalPartitions + "!", thisId));
                        if (GlobalVariables.StopLoops) { break; }
                    }
                }
                // nonspecific search
                else if (SearchParameters.SearchType == SearchType.NonSpecific)
                {
                    SpectralMatch[][] fileSpecificPsmsSeparatedByFdrCategory = new SpectralMatch[numFdrCategories][]; //generate an array of all possible locals
                    for (int i = 0; i < numFdrCategories; i++) //only add if we're using for FDR, else ignore it as null.
                    {
                        fileSpecificPsmsSeparatedByFdrCategory[i] = new SpectralMatch[arrayOfMs2ScansSortedByMass.Length];
                    }

                    //create params for N, C, or both if semi
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

                    //Compress array of deconvoluted ms2 scans to avoid searching the same ms2 multiple times while still identifying coisolated peptides
                    List<int>[] coisolationIndex = new List<int>[] { new List<int>() };
                    if (arrayOfMs2ScansSortedByMass.Length != 0)
                    {
                        int maxScanNumber = arrayOfMs2ScansSortedByMass.Max(x => x.OneBasedScanNumber);
                        coisolationIndex = new List<int>[maxScanNumber + 1];
                        for (int i = 0; i < arrayOfMs2ScansSortedByMass.Length; i++)
                        {
                            int scanNumber = arrayOfMs2ScansSortedByMass[i].OneBasedScanNumber;
                            if (coisolationIndex[scanNumber] == null)
                            {
                                coisolationIndex[scanNumber] = new List<int> { i };
                            }
                            else
                            {
                                coisolationIndex[scanNumber].Add(i);
                            }
                        }
                        coisolationIndex = coisolationIndex.Where(x => x != null).ToArray();
                    }

                    //foreach terminus we're going to look at
                    foreach (CommonParameters paramToUse in paramsToUse)
                    {
                        var proteinList = bioPolymerList.Cast<Protein>().ToList();
                        // scoped to indexing/searching only; paramToUse still carries the configured terminus
                        // to the spectral-library step below
                        CommonParameters indexParams = RaisePartitionsToFitMemory(proteinList, paramToUse, fixedModifications,
                            variableModifications, SearchParameters.SilacLabels, SearchParameters.StartTurnoverLabel,
                            SearchParameters.EndTurnoverLabel, SearchParameters.MaxFragmentSize, ref decidedPartitions);

                        //foreach database partition
                        for (int currentPartition = 0; currentPartition < indexParams.TotalPartitions; currentPartition++)
                        {
                            List<PeptideWithSetModifications> peptideIndex = null;

                            List<Protein> proteinListSubset = proteinList.GetRange(currentPartition * proteinList.Count / indexParams.TotalPartitions,
                                ((currentPartition + 1) * proteinList.Count / indexParams.TotalPartitions) - (currentPartition * proteinList.Count / indexParams.TotalPartitions))
                                .ToList(); // assume that only proteins are used in non-specific search

                            FragmentIndex fragmentIndex = null;
                            List<int>[] precursorIndex = null;

                            Status("Getting fragment dictionary...", new List<string> { taskId });
                            var indexEngine = new IndexingEngine(proteinListSubset, variableModifications, fixedModifications, SearchParameters.SilacLabels,
                                SearchParameters.StartTurnoverLabel, SearchParameters.EndTurnoverLabel, currentPartition, SearchParameters.DecoyType, indexParams, FileSpecificParameters,
                                SearchParameters.MaxFragmentSize, true, dbFilenameList.Select(p => new FileInfo(p.FilePath)).ToList(), SearchParameters.TCAmbiguity, new List<string> { taskId });
                            lock (indexLock)
                            {
                                GenerateIndexes(indexEngine, dbFilenameList, ref peptideIndex, ref fragmentIndex, ref precursorIndex, proteinList, taskId);
                            }

                            Status("Searching files...", taskId);

                            new NonSpecificEnzymeSearchEngine(fileSpecificPsmsSeparatedByFdrCategory, arrayOfMs2ScansSortedByMass, coisolationIndex, peptideIndex, fragmentIndex,
                                precursorIndex, currentPartition, indexParams, this.FileSpecificParameters, variableModifications, massDiffAcceptor,
                                SearchParameters.MaximumMassThatFragmentIonScoreIsDoubled, thisId).Run();

                            ReportProgress(new ProgressEventArgs(100, "Done with search " + (currentPartition + 1) + "/" + indexParams.TotalPartitions + "!", thisId));
                            if (GlobalVariables.StopLoops) { break; }
                        }

                        // Calculate spectral angles for non-specific search after each terminus iteration
                        // This must happen before FDR analysis since spectral angles are used in PEP calculation
                        // We use paramToUse (not combinedParams) because the FragmentationTerminus is used when 
                        // generating theoretical fragments for decoy spectra fallback in CalculateSpectralAngles
                        if (spectralLibrary != null)
                        {
                            Status("Calculating spectral library similarity...", thisId);
                            foreach (var categoryPsms in fileSpecificPsmsSeparatedByFdrCategory)
                            {
                                if (categoryPsms != null)
                                {
                                    SpectralLibrarySearchFunction.CalculateSpectralAngles(spectralLibrary, categoryPsms, arrayOfMs2ScansSortedByMass, paramToUse);
                                }
                            }
                            ReportProgress(new ProgressEventArgs(100, "Done with spectral library similarity!", thisId));
                        }
                    }

                    lock (psmLock)
                    {
                        for (int i = 0; i < allCategorySpecificPsms.Length; i++)
                        {
                            if (allCategorySpecificPsms[i] != null)
                            {
                                allCategorySpecificPsms[i].AddRange(fileSpecificPsmsSeparatedByFdrCategory[i]);
                            }
                        }
                    }
                }
                // classic search
                else
                {
                    Status("Starting search...", thisId);
                    var newClassicSearchEngine = new ClassicSearchEngine(fileSpecificPsms, arrayOfMs2ScansSortedByMass, variableModifications, fixedModifications, SearchParameters.SilacLabels,
                       SearchParameters.StartTurnoverLabel, SearchParameters.EndTurnoverLabel, bioPolymerList, massDiffAcceptor, combinedParams, this.FileSpecificParameters, spectralLibrary, thisId, SearchParameters.WriteSpectralLibrary, SearchParameters.WriteDigestionProductCountFile);
                    var result = newClassicSearchEngine.Run();

                    // The same proteins (all of them) get digested with each classic search engine, therefor we only need to calculate this for the first file that runs
                    if (!collectedDigestionInformation)
                    {
                        collectedDigestionInformation = true;
                        digestionCountDictionary = (result.MyEngine as ClassicSearchEngine).DigestionCountDictionary;
                    }

                    ReportProgress(new ProgressEventArgs(100, "Done with search!", thisId));
                }

                //look for internal fragments
                if (SearchParameters.MinAllowedInternalFragmentLength != 0)
                {
                    MatchInternalFragmentIons(fileSpecificPsms, arrayOfMs2ScansSortedByMass, combinedParams, SearchParameters.MinAllowedInternalFragmentLength);
                }

                // calculate/set spectral angles if there is a spectral library being used
                // Note: Non-specific search handles spectral angle calculation separately (before adding to allCategorySpecificPsms)
                if (spectralLibrary != null && SearchParameters.SearchType != SearchType.NonSpecific)
                {
                    Status("Calculating spectral library similarity...", thisId);
                    SpectralLibrarySearchFunction.CalculateSpectralAngles(spectralLibrary, fileSpecificPsms, arrayOfMs2ScansSortedByMass, combinedParams);
                    ReportProgress(new ProgressEventArgs(100, "Done with search!", thisId));
                }

                lock (psmLock)
                {
                    allPsms.AddRange(fileSpecificPsms);
                }

                completedFiles++;
                FinishedDataFile(origDataFile, new List<string> { taskId, "Individual Spectra Files", origDataFile });
                ReportProgress(new ProgressEventArgs(completedFiles / currentRawFileList.Count, "Searching...", new List<string> { taskId, "Individual Spectra Files" }));
            }

            if (spectralLibrary != null && SearchParameters.UpdateSpectralLibrary == false)
            {
                spectralLibrary.CloseConnections();
            }

            ReportProgress(new ProgressEventArgs(100, "Done with all searches!", new List<string> { taskId, "Individual Spectra Files" }));

            //resolve category specific fdrs (for speedy semi and nonspecific
            if (SearchParameters.SearchType == SearchType.NonSpecific)
            {
                allPsms = NonSpecificEnzymeSearchEngine.ResolveFdrCategorySpecificPsms(allCategorySpecificPsms, numNotches, taskId, CommonParameters, FileSpecificParameters);
            }

            // Finish writing prose settings that depended on files being loaded in
            ProseCreatedWhileRunning.Append("The combined search database contained " + bioPolymerList.Count(p => !p.IsDecoy)
                + $" non-decoy {GlobalVariables.AnalyteType.GetBioPolymerLabel().ToLower()} entries including " + bioPolymerList.Count(p => p.IsContaminant) + " contaminant sequences. ");

            PostSearchAnalysisParameters parameters = new PostSearchAnalysisParameters
            {
                SearchTaskResults = MyTaskResults,
                SearchTaskId = taskId,
                SearchParameters = SearchParameters,
                BioPolymerList = bioPolymerList,
                AllSpectralMatches = allPsms,
                VariableModifications = variableModifications,
                FixedModifications = fixedModifications,
                ListOfDigestionParams = [.. fileSpecificCommonParams.Select(p => p.DigestionParams)],
                CurrentRawFileList = currentRawFileList,
                MyFileManager = myFileManager,
                NumNotches = numNotches,
                OutputFolder = OutputFolder,
                IndividualResultsOutputFolder = Path.Combine(OutputFolder, "Individual File Results"),
                FlashLfqResults = flashLfqResults,
                FileSettingsList = fileSettingsList,
                NumMs2SpectraPerFile = numMs2SpectraPerFile,
                DatabaseFilenameList = dbFilenameList,
                SpectralLibrary = spectralLibrary
            };
            PostSearchAnalysisTask postProcessing = new PostSearchAnalysisTask
            {
                Parameters = parameters,
                FileSpecificParameters = this.FileSpecificParameters,
                CommonParameters = CommonParameters,
                DigestionCountDictionary = digestionCountDictionary
            };
            return postProcessing.Run();
        }

        private static MassDiffAcceptor ParseSearchMode(string text)
        {
            MassDiffAcceptor massDiffAcceptor = null;

            try
            {
                var split = text.Split(' ');

                switch (split[1])
                {
                    case "dot":
                        double[] massShifts = split[4].Split(',').Select(p => double.Parse(p, CultureInfo.InvariantCulture)).ToArray();
                        string newString = split[2].Replace("�", "");
                        double toleranceValue = double.Parse(newString, CultureInfo.InvariantCulture);
                        if (split[3].ToUpperInvariant().Equals("PPM"))
                        {
                            massDiffAcceptor = new DotMassDiffAcceptor(split[0], massShifts, new PpmTolerance(toleranceValue));
                        }
                        else if (split[3].ToUpperInvariant().Equals("DA"))
                        {
                            massDiffAcceptor = new DotMassDiffAcceptor(split[0], massShifts, new AbsoluteTolerance(toleranceValue));
                        }

                        break;

                    case "interval":
                        IEnumerable<DoubleRange> doubleRanges = Array.ConvertAll(split[2].Split(';'), b => new DoubleRange(double.Parse(b.Trim(new char[] { '[', ']' }).Split(',')[0],
                            CultureInfo.InvariantCulture), double.Parse(b.Trim(new char[] { '[', ']' }).Split(',')[1], CultureInfo.InvariantCulture)));
                        massDiffAcceptor = new IntervalMassDiffAcceptor(split[0], doubleRanges);
                        break;

                    case "OpenSearch":
                        massDiffAcceptor = new OpenSearchMode();
                        break;

                    case "daltonsAroundZero":
                        massDiffAcceptor = new SingleAbsoluteAroundZeroSearchMode(double.Parse(split[2], CultureInfo.InvariantCulture));
                        break;

                    case "ppmAroundZero":
                        massDiffAcceptor = new SinglePpmAroundZeroSearchMode(double.Parse(split[2], CultureInfo.InvariantCulture));
                        break;

                    default:
                        throw new MetaMorpheusException("Unrecognized search mode type: " + split[1]);
                }
            }
            catch (Exception e)
            {
                throw new MetaMorpheusException("Could not parse search mode string: " + e.Message);
            }

            return massDiffAcceptor;
        }

        public static void MatchInternalFragmentIons(SpectralMatch[] fileSpecificPsms, Ms2ScanWithSpecificMass[] arrayOfMs2ScansSortedByMass, CommonParameters combinedParams, int minInternalFragmentLength)
        {
            //for each PSM with an ID
            for (int index = 0; index < fileSpecificPsms.Length; index++)
            {
                SpectralMatch psm = fileSpecificPsms[index];
                if (psm != null && psm.BestMatchingBioPolymersWithSetMods.Count() > 0)
                {
                    //Get the scan
                    Ms2ScanWithSpecificMass scanForThisPsm = arrayOfMs2ScansSortedByMass[index];
                    DissociationType dissociationType = combinedParams.DissociationType == DissociationType.Autodetect ?
                    scanForThisPsm.TheScan.DissociationType.Value : combinedParams.DissociationType;

                    //Get the theoretical peptides
                    List<int> notches = new List<int>();
                    var ambiguousPeptides = psm.BestMatchingBioPolymersWithSetMods.ToList();

                    //get matched ions for each peptide
                    List<List<MatchedFragmentIon>> matchedIonsForAllAmbiguousPeptides = new List<List<MatchedFragmentIon>>();
                    List<Product> internalFragments = new List<Product>();
                    foreach (IBioPolymerWithSetMods peptide in ambiguousPeptides.Select(p => p.SpecificBioPolymer))
                    {
                        internalFragments.Clear();
                        peptide.FragmentInternally(dissociationType, minInternalFragmentLength, internalFragments);
                        //TODO: currently, internal and terminal ions can match to the same observed peaks (much like how b- and y-ions can match to the same peaks). Investigate if we should change that...                        
                        matchedIonsForAllAmbiguousPeptides.Add(MetaMorpheusEngine.MatchFragmentIons(scanForThisPsm, internalFragments, combinedParams));
                    }

                    //Find the max number of matched ions
                    int maxNumMatchedIons = matchedIonsForAllAmbiguousPeptides.Max(x => x.Count);

                    //remove peptides if they have fewer than max-1 matched ions, thus requiring at least two internal ions to disambiguate an ID
                    //if not removed, then add the matched internal ions
                    HashSet<IBioPolymerWithSetMods> PeptidesToMatchingInternalFragments = new HashSet<IBioPolymerWithSetMods>();
                    for (int peptideIndex = 0; peptideIndex < ambiguousPeptides.Count; peptideIndex++)
                    {
                        var thisPeptide = ambiguousPeptides[peptideIndex];
                        //if we should remove the theoretical, remove it
                        if (matchedIonsForAllAmbiguousPeptides[peptideIndex].Count + 1 < maxNumMatchedIons)
                        {
                            psm.RemoveThisAmbiguousPeptide(thisPeptide);
                        }
                        // otherwise add the matched internal ions to the total ions
                        else
                        {
                            IBioPolymerWithSetMods currentPwsm = thisPeptide.SpecificBioPolymer;
                            //check that we haven't already added the matched ions for this peptide
                            if (!PeptidesToMatchingInternalFragments.Contains(currentPwsm))
                            {
                                PeptidesToMatchingInternalFragments.Add(currentPwsm); //record that we've seen this peptide
                                thisPeptide.MatchedIons.AddRange(matchedIonsForAllAmbiguousPeptides[peptideIndex]); //add the matched ions
                            }
                        }
                    }
                }
            }
        }
    }
}
