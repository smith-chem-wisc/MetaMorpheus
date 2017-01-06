using System;
using System.Collections.Generic;

namespace MetaMorpheus
{
    public class DatabaseSearcher
    {
        //private IEnumerable<TandemMassSpectra> datas;
        private int minimumAssumedPrecursorChargeState;
        private int maximumAssumedPrecursorChargeState;
        private int maximumNumberOfPeaks;
        private bool assignChargeStates;
        private bool deisotope;
        private bool onTheFlyDecoys;
        private Protease protease;
        private int maximumMissedCleavages;
        private InitiatorMethionineBehavior initiatorMethionineBehavior;
        private int maximumVariableModificationIsoforms;
        //private MassTolerance precursorMassTolerance;
        private IEnumerable<double> acceptedPrecursorMassErrors;
        //private MassTolerance productMassTolerance;
        private double maximumFalseDiscoveryRate;
        private string outputFolder;
        private bool bIons;
        private bool yIons;
        private IEnumerable<Protein> allProteins;
        private IEnumerable<MorpheusModification> variableModifications;
        private IEnumerable<MorpheusModification> fixedModifications;
        private string extraLogStuff;
        private int datasCount;
        private int max_mods_for_peptide;

        public event EventHandler<string> outputHandler;

        public event EventHandler<int> progressHandler;

        public event EventHandler<string> labelStatusHandler;

        public event EventHandler<string> finishedFileHandler;

        public virtual void OnFileFinished(string e)
        {
            finishedFileHandler?.Invoke(this, e);
        }

        public virtual void OnOutput(string e)
        {
            outputHandler?.Invoke(this, e);
        }

        public virtual void OnProgress(int e)
        {
            progressHandler?.Invoke(this, e);
        }

        public virtual void OnStatusUpdate(string e)
        {
            labelStatusHandler?.Invoke(this, e);
        }

        //public DatabaseSearcher(IEnumerable<TandemMassSpectra> datas,
        //    int minimumAssumedPrecursorChargeState, int maximumAssumedPrecursorChargeState,
        //    int maximumNumberOfPeaks,
        //    bool assignChargeStates, bool deisotope, bool onTheFlyDecoys,
        //    Protease protease, int maximumMissedCleavages, InitiatorMethionineBehavior initiatorMethionineBehavior,
        //    int maximumVariableModificationIsoforms, int max_mods_for_peptide,
        //    MassTolerance precursorMassTolerance,
        //    IEnumerable<double> acceptedPrecursorMassErrors,
        //    MassTolerance productMassTolerance,
        //    double maximumFalseDiscoveryRate,
        //    string outputFolder,
        //    bool bIons,
        //    bool yIons,
        //    IEnumerable<Protein> allProteins, IEnumerable<MorpheusModification> variableModifications, IEnumerable<MorpheusModification> fixedModifications,
        //    string extraLogStuff,
        //    int datasCount)
        //{
        //    this.datas = datas;
        //    this.assignChargeStates = assignChargeStates;
        //    this.deisotope = deisotope;
        //    this.onTheFlyDecoys = onTheFlyDecoys;
        //    this.protease = protease;
        //    this.maximumMissedCleavages = maximumMissedCleavages;
        //    this.initiatorMethionineBehavior = initiatorMethionineBehavior;
        //    this.minimumAssumedPrecursorChargeState = minimumAssumedPrecursorChargeState;
        //    this.maximumAssumedPrecursorChargeState = maximumAssumedPrecursorChargeState;
        //    this.maximumNumberOfPeaks = maximumNumberOfPeaks;
        //    this.precursorMassTolerance = precursorMassTolerance;
        //    this.acceptedPrecursorMassErrors = acceptedPrecursorMassErrors;
        //    this.productMassTolerance = productMassTolerance;
        //    this.maximumFalseDiscoveryRate = maximumFalseDiscoveryRate;
        //    this.outputFolder = outputFolder;
        //    this.maximumVariableModificationIsoforms = maximumVariableModificationIsoforms;
        //    this.bIons = bIons;
        //    this.yIons = yIons;
        //    this.allProteins = allProteins;
        //    this.variableModifications = variableModifications;
        //    this.fixedModifications = fixedModifications;
        //    this.extraLogStuff = extraLogStuff;
        //    this.datasCount = datasCount;
        //    this.max_mods_for_peptide = max_mods_for_peptide;
        //}

        public void DoSearch()
        {
            //StreamWriter overall_log = null;
            //StreamWriter summary = null;
            //StreamWriter log = null;

            //DateTime overall_start = DateTime.Now;

            //OnStatusUpdate("Initializing...");

            //OnProgress(0);

            //outputFolder = Path.GetFullPath(outputFolder);

            //int total_proteins = allProteins.Count();

            //OnOutput("total_proteins: " + total_proteins);

            //int num_target_peptides = 0;
            //int num_decoy_peptides = 0;
            //double decoys_over_targets_peptide_ratio = double.NaN;

            //int total_spectra = 0;
            //List<PeptideSpectrumMatch> aggregate_psms = null;

            //if (datasCount > 1)
            //{
            //    overall_log = new StreamWriter(Path.Combine(outputFolder, "log.txt"));
            //    overall_log.AutoFlush = true;

            //    overall_log.WriteLine("Custom MetaMorpheus LOG");
            //    overall_log.WriteLine();

            //    overall_log.WriteLine("PARAMETERS");
            //    string data_filepaths = null;
            //    foreach (var data_filepath in datas)
            //    {
            //        data_filepaths += data_filepath.filename.ToString() + ", ";
            //    }
            //    //data_filepaths = data_filepaths.Substring(0, data_filepaths.Length - 2);
            //    //overall_log.WriteLine("Input Data Files: " + data_filepaths);
            //    overall_log.WriteLine("Unknown Precursor Charge State Range: " + minimumAssumedPrecursorChargeState.ToString("+0;-0;0") + ".." + maximumAssumedPrecursorChargeState.ToString("+0;-0;0"));
            //    overall_log.WriteLine("Maximum Number of MS/MS Peaks: " + (maximumNumberOfPeaks >= 0 ? maximumNumberOfPeaks.ToString() : "disabled"));
            //    overall_log.WriteLine("Assign Charge States: " + assignChargeStates.ToString().ToLower());
            //    overall_log.WriteLine("De-isotope: " + deisotope.ToString().ToLower());
            //    overall_log.WriteLine("Create Target–Decoy Database On The Fly: " + onTheFlyDecoys.ToString().ToLower());
            //    overall_log.WriteLine("Protease: " + protease.ToString());
            //    overall_log.WriteLine("Maximum Missed Cleavages: " + maximumMissedCleavages.ToString());
            //    overall_log.WriteLine("Initiator Methionine Behavior: " + initiatorMethionineBehavior.ToString().ToLower());
            //    overall_log.WriteLine("Maximum Variable Modification Isoforms Per Peptide: " + maximumVariableModificationIsoforms.ToString());
            //    overall_log.WriteLine("Precursor Mass Tolerance: ±" + precursorMassTolerance.Value.ToString(CultureInfo.InvariantCulture) + ' ' + precursorMassTolerance.Units.ToString() + " (" + ')');
            //    overall_log.WriteLine("Accepted Precursor Mass Errors: " + string.Join("; ", acceptedPrecursorMassErrors) + " Da");
            //    overall_log.WriteLine("Product Mass Tolerance: ±" + productMassTolerance.Value.ToString(CultureInfo.InvariantCulture) + ' ' + productMassTolerance.Units.ToString() + " (" + ')');
            //    overall_log.WriteLine("Maximum False Discovery Rate: " + (maximumFalseDiscoveryRate * 100).ToString(CultureInfo.InvariantCulture) + '%');
            //    overall_log.WriteLine("Output Folder: " + outputFolder.ToString());
            //    overall_log.WriteLine(extraLogStuff);
            //    overall_log.WriteLine();

            //    overall_log.WriteLine("RESULTS");

            //    overall_log.WriteLine(total_proteins.ToString("N0") + " target proteins");

            //    aggregate_psms = new List<PeptideSpectrumMatch>();
            //}

            //summary = new StreamWriter(Path.Combine(outputFolder, "summary.tsv"));
            //summary.AutoFlush = true;

            //summary.WriteLine("Dataset\tProteins\tMS/MS Spectra\tPSM MetaMorpheus Score Threshold\tTarget PSMs\tDecoy PSMs\tPSM FDR (%)\tUnique Peptide MetaMorpheus Score Threshold\tUnique Target Peptides\tUnique Decoy Peptides\tUnique Peptide FDR (%)\tProtein Group Summed MetaMorpheus Score Threshold\tTarget Protein Groups\tDecoy Protein Groups\tProtein Group FDR (%)");

            //foreach (TandemMassSpectra spectra in datas)
            //{
            //    GC.Collect();

            //    DateTime start = DateTime.Now;

            //    OnOutput(spectra.filename);

            //    OnProgress(0);

            //    log = new StreamWriter(Path.Combine(outputFolder, Path.GetFileNameWithoutExtension(spectra.filename) + ".log.txt"));
            //    log.AutoFlush = true;

            //    log.WriteLine("Custom MetaMorpheus LOG");
            //    log.WriteLine();

            //    log.WriteLine("PARAMETERS");
            //    log.WriteLine("Input Data File: " + spectra.filename);
            //    log.WriteLine("Unknown Precursor Charge State Range: " + minimumAssumedPrecursorChargeState.ToString("+0;-0;0") + ".." + maximumAssumedPrecursorChargeState.ToString("+0;-0;0"));
            //    log.WriteLine("Maximum Number of MS/MS Peaks: " + (maximumNumberOfPeaks >= 0 ? maximumNumberOfPeaks.ToString() : "disabled"));
            //    log.WriteLine("Assign Charge States: " + assignChargeStates.ToString().ToLower());
            //    log.WriteLine("De-isotope: " + deisotope.ToString().ToLower());
            //    log.WriteLine("Create Target–Decoy Database On The Fly: " + onTheFlyDecoys.ToString().ToLower());
            //    log.WriteLine("Protease: " + protease.ToString());
            //    log.WriteLine("Maximum Missed Cleavages: " + maximumMissedCleavages.ToString());
            //    log.WriteLine("Initiator Methionine Behavior: " + initiatorMethionineBehavior.ToString().ToLower());
            //    log.WriteLine("Maximum Variable Modification Isoforms Per Peptide: " + maximumVariableModificationIsoforms.ToString());
            //    log.WriteLine("Precursor Mass Tolerance: ±" + precursorMassTolerance.Value.ToString(CultureInfo.InvariantCulture) + ' ' + precursorMassTolerance.Units.ToString() + " (" + ')');
            //    log.WriteLine("Accepted Precursor Mass Errors: " + string.Join("; ", acceptedPrecursorMassErrors) + " Da");
            //    log.WriteLine("Product Mass Tolerance: ±" + productMassTolerance.Value.ToString(CultureInfo.InvariantCulture) + ' ' + productMassTolerance.Units.ToString() + " (" + ')');
            //    log.WriteLine("Maximum False Discovery Rate: " + (maximumFalseDiscoveryRate * 100).ToString(CultureInfo.InvariantCulture) + '%');
            //    log.WriteLine("Output Folder: " + outputFolder.ToString());
            //    log.WriteLine(extraLogStuff);
            //    log.WriteLine();

            //    log.WriteLine("RESULTS");

            //    log.WriteLine(total_proteins.ToString("N0") + " target proteins");

            //    OnStatusUpdate("Extracting and preprocessing MS/MS spectra...");

            //    OnProgress(0);

            //    spectra.UpdateProgress += progressHandler;

            //    spectra.Load(spectra.filename, minimumAssumedPrecursorChargeState, maximumAssumedPrecursorChargeState,
            //         maximumNumberOfPeaks,
            //         assignChargeStates, deisotope, productMassTolerance);

            //    spectra.UpdateProgress -= progressHandler;

            //    if (datasCount > 1)
            //    {
            //        total_spectra += spectra.Count;
            //    }

            //    OnStatusUpdate("Searching MS/MS spectra...");
            //    OnProgress(0);

            //    PeptideSpectrumMatch[] psms = null;
            //    if (spectra.Count > 0)
            //    {
            //        int max_spectrum_number = 0;
            //        foreach (TandemMassSpectrum spectrum in spectra)
            //        {
            //            if (spectrum.SpectrumNumber > max_spectrum_number)
            //            {
            //                max_spectrum_number = spectrum.SpectrumNumber;
            //            }
            //        }

            //        psms = new PeptideSpectrumMatch[max_spectrum_number];

            //        spectra.Sort(TandemMassSpectrum.AscendingPrecursorMassComparison);
            //    }

            //    Dictionary<string, bool> peptides_observed = null;
            //    peptides_observed = new Dictionary<string, bool>();

            //    HashSet<string> level3_observed_target = new HashSet<string>();
            //    HashSet<string> level3_observed_decoy = new HashSet<string>();

            //    object progress_lock = new object();
            //    int proteins = 0;
            //    //int old_progress = 0;
            //    ParallelOptions parallel_options = new ParallelOptions();
            //    //parallel_options.MaxDegreeOfParallelism = 1;

            //    var productTypes = new List<ProductType>();
            //    if (bIons)
            //        productTypes.Add(ProductType.b);
            //    if (yIons)
            //        productTypes.Add(ProductType.y);

            //    //Parallel.ForEach(allProteins, parallel_options, protein =>
            //    //{
            //    //    foreach (PeptideWithPossibleModifications bigFragment in protein.GetBigFragments())
            //    //    {
            //    //        foreach (PeptideWithPossibleModifications peptide in bigFragment.Digest(protease, maximumMissedCleavages, initiatorMethionineBehavior))
            //    //        {
            //    //            if (peptide.Length == 1)
            //    //                continue;

            //    //            if (peptide.OneBasedPossibleLocalizedModifications.Count == 0)
            //    //            {
            //    //                lock (level3_observed_target) lock (level3_observed_decoy)
            //    //                    {
            //    //                        var hc = peptide.BaseSequence;
            //    //                        var decoyObserved = level3_observed_decoy.Contains(hc);
            //    //                        // Seen as decoy, keep going (even if target)
            //    //                        if (decoyObserved)
            //    //                            continue;
            //    //                        // Seen as target and is target, keep going
            //    //                        if (level3_observed_target.Contains(hc) && !protein.isDecoy)
            //    //                            continue;
            //    //                        // Decoy and not seen as decoy. Go ahead!
            //    //                        if (!decoyObserved && protein.isDecoy)
            //    //                            level3_observed_decoy.Add(hc);
            //    //                        // Target and not seen. Go ahead!
            //    //                        if (!protein.isDecoy)
            //    //                            level3_observed_target.Add(hc);
            //    //                    }
            //    //            }

            //    //            //foreach (PeptideWithSetModifications modified_peptide in peptide.GetPeptideWithSetModifications(variableModifications, maximumVariableModificationIsoforms, max_mods_for_peptide))
            //    //            //{
            //    //            //    foreach (Tuple<TandemMassSpectrum, int> theTuple in spectra.GetTandemMassSpectraInMassRanges(modified_peptide.MonoisotopicMass, acceptedPrecursorMassErrors, precursorMassTolerance))
            //    //            //    {
            //    //            //        PeptideSpectrumMatch psm = new PeptideSpectrumMatch(theTuple.Item1, modified_peptide, productMassTolerance, psms[theTuple.Item1.SpectrumNumber - 1] == null ? 1 : psms[theTuple.Item1.SpectrumNumber - 1].MatchingProducts, productTypes, theTuple.Item2);
            //    //            //        lock (psms)
            //    //            //        {
            //    //            //            PeptideSpectrumMatch current_best_psm = psms[theTuple.Item1.SpectrumNumber - 1];
            //    //            //            if (current_best_psm == null || PeptideSpectrumMatch.DescendingMetaMorpheusScoreComparison(psm, current_best_psm) < 0)
            //    //            //                psms[theTuple.Item1.SpectrumNumber - 1] = psm;
            //    //            //        }
            //    //            //    }
            //    //            //}
            //    //        }
            //    //    }

            //    //    lock (progress_lock)
            //    //    {
            //    //        proteins++;
            //    //        int new_progress = (int)((double)proteins / (total_proteins) * 100);
            //    //        if (new_progress > old_progress)
            //    //        {
            //    //            OnProgress(new_progress);
            //    //            old_progress = new_progress;
            //    //        }
            //    //    }
            //    //});

            //    OnOutput("proteins: " + proteins);

            //    OnStatusUpdate("Performing post-search analyses...");
            //    OnProgress(0);

            //    log.WriteLine((num_target_peptides + num_decoy_peptides).ToString("N0") + " total (" + num_target_peptides.ToString("N0") + " target + " + num_decoy_peptides.ToString("N0") + " decoy) non-unique peptides");
            //    decoys_over_targets_peptide_ratio = (double)num_decoy_peptides / num_target_peptides;

            //    log.WriteLine(spectra.Count.ToString("N0") + " MS/MS spectra");

            //    List<PeptideSpectrumMatch> psms_no_nulls;
            //    if (psms != null)
            //    {
            //        psms_no_nulls = new List<PeptideSpectrumMatch>(psms.Length);
            //        foreach (PeptideSpectrumMatch psm in psms)
            //        {
            //            if (psm != null)
            //            {
            //                psms_no_nulls.Add(psm);
            //            }
            //        }

            //        if (datasCount > 1)
            //        {
            //            aggregate_psms.AddRange(psms_no_nulls);
            //        }
            //    }
            //    else
            //    {
            //        psms_no_nulls = new List<PeptideSpectrumMatch>(0);
            //    }

            //    List<PeptideSpectrumMatch> sorted_psms = new List<PeptideSpectrumMatch>(psms_no_nulls);
            //    sorted_psms.Sort(PeptideSpectrumMatch.DescendingMetaMorpheusScoreComparison);

            //    foreach (var psm in sorted_psms)
            //    {
            //        TandemMassSpectrum s = spectra[psm.SpectrumIndexHere];
            //        psm.ComputeIonMatchesAndCounts(s);
            //        psm.GetLocalizedScores(s);
            //    }

            //    IEnumerable<IdentificationWithFalseDiscoveryRate<PeptideSpectrumMatch>> psms_with_fdr = FalseDiscoveryRate.DoFalseDiscoveryRateAnalysis(sorted_psms);
            //    Exporters.WriteToTabDelimitedTextFile(psms_with_fdr, Path.Combine(outputFolder, Path.GetFileNameWithoutExtension(spectra.filename) + ".psmtsv"));
            //    OnFileFinished(Path.Combine(outputFolder, Path.GetFileNameWithoutExtension(spectra.filename) + ".psmtsv"));

            //    //Exporters.MyAnalysis(psms_with_fdr, Path.Combine(outputFolder, Path.GetFileNameWithoutExtension(spectra.filename) + ".mytsv"));
            //    //OnFileFinished(Path.Combine(outputFolder, Path.GetFileNameWithoutExtension(spectra.filename) + ".mytsv"));

            //    double psm_score_threshold = double.NegativeInfinity;
            //    int target_psms = sorted_psms.Count;
            //    int decoy_psms = 0;
            //    double psm_fdr = double.NaN;
            //    if (decoys_over_targets_peptide_ratio == 0.0)
            //    {
            //        log.WriteLine(sorted_psms.Count.ToString("N0") + " PSMs (unknown FDR)");
            //    }
            //    else
            //    {
            //        FalseDiscoveryRate.DetermineMaximumIdentifications(psms_with_fdr, false, maximumFalseDiscoveryRate, out psm_score_threshold, out target_psms, out decoy_psms, out psm_fdr);
            //        log.WriteLine(target_psms.ToString("N0") + " target (" + decoy_psms.ToString("N0") + " decoy) PSMs at " + psm_fdr.ToString("0.000%") + " PSM FDR (" + psm_score_threshold.ToString("0.000") + " MetaMorpheus score threshold)");
            //    }

            //    Dictionary<string, PeptideSpectrumMatch> peptides = new Dictionary<string, PeptideSpectrumMatch>();

            //    List<PeptideSpectrumMatch> sorted_peptides = new List<PeptideSpectrumMatch>(peptides.Values);

            //    double peptide_score_threshold = double.NegativeInfinity;
            //    int target_peptides = sorted_peptides.Count;
            //    int decoy_peptides = 0;
            //    double peptide_fdr = double.NaN;
            //    if (decoys_over_targets_peptide_ratio == 0.0)
            //    {
            //        log.WriteLine(sorted_peptides.Count.ToString("N0") + " unique peptides (unknown FDR)");
            //    }
            //    else
            //    {
            //        log.WriteLine(target_peptides.ToString("N0") + " unique target (" + decoy_peptides.ToString("N0") + " decoy) peptides at " + peptide_fdr.ToString("0.000%") + " unique peptide FDR (" + peptide_score_threshold.ToString("0.000") + " MetaMorpheus score threshold)");
            //    }

            //    Exporters.WriteMZIdentMLFile(Path.Combine(outputFolder, Path.GetFileNameWithoutExtension(spectra.filename) + ".mzid"),
            //        new TandemMassSpectra[] { spectra },
            //        minimumAssumedPrecursorChargeState, maximumAssumedPrecursorChargeState,
            //        maximumNumberOfPeaks,
            //        assignChargeStates, deisotope, onTheFlyDecoys, total_proteins,
            //        protease, maximumMissedCleavages, initiatorMethionineBehavior,
            //        maximumVariableModificationIsoforms,
            //        precursorMassTolerance,
            //        acceptedPrecursorMassErrors,
            //        productMassTolerance,
            //        maximumFalseDiscoveryRate,
            //        outputFolder,
            //        psms_with_fdr, allProteins);

            //    OnFileFinished(Path.Combine(outputFolder, Path.GetFileNameWithoutExtension(spectra.filename) + ".mzid"));

            //    DateTime stop = DateTime.Now;
            //    log.WriteLine((stop - start).TotalMinutes.ToString("0.00") + " minutes to analyze");

            //    log.Close();

            //    summary.Write(spectra.filename + '\t');
            //    summary.Write(proteins.ToString() + '\t');
            //    summary.Write(spectra.Count.ToString() + '\t');
            //    summary.Write(psm_score_threshold.ToString("0.000") + '\t');
            //    summary.Write(target_psms.ToString() + '\t');
            //    summary.Write(decoy_psms.ToString() + '\t');
            //    summary.Write(psm_fdr.ToString("0.000%") + '\t');
            //    summary.Write(peptide_score_threshold.ToString("0.000") + '\t');
            //    summary.Write(target_peptides.ToString() + '\t');
            //    summary.Write(decoy_peptides.ToString() + '\t');
            //    summary.Write(peptide_fdr.ToString("0.000%") + '\t');
            //    summary.WriteLine();

            //    OnFileFinished(spectra.filename);
            //}

            //if (datasCount > 1 && overall_log != null)
            //{
            //    OnStatusUpdate("Performing aggregate post-search analyses...");
            //    OnProgress(0);

            //    overall_log.WriteLine((num_target_peptides + num_decoy_peptides).ToString("N0") + " total (" + num_target_peptides.ToString("N0") + " target + " + num_decoy_peptides.ToString("N0") + " decoy) non-unique peptides");

            //    HashSet<string> prefixes = new HashSet<string>();
            //    prefixes.Add("aggregate");

            //    overall_log.WriteLine(total_spectra.ToString("N0") + " MS/MS spectra");

            //    aggregate_psms.Sort(PeptideSpectrumMatch.DescendingMetaMorpheusScoreComparison);

            //    IEnumerable<IdentificationWithFalseDiscoveryRate<PeptideSpectrumMatch>> aggregate_psms_with_fdr = FalseDiscoveryRate.DoFalseDiscoveryRateAnalysis(aggregate_psms);
            //    Exporters.WriteToTabDelimitedTextFile(aggregate_psms_with_fdr, Path.Combine(outputFolder, "aggregate.psmtsv"));
            //    OnFileFinished(Path.Combine(outputFolder, "aggregate.psmtsv"));

            //    //Exporters.MyAnalysis(aggregate_psms_with_fdr, Path.Combine(outputFolder, "aggregate.mytsv"));
            //    //OnFileFinished(Path.Combine(outputFolder, "aggregate.mytsv"));

            //    double aggregate_psm_score_threshold;
            //    int aggregate_target_psms;
            //    int aggregate_decoy_psms;
            //    double aggregate_psm_fdr;
            //    FalseDiscoveryRate.DetermineMaximumIdentifications(aggregate_psms_with_fdr, false, maximumFalseDiscoveryRate, out aggregate_psm_score_threshold, out aggregate_target_psms, out aggregate_decoy_psms, out aggregate_psm_fdr);
            //    overall_log.WriteLine(aggregate_target_psms.ToString("N0") + " target (" + aggregate_decoy_psms.ToString("N0") + " decoy) aggregate PSMs at " + aggregate_psm_fdr.ToString("0.000%") + " PSM FDR (" + aggregate_psm_score_threshold.ToString("0.000") + " MetaMorpheus score threshold)");

            //    Dictionary<string, PeptideSpectrumMatch> aggregate_peptides = new Dictionary<string, PeptideSpectrumMatch>();

            //    Exporters.WriteMZIdentMLFile(Path.Combine(outputFolder, "aggregate.mzid"),
            //                        datas,
            //                        minimumAssumedPrecursorChargeState, maximumAssumedPrecursorChargeState,
            //                        maximumNumberOfPeaks,
            //                        assignChargeStates, deisotope, onTheFlyDecoys, total_proteins,
            //                        protease, maximumMissedCleavages, initiatorMethionineBehavior,
            //                        maximumVariableModificationIsoforms,
            //                        precursorMassTolerance,
            //                        acceptedPrecursorMassErrors,
            //                        productMassTolerance,
            //                        maximumFalseDiscoveryRate,
            //                        outputFolder,
            //                        aggregate_psms_with_fdr, allProteins);

            //    OnFileFinished(Path.Combine(outputFolder, "aggregate.mzid"));

            //    DateTime overall_stop = DateTime.Now;
            //    overall_log.WriteLine((overall_stop - overall_start).TotalMinutes.ToString("0.00") + " minutes to analyze");

            //    overall_log.Close();

            //    summary.Write("AGGREGATE" + '\t');
            //    summary.Write(total_proteins.ToString() + '\t');
            //    summary.Write(total_spectra.ToString() + '\t');
            //    summary.Write(aggregate_psm_score_threshold.ToString("0.000") + '\t');
            //    summary.Write(aggregate_target_psms.ToString() + '\t');
            //    summary.Write(aggregate_decoy_psms.ToString() + '\t');
            //    summary.Write(aggregate_psm_fdr.ToString("0.000%") + '\t');
            //    summary.WriteLine();
            //}

            //summary.Close();
        }
    }
}