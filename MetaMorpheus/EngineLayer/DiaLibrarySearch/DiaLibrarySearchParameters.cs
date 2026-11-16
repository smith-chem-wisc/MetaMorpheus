#nullable enable
using System;
using System.Collections.Generic;

namespace EngineLayer.DiaLibrarySearch;

/// <param name="FragmentTolerancePpm">How far an observed peak may sit from a library fragment and still be it.</param>
/// <param name="IrtHalfWindow">
/// How far, in iRT, a precursor's apex may fall from its library iRT. Wide enough here to absorb a run whose true
/// RT-to-iRT curve bows away from the linear map it was searched with.
/// </param>
/// <param name="TopFragmentCount">How many of a library entry's most intense fragments are read and scored.</param>
/// <param name="ApexHalfWidthScans">Scans on either side of an apex over which fragment co-elution is measured.</param>
/// <param name="PrecursorSampleStride">
/// Score only precursors whose library index is a multiple of this; 1 scores them all. Targets and decoys are sampled
/// alike, so a sampled search's q-values stay fair. Calibration's first pass uses it to spend less time.
/// </param>
/// <param name="MaxApexCandidates">
/// Candidate peaks scored per precursor. The classifier's score picks among them, so a precursor whose true peak scores
/// second to interference on the apex score can still be found. 30 by default: over 10, two-seed means gained 0.52% (PXD022589 HF-X) and 2.43% (PXD005573 1 h) at a matched paired
/// entrapment FDP of 1%, and targets over the curve gained alike on an E. coli arm with permuted decoys, with its false hits within counting noise. 20 was level within noise.
/// Held-out candidate rows grow by about 63%.
/// </param>
/// <param name="ClassifierTrainingQValue">
/// Targets at or below this q-value train the classifier as positives: 0.01, as Percolator. A looser cutoff (0.15)
/// once helped on a library of DIA-NN's own IDs, but on a whole-proteome library its positives are largely false targets and the network ranked hundreds of decoys among the best targets (4,896 against 20,966 precursors at 1%). With interference removal, 0.01 is as good or better on both libraries (results/entrapment/README.md).
/// </param>
/// <param name="InterferenceExplainedFragments">
/// Interference removal (<see cref="DiaInterferenceRemoval"/>): a match is dropped when a better, co-eluting match explains
/// at least this many of its six most intense library fragments. 0 turns removal off. 3 by default: against 4, two-seed
/// means gained 0.75% (HF-X) and 1.21% (held-out); 6 lost 7.4% / 7.7%. Calibration keeps 4
/// (<see cref="DiaIrtSelfCalibration.CalibrationInterferenceExplainedFragments"/>).
/// </param>
/// <param name="ClassifierNetworkMembers">
/// Networks in the classifier's ensemble, trained 10 epochs each; calibration trains as many. 24 by default, twice DIA-NN's
/// 12: at the 64-32-16 network, two-seed means at a matched entrapment FDP of 1% gained 0.12% (HF-X) and 1.60% (PXD005573
/// 1 h) over 12. 12 had gained 1.6% (1 h) and 3.9% (HF-X) over 5 on the whole-proteome library, at about 1.5-2x the search
/// time.
/// </param>
/// <param name="ClassifierNetworkEpochs">Training epochs per network (DIA-NN uses 1).</param>
/// <param name="ClassifierNetworkTrainingSample">
/// Which rows train the network when a fold has more than <see cref="MaxNetworkTrainingRows"/>: a random sample, or the
/// highest-ranked targets and decoys (half each), as DIA-NN trains after removing low-confidence identifications, or the
/// highest-ranked rows whatever their label (one score threshold). Half-and-half over random: +9.4% (PXD022589 HF-X) and
/// +9.9% (PXD005573 1 h) at a matched entrapment FDP of 1%. A random sample holds about one real target in twenty on a
/// whole-proteome library. Pooled by default: the half-and-half cutoffs leave a faint band trained on decoys only; pooled at
/// 180,000 rows gained 3.58% (HF-X) and 2.23% (1 h) in two-seed means over half-and-half at 250,000.
/// </param>
/// <param name="SiblingTopCandidateOnly">
/// Sibling support (SiblingCoElution, SiblingApexDeltaMinutes) from each other charge state's top candidate peak only,
/// rather than from any of its candidates near this apex, where a decoy gets several chances for noise to line up. On by
/// default: at a matched entrapment FDP of 1%, HF-X 57,200 -> 59,190 (+3.5%), PXD005573 1 h unchanged (39,280 -> 39,170).
/// </param>
/// <param name="MinimumApexFragments">
/// A cheap gate before the full features, as DIA-NN's: a scan can be a candidate apex only if at least this many of the six
/// most intense library fragments are seen there. 0 (default) turns it off.
/// </param>
/// <param name="FragmentApexCandidate">
/// One more candidate apex by a different rule: the apex of the best single fragment's smoothed trace over the whole window,
/// when no candidate already lies within <see cref="ApexHalfWidthScans"/> of it.
/// </param>
/// <param name="ClassifierNormalizationGroups">
/// When set, each classifier fold normalises its scores on a random sample of this many training groups instead of scoring
/// every training row (most of the classifier's scoring time). Null scores every row. 300,000 by default: IDs unchanged
/// within 0.2% at a matched paired entrapment FDP of 1%, and the classifier about a sixth faster.
/// </param>
/// <param name="MostIntenseFragmentPeak">
/// Read each fragment from the most intense peak within tolerance (DIA-NN 1.8) rather than the nearest one. On by default:
/// at a matched paired entrapment FDP of 1%, +0.2% (HF-X) and +1.4% (PXD005573 1 h).
/// </param>
/// <param name="DiaNnPeakFinding">
/// Candidate apexes as DIA-NN 1.8 picks them (<see cref="DiaLibrarySearchEngine.DiaNnCandidateApexes"/>) rather than by the
/// apex score, still at most <see cref="MaxApexCandidates"/>. On by default with 6 candidates: neutral at 3, +1.3% (HF-X) and
/// +2.8% (PXD005573 1 h) at 6, at a matched paired entrapment FDP of 1%; with <see cref="DiaNnScores"/>, +1.4% and +3.6%.
/// </param>
/// <param name="DiaNnPeakFindingMs1">With <see cref="DiaNnPeakFinding"/>, add the MS1 trace's correlation with the best fragment to the candidate score, as DIA-NN does.</param>
/// <param name="DiaNnSignalShare">DIA-NN's pSig: the six core fragments' shares of their window signal, as SignalShare1-6. On by default (three seeds: +0.9% / +0.6% at matched FDP).</param>
/// <param name="ClassifierLinearIterations">
/// Rankings of the training rows before the network: the first by the best single feature, each later one by a refit
/// linear discriminant, which also picks each precursor's top candidate peak for the network (DIA-NN refits about 8 times).
/// </param>
/// <param name="FragmentMinimumMz">With <see cref="DiaNnFragmentFilter"/>, the lowest fragment m/z scored (DIA-NN: 200).</param>
/// <param name="FragmentMinimumResidues">With <see cref="DiaNnFragmentFilter"/>, the fewest residues a scored fragment spans (DIA-NN: 3).</param>
/// <param name="FragmentMaximumCharge">With <see cref="DiaNnFragmentFilter"/>, the highest fragment charge scored; 0 (default) sets no limit.</param>
/// <param name="FragmentChargeBelowPrecursor">With <see cref="DiaNnFragmentFilter"/>, score only fragments of lower charge than their precursor (a 2+ precursor's 1+ fragments, a 3+ precursor's 1+ and 2+).</param>
/// <param name="Ms1Offset">The run's MS1 offset over retention time (fitted in calibration); every MS1 lookup is shifted by it. Null shifts nothing.</param>
/// <param name="Ms1PeakFeatures">Features: the M0 mass error over the co-elution window (intensity-weighted) and the M+2 isotope trace's co-elution.</param>
/// <param name="ClassifierNetworkPositiveQValue">With the confident training sample, positives are only targets passing this q-value in the linear ranking (with as many top decoys); null keeps the top half-cap of targets.</param>
/// <param name="ClassifierNetworkLayers">The classifier network's hidden layers, input side first; null is <see cref="DefaultClassifierNetworkLayers"/> (64-32-16).</param>
/// <param name="SqrtCoElution">A feature: the core fragments' co-elution on square-root traces, which damp the apex so the peak's flanks weigh more.</param>
/// <param name="Ms1ToleranceSpreadMultiple">When set, the MS1 tolerance follows the run's own mass accuracy: this multiple of the calibrated offset's residual spread (<see cref="EffectiveMs1TolerancePpm"/>).</param>
/// <param name="Ms1PeakEnvelopeTight">A feature: the <see cref="Ms1PeakEnvelope"/> cosine again from MS1 peaks within 0.6x the MS1 tolerance only.</param>
/// <param name="Ms1Intensity">Features: log10(1 + M0 summed over the peak's MS1 scans), and that less log10(1 + the fragments' summed peak area).</param>
/// <param name="Ms1PeakSignalFraction">A feature: the share of the extraction window's MS1 M0 signal inside the peak (the MS1 counterpart of PeakSignalFraction).</param>
/// <param name="Ms1PeakPoints">A feature: the share of the peak's distinct MS1 scans in which an M0 is found.</param>
/// <param name="InterferenceRtToleranceCycles">
/// How far apart, in cycles, two apexes may be and still co-elute for interference removal; null for
/// <see cref="ApexHalfWidthScans"/> + 1.
/// </param>
/// <param name="ExtraFragmentFallback">
/// A precursor with no fragments beyond the top N gets the main fragments' co-elution, matched fraction and weighted
/// co-elution as its extra features, instead of zeros that read as extras not co-eluting.
/// </param>
/// <param name="CoElutionHalfWidthScans">
/// Scans on either side of an apex over which its co-elution features are measured; null for
/// <see cref="ApexHalfWidthScans"/>, which still finds and spaces the candidate apexes.
/// </param>
/// <param name="ClassifierNetworkTargetFraction">
/// The confident network sample's share of targets (the top targets by the linear model), the rest the top decoys; one half
/// by default.
/// </param>
/// <param name="Ms1CorrelationHalfWidthScans">
/// Scans on either side of the apex over which the MS1 M0 and M+1 traces are correlated with the fragment profile, as
/// DIA-NN correlates MS1 over its whole scan window; null for the co-elution window.
/// </param>
/// <param name="LibraryIrtFeature">A feature: the library iRT (DIA-NN's pRT), so the RT error can be weighed by where on the gradient it falls.</param>
/// <param name="Ms1PeakEnvelope">A feature: the MS1 M-1 to M3 envelope summed over the co-elution window's MS1 scans, against the expected pattern with nothing at M-1.</param>
/// <param name="Ms1TightCorrelation">A feature: MS1 M0 co-elution from peaks within half the MS1 tolerance, as DIA-NN scores MS1 at several tolerances.</param>
/// <param name="Ms1TolerancePpm">
/// MS1 tolerance for the precursor's traces, isotope envelope and mass error, around the offset-corrected m/z
/// (<see cref="Ms1Offset"/>). 5 ppm by default, with the calibrated offset applied: against 10 ppm without it, two-seed means
/// at a matched paired entrapment FDP of 1% gain 1.6% (HF-X) and 1.3% (PXD005573 1 h). 10 ppm against 20 gained 1.1% and
/// 1.6%. Without the offset, 5 ppm clipped real precursors as the run's error drifted up to 3 ppm off zero.
/// </param>
/// <param name="DiaNnFragmentFilter">
/// Score only fragments DIA-NN would (at least 3 residues, 200-1800 m/z), when at least 3 of a precursor's qualify. On by
/// default: two-seed means at a matched paired entrapment FDP of 1%, HF-X +0.8%, PXD005573 1 h +3.9%.
/// </param>
/// <param name="MaxToleranceCoElution">
/// A feature as DIA-NN scores co-elution: each core fragment's best correlation with the profile among its traces at the full
/// tolerance, 0.45x and 0.2x (peaks beyond the fraction dropped), averaged. On by default in the main search, with <see cref="Ms1TightCorrelation"/>, <see cref="Ms1PeakFeatures"/>,
/// <see cref="SqrtCoElution"/>, <see cref="Ms1PeakEnvelopeTight"/>, <see cref="Ms1PeakSignalFraction"/> and <see cref="Ms1PeakPoints"/>. Each lost alone under earlier classifiers; together, under the
/// tree blend at 30 candidates, two-seed means gained 0.91% (PXD022589 HF-X) and 1.01% (PXD005573 1 h) at a matched paired entrapment FDP of 1%. Calibration keeps them off.
/// </param>
/// <param name="ExtraFragmentCount">
/// Library fragments after the top <see cref="TopFragmentCount"/> (by library intensity) read as their own features
/// (ExtraCoElution, ExtraMatchedFraction, ExtraWeightedCoElution), as DIA-NN adds its remaining fragments: averaged into the
/// core scores they dilute them, since most are faint. 12 by default (fragments 13-24): +0.9% HF-X, +0.5% PXD005573 1 h
/// (two-seed means, matched paired entrapment FDP 1%). 0 reads none.
/// </param>
/// <param name="DiaNnScores">
/// DIA-NN 1.8 scores we otherwise lack: MinCorr, NFCorr and ShadowCorr (see <see cref="DiaPrecursorMatch.FeatureNames"/>).
/// Off, they are 0. On (default), seven more traces are read per precursor: +1.0% (HF-X) and +1.4% (PXD005573 1 h) at a
/// matched paired entrapment FDP of 1%.
/// </param>
/// <param name="ClassifierNetworkPasses">
/// Network training passes: the second learns from the candidate peaks the first network picked, as DIA-NN trains twice.
/// Two by default: with the pooled training sample, two-seed means gained 2.82% (PXD022589 HF-X) and 2.00% (PXD005573 1 h)
/// at a matched entrapment FDP of 1% over one pass, as the first network ranks real faint targets into the sample.
/// </param>
/// <param name="ClassifierSeed">
/// The classifier's random seed. 0 is the default; other values give the run-to-run noise a result must exceed.
/// </param>
/// <param name="InterferenceSameMzOnly">
/// As DIA-NN, interference removal drops a match only for a better match with the same precursor m/z (or with it as the +1
/// isotope); false lets any co-eluting match in the same or a neighbouring window explain it. On by default: with the
/// 3-fragment rule, two-seed means gained 0.61% (HF-X) and 0.85% (held-out), with calibration pairing the same way.
/// </param>
/// <param name="MaxNetworkTrainingRows">
/// When set, the classifier's network trains on a random subsample of at most this many rows per fold (every row is still
/// scored). Null trains on all. 250,000, about DIA-NN's training size, took the whole-proteome classifier from 10:10 to
/// 7:03. Default 120,000 with the pooled sample: a larger pooled sample admits decoy-like false targets (375,000 lost 5.3%
/// on PXD005573 1 h). Under two network passes, 120,000 gained 0.95% (HF-X) and 1.06% (held-out) over 180,000. Smaller
/// libraries never reach it.
/// </param>
/// <param name="ClassifierTreeBlend">
/// Blend the classifier's order with gradient-boosted trees' before q-values (<see cref="DiaTreeBlend"/>). On by default: two-seed
/// means gained 2.92% (PXD022589 HF-X) and 2.78% (PXD005573 1 h) at a matched paired entrapment FDP of 1%, for about a minute on
/// HF-X with the machine to itself; competing work slowed it to nearly 8.
/// </param>
/// <param name="CandidatePeaksFromSummedCore">
/// With <see cref="DiaNnPeakFinding"/>, also consider scans where the smoothed sum of the core fragments peaks, in the gaps the
/// best-fragment rule leaves. When a peak's two top scans are close, each one's best fragment can peak on the other, and the whole peak is
/// lost however well its fragments co-elute. A replay found this stopped 482 of 963 (HF-X) and 1,361 of 1,711 (PXD005573 1 h) of the DIA-NN
/// peaks that never became a candidate. These compete on score for the <see cref="MaxApexCandidates"/> places, but no candidate moves: a first
/// form that let them moved found precursors off their peaks and lost more than it rescued. On by default: two-seed means gained 0.92%
/// (PXD022589 HF-X) and 1.07% (PXD005573 1 h) at a matched paired entrapment FDP of 1%. Below 0.5% q it loses some targets to decoys that now
/// score near the very top (PXD005573 1 h at 0.25% q: 3.1-8.0% fewer).
/// </param>
/// <param name="InterferenceOwnIonFloor">
/// When above 0, interference removal also drops a winner whose own ions, those not on its rival's list, co-elute with the ions the two share below this library-weighted
/// correlation (<see cref="DiaInterferenceRemoval.OwnIonCoElution"/>). Off by default.
/// </param>
/// <param name="InterferenceOwnIonShareCeiling">
/// With <see cref="InterferenceOwnIonFloor"/> on, only a winner whose own ions carry no more than this share of its library intensity is judged (<see cref="DiaInterferenceRemoval.OwnIonShare"/>):
/// one with almost nothing of its own cannot be told from its rival, and one with much of its own that reads poorly is usually real. 1 judges every winner.
/// </param>
/// <param name="NeighbourEvidence">
/// After the tree blend, measure how much of each near-threshold precursor's fragment signal confident precursors eluting with it already explain, and re-rank on the score and that
/// evidence (<see cref="DiaNeighbourEvidence"/>). Rows at the very top keep their places. On by default: two-seed means gained 1.91% (PXD022589 HF-X) and 1.42% (PXD005573 1 h) at a matched paired
/// entrapment FDP of 1%, and on a species arm that cannot reward removing near-isomers, 2.37% and 2.63%. DIA-NN's own list rises at fixed q. Calibration keeps it off.
/// </param>
/// <param name="FragmentIntensities">
/// When set, the library's fragment intensities are read through this table, in choosing a precursor's strongest fragments as well as in scoring them (<see cref="DiaFragmentIntensityModel"/>).
/// Calibration fits one (<see cref="DiaIrtCalibration.FragmentIntensities"/>), and applying a calibration sets it. Two-seed means with it gained 3.35% (PXD022589 HF-X) and 5.35% (PXD005573 1 h) at a matched
/// paired entrapment FDP of 1%, and on a species arm 2.36% and 6.36% at seed 0, with E. coli hits rising in proportion. Null reads the library as it is.
/// </param>
/// <param name="IrtComposition">
/// When set, a precursor is expected this far from its library iRT, in the window read as well as in the error scored (<see cref="DiaIrtCompositionModel"/>). Calibration fits one
/// (<see cref="DiaIrtCalibration.IrtComposition"/>), and applying a calibration sets it. Two-seed means with it gained 2.22% (PXD022589 HF-X) and 4.45% (PXD005573 1 h) at a matched paired
/// entrapment FDP of 1%. On a species arm it gained 3.46% on HF-X and lost 0.92% on the 1 h run (two seeds), and gained 3.9-5.9% on four further PXD005573 HeLa runs. Null expects it at its
/// library value.
/// </param>
/// <param name="ClassifierModel">The model the rescorer fits in each fold: a linear discriminant, or a small neural-network ensemble (DIA-NN's approach).</param>
public sealed record DiaLibrarySearchParameters(double FragmentTolerancePpm = 20, double IrtHalfWindow = 20,
    int TopFragmentCount = 12, int ApexHalfWidthScans = 3, int PrecursorSampleStride = 1, int MaxApexCandidates = 30, double ClassifierTrainingQValue = 0.01,
    StatisticalModels.RescoreModel ClassifierModel = StatisticalModels.RescoreModel.NeuralNetworkEnsemble,
    int InterferenceExplainedFragments = 3,
    int? MaxNetworkTrainingRows = 120_000,
    bool InterferenceSameMzOnly = true,
    int ClassifierSeed = 0,
    int ClassifierNetworkMembers = 24,
    int ClassifierNetworkEpochs = StatisticalModels.TargetDecoyRescorer.NetworkEpochs,
    int ClassifierNetworkPasses = 2,
    StatisticalModels.NetworkTrainingSample ClassifierNetworkTrainingSample = StatisticalModels.NetworkTrainingSample.ConfidentPooled,
    bool SiblingTopCandidateOnly = true,
    int MinimumApexFragments = 0,
    bool FragmentApexCandidate = false,
    int? ClassifierNormalizationGroups = 300_000,
    bool MostIntenseFragmentPeak = true,
    bool DiaNnPeakFinding = true,
    bool DiaNnScores = true,
    bool DiaNnPeakFindingMs1 = false,
    bool DiaNnSignalShare = true,
    int ClassifierLinearIterations = 3,
    int ExtraFragmentCount = 12,
    bool MaxToleranceCoElution = true,
    bool DiaNnFragmentFilter = true,
    double FragmentMinimumMz = 200,
    int FragmentMinimumResidues = 3,
    int FragmentMaximumCharge = 0,
    bool FragmentChargeBelowPrecursor = false,
    double Ms1TolerancePpm = 5,
    bool Ms1TightCorrelation = true,
    int[]? ClassifierNetworkLayers = null,
    double? ClassifierNetworkPositiveQValue = null,
    bool Ms1PeakFeatures = true,
    Ms1OffsetModel? Ms1Offset = null,
    bool SqrtCoElution = true,
    bool Ms1PeakEnvelope = true,
    double? Ms1ToleranceSpreadMultiple = null,
    bool Ms1PeakEnvelopeTight = true,
    bool Ms1Intensity = true,
    bool Ms1PeakSignalFraction = true,
    bool Ms1PeakPoints = true,
    double? InterferenceRtToleranceCycles = null,
    bool ExtraFragmentFallback = false,
    int? CoElutionHalfWidthScans = null,
    bool LibraryIrtFeature = false,
    int? Ms1CorrelationHalfWidthScans = null,
    double ClassifierNetworkTargetFraction = 0.5,
    bool ClassifierTreeBlend = true,
    bool CandidatePeaksFromSummedCore = true,
    double InterferenceOwnIonFloor = 0,
    double InterferenceOwnIonShareCeiling = 1,
    bool NeighbourEvidence = true,
    DiaFragmentIntensityModel? FragmentIntensities = null,
    DiaIrtCompositionModel? IrtComposition = null)
{
    /// <summary>
    /// The MS1 tolerance a search uses: <see cref="Ms1ToleranceSpreadMultiple"/> times the calibrated offset's residual
    /// spread (<see cref="Ms1OffsetModel.ResidualSpreadPpm"/>) when both are known, else <see cref="Ms1TolerancePpm"/>.
    /// </summary>
    public double EffectiveMs1TolerancePpm =>
        Ms1ToleranceSpreadMultiple is { } multiple && Ms1Offset is { ResidualSpreadPpm: var spread } && double.IsFinite(spread) && spread > 0
            ? multiple * spread
            : Ms1TolerancePpm;

    /// <summary>
    /// The classifier network's hidden layers when <see cref="ClassifierNetworkLayers"/> is null: 64-32-16, wider than the
    /// rescorer's 25-20-15-10-5 (DIA-NN 2020's). With the pooled sample at 120,000 rows under two passes, two-seed means at a
    /// matched entrapment FDP of 1% gained 1.21% (HF-X) and 1.09% (held-out). Under the half-and-half sample the same network
    /// gained 0.9% and lost 0.3%.
    /// </summary>
    public static IReadOnlyList<int> DefaultClassifierNetworkLayers { get; } = [64, 32, 16];

    /// <summary>The hidden layers the classifier network trains: <see cref="ClassifierNetworkLayers"/>, else <see cref="DefaultClassifierNetworkLayers"/>.</summary>
    public IReadOnlyList<int> EffectiveClassifierNetworkLayers => ClassifierNetworkLayers ?? DefaultClassifierNetworkLayers;

    private readonly int _minimumApexFragments = MinimumApexFragments is >= 0 and <= 6 ? MinimumApexFragments
        : throw new ArgumentOutOfRangeException(nameof(MinimumApexFragments), MinimumApexFragments, "The gate counts the six most intense fragments, so 0 to 6.");

    /// <exception cref="ArgumentOutOfRangeException">Outside 0 to 6.</exception>
    public int MinimumApexFragments
    {
        get => _minimumApexFragments;
        init => _minimumApexFragments = value is >= 0 and <= 6 ? value
            : throw new ArgumentOutOfRangeException(nameof(MinimumApexFragments), value, "The gate counts the six most intense fragments, so 0 to 6.");
    }

    private readonly int _extraFragmentCount = ExtraFragmentCount >= 0 ? ExtraFragmentCount
        : throw new ArgumentOutOfRangeException(nameof(ExtraFragmentCount), ExtraFragmentCount, "The extra fragment count cannot be negative.");

    /// <exception cref="ArgumentOutOfRangeException">Negative.</exception>
    public int ExtraFragmentCount
    {
        get => _extraFragmentCount;
        init => _extraFragmentCount = value >= 0 ? value
            : throw new ArgumentOutOfRangeException(nameof(ExtraFragmentCount), value, "The extra fragment count cannot be negative.");
    }

    private readonly int _classifierLinearIterations = ClassifierLinearIterations >= 1 ? ClassifierLinearIterations
        : throw new ArgumentOutOfRangeException(nameof(ClassifierLinearIterations), ClassifierLinearIterations, "The linear model needs at least one iteration.");

    /// <exception cref="ArgumentOutOfRangeException">Fewer than one iteration.</exception>
    public int ClassifierLinearIterations
    {
        get => _classifierLinearIterations;
        init => _classifierLinearIterations = value >= 1 ? value
            : throw new ArgumentOutOfRangeException(nameof(ClassifierLinearIterations), value, "The linear model needs at least one iteration.");
    }

    private readonly int _classifierNetworkPasses = ClassifierNetworkPasses >= 1 ? ClassifierNetworkPasses
        : throw new ArgumentOutOfRangeException(nameof(ClassifierNetworkPasses), ClassifierNetworkPasses, "The network needs at least one training pass.");

    /// <exception cref="ArgumentOutOfRangeException">Fewer than one pass.</exception>
    public int ClassifierNetworkPasses
    {
        get => _classifierNetworkPasses;
        init => _classifierNetworkPasses = value >= 1 ? value
            : throw new ArgumentOutOfRangeException(nameof(ClassifierNetworkPasses), value, "The network needs at least one training pass.");
    }

    private readonly int _precursorSampleStride = Positive(PrecursorSampleStride);

    /// <exception cref="ArgumentOutOfRangeException">The stride is less than 1.</exception>
    public int PrecursorSampleStride
    {
        get => _precursorSampleStride;
        init => _precursorSampleStride = Positive(value);
    }

    private readonly int _precursorSampleOffset;

    /// <summary>
    /// Which of the stride's samples is scored: precursors whose library index leaves this remainder. 0 by default. A sampled
    /// search with another offset draws a disjoint sample, so calibration can be repeated to see how much its sample matters.
    /// A search at stride 1 scores every precursor whatever the offset.
    /// </summary>
    /// <exception cref="ArgumentOutOfRangeException">The offset is negative or not less than the stride.</exception>
    public int PrecursorSampleOffset
    {
        get => _precursorSampleOffset;
        init => _precursorSampleOffset = value >= 0 && value < PrecursorSampleStride ? value
            : throw new ArgumentOutOfRangeException(nameof(PrecursorSampleOffset), value, "The sample offset must be at least 0 and less than the stride.");
    }

    private static int Positive(int stride) => stride >= 1 ? stride
        : throw new ArgumentOutOfRangeException(nameof(PrecursorSampleStride), stride, "The precursor sample stride must be at least 1.");
}