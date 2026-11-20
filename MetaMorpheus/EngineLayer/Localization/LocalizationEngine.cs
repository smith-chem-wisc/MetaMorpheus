using MassSpectrometry;
using Omics;
using Omics.Fragmentation;
using Omics.Modifications;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Threading;
using System.Threading.Tasks;

namespace EngineLayer.Localization
{
    /// <summary>
    /// Collects, for each match, the spectral evidence for each place its modifications could be, from which
    /// localization probabilities are computed once every spectra file has been read.
    ///
    /// Two kinds of placement are tried:
    /// - an unidentified mass difference (observed minus theoretical precursor mass, from an open, notch or custom
    ///   search) at each residue, scored into <see cref="SpectralMatch.LocalizedScores"/>;
    /// - every arrangement of a match's identified modifications over the residues each could occupy, for matches whose
    ///   peptide explains the precursor mass.
    /// A match whose peptide explains the precursor mass gets the second when it is enabled, as its mass difference is
    /// nothing to place. The engine should not be run for exact mass or missed monoisotopic mass-difference localization.
    /// </summary>
    public class LocalizationEngine : MetaMorpheusEngine
    {
        private readonly IEnumerable<SpectralMatch> AllResultingIdentifications;
        private readonly MsDataFile MyMsDataFile;
        private readonly bool LocalizeMassDifferences;
        private readonly IReadOnlyCollection<Modification> PlaceableModifications;
        private readonly int MaxPlacements;
        private readonly LocalizationSiteRule SiteRule;

        /// <summary>
        /// Localizes the mass difference of each match at each residue.
        /// </summary>
        public LocalizationEngine(IEnumerable<SpectralMatch> allResultingIdentifications, MsDataFile myMsDataFile, CommonParameters commonParameters, List<(string fileName, CommonParameters fileSpecificParameters)> fileSpecificParameters, List<string> nestedIds)
            : this(allResultingIdentifications, myMsDataFile, commonParameters, fileSpecificParameters, nestedIds, localizeMassDifferences: true, placeableModifications: null)
        {
        }

        /// <param name="localizeMassDifferences">Whether to place each match's mass difference at each residue, filling
        /// <see cref="SpectralMatch.LocalizedScores"/>.</param>
        /// <param name="placeableModifications">When not null, collect evidence for every arrangement of each match's
        /// identified modifications. These are the modifications one may be swapped for to reach a residue its own
        /// motif does not allow, typically the search's variable modifications. Under
        /// <see cref="LocalizationSiteRule.AnyResidue"/> arrangements are collected whether or not it is given.</param>
        /// <param name="maxPlacements">The most arrangements enumerated for one match; a match with more gets no
        /// placement evidence.</param>
        /// <param name="siteRule">Which residues an identified modification may be moved to: those its motif allows
        /// (the default), or any residue, for experiments where a modification's chemistry is treated as unknown.</param>
        public LocalizationEngine(IEnumerable<SpectralMatch> allResultingIdentifications, MsDataFile myMsDataFile, CommonParameters commonParameters,
            List<(string fileName, CommonParameters fileSpecificParameters)> fileSpecificParameters, List<string> nestedIds,
            bool localizeMassDifferences, IReadOnlyCollection<Modification> placeableModifications, int maxPlacements = PositionalIsoforms.DefaultMaxPlacements,
            LocalizationSiteRule siteRule = LocalizationSiteRule.Motif)
            : base(commonParameters, fileSpecificParameters, nestedIds)
        {
            AllResultingIdentifications = allResultingIdentifications;
            MyMsDataFile = myMsDataFile;
            LocalizeMassDifferences = localizeMassDifferences;
            PlaceableModifications = placeableModifications;
            MaxPlacements = maxPlacements;
            SiteRule = siteRule;
        }

        /// <summary>
        /// Whether to keep, on each match, the <see cref="LocalizationEvidence"/> of its mass-difference placements, for
        /// <see cref="LocalizationProbabilityEngine"/>. On by default; a search that computes no probabilities turns it
        /// off to save the memory. Placing identified modifications collects evidence and nothing else, so it is kept
        /// whatever this says.
        /// </summary>
        public bool KeepMassDifferenceEvidence { get; init; } = true;

        /// <summary>
        /// Whether ion-trap scans are matched at <see cref="CommonParameters.ProductMassTolerance_LowRes"/> rather than
        /// <see cref="CommonParameters.ProductMassTolerance"/>, which is calibrated on Orbitrap fragments: a method that
        /// alternates ion-trap and Orbitrap MS2 in one file otherwise matches the ion-trap scans at a ppm tolerance far
        /// tighter than their real error, missing many of their true fragments.
        /// </summary>
        public bool UseLowResToleranceForIonTraps { get; init; }

        /// <summary>
        /// Whether a fragment and its neutral-loss forms count as separate pieces of localization evidence. Off by
        /// default: they report one backbone cleavage, so counting each would reward a placement twice for one cleavage
        /// (and charge it twice when the cleavage is missing). Only the evidence is affected; scores are not.
        /// </summary>
        public bool CountNeutralLossesSeparately { get; init; }

        /// <summary>
        /// Whether ion-trap scans also count a fragment and its neutral-loss forms once per cleavage. Off by default:
        /// at the wide ion-trap tolerance a chance peak often matches one form of a cleavage but seldom both, so seeing
        /// both is corroboration worth keeping (on PXD070496 collapsing them made ion-trap probabilities overconfident,
        /// while it helped Orbitrap ones).
        /// </summary>
        public bool CollapseIonTrapNeutralLosses { get; init; }

        private bool PerCleavage(MsDataScan scan) => !CountNeutralLossesSeparately && (CollapseIonTrapNeutralLosses || !IsIonTrap(scan));

        private static bool IsIonTrap(MsDataScan scan) => scan.MzAnalyzer is MZAnalyzerType.IonTrap2D or MZAnalyzerType.IonTrap3D;

        protected override MetaMorpheusEngineResults RunSpecific()
        {
            // any-residue placement swaps nothing, so it needs no list of modifications to swap for
            bool placeModifications = PlaceableModifications != null || SiteRule == LocalizationSiteRule.AnyResidue;

            // Mass differences are localized only for matches with one full sequence. Modifications are also placed for a
            // match ambiguous only in where its modifications sit, as every candidate is one of its arrangements.
            SpectralMatch[] psms = AllResultingIdentifications
                .Where(b => b.FullSequence != null || (placeModifications && b.GetSiteFreeIdentity() != null))
                .ToArray();

            if (psms.Length == 0)
                return new LocalizationEngineResults(this);

            var fixedModifications = new HashSet<(string, string)>(CommonParameters.ListOfModsFixed);
            bool IsFixed(Modification mod) => fixedModifications.Contains((mod.ModificationType, mod.IdWithMotif));

            double psmsSearched = 0;
            int oldPercentProgress = 0;
            int matchesWithTooManyPlacements = 0;

            ReportProgress(new ProgressEventArgs(oldPercentProgress, "Localizing mass-differences... ", NestedIds));

            Parallel.ForEach(Partitioner.Create(0, psms.Length),
                new ParallelOptions { MaxDegreeOfParallelism = CommonParameters.MaxThreadsToUsePerFile },
                (range, loopState) =>
                {
                    List<Product> products = new List<Product>();

                    for (int i = range.Item1; i < range.Item2; i++)
                    {
                        SpectralMatch psm = psms[i];

                        // Stop loop if canceled
                        if (GlobalVariables.StopLoops) { break; }

                        MsDataScan scan = MyMsDataFile.GetOneBasedScan(psm.ScanNumber);
                        Ms2ScanWithSpecificMass scanWithSpecificMass = new Ms2ScanWithSpecificMass(scan, psm.ScanPrecursorMonoisotopicPeakMz, psm.ScanPrecursorCharge, psm.FullFilePath, CommonParameters);
                        IBioPolymerWithSetMods peptide = psm.BestMatchingBioPolymersWithSetMods.First().SpecificBioPolymer;
                        // A most-abundant search admits matches whose deconvoluted monoisotopic peak is off by
                        // whole isotopologues; that offset must come out before the difference is localized,
                        // or this engine localizes a ~1-2 Da modification that does not exist.
                        double observedMass = psm.GetObservedMonoisotopicMass(peptide.MonoisotopicMass, CommonParameters);
                        double massDifference = observedMass - peptide.MonoisotopicMass;

                        // the evidence for each placement is kept, relative to the base peak, for localization probabilities
                        double basePeakIntensity = scan.MassSpectrum.YofPeakWithHighestY ?? 0;
                        bool precursorExplained = CommonParameters.PrecursorMassTolerance.Within(observedMass, peptide.MonoisotopicMass);
                        bool placeThisMatchesModifications = placeModifications && precursorExplained;

                        var observable = new ObservableRange(scan, scanWithSpecificMass, UseLowResToleranceForIonTraps && IsIonTrap(scan), PerCleavage(scan));

                        if (LocalizeMassDifferences && psm.FullSequence != null)
                        {
                            bool keepEvidence = KeepMassDifferenceEvidence && basePeakIntensity > 0 && !placeThisMatchesModifications;
                            LocalizeMassDifference(psm, peptide, massDifference, scan, scanWithSpecificMass, products,
                                keepEvidence ? new LocalizationEvidenceBuilder(basePeakIntensity, ScanTypeOf(scan), PerCleavage(scan)) : null, precursorExplained, observable);
                        }

                        if (placeThisMatchesModifications && basePeakIntensity > 0)
                        {
                            if (!PlaceModifications(psm, peptide, massDifference, scan, scanWithSpecificMass, products, basePeakIntensity, IsFixed, observable))
                            {
                                Interlocked.Increment(ref matchesWithTooManyPlacements);
                            }
                        }

                        // report search progress (PSMs localized so far out of total)
                        psmsSearched++;
                        var percentProgress = (int)((psmsSearched / psms.Length) * 100);

                        if (percentProgress > oldPercentProgress)
                        {
                            oldPercentProgress = percentProgress;
                            ReportProgress(new ProgressEventArgs(percentProgress, "Localizing mass-differences... ", NestedIds));
                        }
                    }
                });

            return new LocalizationEngineResults(this, matchesWithTooManyPlacements);
        }

        /// <summary>
        /// Scores the mass difference at each residue into <see cref="SpectralMatch.LocalizedScores"/>, and when
        /// <paramref name="evidence"/> is given, keeps what each placement matched. The first and last residues also
        /// stand for the peptide termini, which fragment identically. When the peptide already explains the precursor
        /// there is no mass difference to place: the evidence is then the unmodified peptide alone.
        /// </summary>
        private void LocalizeMassDifference(SpectralMatch psm, IBioPolymerWithSetMods peptide, double massDifference, MsDataScan scan,
            Ms2ScanWithSpecificMass scanWithSpecificMass, List<Product> products, LocalizationEvidenceBuilder evidence, bool precursorExplained,
            ObservableRange range)
        {
            if (evidence != null && precursorExplained)
            {
                peptide.Fragment(CommonParameters.DissociationType, CommonParameters.DigestionParams.FragmentationTerminus, products, CommonParameters.FragmentationParameters);
                var ions = MatchFragmentIons(scanWithSpecificMass, products, CommonParameters, isLowRes: range.IsLowRes);
                evidence.AddCandidate(new LocalizationPlacement(System.Array.Empty<LocalizationSite>()), CalculatePeptideScore(scan, ions), ions, range.CountObservable(products));
                range.RecordUnmatched(evidence, scanWithSpecificMass, CommonParameters);
                psm.LocalizationEvidence = evidence.Build(massDifference, peptide.Length, IsProteinNTerminal(peptide), IsProteinCTerminal(peptide));
                evidence = null;
            }

            var localizedScores = new List<double>();
            for (int r = 0; r < peptide.Length; r++)
            {
                // create new PeptideWithSetMods with unidentified mass difference at the given residue
                var peptideWithLocalizedMassDiff = peptide.Localize(r, massDifference);

                // this is the list of theoretical products for this peptide with mass-difference on this residue
                peptideWithLocalizedMassDiff.Fragment(CommonParameters.DissociationType, CommonParameters.DigestionParams.FragmentationTerminus, products, CommonParameters.FragmentationParameters);

                var matchedIons = MatchFragmentIons(scanWithSpecificMass, products, CommonParameters, isLowRes: range.IsLowRes);

                // score when the mass-diff is on this residue
                double localizedScore = CalculatePeptideScore(scan, matchedIons);

                localizedScores.Add(localizedScore);
                evidence?.AddCandidate(new LocalizationPlacement(new[] { new LocalizationSite(r, peptide.BaseSequence[r]) }), localizedScore, matchedIons,
                    range.CountObservable(products));
            }

            psm.LocalizedScores = localizedScores;
            if (evidence != null)
            {
                range.RecordUnmatched(evidence, scanWithSpecificMass, CommonParameters);
                psm.LocalizationEvidence = evidence.Build(massDifference, peptide.Length, IsProteinNTerminal(peptide), IsProteinCTerminal(peptide));
            }
        }

        /// <summary>
        /// Keeps what each arrangement of the match's identified modifications matched. A match with nothing to move has
        /// one arrangement, kept so that its matched peaks are there to fit a model of matched-peak intensities from.
        /// </summary>
        /// <returns>False when the match has more arrangements than are enumerated, leaving it without evidence.</returns>
        private bool PlaceModifications(SpectralMatch psm, IBioPolymerWithSetMods peptide, double massDifference, MsDataScan scan,
            Ms2ScanWithSpecificMass scanWithSpecificMass, List<Product> products, double basePeakIntensity,
            System.Func<Modification, bool> isFixed, ObservableRange range)
        {
            if (!PositionalIsoforms.TryEnumerate(peptide, PlaceableModifications, isFixed, MaxPlacements, out var isoforms, SiteRule))
            {
                return false;
            }

            var evidence = new LocalizationEvidenceBuilder(basePeakIntensity, ScanTypeOf(scan), PerCleavage(scan));
            // the observable fragments every arrangement predicts alike: whichever is right, these are predicted
            Dictionary<(ProductType, long), int> shared = null;
            foreach (PositionalIsoform isoform in isoforms)
            {
                IBioPolymerWithSetMods arranged = PositionalIsoforms.WithModifications(peptide, isoform.AllModsOneIsNterminus);
                if (arranged == null)
                {
                    return true; // a kind of biopolymer that cannot be rebuilt; nothing to place
                }

                arranged.Fragment(CommonParameters.DissociationType, CommonParameters.DigestionParams.FragmentationTerminus, products, CommonParameters.FragmentationParameters);
                var matchedIons = MatchFragmentIons(scanWithSpecificMass, products, CommonParameters, isLowRes: range.IsLowRes);
                evidence.AddCandidate(isoform.Placement, CalculatePeptideScore(scan, matchedIons), matchedIons, range.CountObservable(products));

                var masses = range.ObservableMassCounts(products);
                if (shared == null)
                {
                    shared = masses;
                }
                else
                {
                    foreach (var key in shared.Keys.ToList())
                    {
                        shared[key] = System.Math.Min(shared[key], masses.GetValueOrDefault(key));
                    }
                }
            }

            evidence.SetSharedPredictedFragments(shared?.Values.Sum() ?? 0);
            range.RecordUnmatched(evidence, scanWithSpecificMass, CommonParameters);
            psm.LocalizationEvidence = evidence.Build(massDifference, peptide.Length, IsProteinNTerminal(peptide), IsProteinCTerminal(peptide));
            return true;
        }

        /// <summary>
        /// A fragment's intact neutral mass (to 0.1 mDa): a fragment and its neutral-loss forms, which report one backbone
        /// cleavage, share it.
        /// </summary>
        internal static long CleavageKey(Product product) => (long)System.Math.Round((product.NeutralMass + product.NeutralLoss) * 1e4);

        /// <summary>The scan's analyzer and activation, e.g. "IonTrap2D CID": scans alike in both share intensity statistics.</summary>
        internal static string ScanTypeOf(MsDataScan scan) => $"{scan.MzAnalyzer} {scan.DissociationType?.ToString() ?? "Unknown"}";

        /// <summary>
        /// The fragments a scan could show: those with an m/z inside the scan window at some charge from 1 to the
        /// precursor's (an m/z falls as the charge rises, so the m/z at charge 1 must reach the window's low end and the
        /// m/z at the precursor's charge must not pass its high end), and no heavier than the precursor.
        /// </summary>
        internal readonly struct ObservableRange
        {
            private readonly double _lowestMz;
            private readonly double _highestMz;
            private readonly int _charge;
            private readonly int _unitCharge;
            private readonly double _precursorMass;
            private readonly bool _perCleavage;

            /// <param name="isLowRes">Whether the scan's fragments are matched at the low-resolution tolerance
            /// (<see cref="CommonParameters.ProductMassTolerance_LowRes"/>).</param>
            /// <param name="perCleavage">Whether fragments are counted per backbone cleavage, a fragment and its
            /// neutral-loss forms once.</param>
            public ObservableRange(MsDataScan scan, Ms2ScanWithSpecificMass scanWithSpecificMass, bool isLowRes = false, bool perCleavage = false)
            {
                IsLowRes = isLowRes;
                _perCleavage = perCleavage;
                MzLibUtil.MzRange window = scan.ScanWindowRange;
                if (window == null || !(window.Maximum > window.Minimum))
                {
                    window = scan.MassSpectrum.Range;
                }
                _lowestMz = window?.Minimum ?? 0;
                _highestMz = window?.Maximum ?? double.MaxValue;
                _charge = scanWithSpecificMass.PrecursorCharge == 0 ? 1 : scanWithSpecificMass.PrecursorCharge;
                _unitCharge = System.Math.Sign(_charge);
                _precursorMass = scanWithSpecificMass.PrecursorMass;
            }

            /// <summary>Whether the scan's fragments are matched at the low-resolution tolerance.</summary>
            public bool IsLowRes { get; }

            public double LowestMass => Chemistry.ClassExtensions.ToMass(_lowestMz, _unitCharge);
            public double HighestMass => System.Math.Min(_precursorMass, Chemistry.ClassExtensions.ToMass(_highestMz, _charge));

            /// <summary>
            /// The number of products (diagnostic and precursor ions aside) the scan could show; per cleavage, the number
            /// of backbone cleavages with at least one form (intact or neutral-loss) it could show.
            /// </summary>
            public int CountObservable(List<Product> products) => ObservableMassCounts(products).Values.Sum();

            /// <summary>
            /// The products <see cref="CountObservable"/> counts, keyed by fragment type and neutral mass (to 0.1 mDa)
            /// with how many share the key, so that the fragments every candidate predicts alike can be counted. Per
            /// cleavage the mass is the intact fragment's, so a fragment and its neutral-loss forms share one key, once.
            /// </summary>
            public Dictionary<(ProductType, long), int> ObservableMassCounts(List<Product> products)
            {
                double lowest = LowestMass, highest = HighestMass;
                var counts = new Dictionary<(ProductType, long), int>(products.Count);
                foreach (Product product in products)
                {
                    double mass = product.NeutralMass;
                    if (product.ProductType is ProductType.D or ProductType.M || double.IsNaN(mass) || mass < lowest || mass > highest)
                    {
                        continue;
                    }
                    if (_perCleavage)
                    {
                        counts[(product.ProductType, CleavageKey(product))] = 1;
                    }
                    else
                    {
                        var key = (product.ProductType, (long)System.Math.Round(mass * 1e4));
                        counts[key] = counts.GetValueOrDefault(key) + 1;
                    }
                }
                return counts;
            }

            /// <summary>
            /// Records the envelopes no candidate matched, for the chance-peak model. Scans matched peak by peak
            /// (XCorr-processed) or not deconvoluted have none recorded.
            /// </summary>
            public void RecordUnmatched(LocalizationEvidenceBuilder evidence, Ms2ScanWithSpecificMass scan, CommonParameters commonParameters)
            {
                if (scan.ExperimentalFragments == null || scan.TheScan.MassSpectrum.XcorrProcessed)
                {
                    return;
                }
                evidence.RecordUnmatchedEnvelopes(scan.ExperimentalFragments,
                    IsLowRes ? commonParameters.ProductMassTolerance_LowRes : commonParameters.ProductMassTolerance, LowestMass, HighestMass, _charge);
            }
        }

        /// <summary>
        /// Whether the peptide begins its protein, counting a peptide that begins after a removed initiator methionine.
        /// </summary>
        internal static bool IsProteinNTerminal(IBioPolymerWithSetMods peptide)
        {
            return peptide.OneBasedStartResidue == 1
                || (peptide.OneBasedStartResidue == 2 && peptide.Parent?.BaseSequence is { Length: > 0 } sequence && sequence[0] == 'M');
        }

        /// <summary>
        /// Whether the peptide ends its protein.
        /// </summary>
        internal static bool IsProteinCTerminal(IBioPolymerWithSetMods peptide)
        {
            return peptide.Parent != null && peptide.OneBasedEndResidue == peptide.Parent.Length;
        }
    }
}
