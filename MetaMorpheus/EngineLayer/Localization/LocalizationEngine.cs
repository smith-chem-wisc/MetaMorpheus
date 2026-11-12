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
        /// motif does not allow, typically the search's variable modifications.</param>
        /// <param name="maxPlacements">The most arrangements enumerated for one match; a match with more gets no
        /// placement evidence.</param>
        public LocalizationEngine(IEnumerable<SpectralMatch> allResultingIdentifications, MsDataFile myMsDataFile, CommonParameters commonParameters,
            List<(string fileName, CommonParameters fileSpecificParameters)> fileSpecificParameters, List<string> nestedIds,
            bool localizeMassDifferences, IReadOnlyCollection<Modification> placeableModifications, int maxPlacements = PositionalIsoforms.DefaultMaxPlacements)
            : base(commonParameters, fileSpecificParameters, nestedIds)
        {
            AllResultingIdentifications = allResultingIdentifications;
            MyMsDataFile = myMsDataFile;
            LocalizeMassDifferences = localizeMassDifferences;
            PlaceableModifications = placeableModifications;
            MaxPlacements = maxPlacements;
        }

        protected override MetaMorpheusEngineResults RunSpecific()
        {
            bool placeModifications = PlaceableModifications != null;

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
                        bool placeThisMatchesModifications = placeModifications && CommonParameters.PrecursorMassTolerance.Within(observedMass, peptide.MonoisotopicMass);

                        if (LocalizeMassDifferences && psm.FullSequence != null)
                        {
                            bool keepEvidence = basePeakIntensity > 0 && !placeThisMatchesModifications;
                            LocalizeMassDifference(psm, peptide, massDifference, scan, scanWithSpecificMass, products,
                                keepEvidence ? new LocalizationEvidenceBuilder(basePeakIntensity) : null);
                        }

                        if (placeThisMatchesModifications && basePeakIntensity > 0)
                        {
                            if (!PlaceModifications(psm, peptide, massDifference, scan, scanWithSpecificMass, products, basePeakIntensity, IsFixed))
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
        /// stand for the peptide termini, which fragment identically.
        /// </summary>
        private void LocalizeMassDifference(SpectralMatch psm, IBioPolymerWithSetMods peptide, double massDifference, MsDataScan scan,
            Ms2ScanWithSpecificMass scanWithSpecificMass, List<Product> products, LocalizationEvidenceBuilder evidence)
        {
            var localizedScores = new List<double>();
            for (int r = 0; r < peptide.Length; r++)
            {
                // create new PeptideWithSetMods with unidentified mass difference at the given residue
                var peptideWithLocalizedMassDiff = peptide.Localize(r, massDifference);

                // this is the list of theoretical products for this peptide with mass-difference on this residue
                peptideWithLocalizedMassDiff.Fragment(CommonParameters.DissociationType, CommonParameters.DigestionParams.FragmentationTerminus, products, CommonParameters.FragmentationParameters);

                var matchedIons = MatchFragmentIons(scanWithSpecificMass, products, CommonParameters);

                // score when the mass-diff is on this residue
                double localizedScore = CalculatePeptideScore(scan, matchedIons);

                localizedScores.Add(localizedScore);
                evidence?.AddCandidate(new LocalizationPlacement(new[] { new LocalizationSite(r, peptide.BaseSequence[r]) }), localizedScore, matchedIons);
            }

            psm.LocalizedScores = localizedScores;
            if (evidence != null)
            {
                psm.LocalizationEvidence = evidence.Build(massDifference, IsProteinNTerminal(peptide), IsProteinCTerminal(peptide));
            }
        }

        /// <summary>
        /// Keeps what each arrangement of the match's identified modifications matched. A match with nothing to move has
        /// one arrangement, kept so that its matched peaks are there to fit a model of matched-peak intensities from.
        /// </summary>
        /// <returns>False when the match has more arrangements than are enumerated, leaving it without evidence.</returns>
        private bool PlaceModifications(SpectralMatch psm, IBioPolymerWithSetMods peptide, double massDifference, MsDataScan scan,
            Ms2ScanWithSpecificMass scanWithSpecificMass, List<Product> products, double basePeakIntensity,
            System.Func<Modification, bool> isFixed)
        {
            if (!PositionalIsoforms.TryEnumerate(peptide, PlaceableModifications, isFixed, MaxPlacements, out var isoforms))
            {
                return false;
            }

            var evidence = new LocalizationEvidenceBuilder(basePeakIntensity);
            foreach (PositionalIsoform isoform in isoforms)
            {
                IBioPolymerWithSetMods arranged = PositionalIsoforms.WithModifications(peptide, isoform.AllModsOneIsNterminus);
                if (arranged == null)
                {
                    return true; // a kind of biopolymer that cannot be rebuilt; nothing to place
                }

                arranged.Fragment(CommonParameters.DissociationType, CommonParameters.DigestionParams.FragmentationTerminus, products, CommonParameters.FragmentationParameters);
                var matchedIons = MatchFragmentIons(scanWithSpecificMass, products, CommonParameters);
                evidence.AddCandidate(isoform.Placement, CalculatePeptideScore(scan, matchedIons), matchedIons);
            }

            psm.LocalizationEvidence = evidence.Build(massDifference, IsProteinNTerminal(peptide), IsProteinCTerminal(peptide));
            return true;
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
