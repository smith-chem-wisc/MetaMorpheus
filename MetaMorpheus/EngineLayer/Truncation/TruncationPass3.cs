using System.Collections.Generic;
using System.Linq;
using MassSpectrometry;
using Omics.Fragmentation;
using Proteomics.ProteolyticDigestion;

namespace EngineLayer.Truncation
{
    /// <summary>
    /// A Pass 3 truncation PSM: the standard-scored <see cref="SpectralMatch"/> for a deduced truncated
    /// form, plus the <c>Truncation Product Type</c> label (#13) and the parent protein accession(s)
    /// (pipe-separated once duplicate truncations are collapsed, #10).
    /// </summary>
    public class TruncationPsm
    {
        public SpectralMatch SpectralMatch { get; init; }
        public PeptideWithSetModifications TruncatedForm { get; init; }
        public string TruncationProductType { get; init; }
        public string ProteinAccessions { get; init; }
        public int ScanIndex { get; init; }
        public double Score => SpectralMatch.Score;
    }

    /// <summary>
    /// Pass 3 (01_Architecture.md decisions #9-#13): for each Pass 2 winner, chop the terminus OPPOSITE
    /// the winning ion series (a winning N-terminal-ion series means the C-terminus was lost) until the
    /// precursor mass matches, then score the deduced truncated form against the scan with standard
    /// both-series MetaMorpheus scoring and emit a <see cref="TruncationPsm"/>. Duplicate truncations
    /// (same truncated FullSequence + scan) collapse to one PSM, preserving protein-accession ambiguity.
    /// </summary>
    public static class TruncationPass3
    {
        public const string NTerminalTruncation = "N-terminal truncation";
        public const string CTerminalTruncation = "C-terminal truncation";

        /// <summary>
        /// Canonical initiator-methionine excision (#13) — a standard N-terminal form, NOT a truncation.
        /// Distinguished so a clean 1-Met chop is not miscounted as an N-terminal truncation (Chen 2019:
        /// ~17% of reported "truncations" were plain NME, another ~13% NME + N-acetyl).
        /// </summary>
        public const string NTerminalMetExcision = "N-terminal Met excision";

        /// <summary>Initiator-Met excision plus N-acetylation on the new N-terminus (#13).</summary>
        public const string NTerminalMetExcisionPlusAcetyl = "N-terminal Met excision + acetylation";

        /// <summary>Description-column label for a Pass 1 intact match carried into the pooled output (#4a).</summary>
        public const string FullLength = "full-length";

        /// <summary>The terminus chopped is the opposite of the matched ion series.</summary>
        public static FragmentationTerminus TerminusToChop(FragmentationTerminus winningSeries) =>
            winningSeries == FragmentationTerminus.N ? FragmentationTerminus.C : FragmentationTerminus.N;

        /// <summary>
        /// Chops, scores, and builds a Pass 3 PSM for one Pass 2 winner, or returns null if no chop
        /// reaches the precursor mass (the scan is then dropped, #9 — no mass-shift absorption).
        /// </summary>
        public static TruncationPsm ScoreTruncation(TruncationParentSelection winner, CommonParameters commonParameters,
            MassDiffAcceptor chopAcceptor, TruncationTimings timings = null)
        {
            FragmentationTerminus terminusToChop = TerminusToChop(winner.WinningSeries);

            var stopwatch = timings != null ? System.Diagnostics.Stopwatch.StartNew() : null;
            ChopResult chop = ProteoformChopper.ChopUntilMassMatches(
                winner.WinningParent.Proteoform, terminusToChop, winner.Scan.PrecursorMass, chopAcceptor);
            if (timings != null)
            {
                timings.ChoppingSeconds += stopwatch.Elapsed.TotalSeconds;
                stopwatch.Restart();
            }
            if (chop == null)
            {
                return null;
            }

            // Pass 3 uses BOTH ion series (the single-series scoring of Pass 2 was only for parent selection, #12).
            var products = new List<Product>();
            DissociationType dissociationType = TruncationSearchEngine.ResolveDissociationType(commonParameters, winner.Scan);
            chop.TruncatedForm.Fragment(dissociationType, FragmentationTerminus.Both, products, commonParameters.FragmentationParameters);
            List<MatchedFragmentIon> matchedIons = MetaMorpheusEngine.MatchFragmentIons(winner.Scan, products, commonParameters);
            double score = MetaMorpheusEngine.CalculatePeptideScore(winner.Scan.TheScan, matchedIons);
            if (timings != null)
            {
                timings.ScoringSeconds += stopwatch.Elapsed.TotalSeconds;
            }

            var psm = new PeptideSpectralMatch(chop.TruncatedForm, chop.Notch, score, winner.ScanIndex, winner.Scan, commonParameters, matchedIons);

            return new TruncationPsm
            {
                SpectralMatch = psm,
                TruncatedForm = chop.TruncatedForm,
                TruncationProductType = LabelFor(chop.TruncatedForm),
                ProteinAccessions = winner.WinningParent.ProteinAccession,
                ScanIndex = winner.ScanIndex
            };
        }

        /// <summary>
        /// Reads the truncation-type label from the chopped form's Description prefix (set by the chopper),
        /// so the PSM's product type stays consistent with what was written to the Description column —
        /// including the NME / NME+acetyl classifications (#13). Longer prefixes are tested first.
        /// </summary>
        private static string LabelFor(PeptideWithSetModifications truncatedForm)
        {
            string description = truncatedForm.Description ?? string.Empty;
            if (description.StartsWith(NTerminalMetExcisionPlusAcetyl)) return NTerminalMetExcisionPlusAcetyl;
            if (description.StartsWith(NTerminalMetExcision)) return NTerminalMetExcision;
            if (description.StartsWith(NTerminalTruncation)) return NTerminalTruncation;
            return CTerminalTruncation;
        }

        /// <summary>
        /// Carries a Pass 1 intact match into the pooled truncation output as a full-length form (#4a/#18):
        /// rebuilds the matched proteoform with a "full-length" Description (so it rides the standard
        /// Description column, #13) and re-scores it against the scan at notch 0 with both-series scoring.
        /// The base/full sequence is unchanged, so proteoform dedup still groups it correctly.
        /// </summary>
        public static SpectralMatch InheritAsFullLength(SpectralMatch pass1Psm, Ms2ScanWithSpecificMass scan,
            int scanIndex, CommonParameters commonParameters)
        {
            var source = (PeptideWithSetModifications)pass1Psm.BestMatchingBioPolymersWithSetMods.First().SpecificBioPolymer;

            // Description has no public setter; rebuild via the standard ctor with the full-length label.
            var fullLength = new PeptideWithSetModifications(source.Protein, source.DigestionParams,
                source.OneBasedStartResidueInProtein, source.OneBasedEndResidueInProtein,
                source.CleavageSpecificityForFdrCategory, FullLength, source.MissedCleavages,
                source.AllModsOneIsNterminus, source.NumFixedMods);

            var products = new List<Product>();
            DissociationType dissociationType = TruncationSearchEngine.ResolveDissociationType(commonParameters, scan);
            fullLength.Fragment(dissociationType, FragmentationTerminus.Both, products, commonParameters.FragmentationParameters);
            List<MatchedFragmentIon> matchedIons = MetaMorpheusEngine.MatchFragmentIons(scan, products, commonParameters);
            double score = MetaMorpheusEngine.CalculatePeptideScore(scan.TheScan, matchedIons);

            return new PeptideSpectralMatch(fullLength, 0, score, scanIndex, scan, commonParameters, matchedIons);
        }

        /// <summary>
        /// Collapses PSMs that chopped to the same truncated FullSequence for the same scan into one
        /// <see cref="SpectralMatch"/> (#10): the alternates are added to the representative as ambiguous
        /// matches (equal score), so the standard writer pipe-merges their parent protein accessions
        /// natively once <see cref="SpectralMatch.ResolveAllAmbiguities"/> runs.
        /// </summary>
        public static List<SpectralMatch> CollapseDuplicateTruncations(IEnumerable<TruncationPsm> psms)
        {
            var collapsed = new List<SpectralMatch>();

            foreach (var group in psms.GroupBy(p => (p.TruncatedForm.FullSequence, p.ScanIndex)))
            {
                List<TruncationPsm> members = group.ToList();
                SpectralMatch representative = members[0].SpectralMatch;

                for (int i = 1; i < members.Count; i++)
                {
                    representative.AddOrReplace(members[i].TruncatedForm, representative.Score,
                        members[i].SpectralMatch.Notch ?? 0, reportAllAmbiguity: true, new List<MatchedFragmentIon>());
                }

                collapsed.Add(representative);
            }

            return collapsed;
        }
    }
}
