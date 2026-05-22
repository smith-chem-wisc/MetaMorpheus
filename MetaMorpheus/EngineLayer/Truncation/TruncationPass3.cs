using System.Collections.Generic;
using System.Linq;
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

        /// <summary>The terminus chopped is the opposite of the matched ion series.</summary>
        public static FragmentationTerminus TerminusToChop(FragmentationTerminus winningSeries) =>
            winningSeries == FragmentationTerminus.N ? FragmentationTerminus.C : FragmentationTerminus.N;

        /// <summary>
        /// Chops, scores, and builds a Pass 3 PSM for one Pass 2 winner, or returns null if no chop
        /// reaches the precursor mass (the scan is then dropped, #9 — no mass-shift absorption).
        /// </summary>
        public static TruncationPsm ScoreTruncation(TruncationParentSelection winner, CommonParameters commonParameters, MassDiffAcceptor chopAcceptor)
        {
            FragmentationTerminus terminusToChop = TerminusToChop(winner.WinningSeries);

            ChopResult chop = ProteoformChopper.ChopUntilMassMatches(
                winner.WinningParent.Proteoform, terminusToChop, winner.Scan.PrecursorMass, chopAcceptor);
            if (chop == null)
            {
                return null;
            }

            // Pass 3 uses BOTH ion series (the single-series scoring of Pass 2 was only for parent selection, #12).
            var products = new List<Product>();
            chop.TruncatedForm.Fragment(commonParameters.DissociationType, FragmentationTerminus.Both, products, commonParameters.FragmentationParameters);
            List<MatchedFragmentIon> matchedIons = MetaMorpheusEngine.MatchFragmentIons(winner.Scan, products, commonParameters);
            double score = MetaMorpheusEngine.CalculatePeptideScore(winner.Scan.TheScan, matchedIons);

            var psm = new PeptideSpectralMatch(chop.TruncatedForm, chop.Notch, score, winner.ScanIndex, winner.Scan, commonParameters, matchedIons);

            return new TruncationPsm
            {
                SpectralMatch = psm,
                TruncatedForm = chop.TruncatedForm,
                TruncationProductType = terminusToChop == FragmentationTerminus.N ? NTerminalTruncation : CTerminalTruncation,
                ProteinAccessions = winner.WinningParent.ProteinAccession,
                ScanIndex = winner.ScanIndex
            };
        }

        /// <summary>
        /// Collapses PSMs that chopped to the same truncated FullSequence for the same scan into one PSM
        /// (#10), merging their parent protein accessions into a pipe-separated, MetaMorpheus-style list.
        /// </summary>
        public static List<TruncationPsm> CollapseDuplicateTruncations(IEnumerable<TruncationPsm> psms)
        {
            return psms
                .GroupBy(p => (p.TruncatedForm.FullSequence, p.ScanIndex))
                .Select(group =>
                {
                    TruncationPsm representative = group.First();
                    string mergedAccessions = string.Join("|",
                        group.Select(p => p.ProteinAccessions).Distinct().OrderBy(a => a));

                    return new TruncationPsm
                    {
                        SpectralMatch = representative.SpectralMatch,
                        TruncatedForm = representative.TruncatedForm,
                        TruncationProductType = representative.TruncationProductType,
                        ProteinAccessions = mergedAccessions,
                        ScanIndex = representative.ScanIndex
                    };
                })
                .ToList();
        }
    }
}
