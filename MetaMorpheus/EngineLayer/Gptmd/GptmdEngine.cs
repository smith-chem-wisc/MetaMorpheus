using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Linq;
using Omics;
using Omics.BioPolymer;
using Omics.Modifications;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System.Threading.Tasks;
using System.Collections.Concurrent;
using MassSpectrometry;
using Omics.Fragmentation;
using System.Threading;

namespace EngineLayer.Gptmd
{
    public class GptmdEngine : MetaMorpheusEngine
    {
        // It is assumed that there is a one-to-one correspondence between the psms in AllIdentifications and the scans in AllScans
        private readonly List<SpectralMatch> AllIdentifications;
        private readonly IEnumerable<Tuple<double, double>> Combos;
        private readonly List<Modification> GptmdModifications;
        private readonly Dictionary<string, Tolerance> FilePathToPrecursorMassTolerance; // this exists because of file-specific tolerances
        //The ScoreTolerance property is used to differentiatie when a PTM candidate is added to a peptide. We check the score at each position and then add that mod where the score is highest.
        private readonly double ScoreTolerance = 0.1;
        private static readonly double QValueNotchThreshold = 0.05;
        public Dictionary<string, HashSet<Tuple<int, Modification>>> ModDictionary { get; init; }
        private readonly List<IGptmdFilter> Filters;
        private readonly bool AddAlternativeSitesForDatabaseMods;
        private readonly List<Modification> FixedModifications;

        public GptmdEngine(
            List<SpectralMatch> allIdentifications,
            List<Modification> gptmdModifications,
            IEnumerable<Tuple<double, double>> combos,
            Dictionary<string, Tolerance> filePathToPrecursorMassTolerance,
            CommonParameters commonParameters,
            List<(string fileName, CommonParameters fileSpecificParameters)> fileSpecificParameters,
            List<string> nestedIds,
            Dictionary<string, HashSet<Tuple<int, Modification>>> modDictionary,
            List<IGptmdFilter> filters = null,
            bool addAlternativeSitesForDatabaseMods = false,
            List<Modification> fixedModifications = null)
            : base(commonParameters, fileSpecificParameters, nestedIds)
        {
            AllIdentifications = allIdentifications;
            GptmdModifications = gptmdModifications;
            Combos = combos;
            FilePathToPrecursorMassTolerance = filePathToPrecursorMassTolerance;
            ModDictionary = modDictionary ?? new Dictionary<string, HashSet<Tuple<int, Modification>>>();
            Filters = filters ?? new List<IGptmdFilter>();
            AddAlternativeSitesForDatabaseMods = addAlternativeSitesForDatabaseMods;
            FixedModifications = fixedModifications ?? new List<Modification>();
        }

        public static bool ModFits(Modification attemptToLocalize, IBioPolymer protein, int peptideOneBasedIndex, int peptideLength, int proteinOneBasedIndex)
        {
            var motif = attemptToLocalize.Target;
            var hehe = motif.ToString().IndexOf(motif.ToString().First(b => char.IsUpper(b)));
            var proteinToMotifOffset = proteinOneBasedIndex - hehe - 1;
            var indexUp = 0;
            // Look up starting at and including the capital letter

            while (indexUp < motif.ToString().Length)
            {
                if (indexUp + proteinToMotifOffset < 0 || indexUp + proteinToMotifOffset >= protein.Length || (!char.ToUpper(motif.ToString()[indexUp]).Equals('X') && !char.ToUpper(motif.ToString()[indexUp]).Equals(protein.BaseSequence[indexUp + proteinToMotifOffset])))
                    return false;
                indexUp++;
            }
            // if a UniProt mod already exists at this location with the same mass, don't annotate the GPTMD mod
            if (protein.OneBasedPossibleLocalizedModifications.TryGetValue(proteinOneBasedIndex, out List<Modification> modsAtThisLocation)
                && modsAtThisLocation.Any(m => m.ModificationType == "UniProt" && Math.Abs(m.MonoisotopicMass.Value - attemptToLocalize.MonoisotopicMass.Value) < 0.005))
            {
                return false;
            }

            if (attemptToLocalize.LocationRestriction == "Anywhere.")
                return true;
            if (attemptToLocalize.LocationRestriction is "N-terminal." or "5'-terminal." && (proteinOneBasedIndex <= 2))
                return true;
            if (attemptToLocalize.LocationRestriction is "Peptide N-terminal." or "Oligo 5'-terminal." && peptideOneBasedIndex == 1)
                return true;
            if (attemptToLocalize.LocationRestriction is "Peptide C-terminal." or "Oligo 3'-terminal." && peptideOneBasedIndex == peptideLength)
                return true;
            if (attemptToLocalize.LocationRestriction is "C-terminal." or "3'-terminal." && proteinOneBasedIndex == protein.Length)
                return true;
            return false;
        }

        protected override MetaMorpheusEngineResults RunSpecific()
        {
            var modDict = new ConcurrentDictionary<string, ConcurrentBag<Tuple<int, Modification>>>();
            int modsAdded = 0;

            int maxThreadsPerFile = CommonParameters.MaxThreadsToUsePerFile;

            // Keep psms if they are below q value notch threshold, or notch ambiguous results with hypothesis below threshold. 
            var psms = AllIdentifications.Where(b => b.FdrInfo.QValueNotch <= QValueNotchThreshold || (Math.Abs(b.FdrInfo.QValueNotch - 2) < 0.001 && b.BestMatchingBioPolymersWithSetMods.Any(p => p.QValueNotch.HasValue && p.QValueNotch <= QValueNotchThreshold))).ToList();
            if (psms.Any() == false)
            {
                return new GptmdResults(this, ModDictionary, 0);
            }
            Parallel.ForEach(Partitioner.Create(0, psms.Count), new ParallelOptions() { MaxDegreeOfParallelism = maxThreadsPerFile }, (range) =>
            {
                for (int i = range.Item1; i < range.Item2; i++)
                {
                    // Extract all necessary information from the PSM
                    var psm = psms[i];
                    var dissociationType = CommonParameters.DissociationType == DissociationType.Autodetect ?
                        psms[i].Ms2Scan.DissociationType.Value : CommonParameters.DissociationType;
                    var scan = psm.Ms2Scan;
                    var precursorMass = psm.ScanPrecursorMass;
                    var precursorCharge = psm.ScanPrecursorCharge;
                    var fileName = psm.FullFilePath;
                    var originalScore = psm.Score;
                    // Built on first use and shared by both passes: its constructor deconvolutes the scan, and most PSMs
                    // never have a candidate to score.
                    var ms2ScanWithSpecificMass = new Lazy<Ms2ScanWithSpecificMass>(
                        () => new Ms2ScanWithSpecificMass(scan, precursorMass, precursorCharge, fileName, CommonParameters), isThreadSafe: false);
                    var peptideTheorProducts = new List<Product>();
                    List<(int site, Modification mod, string proteinAccession)> bestMatches = [];

                    foreach (var hypothesis in psm.BestMatchingBioPolymersWithSetMods.Select(v => v))
                    {
                        // skip hypotheses that are above the q value notch threshold - will be null unless psm is notch ambiguous
                        if (hypothesis.QValueNotch.HasValue && hypothesis.QValueNotch > QValueNotchThreshold)
                            continue;

                        bestMatches.Clear();
                        var pepWithSetMods = hypothesis.SpecificBioPolymer;

                        // Each candidate mass (peptide + modification) is tested against the observation the way
                        // this search matches: monoisotopic-to-monoisotopic by default, apex-to-apex within the
                        // isotopologue apex tolerance in most-abundant mode. Testing per candidate - rather than
                        // correcting the precursor mass once - is what keeps a modification mass from being
                        // confused with a whole-isotopologue apex misprediction: both land in the same residual
                        // and cannot be separated after the fact.
                        var possibleModifications = GetPossibleMods(psm, GptmdModifications, Combos,
                            FilePathToPrecursorMassTolerance[fileName], pepWithSetMods, CommonParameters);

                        double bestScore = 0;
                        CleavageBlockingPolicy? cleavageBlocking = null; // built on first use, once per peptide

                        foreach (var mod in possibleModifications)
                        {
                            if (!mod.MonoisotopicMass.HasValue)
                                continue;

                            for (int pepSeqIndex = 0; pepSeqIndex < pepWithSetMods.Length; pepSeqIndex++)
                            {
                                int indexInProtein = pepWithSetMods.OneBasedStartResidue + pepSeqIndex;
                                if (!ModFits(mod, pepWithSetMods.Parent, pepSeqIndex + 1, pepWithSetMods.Length, indexInProtein))
                                    continue;
                                if (pepSeqIndex == pepWithSetMods.Length - 1
                                    && (cleavageBlocking ??= CleavageBlockingFor(pepWithSetMods, GptmdModifications)).Blocks(mod, indexInProtein))
                                    continue;

                                var newPep = pepWithSetMods.Localize(pepSeqIndex, mod.MonoisotopicMass.Value);
                                double score = ScoreCandidate(newPep, dissociationType, ms2ScanWithSpecificMass.Value, peptideTheorProducts, out var matchedIons);

                                int modSiteInPeptideResidue = pepSeqIndex + 1;
                                if (Filters.Any(f => !f.Passes(newPep, psm, score, originalScore, matchedIons, modSiteInPeptideResidue, pepWithSetMods.Length, mod)))
                                    continue;

                                if (score < bestScore - ScoreTolerance)
                                    continue;

                                var (adjustedSite, accession) = ResolveSiteInDatabase(pepWithSetMods, indexInProtein);

                                if (score > bestScore + ScoreTolerance) // new high score, reset list
                                {
                                    bestScore = score;
                                    bestMatches.Clear();
                                    bestMatches.Add((adjustedSite, mod, accession));
                                }
                                else if (Math.Abs(score - bestScore) <= ScoreTolerance)
                                {
                                    bestMatches.Add((adjustedSite, mod, accession));
                                }
                            }
                        }

                        foreach (var match in bestMatches)
                        {
                            AddIndexedMod(modDict, match.proteinAccession, new Tuple<int, Modification>(match.site, match.mod));
                            Interlocked.Increment(ref modsAdded);
                        }

                        if (!AddAlternativeSitesForDatabaseMods)
                            continue;

                        foreach (var match in FindAlternativeSitesForDatabaseMods(psm, pepWithSetMods, FilePathToPrecursorMassTolerance[fileName], dissociationType,
                                     ms2ScanWithSpecificMass, peptideTheorProducts))
                        {
                            AddIndexedMod(modDict, match.proteinAccession, new Tuple<int, Modification>(match.site, match.mod));
                            Interlocked.Increment(ref modsAdded);
                        }
                    }
                }
            });

            //UpdateModDictionary(modDict);
            foreach(var kvp in modDict)
                ModDictionary.MergeOrCreate(kvp.Key, new HashSet<Tuple<int, Modification>>(kvp.Value));
            return new GptmdResults(this, ModDictionary, modsAdded);
        }

        /// <summary>
        /// Other sites for the modifications the database put on this peptide.
        /// <para>
        /// GPTMD otherwise looks only at precursor mass a peptide leaves unexplained. A spectrum matched to a database
        /// modification (a UniProt phosphoserine, say) leaves none, so that modification is never tried at another
        /// residue, and the search can only report the annotated one.
        /// </para>
        /// <para>
        /// Here, when the precursor matches the peptide with its database modifications as they are, each one is moved
        /// (not duplicated) to every other residue where a GPTMD modification of the same chemistry fits. Each move is
        /// scored and filtered as an ordinary GPTMD candidate is. A move to the same amino acid is written as the
        /// database's own modification, so the two placements carry the same neutral losses and diagnostic ions in the
        /// search; a move to another amino acid is written as the GPTMD modification that fits there.
        /// </para>
        /// <para>
        /// The bar is the best placement the database already offers: the original, and any other residue of the
        /// peptide where the database has the same chemistry (UniProt's, or an earlier round's), all scored the same
        /// way. Every move within <see cref="ScoreTolerance"/> of the best score, the bar included, is returned: a
        /// spectrum that cannot tell two residues of one amino acid apart puts both in the database, so the search can
        /// report the site as ambiguous. A move to another amino acid must beat the bar by more than the tolerance: a tie
        /// is no evidence for it, and the same formula on another residue is often another modification
        /// (4-hydroxyproline and Hydroxylation on N). Offered placements are never added again, so a second round over
        /// the database this one writes adds nothing new from the same spectra, whatever names its sites carry.
        /// </para>
        /// <para>
        /// A modification on a residue a fixed modification targets is not moved: the search puts the fixed modification
        /// back there, so the moved form would no longer match this precursor.
        /// </para>
        /// </summary>
        private List<(int site, Modification mod, string proteinAccession)> FindAlternativeSitesForDatabaseMods(SpectralMatch psm,
            IBioPolymerWithSetMods peptide, Tolerance precursorTolerance, DissociationType dissociationType,
            Lazy<Ms2ScanWithSpecificMass> scan, List<Product> theoreticalProducts)
        {
            var alternatives = new List<(int site, Modification mod, string proteinAccession)>();
            if (!psm.MatchesCandidateMass(peptide.MonoisotopicMass, precursorTolerance, CommonParameters))
                return alternatives;

            foreach (var (key, databaseMod) in peptide.AllModsOneIsNterminus)
            {
                // zero-based residue; modifications on the termini (keys 1 and Length + 2) are not moved
                int fromIndex = key - 2;
                if (fromIndex < 0 || fromIndex >= peptide.Length || !databaseMod.MonoisotopicMass.HasValue
                    || !IsDatabaseMod(peptide, fromIndex, databaseMod)
                    || FixedModificationFits(peptide, fromIndex))
                    continue;

                var sameChemistry = GptmdModifications.Where(m => m.ValidModification == true && IsSameChemistry(m, databaseMod)).ToList();
                if (sameChemistry.Count == 0)
                    continue;

                double mass = databaseMod.MonoisotopicMass.Value;

                // the original placement, scored the way every GPTMD candidate is (as a mass-only modification)
                double originalPlacementScore = ScoreCandidate(peptide.Localize(fromIndex, 0), dissociationType, scan.Value, theoreticalProducts, out _);
                double bestOfferedScore = originalPlacementScore;
                var moves = new List<(int site, Modification mod, string proteinAccession, double score, bool otherAminoAcid)>();

                for (int toIndex = 0; toIndex < peptide.Length; toIndex++)
                {
                    // a residue that already carries a modification cannot take a second one in the search
                    if (toIndex == fromIndex || peptide.AllModsOneIsNterminus.ContainsKey(toIndex + 2))
                        continue;

                    int indexInProtein = peptide.OneBasedStartResidue + toIndex;
                    bool otherAminoAcid = peptide.BaseSequence[toIndex] != peptide.BaseSequence[fromIndex];
                    bool offered = HasSameChemistryAt(peptide.Parent, indexInProtein, databaseMod);
                    // On the same amino acid the database's own modification is written, so the search scores both
                    // placements with the same fragments: UniProt's Phosphoserine carries no neutral loss, GPTMD's
                    // Phosphorylation on S carries H3PO4, and those losses would outvote the ions that place the site.
                    var fitting = offered ? new List<Modification>()
                        : !otherAminoAcid && ModFits(databaseMod, peptide.Parent, toIndex + 1, peptide.Length, indexInProtein) ? new List<Modification> { databaseMod }
                        : sameChemistry.Where(m => ModFits(m, peptide.Parent, toIndex + 1, peptide.Length, indexInProtein)).ToList();
                    if (!offered && fitting.Count == 0)
                        continue;

                    // the moved form is the same peptidoform whichever same-chemistry modification is written
                    var moved = peptide.Localize(fromIndex, -mass).Localize(toIndex, mass);
                    double score = ScoreCandidate(moved, dissociationType, scan.Value, theoreticalProducts, out var matchedIons);
                    if (offered)
                    {
                        bestOfferedScore = Math.Max(bestOfferedScore, score);
                        continue;
                    }

                    var (site, accession) = ResolveSiteInDatabase(peptide, indexInProtein);
                    foreach (var mod in fitting.Where(m => !Filters.Any(f => !f.Passes(moved, psm, score, originalPlacementScore, matchedIons, toIndex + 1, peptide.Length, m))))
                        moves.Add((site, mod, accession, score, otherAminoAcid));
                }

                var eligible = moves.Where(m => !m.otherAminoAcid || m.score > bestOfferedScore + ScoreTolerance)
                    .Select(m => ((m.site, m.mod, m.proteinAccession), m.score));
                alternatives.AddRange(WithinToleranceOfBest(eligible, bestOfferedScore, ScoreTolerance));
            }

            return alternatives;
        }

        /// <summary>
        /// Whether the database, rather than a fixed or variable modification, put this modification at this residue.
        /// <see cref="Modification.Equals(object)"/> compares IdWithMotif, OriginalId, ModificationType and mass, so a
        /// copy read from a database's own modification list matches too.
        /// </summary>
        private static bool IsDatabaseMod(IBioPolymerWithSetMods peptide, int zeroBasedResidue, Modification mod)
        {
            return peptide.Parent.OneBasedPossibleLocalizedModifications.TryGetValue(peptide.OneBasedStartResidue + zeroBasedResidue, out var modsHere)
                && modsHere.Contains(mod);
        }

        /// <summary>
        /// Whether one of the search's fixed modifications targets this residue of the peptide.
        /// </summary>
        private bool FixedModificationFits(IBioPolymerWithSetMods peptide, int zeroBasedResidue)
        {
            return FixedModifications.Any(f => f.LocationRestriction == "Anywhere."
                && ModificationLocalization.ModFits(f, peptide.Parent.BaseSequence, zeroBasedResidue + 1, peptide.Length, peptide.OneBasedStartResidue + zeroBasedResidue));
        }

        /// <summary>
        /// Whether the database already offers a modification of this chemistry at this residue, from UniProt, another
        /// annotation, or an earlier round of GPTMD.
        /// </summary>
        private static bool HasSameChemistryAt(IBioPolymer protein, int oneBasedResidue, Modification mod)
        {
            return protein.OneBasedPossibleLocalizedModifications.TryGetValue(oneBasedResidue, out var modsHere)
                && modsHere.Any(m => IsSameChemistry(m, mod));
        }

        /// <summary>
        /// Same chemical formula, or, when either has none, the same mass within the 0.005 Da that <see cref="ModFits"/>
        /// uses. UniProt's "Phosphoserine" and GPTMD's "Phosphorylation on T" are the same chemistry.
        /// </summary>
        private static bool IsSameChemistry(Modification a, Modification b)
        {
            if (a.ChemicalFormula != null && b.ChemicalFormula != null)
                return a.ChemicalFormula.Equals(b.ChemicalFormula);
            return a.MonoisotopicMass.HasValue && b.MonoisotopicMass.HasValue
                && Math.Abs(a.MonoisotopicMass.Value - b.MonoisotopicMass.Value) < 0.005;
        }

        private double ScoreCandidate(IBioPolymerWithSetMods candidate, DissociationType dissociationType, Ms2ScanWithSpecificMass scan,
            List<Product> theoreticalProducts, out List<MatchedFragmentIon> matchedIons)
        {
            theoreticalProducts.Clear();
            candidate.Fragment(dissociationType, CommonParameters.DigestionParams.FragmentationTerminus, theoreticalProducts, CommonParameters.FragmentationParameters);
            matchedIons = MatchFragmentIons(scan, theoreticalProducts, CommonParameters, matchAllCharges: false);
            return CalculatePeptideScore(scan.TheScan, matchedIons);
        }

        /// <summary>
        /// The accession and one-based site under which a modification found at <paramref name="indexInProtein"/> (in the
        /// peptide's protein, its variants applied) is written, so that read back from the GPTMD database it lands on the
        /// residue it was found on.
        /// <para>
        /// A sequence variant carries its own modifications (<see cref="SequenceVariation.OneBasedModifications"/>, keyed
        /// by position in the protein with that variant applied), which ProteinDbWriter writes from the entry under the
        /// single-variant accession. A site goes there when the variant edits its residue, or when the peptide spans that
        /// variant's edited residues (for a deletion, the junction it leaves) and no other variant's: the site was then
        /// seen only with the variant. Any other
        /// site goes to the consensus protein at its consensus position. A peptide spanning two variants' edits needs
        /// both, which neither variant's modifications can say, so its sites off the edits go to the consensus.
        /// </para>
        /// </summary>
        private static (int site, string accession) ResolveSiteInDatabase(IBioPolymerWithSetMods peptide, int indexInProtein)
        {
            var protein = peptide.Parent;
            if (protein.AppliedSequenceVariations is not { Count: > 0 })
                return (indexInProtein, protein.Accession);

            int consensusSite = VariantApplication.RestoreModificationIndex(protein, indexInProtein);
            var variant = protein.AppliedSequenceVariations.FirstOrDefault(v => VariantApplication.IsSequenceVariantModification(v, indexInProtein));
            if (variant == null)
            {
                var spanned = protein.AppliedSequenceVariations.Where(v => PeptideHoldsVariant(peptide, v)).ToList();
                if (spanned.Count != 1)
                    return (consensusSite, protein.ConsensusVariant.Accession);
                variant = spanned[0];
            }

            // the variant's own position: every other applied variant's length change undone, this one's kept
            int siteInVariant = indexInProtein > LastEditedResidue(variant)
                ? consensusSite + variant.VariantSequence.Length - variant.OriginalSequence.Length
                : consensusSite;

            // the accession the writer looks for is built from the variant as the consensus protein lists it
            int consensusBegin = VariantApplication.RestoreModificationIndex(protein, variant.OneBasedBeginPosition);
            var listed = protein.ConsensusVariant.SequenceVariations.FirstOrDefault(v => v.OneBasedBeginPosition == consensusBegin
                && v.OriginalSequence == variant.OriginalSequence && v.VariantSequence == variant.VariantSequence) ?? variant;
            return (siteInVariant, VariantApplication.GetAccession(protein, new[] { listed }));
        }

        /// <summary>
        /// The positions, in the variant-applied protein, of the residues an applied variant edits (not those an anchored
        /// indel keeps).
        /// </summary>
        private static IEnumerable<int> EditedResidues(SequenceVariation appliedVariant)
        {
            return Enumerable.Range(appliedVariant.OneBasedBeginPosition, appliedVariant.OneBasedEndPosition - appliedVariant.OneBasedBeginPosition + 1)
                .Where(r => VariantApplication.IsSequenceVariantModification(appliedVariant, r));
        }

        /// <summary>
        /// The last residue, in the variant-applied protein, an applied variant edits. A deletion that keeps every
        /// residue it leaves (ID -> I) edits none, and for it this is the residue just before the deleted ones were.
        /// Residues after it are downstream of the variant.
        /// </summary>
        private static int LastEditedResidue(SequenceVariation appliedVariant)
        {
            var edited = EditedResidues(appliedVariant).ToList();
            if (edited.Count > 0)
                return edited.Max();

            // mzLib counts a deletion's kept residues from its start first (ID -> I keeps the I, DI -> I keeps none there)
            int keptPrefix = 0;
            while (keptPrefix < appliedVariant.VariantSequence.Length
                   && appliedVariant.OriginalSequence[keptPrefix] == appliedVariant.VariantSequence[keptPrefix])
                keptPrefix++;
            return appliedVariant.OneBasedBeginPosition + keptPrefix - 1;
        }

        /// <summary>
        /// Whether the peptide's sequence exists only with this applied variant: it covers a residue the variant edits,
        /// or, for a deletion, the residues either side of where the deleted ones were (or the protein end they reach).
        /// </summary>
        private static bool PeptideHoldsVariant(IBioPolymerWithSetMods peptide, SequenceVariation appliedVariant)
        {
            var edited = EditedResidues(appliedVariant).ToList();
            if (edited.Count > 0)
                return edited.Any(r => peptide.OneBasedStartResidue <= r && r <= peptide.OneBasedEndResidue);

            int beforeJunction = LastEditedResidue(appliedVariant);
            return peptide.OneBasedStartResidue <= Math.Max(beforeJunction, 1)
                && Math.Min(beforeJunction + 1, peptide.Parent.Length) <= peptide.OneBasedEndResidue;
        }

        /// <summary>
        /// The candidates scoring within <paramref name="tolerance"/> of the best score, the best taken over the candidates
        /// and <paramref name="baselineScore"/> together (the baseline itself is not returned). Every candidate is scored
        /// before any is dropped, so the order they come in does not change which are kept.
        /// </summary>
        public static List<T> WithinToleranceOfBest<T>(IEnumerable<(T Candidate, double Score)> scored, double baselineScore, double tolerance)
        {
            var all = scored.ToList();
            double best = all.Select(s => s.Score).Append(baselineScore).Max();
            return all.Where(s => s.Score >= best - tolerance).Select(s => s.Candidate).ToList();
        }

        /// <summary>
        /// The cleavage-blocking policy the search's digestion applied to <paramref name="peptide"/>. Its
        /// Blocks(mod, residue) on the peptide's last residue is true when the modification there would block the
        /// cleavage that made the peptide's C-terminus -- trypsin does not cut after an acetylated lysine, so a tryptic
        /// peptide cannot end in one. GPTMD then must not localize it there, and tries the peptide's other residues instead.
        /// </summary>
        /// <remarks>
        /// The search's digestion drops exactly these peptidoforms when RespectCleavageBlockingModifications is on, and
        /// this asks mzLib's CleavageBlockingPolicy the same question, read from the digestion settings the peptide was
        /// produced under (the file-specific ones). So GPTMD never writes a site the search would refuse to identify,
        /// and it is inert wherever the search is: flag off, a semi- or non-specific search, the protein's own
        /// C-terminus, a protease that does not cut after the residue, a site a wildcard motif also cuts, and RNA.
        ///
        /// One policy serves every candidate on the peptide. Built from all GPTMD modifications, it is active whenever
        /// a policy built from the one candidate would be, and Blocks still asks whether that candidate blocks.
        /// </remarks>
        private static CleavageBlockingPolicy CleavageBlockingFor(IBioPolymerWithSetMods peptide, List<Modification> gptmdModifications) =>
            peptide.DigestionParams is DigestionParams digestionParams && peptide.Parent is Protein protein
                ? CleavageBlockingPolicy.For(digestionParams, gptmdModifications, protein)
                : default;

        private static void AddIndexedMod(ConcurrentDictionary<string, ConcurrentBag<Tuple<int, Modification>>> modDict, string proteinAccession, Tuple<int, Modification> indexedMod)
        {
            modDict.AddOrUpdate(proteinAccession,
                new ConcurrentBag<Tuple<int, Modification>> { indexedMod },
                (key, existingBag) =>
                {
                    existingBag.Add(indexedMod);
                    return existingBag;
                });
        }

        /// <summary>
        /// The modifications whose mass the observed precursor supports for this candidate.
        /// <para>
        /// Each candidate mass (peptide + already-accumulated shift + modification) is tested against the
        /// observation via <see cref="PrecursorMassExtensions.MatchesCandidateMass"/>, which applies the
        /// current search's matching rule: monoisotopic-to-monoisotopic by default, apex-to-apex within the
        /// isotopologue apex tolerance in most-abundant mode. Previously this compared a single precursor
        /// mass against each candidate, which in most-abundant mode let a whole-isotopologue apex
        /// misprediction masquerade as a ~1-2 Da modification.
        /// </para>
        /// <paramref name="candidateMassOffset"/> accumulates the shift already committed to by an enclosing
        /// combo, replacing the old "subtract it from the target mass" recursion — equivalent arithmetic,
        /// but it keeps the OBSERVED mass untouched so the apex test always sees what was really measured.
        /// </summary>
        private static IEnumerable<Modification> GetPossibleMods(SpectralMatch psm, IEnumerable<Modification> allMods, IEnumerable<Tuple<double, double>> combos, Tolerance precursorTolerance, IBioPolymerWithSetMods peptideWithSetModifications, CommonParameters commonParameters, double candidateMassOffset = 0)
        {
            double peptideMass = peptideWithSetModifications.MonoisotopicMass + candidateMassOffset;

            foreach (var Mod in allMods.Where(b => b.ValidModification == true))
            {
                //TODO: not necessarily here. I think we're creating ambiguity. If we're going to add a gptmd mod to a peptide that already has that mod, then we need info
                // to suggest that it is at a postion other than that in the database. could be presence of frag for unmodified or presence of frag with modified at alternative location.
                if (psm.MatchesCandidateMass(peptideMass + (double)Mod.MonoisotopicMass, precursorTolerance, commonParameters))
                    yield return Mod;
                foreach (var modOnPsm in peptideWithSetModifications.AllModsOneIsNterminus.Values.Where(b => b.ValidModification == true))
                    if (modOnPsm.Target.Equals(Mod.Target))
                    {
                        if (psm.MatchesCandidateMass(peptideMass + (double)Mod.MonoisotopicMass - (double)modOnPsm.MonoisotopicMass, precursorTolerance, commonParameters))
                            yield return Mod;
                    }
            }

            foreach (var combo in combos)
            {
                var m1 = combo.Item1;
                var m2 = combo.Item2;
                var combined = m1 + m2;
                if (psm.MatchesCandidateMass(peptideMass + combined, precursorTolerance, commonParameters))
                {
                    foreach (var mod in GetPossibleMods(psm, allMods, combos, precursorTolerance, peptideWithSetModifications, commonParameters, candidateMassOffset + m1))
                        yield return mod;
                }
            }
        }
    }
}
