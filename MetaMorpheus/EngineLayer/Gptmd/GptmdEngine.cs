using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Linq;
using Omics;
using Omics.Modifications;
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
        public Dictionary<string, HashSet<Tuple<int, Modification>>> ModDictionary { get; init; }
        private readonly List<IGptmdFilter> Filters;

        public GptmdEngine(
            List<SpectralMatch> allIdentifications, 
            List<Modification> gptmdModifications, 
            IEnumerable<Tuple<double, double>> combos, 
            Dictionary<string, Tolerance> filePathToPrecursorMassTolerance, 
            CommonParameters commonParameters, 
            List<(string fileName, CommonParameters fileSpecificParameters)> fileSpecificParameters, 
            List<string> nestedIds,
            Dictionary<string, HashSet<Tuple<int, Modification>>> modDictionary,
            List<IGptmdFilter> filters = null) 
            : base(commonParameters, fileSpecificParameters, nestedIds)
        {
            AllIdentifications = allIdentifications;
            GptmdModifications = gptmdModifications;
            Combos = combos;
            FilePathToPrecursorMassTolerance = filePathToPrecursorMassTolerance;
            ModDictionary = modDictionary ?? new Dictionary<string, HashSet<Tuple<int, Modification>>>();
            Filters = filters ?? new List<IGptmdFilter>();
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
            var psms = AllIdentifications.Where(b => b.FdrInfo.QValueNotch <= 0.05 && !b.IsDecoy).ToList();
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
                    Ms2ScanWithSpecificMass ms2ScanWithSpecificMass = null;
                    var peptideTheorProducts = new List<Product>();
                    List<(int site, Modification mod, string proteinAccession)> bestMatches = [];

                    foreach (var pepWithSetMods in psm.BestMatchingBioPolymersWithSetMods.Select(v => v.SpecificBioPolymer))
                    {
                        bestMatches.Clear();
                        var isVariantProtein = pepWithSetMods.Parent != pepWithSetMods.Parent.ConsensusVariant;
                        var possibleModifications = GetPossibleMods(precursorMass, GptmdModifications, Combos,
                            FilePathToPrecursorMassTolerance[fileName], pepWithSetMods);

                        double bestScore = 0;

                        foreach (var mod in possibleModifications)
                        {
                            if (!mod.MonoisotopicMass.HasValue)
                                continue;

                            for (int pepSeqIndex = 0; pepSeqIndex < pepWithSetMods.Length; pepSeqIndex++)
                            {
                                int indexInProtein = pepWithSetMods.OneBasedStartResidue + pepSeqIndex;
                                if (!ModFits(mod, pepWithSetMods.Parent, pepSeqIndex + 1, pepWithSetMods.Length, indexInProtein))
                                    continue;

                                var newPep = pepWithSetMods.Localize(pepSeqIndex, mod.MonoisotopicMass.Value);
                                peptideTheorProducts.Clear();
                                newPep.Fragment(dissociationType, CommonParameters.DigestionParams.FragmentationTerminus, peptideTheorProducts);

                                ms2ScanWithSpecificMass ??= new Ms2ScanWithSpecificMass(scan, precursorMass, precursorCharge, fileName, CommonParameters);
                                var matchedIons = MatchFragmentIons(ms2ScanWithSpecificMass, peptideTheorProducts, CommonParameters, matchAllCharges: false);
                                double score = CalculatePeptideScore(scan, matchedIons);

                                int modSiteInProteinIndex = pepWithSetMods.OneBasedStartResidue + pepSeqIndex + 1;
                                int modSiteInPeptideIndex = pepSeqIndex + 2; // plus 2 is to translate from zero based string array index to OneBasedModification index
                                if (!Filters.All(f => f.Passes(newPep, psm, score, originalScore, matchedIons, modSiteInPeptideIndex, pepWithSetMods.Length, mod)))
                                    continue;

                                if (score < bestScore - ScoreTolerance)
                                    continue;

                                // resolve variant protein location
                                string accession;
                                int adjustedSite = modSiteInProteinIndex;
                                if (!isVariantProtein)
                                {
                                    accession = pepWithSetMods.Parent.Accession;
                                }
                                else
                                {
                                    accession = null;
                                    int offset = 0;
                                    foreach (var variant in pepWithSetMods.Parent.AppliedSequenceVariations.OrderBy(v => v.OneBasedBeginPosition))
                                    {
                                        bool modIsBeforeVariant = indexInProtein < variant.OneBasedBeginPosition + offset;
                                        bool modIsOnVariant = variant.OneBasedBeginPosition + offset <= indexInProtein &&
                                                              indexInProtein <= variant.OneBasedEndPosition + offset;

                                        if (modIsOnVariant)
                                        {
                                            accession = pepWithSetMods.Parent.Accession;
                                            break;
                                        }

                                        if (modIsBeforeVariant)
                                        {
                                            accession = pepWithSetMods.Parent.ConsensusVariant.Accession;
                                            adjustedSite = indexInProtein - offset;
                                            break;
                                        }

                                        offset += variant.VariantSequence.Length - variant.OriginalSequence.Length;
                                    }

                                    if (accession == null)
                                    {
                                        accession = pepWithSetMods.Parent.ConsensusVariant.Accession;
                                        adjustedSite = indexInProtein - offset;
                                    }
                                }

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
                    }
                }
            });

            //UpdateModDictionary(modDict);
            foreach(var kvp in modDict)
                ModDictionary.MergeOrCreate(kvp.Key, new HashSet<Tuple<int, Modification>>(kvp.Value));
            return new GptmdResults(this, ModDictionary, modsAdded);
        }

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

        private static IEnumerable<Modification> GetPossibleMods(double totalMassToGetTo, IEnumerable<Modification> allMods, IEnumerable<Tuple<double, double>> combos, Tolerance precursorTolerance, IBioPolymerWithSetMods peptideWithSetModifications)
        {
            foreach (var Mod in allMods.Where(b => b.ValidModification == true))
            {
                //TODO: not necessarily here. I think we're creating ambiguity. If we're going to add a gptmd mod to a peptide that already has that mod, then we need info
                // to suggest that it is at a postion other than that in the database. could be presence of frag for unmodified or presence of frag with modified at alternative location.
                if (precursorTolerance.Within(totalMassToGetTo, peptideWithSetModifications.MonoisotopicMass + (double)Mod.MonoisotopicMass))
                    yield return Mod;
                foreach (var modOnPsm in peptideWithSetModifications.AllModsOneIsNterminus.Values.Where(b => b.ValidModification == true))
                    if (modOnPsm.Target.Equals(Mod.Target))
                    {
                        if (precursorTolerance.Within(totalMassToGetTo, peptideWithSetModifications.MonoisotopicMass + (double)Mod.MonoisotopicMass - (double)modOnPsm.MonoisotopicMass))
                            yield return Mod;
                    }
            }

            foreach (var combo in combos)
            {
                var m1 = combo.Item1;
                var m2 = combo.Item2;
                var combined = m1 + m2;
                if (precursorTolerance.Within(totalMassToGetTo, peptideWithSetModifications.MonoisotopicMass + combined))
                {
                    foreach (var mod in GetPossibleMods(totalMassToGetTo - m1, allMods, combos, precursorTolerance, peptideWithSetModifications))
                        yield return mod;
                }
            }
        }
    }
}