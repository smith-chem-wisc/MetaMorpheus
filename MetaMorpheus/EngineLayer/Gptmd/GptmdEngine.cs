using MzLibUtil;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Linq;
using Omics;
using Omics.Modifications;
using System.Threading.Tasks;
using System.Collections.Concurrent;
using MassSpectrometry;
using Omics.Fragmentation;

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

        public GptmdEngine(
            List<SpectralMatch> allIdentifications, 
            List<Modification> gptmdModifications, 
            IEnumerable<Tuple<double, double>> combos, 
            Dictionary<string, Tolerance> filePathToPrecursorMassTolerance, 
            CommonParameters commonParameters, 
            List<(string fileName, CommonParameters fileSpecificParameters)> fileSpecificParameters, 
            List<string> nestedIds) 
            : base(commonParameters, fileSpecificParameters, nestedIds)
        {
            AllIdentifications = allIdentifications;
            GptmdModifications = gptmdModifications;
            Combos = combos;
            FilePathToPrecursorMassTolerance = filePathToPrecursorMassTolerance;
            ModDictionary = new Dictionary<string, HashSet<Tuple<int, Modification>>>();
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
                    foreach (var pepWithSetMods in psms[i].BestMatchingBioPolymersWithSetMods.Select(v => v.SpecificBioPolymer))
                    {
                        var isVariantProtein = pepWithSetMods.Parent != pepWithSetMods.Parent.ConsensusVariant;
                        var possibleModifications = GetPossibleMods(psms[i].ScanPrecursorMass, GptmdModifications, Combos, FilePathToPrecursorMassTolerance[psms[i].FullFilePath], pepWithSetMods);

                        if (!isVariantProtein)
                        {
                            foreach (var mod in possibleModifications)
                            {
                                List<int> possibleIndices = Enumerable.Range(0, pepWithSetMods.Length).Where(i => ModFits(mod, pepWithSetMods.Parent, i + 1, pepWithSetMods.Length, pepWithSetMods.OneBasedStartResidue + i)).ToList();
                                if (possibleIndices.Any())
                                {
                                    List<IBioPolymerWithSetMods> newPeptides = new();
                                    foreach (int index in possibleIndices)
                                    {
                                        if (mod.MonoisotopicMass.HasValue)
                                        {
                                            newPeptides.Add(pepWithSetMods.Localize(index, mod.MonoisotopicMass.Value));
                                        }
                                    }

                                    if (newPeptides.Any())
                                    {
                                        var scores = new List<double>();
                                        var dissociationType = CommonParameters.DissociationType == DissociationType.Autodetect ?
                                            psms[i].Ms2Scan.DissociationType.Value : CommonParameters.DissociationType;

                                        scores = CalculatePeptideScores(newPeptides, dissociationType, psms[i]);

                                        // If the score is within tolerance of the highest score, add the mod to the peptide
                                        // If the tolerance is too tight, then the number of identifications in subsequent searches will be reduced
                                            
                                        var highScoreIndices = scores.Select((item, index) => new { item, index })
                                            .Where(x => x.item > (scores.Max() - ScoreTolerance))
                                            .Select(x => x.index)
                                            .ToList();

                                        foreach (var index in highScoreIndices)
                                        {
                                            AddIndexedMod(modDict, pepWithSetMods.Parent.Accession, new Tuple<int, Modification>(pepWithSetMods.OneBasedStartResidue + possibleIndices[index], mod));
                                            System.Threading.Interlocked.Increment(ref modsAdded); ;
                                        }
                                    }
                                }
                            }
                        }
                        // if a variant protein, index to variant protein if on variant, or to the original protein if not
                        else
                        {
                            foreach (var mod in possibleModifications)
                            {
                                for (int j = 0; j < pepWithSetMods.Length; j++)
                                {
                                    int indexInProtein = pepWithSetMods.OneBasedStartResidue + j;

                                    if (ModFits(mod, pepWithSetMods.Parent, j + 1, pepWithSetMods.Length, indexInProtein))
                                    {
                                        bool foundSite = false;
                                        int offset = 0;
                                        foreach (var variant in pepWithSetMods.Parent.AppliedSequenceVariations.OrderBy(v => v.OneBasedBeginPosition))
                                        {
                                            bool modIsBeforeVariant = indexInProtein < variant.OneBasedBeginPosition + offset;
                                            bool modIsOnVariant = variant.OneBasedBeginPosition + offset <= indexInProtein && indexInProtein <= variant.OneBasedEndPosition + offset;

                                            // if a variant protein and the mod is on the variant, index to the variant protein sequence
                                            if (modIsOnVariant)
                                            {
                                                AddIndexedMod(modDict, pepWithSetMods.Parent.Accession, new Tuple<int, Modification>(indexInProtein, mod));
                                                foundSite = true;
                                                System.Threading.Interlocked.Increment(ref modsAdded); ;
                                                break;
                                            }

                                            // otherwise back calculate the index to the original protein sequence
                                            if (modIsBeforeVariant)
                                            {
                                                AddIndexedMod(modDict, pepWithSetMods.Parent.ConsensusVariant.Accession, new Tuple<int, Modification>(indexInProtein - offset, mod));
                                                foundSite = true;
                                                System.Threading.Interlocked.Increment(ref modsAdded); ;
                                                break;
                                            }

                                            offset += variant.VariantSequence.Length - variant.OriginalSequence.Length;
                                        }
                                        if (!foundSite)
                                        {
                                            AddIndexedMod(modDict, pepWithSetMods.Parent.ConsensusVariant.Accession, new Tuple<int, Modification>(indexInProtein - offset, mod));
                                            System.Threading.Interlocked.Increment(ref modsAdded); ;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            });

            //UpdateModDictionary(modDict);
            foreach(var kvp in modDict)
                ModDictionary.MergeOrCreate(kvp.Key, new HashSet<Tuple<int, Modification>>(kvp.Value));
            return new GptmdResults(this, ModDictionary, modsAdded);
        }

        private List<double> CalculatePeptideScores(List<IBioPolymerWithSetMods> newPeptides, DissociationType dissociationType, SpectralMatch psm)
        {
            var scores = new List<double>();

            foreach (var peptide in newPeptides)
            {
                var peptideTheorProducts = new List<Product>();
                peptide.Fragment(dissociationType, CommonParameters.DigestionParams.FragmentationTerminus, peptideTheorProducts);

                var scan = psm.Ms2Scan;
                var precursorMass = psm.ScanPrecursorMass;
                var precursorCharge = psm.ScanPrecursorCharge;
                var fileName = psm.FullFilePath;
                List<MatchedFragmentIon> matchedIons = MatchFragmentIons(new Ms2ScanWithSpecificMass(scan, precursorMass, precursorCharge, fileName, CommonParameters), peptideTheorProducts, CommonParameters, matchAllCharges: false);

                scores.Add(CalculatePeptideScore(scan, matchedIons, false));
            }

            return scores;
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