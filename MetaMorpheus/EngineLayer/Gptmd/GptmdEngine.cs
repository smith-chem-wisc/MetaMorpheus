using MzLibUtil;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Linq;
using Omics;
using Omics.Modifications;
using System.Threading.Tasks;
using System.Collections.Concurrent;

namespace EngineLayer.Gptmd
{
    public class GptmdEngine : MetaMorpheusEngine
    {
        private readonly List<SpectralMatch> AllIdentifications;
        private readonly IEnumerable<Tuple<double, double>> Combos;
        private readonly List<Modification> GptmdModifications;
        private readonly Dictionary<string, Tolerance> FilePathToPrecursorMassTolerance; // this exists because of file-specific tolerances

        public GptmdEngine(List<SpectralMatch> allIdentifications, List<Modification> gptmdModifications, IEnumerable<Tuple<double, double>> combos, Dictionary<string, Tolerance> filePathToPrecursorMassTolerance, CommonParameters commonParameters, List<(string fileName, CommonParameters fileSpecificParameters)> fileSpecificParameters, List<string> nestedIds) : base(commonParameters, fileSpecificParameters, nestedIds)
        {
            AllIdentifications = allIdentifications;
            GptmdModifications = gptmdModifications;
            Combos = combos;
            FilePathToPrecursorMassTolerance = filePathToPrecursorMassTolerance;
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
            if (attemptToLocalize.LocationRestriction == "N-terminal." && (proteinOneBasedIndex <= 2))
                return true;
            if (attemptToLocalize.LocationRestriction == "Peptide N-terminal." && peptideOneBasedIndex == 1)
                return true;
            if (attemptToLocalize.LocationRestriction == "Peptide C-terminal." && peptideOneBasedIndex == peptideLength)
                return true;
            if (attemptToLocalize.LocationRestriction == "C-terminal." && proteinOneBasedIndex == protein.Length)
                return true;
            return false;
        }

        protected override MetaMorpheusEngineResults RunSpecific()
        {
            var modDict = new ConcurrentDictionary<string, ConcurrentBag<Tuple<int, Modification>>>();

            int modsAdded = 0;
            //foreach peptide in each psm and for each modification that matches the notch,
            //add that modification to every allowed residue

            int maxThreadsPerFile = CommonParameters.MaxThreadsToUsePerFile;
            int[] threads = Enumerable.Range(0, maxThreadsPerFile).ToArray();
            Parallel.ForEach(threads, (j) => 
            {
                var ids = AllIdentifications.Where(b => b.FdrInfo.QValueNotch <= 0.05 && !b.IsDecoy).ToList();
                for (; j < ids.Count; j += maxThreadsPerFile)
                {
                    // get file-specific precursor tolerance
                    Tolerance precursorMassTolerance = FilePathToPrecursorMassTolerance[ids[j].FullFilePath];

                    // get mods to annotate database with
                    foreach (var pepWithSetMods in ids[j].BestMatchingBioPolymersWithSetMods.Select(v => v.Peptide as PeptideWithSetModifications))
                    {
                        foreach (Modification mod in GetPossibleMods(ids[j].ScanPrecursorMass, GptmdModifications, Combos, precursorMassTolerance, pepWithSetMods))
                        {
                            var isVariantProtein = pepWithSetMods.Parent != pepWithSetMods.Protein.NonVariantProtein;

                            for (int i = 0; i < pepWithSetMods.Length; i++)
                            {
                                int indexInProtein = pepWithSetMods.OneBasedStartResidue + i;

                                if (ModFits(mod, pepWithSetMods.Parent, i + 1, pepWithSetMods.Length, indexInProtein))
                                {
                                    // if not a variant protein, index to base protein sequence
                                    if (!isVariantProtein)
                                    {
                                        AddIndexedMod(modDict, pepWithSetMods.Protein.Accession, new Tuple<int, Modification>(indexInProtein, mod));
                                        modsAdded++;
                                    }

                                    // if a variant protein, index to variant protein if on variant, or to the original protein if not
                                    else
                                    {
                                        bool foundSite = false;
                                        int offset = 0;
                                        foreach (var variant in pepWithSetMods.Protein.AppliedSequenceVariations.OrderBy(v => v.OneBasedBeginPosition))
                                        {
                                            bool modIsBeforeVariant = indexInProtein < variant.OneBasedBeginPosition + offset;
                                            bool modIsOnVariant = variant.OneBasedBeginPosition + offset <= indexInProtein && indexInProtein <= variant.OneBasedEndPosition + offset;

                                            // if a variant protein and the mod is on the variant, index to the variant protein sequence
                                            if (modIsOnVariant)
                                            {
                                                AddIndexedMod(modDict, pepWithSetMods.Protein.Accession, new Tuple<int, Modification>(indexInProtein, mod));
                                                modsAdded++;
                                                foundSite = true;
                                                break;
                                            }

                                            // otherwise back calculate the index to the original protein sequence
                                            if (modIsBeforeVariant)
                                            {
                                                AddIndexedMod(modDict, pepWithSetMods.Protein.NonVariantProtein.Accession, new Tuple<int, Modification>(indexInProtein - offset, mod));
                                                modsAdded++;
                                                foundSite = true;
                                                break;
                                            }

                                            offset += variant.VariantSequence.Length - variant.OriginalSequence.Length;
                                        }
                                        if (!foundSite)
                                        {
                                            AddIndexedMod(modDict, pepWithSetMods.Protein.NonVariantProtein.Accession, new Tuple<int, Modification>(indexInProtein - offset, mod));
                                            modsAdded++;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            });
            // Convert ConcurrentDictionary to Dictionary with HashSet
            var finalModDict = modDict.ToDictionary(
                kvp => kvp.Key,
                kvp => new HashSet<Tuple<int, Modification>>(kvp.Value)
            );
            return new GptmdResults(this, finalModDict, modsAdded);
        }
        private List<double> CalculatePeptideScores(List<PeptideWithSetModifications> newPeptides, DissociationType dissociationType, SpectralMatch psm)
        {
            var scores = new List<double>();

            foreach (var peptide in newPeptides)
            {
                var peptideTheorProducts = new List<Product>();
                peptide.Fragment(dissociationType, CommonParameters.DigestionParams.FragmentationTerminus, peptideTheorProducts);

                var scan = psm.MsDataScan;
                var precursorMass = psm.ScanPrecursorMass;
                var precursorCharge = psm.ScanPrecursorCharge;
                var fileName = psm.FullFilePath;
                List<MatchedFragmentIon> matchedIons = MatchFragmentIons(new Ms2ScanWithSpecificMass(scan, precursorMass, precursorCharge, fileName, CommonParameters), peptideTheorProducts, CommonParameters, matchAllCharges: false);

                scores.Add(CalculatePeptideScore(psm.MsDataScan, matchedIons, false));
            }

            return scores;
        }
        private static void AddIndexedMod(ConcurrentDictionary<string, ConcurrentBag<Tuple<int, Modification>>> modDict, string proteinAccession, Tuple<int, Modification> indexedMod)
        {
            modDict.AddOrUpdate(proteinAccession, new ConcurrentBag<Tuple<int, Modification>> { indexedMod }, (key, existingBag) =>
            {
                existingBag.Add(indexedMod);
                return existingBag;
            });
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

        private static IEnumerable<Modification> GetPossibleMods(double totalMassToGetTo, IEnumerable<Modification> allMods, IEnumerable<Tuple<double, double>> combos, Tolerance precursorTolerance, PeptideWithSetModifications peptideWithSetModifications)
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