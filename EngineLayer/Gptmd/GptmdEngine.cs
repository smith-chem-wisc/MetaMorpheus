using MzLibUtil;
using Proteomics;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Linq;

namespace EngineLayer.Gptmd
{
    public class GptmdEngine : MetaMorpheusEngine
    {
        private readonly List<PeptideSpectralMatch> AllIdentifications;
        private readonly IEnumerable<Tuple<double, double>> Combos;
        private readonly List<ModificationWithMass> GptmdModifications;
        private readonly Dictionary<string, Tolerance> FilePathToPrecursorMassTolerance; // this exists because of file-specific tolerances

        public GptmdEngine(List<PeptideSpectralMatch> allIdentifications, List<ModificationWithMass> gptmdModifications, IEnumerable<Tuple<double, double>> combos, Dictionary<string, Tolerance> filePathToPrecursorMassTolerance, CommonParameters commonParameters, List<string> nestedIds) : base(commonParameters, nestedIds)
        {
            AllIdentifications = allIdentifications;
            GptmdModifications = gptmdModifications;
            Combos = combos;
            FilePathToPrecursorMassTolerance = filePathToPrecursorMassTolerance;
        }

        public static bool ModFits(ModificationWithMass attemptToLocalize, Protein protein, int peptideOneBasedIndex, int peptideLength, int proteinOneBasedIndex)
        {
            var motif = attemptToLocalize.motif;
            // First find the capital letter...
            var hehe = motif.ToString().IndexOf(motif.ToString().First(b => char.IsUpper(b)));

            var proteinToMotifOffset = proteinOneBasedIndex - hehe - 1;
            var indexUp = 0;
            // Look up starting at and including the capital letter
            while (indexUp < motif.ToString().Length)
            {
                if (indexUp + proteinToMotifOffset < 0
                    || indexUp + proteinToMotifOffset >= protein.Length
                    || (!char.ToUpper(motif.ToString()[indexUp]).Equals('X')
                        && !char.ToUpper(motif.ToString()[indexUp]).Equals(protein.BaseSequence[indexUp + proteinToMotifOffset])))
                {
                    return false;
                }
                indexUp++;
            }
            if (attemptToLocalize.terminusLocalization == TerminusLocalization.NProt
                && (proteinOneBasedIndex > 2))
            {
                return false;
            }
            if (attemptToLocalize.terminusLocalization == TerminusLocalization.NPep
                && peptideOneBasedIndex > 1)
            {
                return false;
            }
            if (attemptToLocalize.terminusLocalization == TerminusLocalization.PepC
                && peptideOneBasedIndex < peptideLength)
            {
                return false;
            }
            if (attemptToLocalize.terminusLocalization == TerminusLocalization.ProtC
                && proteinOneBasedIndex < protein.Length)
            {
                return false;
            }
            return true;
        }

        protected override MetaMorpheusEngineResults RunSpecific()
        {
            var Mods = new Dictionary<string, HashSet<Tuple<int, Modification>>>();

            int modsAdded = 0;
            //foreach peptide in each psm and for each modification that matches the notch,
            //add that modification to every allowed residue
            foreach (var psm in AllIdentifications.Where(b => b.FdrInfo.QValueNotch <= 0.05 && !b.IsDecoy))
            {
                // get file-specific precursor tolerance
                Tolerance precursorMassTolerance = FilePathToPrecursorMassTolerance[psm.FullFilePath];

                // get mods to annotate database with
                foreach (var pepWithSetMods in psm.CompactPeptides.SelectMany(b => b.Value.Item2))
                {
                    foreach (ModificationWithMass mod in GetPossibleMods(psm.ScanPrecursorMass, GptmdModifications, Combos, precursorMassTolerance, pepWithSetMods))
                    {
                        var proteinAccession = pepWithSetMods.Protein.Accession;

                        for (int i = 0; i < pepWithSetMods.Length; i++)
                        {
                            int indexInProtein = pepWithSetMods.OneBasedStartResidueInProtein + i;

                            if (ModFits(mod, pepWithSetMods.Protein, i + 1, pepWithSetMods.Length, indexInProtein))
                            {
                                if (!Mods.ContainsKey(proteinAccession))
                                {
                                    Mods[proteinAccession] = new HashSet<Tuple<int, Modification>>();
                                }

                                var modWithIndexInProtein = new Tuple<int, Modification>(indexInProtein, mod);

                                if (!Mods[proteinAccession].Contains(modWithIndexInProtein))
                                {
                                    Mods[proteinAccession].Add(modWithIndexInProtein);
                                    modsAdded++;
                                }
                            }
                        }
                    }
                }
            }

            return new GptmdResults(this, Mods, modsAdded);
        }

        private static IEnumerable<ModificationWithMass> GetPossibleMods(double totalMassToGetTo, IEnumerable<ModificationWithMass> allMods, IEnumerable<Tuple<double, double>> combos, Tolerance precursorTolerance, PeptideWithSetModifications peptideWithSetModifications)
        {
            foreach (var Mod in allMods)
            {
                if (precursorTolerance.Within(totalMassToGetTo, peptideWithSetModifications.MonoisotopicMass + Mod.monoisotopicMass))
                {
                    yield return Mod;
                }
                foreach (var modOnPsm in peptideWithSetModifications.AllModsOneIsNterminus.Values)
                {
                    if (modOnPsm.motif.Equals(Mod.motif))
                    {
                        if (precursorTolerance.Within(totalMassToGetTo, peptideWithSetModifications.MonoisotopicMass + Mod.monoisotopicMass - modOnPsm.monoisotopicMass))
                        {
                            yield return Mod;
                        }
                    }
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
                    {
                        yield return mod;
                    }
                    foreach (var mod in GetPossibleMods(totalMassToGetTo - m2, allMods, combos, precursorTolerance, peptideWithSetModifications))
                    {
                        yield return mod;
                    }
                }
            }
        }
    }
}