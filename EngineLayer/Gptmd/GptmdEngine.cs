using MzLibUtil;
using Proteomics;
using System;
using System.Collections.Generic;
using System.Linq;

namespace EngineLayer.Gptmd
{
    public class GptmdEngine : MetaMorpheusEngine
    {

        #region Private Fields

        private readonly List<Psm> allIdentifications;
        private readonly IEnumerable<Tuple<double, double>> combos;
        private readonly List<ModificationWithMass> gptmdModifications;
        private readonly Tolerance precursorMassTolerance;

        #endregion Private Fields

        #region Public Constructors

        public GptmdEngine(List<Psm> allIdentifications, List<ModificationWithMass> gptmdModifications, IEnumerable<Tuple<double, double>> combos, Tolerance precursorMassTolerance, List<string> nestedIds) : base(nestedIds)
        {
            this.allIdentifications = allIdentifications;
            this.gptmdModifications = gptmdModifications;
            this.combos = combos;
            this.precursorMassTolerance = precursorMassTolerance;
        }

        #endregion Public Constructors

        #region Public Methods

        public static bool ModFits(ModificationWithMass attemptToLocalize, Protein protein, int peptideOneBasedIndex, int peptideLength, int proteinOneBasedIndex)
        {
            var motif = attemptToLocalize.motif.Motif;
            // First find the capital letter...
            var hehe = motif.IndexOf(motif.First(b => char.IsUpper(b)));

            var proteinToMotifOffset = proteinOneBasedIndex - hehe - 1;
            var indexUp = 0;
            // Look up starting at and including the capital letter
            while (indexUp < motif.Length)
            {
                if (indexUp + proteinToMotifOffset < 0 || indexUp + proteinToMotifOffset >= protein.Length || (!char.ToUpper(motif[indexUp]).Equals('X') && !char.ToUpper(motif[indexUp]).Equals(protein.BaseSequence[indexUp + proteinToMotifOffset])))
                    return false;
                indexUp++;
            }
            if (attemptToLocalize.terminusLocalization == ModificationSites.NProt && (proteinOneBasedIndex > 2))
                return false;
            if (attemptToLocalize.terminusLocalization == ModificationSites.NPep && peptideOneBasedIndex > 1)
                return false;
            if (attemptToLocalize.terminusLocalization == ModificationSites.PepC && peptideOneBasedIndex < peptideLength)
                return false;
            if (attemptToLocalize.terminusLocalization == ModificationSites.ProtC && proteinOneBasedIndex < protein.Length)
                return false;
            return true;
        }

        #endregion Public Methods

        #region Protected Methods

        protected override MetaMorpheusEngineResults RunSpecific()
        {
            var Mods = new Dictionary<string, HashSet<Tuple<int, Modification>>>();

            int modsAdded = 0;
            // Look at all confident identifications (with notch q value less than 0.01)
            // Of those only targets (do not add modifications for decoy peptides)
            foreach (var ye in allIdentifications.Where(b => b.FdrInfo.QValueNotch <= 0.01 && !b.IsDecoy))
            {
                var baseSequence = ye.MostProbableProteinInfo.BaseSequence;
                foreach (var peptide in ye.MostProbableProteinInfo.PeptidesWithSetModifications)
                    foreach (ModificationWithMass mod in GetPossibleMods(ye.ScanPrecursorMass, gptmdModifications, combos, precursorMassTolerance, peptide))
                    {
                        var proteinAcession = peptide.Protein.Accession;
                        for (int i = 0; i < baseSequence.Length; i++)
                        {
                            int indexInProtein = peptide.OneBasedStartResidueInProtein + i;
                            if (ModFits(mod, peptide.Protein, i + 1, baseSequence.Length, indexInProtein))
                            {
                                if (!Mods.ContainsKey(proteinAcession))
                                    Mods[proteinAcession] = new HashSet<Tuple<int, Modification>>();
                                var theTuple = new Tuple<int, Modification>(indexInProtein, mod);
                                if (!Mods[proteinAcession].Contains(theTuple))
                                {
                                    Mods[proteinAcession].Add(theTuple);
                                    modsAdded++;
                                }
                            }
                        }
                    }
            }
            return new GptmdResults(this, Mods, modsAdded);
        }

        #endregion Protected Methods

        #region Private Methods

        private static IEnumerable<ModificationWithMass> GetPossibleMods(double totalMassToGetTo, IEnumerable<ModificationWithMass> allMods, IEnumerable<Tuple<double, double>> combos, Tolerance precursorTolerance, PeptideWithSetModifications peptideWithSetModifications)
        {
            foreach (var Mod in allMods)
            {
                if (precursorTolerance.Within(totalMassToGetTo, peptideWithSetModifications.MonoisotopicMass + Mod.monoisotopicMass))
                    yield return Mod;
                foreach (var modOnPsm in peptideWithSetModifications.allModsOneIsNterminus.Values)
                    if (modOnPsm.motif.Motif.Equals(Mod.motif.Motif))
                    {
                        if (precursorTolerance.Within(totalMassToGetTo, peptideWithSetModifications.MonoisotopicMass + Mod.monoisotopicMass - modOnPsm.monoisotopicMass))
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
                    foreach (var mod in GetPossibleMods(totalMassToGetTo - m2, allMods, combos, precursorTolerance, peptideWithSetModifications))
                        yield return mod;
                }
            }
        }

        #endregion Private Methods

    }
}