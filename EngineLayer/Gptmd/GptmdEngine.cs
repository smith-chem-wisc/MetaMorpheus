using MzLibUtil;
using Proteomics;
using System;
using System.Collections.Generic;
using System.Linq;

namespace EngineLayer.Gptmd
{
    public class GptmdEngine : MyEngine
    {

        #region Private Fields

        private const double missedMonoisopePeak = 1.003;
        private readonly List<NewPsmWithFdr> allIdentifications;
        private readonly IEnumerable<Tuple<double, double>> combos;
        private readonly List<ModificationWithMass> gptmdModifications;
        private readonly bool isotopeErrors;
        private readonly Tolerance precursorMassTolerance;

        #endregion Private Fields

        #region Public Constructors

        public GptmdEngine(List<NewPsmWithFdr> allIdentifications, bool isotopeErrors, List<ModificationWithMass> gptmdModifications, IEnumerable<Tuple<double, double>> combos, Tolerance precursorMassTolerance)
        {
            this.allIdentifications = allIdentifications;
            this.isotopeErrors = isotopeErrors;
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
            if (attemptToLocalize.position == ModificationSites.NProt && (proteinOneBasedIndex > 2))
                return false;
            if (attemptToLocalize.position == ModificationSites.NPep && peptideOneBasedIndex > 1)
                return false;
            if (attemptToLocalize.position == ModificationSites.PepC && peptideOneBasedIndex < peptideLength)
                return false;
            if (attemptToLocalize.position == ModificationSites.ProtC && proteinOneBasedIndex < protein.Length)
                return false;
            return true;
        }

        #endregion Public Methods

        #region Protected Methods

        protected override MyResults RunSpecific()
        {
            var Mods = new Dictionary<string, HashSet<Tuple<int, ModificationWithMass>>>();

            int modsAdded = 0;
            foreach (var ye in allIdentifications.Where(b => b.qValueNotch <= 0.01 && !b.IsDecoy))
            {
                var theDict = ye.thisPSM.peptidesWithSetModifications;
                // Only add to non-ambiguous peptides
                if (theDict.Count == 1)
                {
                    var peptide = theDict.First();
                    var baseSequence = ye.thisPSM.BaseSequence;
                    foreach (ModificationWithMass mod in GetMod(ye.thisPSM.ScanPrecursorMass, ye.thisPSM.PeptideMonoisotopicMass, isotopeErrors, gptmdModifications, combos, precursorMassTolerance))
                    {
                        var proteinAcession = peptide.Protein.Accession;
                        for (int i = 0; i < baseSequence.Length; i++)
                        {
                            int indexInProtein = peptide.OneBasedStartResidueInProtein + i;
                            if (ModFits(mod, peptide.Protein, i + 1, baseSequence.Length, indexInProtein))
                            {
                                if (!Mods.ContainsKey(proteinAcession))
                                    Mods[proteinAcession] = new HashSet<Tuple<int, ModificationWithMass>>();
                                var theTuple = new Tuple<int, ModificationWithMass>(indexInProtein, mod);
                                if (!Mods[proteinAcession].Contains(theTuple))
                                {
                                    Mods[proteinAcession].Add(theTuple);
                                    modsAdded++;
                                }
                            }
                        }
                    }
                }
            }
            return new GptmdResults(this, Mods, modsAdded);
        }

        #endregion Protected Methods

        #region Private Methods

        private static IEnumerable<ModificationWithMass> GetMod(double scanPrecursorMass, double peptideMonoisotopicMass, bool isotopeErrors, IEnumerable<ModificationWithMass> allMods, IEnumerable<Tuple<double, double>> combos, Tolerance precursorTolerance)
        {
            foreach (var Mod in allMods)
            {
                foreach (var massObserved in Mod.massesObserved)
                {
                    if (precursorTolerance.Within(scanPrecursorMass, peptideMonoisotopicMass + massObserved))
                        yield return Mod;
                    if (isotopeErrors && precursorTolerance.Within(scanPrecursorMass - missedMonoisopePeak, peptideMonoisotopicMass + massObserved))
                        yield return Mod;
                }
            }

            foreach (var combo in combos)
            {
                var m1 = combo.Item1;
                var m2 = combo.Item2;
                var combined = m1 + m2;
                if (precursorTolerance.Within(scanPrecursorMass, peptideMonoisotopicMass + combined))
                {
                    foreach (var mod in GetMod(scanPrecursorMass, peptideMonoisotopicMass + m1, isotopeErrors, allMods, combos, precursorTolerance))
                        yield return mod;
                    foreach (var mod in GetMod(scanPrecursorMass, peptideMonoisotopicMass + m2, isotopeErrors, allMods, combos, precursorTolerance))
                        yield return mod;
                }
                if (isotopeErrors && precursorTolerance.Within(scanPrecursorMass - missedMonoisopePeak, peptideMonoisotopicMass + combined))
                {
                    foreach (var mod in GetMod(scanPrecursorMass, peptideMonoisotopicMass + m1, isotopeErrors, allMods, combos, precursorTolerance))
                        yield return mod;
                    foreach (var mod in GetMod(scanPrecursorMass, peptideMonoisotopicMass + m2, isotopeErrors, allMods, combos, precursorTolerance))
                        yield return mod;
                }
            }
        }

        #endregion Private Methods

    }
}