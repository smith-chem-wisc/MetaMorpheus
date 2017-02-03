using Spectra;
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
        private readonly List<MetaMorpheusModification> gptmdModifications;
        private readonly bool isotopeErrors;
        private readonly Tolerance precursorMassTolerance;

        #endregion Private Fields

        #region Public Constructors

        public GptmdEngine(List<NewPsmWithFdr> allIdentifications, bool isotopeErrors, List<MetaMorpheusModification> gptmdModifications, IEnumerable<Tuple<double, double>> combos, Tolerance precursorMassTolerance) : base(2)
        {
            this.allIdentifications = allIdentifications;
            this.isotopeErrors = isotopeErrors;
            this.gptmdModifications = gptmdModifications;
            this.combos = combos;
            this.precursorMassTolerance = precursorMassTolerance;
        }

        #endregion Public Constructors

        #region Protected Methods

        protected override MyResults RunSpecific()
        {
            var Mods = new Dictionary<string, HashSet<Tuple<int, string, string>>>();

            int modsAdded = 0;
            foreach (var ye in allIdentifications.Where(b => b.qValueNotch <= 0.01 && !b.IsDecoy))
            {
                var theDict = ye.thisPSM.peptidesWithSetModifications;
                // Only add to non-ambiguous peptides
                if (theDict.Count == 1)
                {
                    var peptide = theDict.First();
                    var baseSequence = ye.thisPSM.BaseSequence;
                    foreach (MetaMorpheusModification mod in GetMod(ye.thisPSM.ScanPrecursorMass, ye.thisPSM.PeptideMonoisotopicMass, isotopeErrors, gptmdModifications, combos, precursorMassTolerance))
                    {
                        int proteinLength = peptide.Protein.Length;
                        var proteinAcession = peptide.Protein.Accession;
                        for (int i = 0; i < baseSequence.Length; i++)
                        {
                            int indexInProtein = peptide.OneBasedStartResidueInProtein + i;

                            if (ModFits(mod, baseSequence[i], i > 0 ? baseSequence[i - 1] : peptide.PreviousAminoAcid, i + 1, baseSequence.Length, indexInProtein, proteinLength))
                            {
                                if (!Mods.ContainsKey(proteinAcession))
                                    Mods[proteinAcession] = new HashSet<Tuple<int, string, string>>();
                                var theTuple = new Tuple<int, string, string>(indexInProtein, mod.NameInXml, mod.Database);
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

        private static bool ModFits(MetaMorpheusModification attemptToLocalize, char v1, char prevAA, int peptideIndex, int peptideLength, int proteinIndex, int proteinLength)
        {
            if (!attemptToLocalize.AminoAcid.Equals('\0') && !attemptToLocalize.AminoAcid.Equals(v1))
                return false;
            if (!attemptToLocalize.PrevAminoAcid.Equals('\0') && !attemptToLocalize.PrevAminoAcid.Equals(prevAA))
                return false;
            if (attemptToLocalize.ThisModificationType == ModificationType.ProteinNTerminus &&
                ((proteinIndex > 2) || (proteinIndex == 2 && prevAA != 'M')))
                return false;
            if (attemptToLocalize.ThisModificationType == ModificationType.PeptideNTerminus && peptideIndex > 1)
                return false;
            if (attemptToLocalize.ThisModificationType == ModificationType.PeptideCTerminus && peptideIndex < peptideLength)
                return false;
            if (attemptToLocalize.ThisModificationType == ModificationType.ProteinCTerminus && proteinIndex < proteinLength)
                return false;
            return true;
        }

        private static IEnumerable<MetaMorpheusModification> GetMod(double scanPrecursorMass, double peptideMonoisotopicMass, bool isotopeErrors, IEnumerable<MetaMorpheusModification> allMods, IEnumerable<Tuple<double, double>> combos, Tolerance precursorTolerance)
        {
            foreach (var Mod in allMods)
            {
                if (precursorTolerance.Within(scanPrecursorMass, peptideMonoisotopicMass + Mod.ObservedMassShift))
                    yield return Mod;
                if (isotopeErrors && precursorTolerance.Within(scanPrecursorMass - missedMonoisopePeak, peptideMonoisotopicMass + Mod.ObservedMassShift))
                    yield return Mod;
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