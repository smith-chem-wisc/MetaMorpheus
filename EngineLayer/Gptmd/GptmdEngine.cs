
using System;
using System.Collections.Generic;
using System.Linq;

namespace EngineLayer.Gptmd
{
    public class GptmdEngine : MyEngine
    {

        #region Private Fields

        private readonly List<NewPsmWithFdr> allResultingIdentifications;
        private readonly IEnumerable<Tuple<double, double>> combos;
        private readonly List<MorpheusModification> gptmdModifications;
        private readonly bool isotopeErrors;
        private readonly double tolInDaltons;

        #endregion Private Fields

        #region Public Constructors

        public GptmdEngine(List<NewPsmWithFdr> allResultingIdentifications, bool isotopeErrors, List<MorpheusModification> gptmdModifications, IEnumerable<Tuple<double, double>> combos, double tolInDaltons) : base(2)
        {
            this.allResultingIdentifications = allResultingIdentifications;
            this.isotopeErrors = isotopeErrors;
            this.gptmdModifications = gptmdModifications;
            this.combos = combos;
            this.tolInDaltons = tolInDaltons;
        }

        #endregion Public Constructors

        #region Protected Methods

        protected override MyResults RunSpecific()
        {
            var Mods = new Dictionary<string, HashSet<Tuple<int, string, string>>>();

            int modsAdded = 0;
            foreach (var ye in allResultingIdentifications.Where(b => b.qValue <= 0.01 && !b.IsDecoy))
            {
                var theDict = ye.thisPSM.peptidesWithSetModifications;
                // Only add to non-ambiguous peptides
                if (theDict.Count == 1)
                {
                    var peptide = theDict.First();
                    var baseSequence = ye.thisPSM.BaseSequence;
                    double massDiff = ye.thisPSM.ScanPrecursorMass - ye.thisPSM.PeptideMonoisotopicMass;
                    foreach (MorpheusModification mod in GetMod(massDiff, isotopeErrors, gptmdModifications, combos, tolInDaltons))
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

        private static bool ModFits(MorpheusModification attemptToLocalize, char v1, char prevAA, int peptideIndex, int peptideLength, int proteinIndex, int proteinLength)
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

        private static IEnumerable<MorpheusModification> GetMod(double massDiff, bool isotopeErrors, IEnumerable<MorpheusModification> allMods, IEnumerable<Tuple<double, double>> combos, double tolInDaltons)
        {
            foreach (var Mod in allMods)
            {
                if (Mod.ObservedMassShift > massDiff - tolInDaltons && Mod.ObservedMassShift < massDiff + tolInDaltons)
                    yield return Mod;
                if (isotopeErrors && Mod.ObservedMassShift > massDiff - tolInDaltons - 1.003 && Mod.ObservedMassShift < massDiff + tolInDaltons - 1.003)
                    yield return Mod;
            }

            foreach (var combo in combos)
            {
                var m1 = combo.Item1;
                var m2 = combo.Item2;
                var combined = m1 + m2;
                if (combined > massDiff - tolInDaltons && combined < massDiff + tolInDaltons)
                {
                    foreach (var mod in GetMod(m1, isotopeErrors, allMods, combos, tolInDaltons))
                        yield return mod;
                    foreach (var mod in GetMod(m2, isotopeErrors, allMods, combos, tolInDaltons))
                        yield return mod;
                }
                if (isotopeErrors && combined > massDiff - tolInDaltons - 1.003 && combined < massDiff + tolInDaltons - 1.003)
                {
                    foreach (var mod in GetMod(m1, isotopeErrors, allMods, combos, tolInDaltons))
                        yield return mod;
                    foreach (var mod in GetMod(m2, isotopeErrors, allMods, combos, tolInDaltons))
                        yield return mod;
                }
            }
        }

        #endregion Private Methods

    }
}