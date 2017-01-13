using OldInternalLogic;
using System;
using System.Collections.Generic;
using System.Linq;

namespace InternalLogicEngineLayer
{
    public class GPTMDEngine : MyEngine
    {
        #region Private Fields

        private List<NewPsmWithFDR>[] allResultingIdentifications;
        private IEnumerable<Tuple<double, double>> combos;
        private Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>> dict;
        private List<MorpheusModification> gptmdModifications;
        private bool isotopeErrors;
        private List<MorpheusModification> localizeableModifications;
        private double tol;
        private List<MorpheusModification> variableModifications;

        #endregion Private Fields

        #region Public Constructors

        public GPTMDEngine(List<NewPsmWithFDR>[] allResultingIdentifications, Dictionary<CompactPeptide, HashSet<PeptideWithSetModifications>> dict, List<MorpheusModification> variableModifications, List<MorpheusModification> localizeableModifications, bool isotopeErrors, List<MorpheusModification> gptmdModifications, IEnumerable<Tuple<double, double>> combos, double tol) : base(2)
        {
            this.dict = dict;
            this.allResultingIdentifications = allResultingIdentifications;
            this.variableModifications = variableModifications;
            this.localizeableModifications = localizeableModifications;
            this.isotopeErrors = isotopeErrors;
            this.gptmdModifications = gptmdModifications;
            this.combos = combos;
            this.tol = tol;
        }

        #endregion Public Constructors

        #region Public Methods

        protected override void ValidateParams()
        {
            if (allResultingIdentifications == null)
                throw new EngineValidationException("allResultingIdentifications cannot be null");
        }

        #endregion Public Methods

        #region Protected Methods

        protected override MyResults RunSpecific()
        {
            Dictionary<string, HashSet<Tuple<int, string>>> Mods = new Dictionary<string, HashSet<Tuple<int, string>>>();

            int modsAdded = 0;
            foreach (var ye in allResultingIdentifications[0].Where(b => b.QValue <= 0.01 && !b.isDecoy))
            {
                var baseSequence = ye.thisPSM.BaseSequence;
                double massDiff = ye.thisPSM.scanPrecursorMass - ye.thisPSM.PeptideMonoisotopicMass;
                foreach (MorpheusModification mod in GetMod(massDiff, isotopeErrors, gptmdModifications, combos, tol))
                {
                    foreach (var peptide in dict[ye.thisPSM.newPsm.GetCompactPeptide(variableModifications, localizeableModifications)])
                    {
                        int proteinLength = peptide.protein.Length;
                        var proteinAcession = peptide.protein.Accession;
                        for (int i = 0; i < baseSequence.Length; i++)
                        {
                            int indexInProtein = peptide.OneBasedStartResidueInProtein + i;

                            if (ModFits(mod, baseSequence[i], i > 0 ? baseSequence[i - 1] : peptide.PreviousAminoAcid, i + 1, baseSequence.Length, indexInProtein, proteinLength))
                            {
                                if (!Mods.ContainsKey(proteinAcession))
                                    Mods[proteinAcession] = new HashSet<Tuple<int, string>>();
                                var theTuple = new Tuple<int, string>(indexInProtein, mod.NameInXML);
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
            return new GPTMDResults(this, Mods, modsAdded);
        }

        #endregion Protected Methods

        #region Private Methods

        private static bool ModFits(MorpheusModification attemptToLocalize, char v1, char prevAA, int peptideIndex, int peptideLength, int proteinIndex, int proteinLength)
        {
            if (!attemptToLocalize.AminoAcid.Equals('\0') && !attemptToLocalize.AminoAcid.Equals(v1))
                return false;
            if (!attemptToLocalize.PrevAminoAcid.Equals('\0') && !attemptToLocalize.PrevAminoAcid.Equals(prevAA))
                return false;
            if (attemptToLocalize.Type == ModificationType.ProteinNTerminus &&
                ((proteinIndex > 2) || (proteinIndex == 2 && prevAA != 'M')))
                return false;
            if (attemptToLocalize.Type == ModificationType.PeptideNTerminus && peptideIndex > 1)
                return false;
            if (attemptToLocalize.Type == ModificationType.PeptideCTerminus && peptideIndex < peptideLength)
                return false;
            if (attemptToLocalize.Type == ModificationType.ProteinCTerminus && proteinIndex < proteinLength)
                return false;
            return true;
        }

        private static IEnumerable<MorpheusModification> GetMod(double massDiff, bool isotopeErrors, IEnumerable<MorpheusModification> allMods, IEnumerable<Tuple<double, double>> combos, double tol)
        {
            foreach (var Mod in allMods)
            {
                if (Mod.MonoisotopicMassShift > massDiff - tol && Mod.MonoisotopicMassShift < massDiff + tol)
                    yield return Mod;
                if (isotopeErrors && Mod.MonoisotopicMassShift > massDiff - tol - 1.003 && Mod.MonoisotopicMassShift < massDiff + tol - 1.003)
                    yield return Mod;
                if (!double.IsNaN(Mod.AlternativeMassShift) && Mod.AlternativeMassShift > massDiff - tol && Mod.AlternativeMassShift < massDiff + tol)
                    yield return Mod;
                if (!double.IsNaN(Mod.AlternativeMassShift) && isotopeErrors && Mod.AlternativeMassShift > massDiff - tol - 1.003 && Mod.AlternativeMassShift < massDiff + tol - 1.003)
                    yield return Mod;
            }

            foreach (var combo in combos)
            {
                var m1 = combo.Item1;
                var m2 = combo.Item2;
                var combined = m1 + m2;
                if (combined > massDiff - tol && combined < massDiff + tol)
                {
                    foreach (var mod in GetMod(m1, isotopeErrors, allMods, combos, tol))
                        yield return mod;
                    foreach (var mod in GetMod(m2, isotopeErrors, allMods, combos, tol))
                        yield return mod;
                }
                if (isotopeErrors && combined > massDiff - tol - 1.003 && combined < massDiff + tol - 1.003)
                {
                    foreach (var mod in GetMod(m1, isotopeErrors, allMods, combos, tol))
                        yield return mod;
                    foreach (var mod in GetMod(m2, isotopeErrors, allMods, combos, tol))
                        yield return mod;
                }
            }
        }

        #endregion Private Methods
    }
}