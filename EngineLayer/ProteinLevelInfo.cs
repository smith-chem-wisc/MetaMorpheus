using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Linq;

namespace EngineLayer
{
    public class ProteinLevelInfo
    {

        #region Public Constructors

        public ProteinLevelInfo(HashSet<PeptideWithSetModifications> hashSet, Tolerance fragmentTolerance, Ms2ScanWithSpecificMass theScan, List<ProductType> lp)
        {
            PeptidesWithSetModifications = hashSet;
            IsDecoy = PeptidesWithSetModifications.Any(bb => bb.Protein.IsDecoy);
            IsContaminant = PeptidesWithSetModifications.Any(bb => bb.Protein.IsContaminant);
            var representative = PeptidesWithSetModifications.First();
            var MatchedIonDictPositiveIsMatch = new Dictionary<ProductType, double[]>();
            foreach (var huh in lp)
            {
                var df = representative.ProductMassesMightHaveDuplicatesAndNaNs(new List<ProductType> { huh });
                Array.Sort(df);
                double[] matchedIonMassesListPositiveIsMatch = new double[df.Length];
                PsmParent.MatchIons(theScan.TheScan, fragmentTolerance, df, matchedIonMassesListPositiveIsMatch);
                MatchedIonDictPositiveIsMatch.Add(huh, matchedIonMassesListPositiveIsMatch);
            }

            var localizedScores = new List<double>();
            for (int indexToLocalize = 0; indexToLocalize < representative.Length; indexToLocalize++)
            {
                PeptideWithSetModifications localizedPeptide = representative.Localize(indexToLocalize, theScan.PrecursorMass - representative.MonoisotopicMass);

                var gg = localizedPeptide.ProductMassesMightHaveDuplicatesAndNaNs(lp);
                Array.Sort(gg);
                double[] matchedIonMassesListPositiveIsMatch = new double[gg.Length];
                var score = PsmParent.MatchIons(theScan.TheScan, fragmentTolerance, gg, matchedIonMassesListPositiveIsMatch);
                localizedScores.Add(score);
            }

            MatchedIonsListPositiveIsMatch = MatchedIonDictPositiveIsMatch;
            LocalizedScores = localizedScores;
            PeptideMonoisotopicMass = representative.MonoisotopicMass;
            FullSequence = representative.Sequence;
            BaseSequence = representative.BaseSequence;
            MissedCleavages = representative.MissedCleavages;
            NumVariableMods = representative.NumMods - representative.numFixedMods;
            SequenceWithChemicalFormulas = representative.SequenceWithChemicalFormulas;
        }

        #endregion Public Constructors

        #region Public Properties

        public HashSet<PeptideWithSetModifications> PeptidesWithSetModifications { get; }
        public Dictionary<ProductType, double[]> MatchedIonsListPositiveIsMatch { get; }
        public List<double> LocalizedScores { get; }
        public string FullSequence { get; }
        public string BaseSequence { get; }
        public int MissedCleavages { get; }
        public double PeptideMonoisotopicMass { get; }
        public int NumVariableMods { get; }
        public string SequenceWithChemicalFormulas { get; }
        public bool IsContaminant { get; }
        public bool IsDecoy { get; }

        #endregion Public Properties

    }
}