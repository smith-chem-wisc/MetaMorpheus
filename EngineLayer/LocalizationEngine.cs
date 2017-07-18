using MassSpectrometry;
using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Linq;

namespace EngineLayer
{
    public class LocalizationEngine : MetaMorpheusEngine
    {

        #region Private Fields

        private readonly IEnumerable<PsmParent> allResultingIdentifications;
        private readonly List<ProductType> lp;
        private readonly IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> myMsDataFile;
        private readonly Tolerance fragmentTolerance;

        #endregion Private Fields

        #region Public Constructors

        public LocalizationEngine(IEnumerable<PsmParent> allResultingIdentifications, List<ProductType> lp, IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> myMsDataFile, Tolerance fragmentTolerance, List<string> nestedIds) : base(nestedIds)
        {
            this.allResultingIdentifications = allResultingIdentifications;
            this.lp = lp;
            this.myMsDataFile = myMsDataFile;
            this.fragmentTolerance = fragmentTolerance;
        }

        #endregion Public Constructors

        #region Protected Methods

        protected override MetaMorpheusEngineResults RunSpecific()
        {
            foreach (var ok in allResultingIdentifications)
            {
                var MatchedIonDictPositiveIsMatch = new Dictionary<ProductType, double[]>();
                var representative = ok.Pli.PeptidesWithSetModifications.First();
                var theScan = myMsDataFile.GetOneBasedScan(ok.ScanNumber);
                foreach (var huh in lp)
                {
                    var ionMasses = representative.CompactPeptide.ProductMassesMightHaveDuplicatesAndNaNs(new List<ProductType> { huh });
                    Array.Sort(ionMasses);
                    double[] matchedIonMassesListPositiveIsMatch = new double[ionMasses.Length];
                    PsmParent.MatchIons(theScan, fragmentTolerance, ionMasses, matchedIonMassesListPositiveIsMatch);
                    MatchedIonDictPositiveIsMatch.Add(huh, matchedIonMassesListPositiveIsMatch);
                }

                var localizedScores = new List<double>();
                for (int indexToLocalize = 0; indexToLocalize < representative.Length; indexToLocalize++)
                {
                    PeptideWithSetModifications localizedPeptide = representative.Localize(indexToLocalize, ok.ScanPrecursorMass - representative.MonoisotopicMass);

                    var gg = localizedPeptide.CompactPeptide.ProductMassesMightHaveDuplicatesAndNaNs(lp);
                    Array.Sort(gg);
                    double[] matchedIonMassesListPositiveIsMatch = new double[gg.Length];
                    var score = PsmParent.MatchIons(theScan, fragmentTolerance, gg, matchedIonMassesListPositiveIsMatch);
                    localizedScores.Add(score);
                }

                ok.LocalizationResults = new LocalizationResults(MatchedIonDictPositiveIsMatch, localizedScores);
            }
            return new LocalizationEngineResults(this);
        }

        #endregion Protected Methods

    }
}