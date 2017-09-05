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

        private readonly IEnumerable<Psm> allResultingIdentifications;
        private readonly List<ProductType> lp;
        private readonly IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> myMsDataFile;
        private readonly Tolerance fragmentTolerance;
        private readonly bool addCompIons;
        private readonly List<DissociationType> dissociationTypes;

        #endregion Private Fields

        #region Public Constructors

        public LocalizationEngine(IEnumerable<Psm> allResultingIdentifications, List<ProductType> lp, IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> myMsDataFile, Tolerance fragmentTolerance, List<string> nestedIds, bool addCompIons) : base(nestedIds)
        {
            this.allResultingIdentifications = allResultingIdentifications;
            this.lp = lp;
            this.myMsDataFile = myMsDataFile;
            this.fragmentTolerance = fragmentTolerance;
            this.addCompIons = addCompIons;
            this.dissociationTypes = DetermineDissociationType(lp);
        }

        #endregion Public Constructors

        #region Protected Methods

        protected override MetaMorpheusEngineResults RunSpecific()
        {
            TerminusType terminusType = ProductTypeToTerminusType.IdentifyTerminusType(lp);

            foreach (var ok in allResultingIdentifications)
            {
                ok.MatchedIonDictOnlyMatches = new Dictionary<ProductType, double[]>();
                ok.ProductMassErrorDa = new Dictionary<ProductType, double[]>();
                ok.ProductMassErrorPpm = new Dictionary<ProductType, double[]>();
                var theScan = myMsDataFile.GetOneBasedScan(ok.ScanNumber);
                double thePrecursorMass = ok.ScanPrecursorMass;
                foreach (var huh in lp)
                {
                    var ionMasses = ok.CompactPeptides.First().Key.ProductMassesMightHaveDuplicatesAndNaNs(new List<ProductType> { huh });
                    Array.Sort(ionMasses);
                    List<double> matchedIonMassesList = new List<double>();
                    List<double> productMassErrorDaList = new List<double>();
                    List<double> productMassErrorPpmList = new List<double>();
                    MatchIons(theScan, fragmentTolerance, ionMasses, matchedIonMassesList, productMassErrorDaList, productMassErrorPpmList, thePrecursorMass, dissociationTypes, addCompIons);
                    double[] matchedIonMassesOnlyMatches = matchedIonMassesList.ToArray();
                    ok.MatchedIonDictOnlyMatches.Add(huh, matchedIonMassesOnlyMatches);
                    ok.ProductMassErrorDa.Add(huh, productMassErrorDaList.ToArray());
                    ok.ProductMassErrorPpm.Add(huh, productMassErrorPpmList.ToArray());
                }
            }

            foreach (var ok in allResultingIdentifications.Where(b => b.NumDifferentCompactPeptides == 1))
            {
                var theScan = myMsDataFile.GetOneBasedScan(ok.ScanNumber);
                double thePrecursorMass = ok.ScanPrecursorMass;

                if (ok.FullSequence == null)
                    continue;

                var representative = ok.CompactPeptides.First().Value.Item2.First();

                var localizedScores = new List<double>();
                for (int indexToLocalize = 0; indexToLocalize < representative.Length; indexToLocalize++)
                {
                    PeptideWithSetModifications localizedPeptide = representative.Localize(indexToLocalize, ok.ScanPrecursorMass - representative.MonoisotopicMass);

                    var gg = localizedPeptide.CompactPeptide(terminusType).ProductMassesMightHaveDuplicatesAndNaNs(lp);
                    Array.Sort(gg);
                    var score = CalculateClassicScore(theScan, fragmentTolerance, gg, thePrecursorMass, dissociationTypes, addCompIons);
                    localizedScores.Add(score);
                }

                ok.LocalizedScores = localizedScores;
            }
            return new LocalizationEngineResults(this);
        }

        #endregion Protected Methods
    }
}