using MassSpectrometry;
using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Linq;

namespace EngineLayer.Localization
{
    public class LocalizationEngine : MetaMorpheusEngine
    {
        #region Private Fields

        private readonly IEnumerable<PeptideSpectralMatch> allResultingIdentifications;
        private readonly List<ProductType> productTypes;
        private readonly MsDataFile myMsDataFile;
        private readonly Tolerance fragmentTolerance;
        private readonly bool addCompIons;
        private readonly List<DissociationType> dissociationTypes;

        #endregion Private Fields

        #region Public Constructors

        public LocalizationEngine(IEnumerable<PeptideSpectralMatch> allResultingIdentifications, List<ProductType> lp, MsDataFile myMsDataFile, Tolerance fragmentTolerance, List<string> nestedIds, bool addCompIons) : base(nestedIds)
        {
            this.allResultingIdentifications = allResultingIdentifications;
            this.productTypes = lp;
            this.myMsDataFile = myMsDataFile;
            this.fragmentTolerance = fragmentTolerance;
            this.addCompIons = addCompIons;
            this.dissociationTypes = DetermineDissociationType(lp);
        }

        #endregion Public Constructors

        #region Protected Methods

        protected override MetaMorpheusEngineResults RunSpecific()
        {
            TerminusType terminusType = ProductTypeMethod.IdentifyTerminusType(productTypes);

            foreach (PeptideSpectralMatch psm in allResultingIdentifications)
            {
                psm.MatchedIonSeriesDict = new Dictionary<ProductType, int[]>();
                psm.MatchedIonMassToChargeRatioDict = new Dictionary<ProductType, double[]>();
                psm.ProductMassErrorDa = new Dictionary<ProductType, double[]>();
                psm.ProductMassErrorPpm = new Dictionary<ProductType, double[]>();
                psm.MatchedIonIntensitiesDict = new Dictionary<ProductType, double[]>();
                var theScan = myMsDataFile.GetOneBasedScan(psm.ScanNumber);
                double thePrecursorMass = psm.ScanPrecursorMass;
                foreach (ProductType productType in productTypes)
                {
                    var sortedTheoreticalProductMasses = psm.CompactPeptides.First().Key.ProductMassesMightHaveDuplicatesAndNaNs(new List<ProductType> { productType });
                    Array.Sort(sortedTheoreticalProductMasses);
                    List<int> matchedIonSeriesList = new List<int>();
                    List<double> matchedIonMassToChargeRatioList = new List<double>();
                    List<double> productMassErrorDaList = new List<double>();
                    List<double> productMassErrorPpmList = new List<double>();
                    List<double> matchedIonIntensityList = new List<double>();

                    //populate the above lists
                    MatchIons(theScan, fragmentTolerance, sortedTheoreticalProductMasses, matchedIonSeriesList, matchedIonMassToChargeRatioList, productMassErrorDaList, productMassErrorPpmList, matchedIonIntensityList, thePrecursorMass, productType, addCompIons);

                    psm.MatchedIonSeriesDict.Add(productType, matchedIonSeriesList.ToArray());
                    psm.MatchedIonMassToChargeRatioDict.Add(productType, matchedIonMassToChargeRatioList.ToArray());
                    psm.ProductMassErrorDa.Add(productType, productMassErrorDaList.ToArray());
                    psm.ProductMassErrorPpm.Add(productType, productMassErrorPpmList.ToArray());
                    psm.MatchedIonIntensitiesDict.Add(productType, matchedIonIntensityList.ToArray());
                }
            }

            foreach (PeptideSpectralMatch psm in allResultingIdentifications.Where(b => b.NumDifferentCompactPeptides == 1))
            {
                var theScan = myMsDataFile.GetOneBasedScan(psm.ScanNumber);
                double thePrecursorMass = psm.ScanPrecursorMass;

                if (psm.FullSequence == null)
                {
                    continue;
                }

                PeptideWithSetModifications representative = psm.CompactPeptides.First().Value.Item2.First();

                var localizedScores = new List<double>();
                for (int indexToLocalize = 0; indexToLocalize < representative.Length; indexToLocalize++)
                {
                    PeptideWithSetModifications localizedPeptide = representative.Localize(indexToLocalize, psm.ScanPrecursorMass - representative.MonoisotopicMass);

                    var gg = localizedPeptide.CompactPeptide(terminusType).ProductMassesMightHaveDuplicatesAndNaNs(productTypes);
                    Array.Sort(gg);
                    var score = CalculatePeptideScore(theScan, fragmentTolerance, gg, thePrecursorMass, dissociationTypes, addCompIons, 0);
                    localizedScores.Add(score);
                }

                psm.LocalizedScores = localizedScores;
            }
            return new LocalizationEngineResults(this);
        }

        #endregion Protected Methods
    }
}