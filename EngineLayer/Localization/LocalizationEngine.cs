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
        private readonly List<ProductType> lp;
        private readonly MsDataFile myMsDataFile;
        private readonly CommonParameters commonParameters;
        private readonly List<DissociationType> dissociationTypes;

        #endregion Private Fields

        #region Public Constructors

        public LocalizationEngine(IEnumerable<PeptideSpectralMatch> allResultingIdentifications, List<ProductType> lp, MsDataFile myMsDataFile, CommonParameters commonParameters, List<string> nestedIds) : base(nestedIds)
        {
            this.allResultingIdentifications = allResultingIdentifications;
            this.lp = lp;
            this.myMsDataFile = myMsDataFile;
            this.commonParameters = commonParameters;
            this.dissociationTypes = DetermineDissociationType(lp);
        }

        #endregion Public Constructors

        #region Protected Methods

        protected override MetaMorpheusEngineResults RunSpecific()
        {
            TerminusType terminusType = ProductTypeMethod.IdentifyTerminusType(lp);

            foreach (var ok in allResultingIdentifications)
            {
                ok.MatchedIonMassesDict = new Dictionary<ProductType, double[]>();
                ok.ProductMassErrorDa = new Dictionary<ProductType, double[]>();
                ok.ProductMassErrorPpm = new Dictionary<ProductType, double[]>();
                ok.MatchedIonIntensitiesDict = new Dictionary<ProductType, double[]>();
                var theScan = myMsDataFile.GetOneBasedScan(ok.ScanNumber);
                double thePrecursorMass = ok.ScanPrecursorMass;
                foreach (var huh in lp)
                {
                    var ionMasses = ok.CompactPeptides.First().Key.ProductMassesMightHaveDuplicatesAndNaNs(new List<ProductType> { huh });
                    Array.Sort(ionMasses);
                    List<double> matchedIonMassesList = new List<double>();
                    List<double> productMassErrorDaList = new List<double>();
                    List<double> productMassErrorPpmList = new List<double>();
                    List<double> matchedIonIntensityList = new List<double>(); 
                    MatchIonsOld(theScan, commonParameters.ProductMassTolerance, ionMasses, matchedIonMassesList, productMassErrorDaList, productMassErrorPpmList, thePrecursorMass, dissociationTypes, commonParameters.AddCompIons, matchedIonIntensityList); 
                    double[] matchedIonMassesOnlyMatches = matchedIonMassesList.ToArray();
                    ok.MatchedIonMassesDict.Add(huh, matchedIonMassesOnlyMatches);
                    ok.ProductMassErrorDa.Add(huh, productMassErrorDaList.ToArray());
                    ok.ProductMassErrorPpm.Add(huh, productMassErrorPpmList.ToArray());
                    ok.MatchedIonIntensitiesDict.Add(huh, matchedIonIntensityList.ToArray()); 
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
                    var score = CalculatePeptideScoreOld(theScan, commonParameters.ProductMassTolerance, gg, thePrecursorMass, dissociationTypes, commonParameters.AddCompIons, 0);
                    localizedScores.Add(score);
                }

                ok.LocalizedScores = localizedScores;
            }
            return new LocalizationEngineResults(this);
        }

        #endregion Protected Methods
    }
}