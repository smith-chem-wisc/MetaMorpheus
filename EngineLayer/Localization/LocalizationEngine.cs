using MassSpectrometry;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Linq;

namespace EngineLayer.Localization
{
    public class LocalizationEngine : MetaMorpheusEngine
    {
        private readonly IEnumerable<PeptideSpectralMatch> AllResultingIdentifications;
        private readonly List<ProductType> ProductTypes;
        private readonly MsDataFile MyMsDataFile;
        private readonly List<DissociationType> DissociationTypes;

        public LocalizationEngine(IEnumerable<PeptideSpectralMatch> allResultingIdentifications, List<ProductType> productTypes, MsDataFile myMsDataFile,
            CommonParameters commonParameters, List<string> nestedIds)
            : base(commonParameters, nestedIds)
        {
            AllResultingIdentifications = allResultingIdentifications;
            ProductTypes = productTypes;
            MyMsDataFile = myMsDataFile;
            DissociationTypes = DetermineDissociationType(productTypes);
        }

        protected override MetaMorpheusEngineResults RunSpecific()
        {
            TerminusType terminusType = ProductTypeMethods.IdentifyTerminusType(ProductTypes);

            foreach (var psm in AllResultingIdentifications)
            {
                if (GlobalVariables.StopLoops) { break; }

                psm.MatchedIonSeriesDict = new Dictionary<ProductType, int[]>();
                psm.ProductMassErrorDa = new Dictionary<ProductType, double[]>();
                psm.ProductMassErrorPpm = new Dictionary<ProductType, double[]>();
                psm.MatchedIonIntensitiesDict = new Dictionary<ProductType, double[]>();
                var theScan = MyMsDataFile.GetOneBasedScan(psm.ScanNumber);
                double thePrecursorMass = psm.ScanPrecursorMass;
                foreach (ProductType productType in ProductTypes)
                {
                    var sortedTheoreticalProductMasses = psm.CompactPeptides.First().Key.ProductMassesMightHaveDuplicatesAndNaNs(new List<ProductType> { productType });
                    Array.Sort(sortedTheoreticalProductMasses);
                    List<int> matchedIonSeriesList = new List<int>();
                    List<double> matchedIonMassToChargeRatioList = new List<double>();
                    List<double> productMassErrorDaList = new List<double>();
                    List<double> productMassErrorPpmList = new List<double>();
                    List<double> matchedIonIntensityList = new List<double>();

                    //populate the above lists
                    MatchIonsOld(theScan, CommonParameters.ProductMassTolerance, sortedTheoreticalProductMasses, matchedIonSeriesList, matchedIonMassToChargeRatioList, productMassErrorDaList, productMassErrorPpmList, matchedIonIntensityList, thePrecursorMass, productType, CommonParameters.AddCompIons);

                    psm.MatchedIonSeriesDict.Add(productType, matchedIonSeriesList.ToArray());
                    psm.MatchedIonMassToChargeRatioDict.Add(productType, matchedIonMassToChargeRatioList.ToArray());
                    psm.ProductMassErrorDa.Add(productType, productMassErrorDaList.ToArray());
                    psm.ProductMassErrorPpm.Add(productType, productMassErrorPpmList.ToArray());
                    psm.MatchedIonIntensitiesDict.Add(productType, matchedIonIntensityList.ToArray());
                }
            }

            foreach (PeptideSpectralMatch psm in AllResultingIdentifications.Where(b => b.NumDifferentCompactPeptides == 1))
            {
                if (GlobalVariables.StopLoops) { break; }

                var theScan = MyMsDataFile.GetOneBasedScan(psm.ScanNumber);
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

                    var gg = localizedPeptide.CompactPeptide(terminusType).ProductMassesMightHaveDuplicatesAndNaNs(ProductTypes);
                    Array.Sort(gg);
                    var score = CalculatePeptideScoreOld(theScan, CommonParameters.ProductMassTolerance, gg, thePrecursorMass, DissociationTypes, CommonParameters.AddCompIons, 0);
                    localizedScores.Add(score);
                }

                psm.LocalizedScores = localizedScores;
            }
            return new LocalizationEngineResults(this);
        }
    }
}