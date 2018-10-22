using Chemistry;
using MassSpectrometry;
using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using Proteomics.Fragmentation;

namespace EngineLayer
{
    public abstract class MetaMorpheusEngine
    {
        protected static readonly Dictionary<DissociationType, double> complementaryIonConversionDictionary = new Dictionary<DissociationType, double>
        {
            { DissociationType.HCD, Constants.ProtonMass },
            { DissociationType.ETD, 2 * Constants.ProtonMass }, //presence of zplusone (zdot) makes this two instead of one
            { DissociationType.CID, Constants.ProtonMass },
            //TODO: refactor such that complementary ions are generated specifically for their complementary pair. 
            //TODO: create a method to auto-determine the conversion
        };

        protected readonly CommonParameters commonParameters;

        protected readonly List<string> nestedIds;
        
        protected MetaMorpheusEngine(CommonParameters commonParameters, List<string> nestedIds)
        {
            this.commonParameters = commonParameters;
            this.nestedIds = nestedIds;
        }

        public static event EventHandler<SingleEngineEventArgs> StartingSingleEngineHander;

        public static event EventHandler<SingleEngineFinishedEventArgs> FinishedSingleEngineHandler;

        public static event EventHandler<StringEventArgs> OutLabelStatusHandler;

        public static event EventHandler<StringEventArgs> WarnHandler;

        public static event EventHandler<ProgressEventArgs> OutProgressHandler;
        
        public static double CalculatePeptideScore(MsDataScan thisScan, List<MatchedFragmentIon> matchedFragmentIons, double maximumMassThatFragmentIonScoreIsDoubled)
        {
            double score = 0;

            foreach (var fragment in matchedFragmentIons)
            {
                double fragmentScore = 1 + (fragment.Intensity / thisScan.TotalIonCurrent);
                score += fragmentScore;

                if (fragment.NeutralTheoreticalProduct.NeutralMass <= maximumMassThatFragmentIonScoreIsDoubled)
                {
                    score += fragmentScore;
                }
            }

            return score;
        }

        public static List<MatchedFragmentIon> MatchFragmentIons(Ms2ScanWithSpecificMass scan, List<Product> theoreticalProducts, CommonParameters commonParameters)
        {
            var matchedFragmentIons = new List<MatchedFragmentIon>();

            // if the spectrum has no peaks
            if (!scan.ExperimentalFragments.Any())
            {
                return matchedFragmentIons;
            }

            // search for ions in the spectrum
            foreach (Product product in theoreticalProducts)
            {
                // unknown fragment mass; this only happens rarely for sequences with unknown amino acids
                if (double.IsNaN(product.NeutralMass))
                {
                    continue;
                }

                // get the closest peak in the spectrum to the theoretical peak
                var closestExperimentalMass = scan.GetClosestExperimentalFragmentMass(product.NeutralMass);
                
                // is the mass error acceptable?
                if (commonParameters.ProductMassTolerance.Within(closestExperimentalMass.monoisotopicMass, product.NeutralMass) && closestExperimentalMass.charge <= scan.PrecursorCharge)
                {
                    matchedFragmentIons.Add(new MatchedFragmentIon(product, closestExperimentalMass.monoisotopicMass.ToMz(closestExperimentalMass.charge), 
                        closestExperimentalMass.peaks.First().intensity, closestExperimentalMass.charge));
                }
            }
            if (commonParameters.AddCompIons)
            {
                double protonMassShift = complementaryIonConversionDictionary[commonParameters.DissociationType].ToMass(1);

                foreach (Product product in theoreticalProducts)
                {
                    // unknown fragment mass or diagnostic ion or precursor; skip those
                    if (double.IsNaN(product.NeutralMass) || product.ProductType == ProductType.D || product.ProductType == ProductType.M)
                    {
                        continue;
                    }
                    
                    double compIonMass = scan.PrecursorMass + protonMassShift - product.NeutralMass;

                    // get the closest peak in the spectrum to the theoretical peak
                    var closestExperimentalMass = scan.GetClosestExperimentalFragmentMass(compIonMass);

                    // is the mass error acceptable?
                    if (commonParameters.ProductMassTolerance.Within(closestExperimentalMass.monoisotopicMass, compIonMass) && closestExperimentalMass.charge <= scan.PrecursorCharge)
                    {
                        matchedFragmentIons.Add(new MatchedFragmentIon(product, closestExperimentalMass.monoisotopicMass.ToMz(closestExperimentalMass.charge),
                            closestExperimentalMass.totalIntensity, closestExperimentalMass.charge));
                    }
                }
            }

            return matchedFragmentIons;
        }

        public MetaMorpheusEngineResults Run()
        {
            StartingSingleEngine();
            var stopWatch = new Stopwatch();
            stopWatch.Start();
            var myResults = RunSpecific();
            stopWatch.Stop();
            myResults.Time = stopWatch.Elapsed;
            FinishedSingleEngine(myResults);
            return myResults;
        }

        public string GetId()
        {
            return string.Join(",", nestedIds);
        }

        protected void Warn(string v)
        {
            WarnHandler?.Invoke(this, new StringEventArgs(v, nestedIds));
        }

        protected void Status(string v)
        {
            OutLabelStatusHandler?.Invoke(this, new StringEventArgs(v, nestedIds));
        }

        protected void ReportProgress(ProgressEventArgs v)
        {
            OutProgressHandler?.Invoke(this, v);
        }

        protected abstract MetaMorpheusEngineResults RunSpecific();

        private void StartingSingleEngine()
        {
            StartingSingleEngineHander?.Invoke(this, new SingleEngineEventArgs(this));
        }

        private void FinishedSingleEngine(MetaMorpheusEngineResults myResults)
        {
            FinishedSingleEngineHandler?.Invoke(this, new SingleEngineFinishedEventArgs(myResults));
        }
    }
}