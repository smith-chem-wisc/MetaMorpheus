using Chemistry;
using MassSpectrometry;
using Proteomics.Fragmentation;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;

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

        public static double CalculatePeptideXcorr(MsDataScan thisScan, List<Product> peptideTheorProducts)
        {
            double xcorr = 0;

            double[] scanMassArray = new double[thisScan.MassSpectrum.XArray.Count()];
            double[] scanIntensityArray = new double[thisScan.MassSpectrum.YArray.Count()];

            scanMassArray = thisScan.MassSpectrum.XArray.ToArray();
            scanIntensityArray = thisScan.MassSpectrum.YArray.ToArray();
            double scanRangeMinMz = 0;
            double scanRangeMaxMz = 1968.999547;

            MsDataFile.XCorrPrePreprocessing(ref scanIntensityArray, ref scanMassArray, scanRangeMinMz, scanRangeMaxMz, thisScan.IsolationMz.Value);


            for (int i = 0; i < scanIntensityArray.Length; i++)
            {
                if (Double.IsNaN(scanIntensityArray[i])||Double.IsInfinity(scanIntensityArray[i]))
                {
                    scanIntensityArray[i] = 0;
                }
            }




            double[] theoreticalFragmentIonMassArray = scanMassArray;
            double[] theoreticalFragmentIonIntensityArray = new double[theoreticalFragmentIonMassArray.Count()];

            for (int z = 1; z < 2; z++)
            //for (int z = 1; z < thisScan.SelectedIonChargeStateGuess; z++)
            {
                foreach (Product theoreticalProduct in peptideTheorProducts)
                {
                    double theoreticalProductMass = Math.Round(theoreticalProduct.NeutralMass.ToMz(z) / 1.0005079, 0) * 1.0005079;
                    int myIndex = FirstArrayIndexOfDouble(theoreticalFragmentIonMassArray, theoreticalProductMass, 0.000001);
                    if(myIndex > -1)
                    {
                        switch (theoreticalProduct.ProductType)
                        {
                            case ProductType.aDegree:
                            case ProductType.aStar:
                            case ProductType.bDegree:
                            case ProductType.bStar:
                            case ProductType.yDegree:
                            case ProductType.yStar:
                                theoreticalFragmentIonIntensityArray[myIndex] += 0.2;
                                break;
                            default:
                                theoreticalFragmentIonIntensityArray[myIndex] += 1;
                                break;

                        }
                    }
                    

                }
            }

            for (int i = 0; i < scanIntensityArray.Length; i++)
            {
                if (Double.IsNaN(scanIntensityArray[i])|| Double.IsNaN(theoreticalFragmentIonIntensityArray[i])   || Double.IsInfinity(scanIntensityArray[i]) || Double.IsInfinity(theoreticalFragmentIonIntensityArray[i]))
                {
                    int bubba = 1;
                }
                xcorr += scanIntensityArray[i] * theoreticalFragmentIonIntensityArray[i];
            }

            return xcorr;
        }

        private static int FirstArrayIndexOfDouble(double[] a, double m, double t)
        {
            for (int i = 0; i < a.Length; i++)
            {
                if (((a[i] - t) < m) && ((a[i] + t) > m))
                    return i;
            }
            return -1;
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