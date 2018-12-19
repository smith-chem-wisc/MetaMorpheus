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

            if (thisScan.MassSpectrum.XcorrProcessed)
            {
                foreach (var fragment in matchedFragmentIons)
                {
                    switch (fragment.NeutralTheoreticalProduct.ProductType)
                    {
                        case ProductType.aDegree:
                        case ProductType.aStar:
                        case ProductType.bDegree:
                        case ProductType.bStar:
                        case ProductType.yDegree:
                        case ProductType.yStar:
                            score += 0.01 * fragment.Intensity;
                            break;

                        default:
                            score += 1 * fragment.Intensity;
                            break;
                    }
                }
            }
            else
            {
                foreach (var fragment in matchedFragmentIons)
                {
                    double fragmentScore = 1 + (fragment.Intensity / thisScan.TotalIonCurrent);
                    score += fragmentScore;

                    if (fragment.NeutralTheoreticalProduct.NeutralMass <= maximumMassThatFragmentIonScoreIsDoubled)
                    {
                        score += fragmentScore;
                    }
                }
            }

            return score;
        }

        public static double CalculatePeptideXcorr(Ms2ScanWithSpecificMass thisScan, List<Product> peptideTheorProducts)
        {
            double xcorr = 0;
            //double discreteMassBin = Constants.C13MinusC12;
            //double discreteMassBin = 1.001123503868;
            double discreteMassBin = 1.0005079;
            //double discreteMassBin = 1.0005079/2;
            //int multiplier = (int)Math.Round(2000 / discreteMassBin, 0);
            //if (!thisScan.TheScan.MassSpectrum.XcorrProcessed)
            //{
            //    thisScan.TheScan.MassSpectrum.XCorrPrePreprocessing(0, multiplier * discreteMassBin, thisScan.TheScan.IsolationMz.Value, 1.5, discreteMassBin, 0.05);
            //}

            SortedDictionary<int, double> massInt = new SortedDictionary<int, double>();

            for (int z = 1; z < 2; z++)
            //for (int z = 1; z < thisScan.SelectedIonChargeStateGuess; z++)
            {
                foreach (Product theoreticalProduct in peptideTheorProducts)
                {
                    double theoreticalProductMass = Math.Round(theoreticalProduct.NeutralMass.ToMz(z) / discreteMassBin, 0) * discreteMassBin;

                    int index = (int)Math.Round(theoreticalProductMass / discreteMassBin, 0);

                    switch (theoreticalProduct.ProductType)
                    {
                        case ProductType.aDegree:
                        case ProductType.aStar:
                        case ProductType.bDegree:
                        case ProductType.bStar:
                        case ProductType.yDegree:
                        case ProductType.yStar:
                            if (massInt.Keys.Contains(index))
                            {
                                massInt[index] += 0.01;
                            }
                            else
                            {
                                massInt.Add(index, 0.01);
                            }
                            break;

                        default:
                            if (massInt.Keys.Contains(index))
                            {
                                massInt[index] += 1;
                            }
                            else
                            {
                                massInt.Add(index, 1);
                            }
                            break;
                    }
                }
            }

            double[] theoreticalFragmentIonMassArray = new double[massInt.Keys.Count()];
            double[] theoreticalFragmentIonIntensityArray = massInt.Values.ToArray();

            int myI = 0;
            foreach (int key in massInt.Keys)
            {
                theoreticalFragmentIonMassArray[myI] = (double)key * discreteMassBin;
                myI++;
            }

            int min = 0;
            double tolerance = discreteMassBin / 2;
            for (int i = 0; i < theoreticalFragmentIonMassArray.Length; i++)
            {
                for (int j = min; j < thisScan.TheScan.MassSpectrum.XArray.Length; j++)
                {
                    if (MassInTolerance(theoreticalFragmentIonMassArray[i], thisScan.TheScan.MassSpectrum.XArray[j], tolerance))
                    {
                        min = j;// no reason to start loop below this mass for next iteration.
                        xcorr += theoreticalFragmentIonIntensityArray[i] * thisScan.TheScan.MassSpectrum.YArray[j];
                        break;
                    }
                    if ((thisScan.TheScan.MassSpectrum.XArray[j] - tolerance) > theoreticalFragmentIonMassArray[i])
                    {
                        break;
                    }
                }
            }

            return xcorr;
        }

        private static bool MassInTolerance(double theoreticalMass, double experimentalMass, double thomsonTolerance)
        {
            if (Math.Abs(theoreticalMass - experimentalMass) < thomsonTolerance)
            {
                return true;
            }
            else
            {
                return false;
            }
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
            double discreteMassBin = 1.0005079;
            int multiplier = (int)Math.Round(2000 / discreteMassBin, 0);

            // if the spectrum has no peaks
            if (!scan.ExperimentalFragments.Any())
            {
                return matchedFragmentIons;
            }

            if (commonParameters.DissociationType == MassSpectrometry.DissociationType.LowCID)
            {
                if (!scan.TheScan.MassSpectrum.XcorrProcessed)
                {
                    scan.TheScan.MassSpectrum.XCorrPrePreprocessing(0, multiplier * discreteMassBin, scan.TheScan.IsolationMz.Value, 1.5, discreteMassBin, 0.05);
                    Ms2ScanWithSpecificMass.GetNeutralExperimentalFragments(scan.TheScan, commonParameters);
                }
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