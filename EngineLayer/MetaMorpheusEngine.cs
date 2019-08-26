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

        protected readonly CommonParameters CommonParameters;

        protected readonly List<string> NestedIds;

        protected MetaMorpheusEngine(CommonParameters commonParameters, List<string> nestedIds)
        {
            CommonParameters = commonParameters;
            NestedIds = nestedIds;
        }

        public static event EventHandler<SingleEngineEventArgs> StartingSingleEngineHander;

        public static event EventHandler<SingleEngineFinishedEventArgs> FinishedSingleEngineHandler;

        public static event EventHandler<StringEventArgs> OutLabelStatusHandler;

        public static event EventHandler<StringEventArgs> WarnHandler;

        public static event EventHandler<ProgressEventArgs> OutProgressHandler;

        public static double CalculatePeptideScore(MsDataScan thisScan, List<MatchedFragmentIon> matchedFragmentIons)
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

            if (scan.TheScan.MassSpectrum.XcorrProcessed)
            {
                foreach (Product product in theoreticalProducts)
                {
                    // unknown fragment mass; this only happens rarely for sequences with unknown amino acids
                    if (double.IsNaN(product.NeutralMass))
                    {
                        continue;
                    }

                    double theoreticalFragmentMz = Math.Round(product.NeutralMass.ToMz(1) / 1.0005079, 0) * 1.0005079;
                    var closestMzIndex = scan.TheScan.MassSpectrum.GetClosestPeakIndex(theoreticalFragmentMz).Value;

                    if (commonParameters.ProductMassTolerance.Within(scan.TheScan.MassSpectrum.XArray[closestMzIndex], theoreticalFragmentMz))
                    {
                        matchedFragmentIons.Add(new MatchedFragmentIon(product, theoreticalFragmentMz, scan.TheScan.MassSpectrum.YArray[closestMzIndex], 1));
                    }
                }

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

        public static List<MatchedFragmentIon> MatchOriginFragmentIons(Ms2ScanWithSpecificMass scan, List<Product> theoreticalProducts, CommonParameters commonParameters)
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

                if (product.ProductType == ProductType.M)
                {
                    for (int i = 1; i <= scan.PrecursorCharge; i++)
                    {

                        var closestExperimentalMz = scan.GetClosestExperimentalFragmentMz(product.NeutralMass.ToMz(i), out double? intensity);

                        if (closestExperimentalMz.HasValue && commonParameters.ProductMassTolerance.Within(closestExperimentalMz.Value, product.NeutralMass.ToMz(i)))
                        {
                            matchedFragmentIons.Add(new MatchedFragmentIon(product, closestExperimentalMz.Value, intensity.Value, i));
                        }
                    }
                }

                else
                {
                    // get the closest peak in the spectrum to the theoretical peak
                    var closestExperimentalMass = scan.GetClosestExperimentalFragmentMass(product.NeutralMass);

                    // is the mass error acceptable?
                    if (commonParameters.ProductMassTolerance.Within(closestExperimentalMass.monoisotopicMass, product.NeutralMass) && closestExperimentalMass.charge <= scan.PrecursorCharge)
                    {
                        matchedFragmentIons.Add(new MatchedFragmentIon(product, closestExperimentalMass.monoisotopicMass.ToMz(closestExperimentalMass.charge),
                            closestExperimentalMass.peaks.First().intensity, closestExperimentalMass.charge));
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
            return string.Join(",", NestedIds);
        }

        protected void Warn(string v)
        {
            WarnHandler?.Invoke(this, new StringEventArgs(v, NestedIds));
        }

        protected void Status(string v)
        {
            OutLabelStatusHandler?.Invoke(this, new StringEventArgs(v, NestedIds));
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