using Chemistry;
using MassSpectrometry;
using MzLibUtil;
using System;
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using Proteomics.Fragmentation;
using ClassExtensions = Chemistry.ClassExtensions;

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

        /// <summary>
        /// This is where deconvoluted MS2 isotopic envelopes are stored for fragment matching. The first int key is the
        /// one-based MS2 scan number, and the second int key is the rounded mass (to the nearest Dalton).
        /// </summary>
        protected Dictionary<int, IsotopicEnvelope[]> DeconvolutedMs2IsotopicEnvelopes;

        protected Dictionary<int, HashSet<double>> DeconvolutedPeakMzs;

        private static GlobalVariables.ValueComparer<IsotopicEnvelope> IsotopicEnvelopeMassComparer;

        protected MetaMorpheusEngine(CommonParameters commonParameters, List<string> nestedIds)
        {
            this.commonParameters = commonParameters;
            this.nestedIds = nestedIds;
            IsotopicEnvelopeMassComparer = new GlobalVariables.ValueComparer<IsotopicEnvelope>(f => f.monoisotopicMass);
        }
        
        public static event EventHandler<SingleEngineEventArgs> StartingSingleEngineHander;

        public static event EventHandler<SingleEngineFinishedEventArgs> FinishedSingleEngineHandler;

        public static event EventHandler<StringEventArgs> OutLabelStatusHandler;

        public static event EventHandler<StringEventArgs> WarnHandler;

        public static event EventHandler<ProgressEventArgs> OutProgressHandler;

        public void DeconvoluteAndStoreMs2(MsDataScan[] ms2Scans)
        {
            DeconvolutedMs2IsotopicEnvelopes = new Dictionary<int, IsotopicEnvelope[]>();
            DeconvolutedPeakMzs = new Dictionary<int, HashSet<double>>();

            double ms2DeconvolutionPpmTolerance = 5.0;
            int minZ = 1;
            int maxZ = 10;

            foreach (MsDataScan scan in ms2Scans)
            {
                if (DeconvolutedMs2IsotopicEnvelopes.ContainsKey(scan.OneBasedScanNumber))
                {
                    continue;
                }

                // deconvolute the scan
                var isotopicEnvelopes = scan.MassSpectrum.Deconvolute(scan.MassSpectrum.Range, minZ, maxZ,
                    ms2DeconvolutionPpmTolerance, commonParameters.DeconvolutionIntensityRatio).OrderBy(p => p.monoisotopicMass).ToArray();

                // store the scan's deconvoluted envelopes
                DeconvolutedMs2IsotopicEnvelopes.Add(scan.OneBasedScanNumber, isotopicEnvelopes);

                // store the scan's deconvoluted envelope peaks
                HashSet<double> deconvolutedMzForThisScan = new HashSet<double>();

                foreach (IsotopicEnvelope isotopicEnvelope in isotopicEnvelopes)
                {
                    foreach (var peak in isotopicEnvelope.peaks)
                    {
                        deconvolutedMzForThisScan.Add(ClassExtensions.RoundedDouble(peak.mz).Value);
                    }
                }

                DeconvolutedPeakMzs.Add(scan.OneBasedScanNumber, deconvolutedMzForThisScan);
            }
        }

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

        public static List<MatchedFragmentIon> MatchFragmentIons(MsDataScan scan, List<Product> theoreticalProducts, CommonParameters commonParameters,
            double precursorMass, int precursorCharge, Dictionary<int, IsotopicEnvelope[]> deconvolutedMs2IsotopicEnvelopes, Dictionary<int, HashSet<double>> deconvolutedPeakMzs)
        {
            var matchedFragmentIons = new List<MatchedFragmentIon>();

            // if the spectrum has no peaks
            if (scan.MassSpectrum.Size == 0)
            {
                return matchedFragmentIons;
            }

            //search for ions in the spectrum
            foreach (Product product in theoreticalProducts)
            {
                // unknown fragment mass; this only happens rarely for sequences with unknown amino acids
                if (double.IsNaN(product.NeutralMass) || product.NeutralMass <= 0)
                {
                    continue;
                }

                // look for higher-charge fragments via deconvolution
                if (commonParameters.DeconvoluteMs2)
                {
                    IsotopicEnvelope bestEnvelope = BinarySearchForDeconvolutedMass(deconvolutedMs2IsotopicEnvelopes[scan.OneBasedScanNumber], product);

                    if (bestEnvelope != null && commonParameters.ProductMassTolerance.Within(bestEnvelope.monoisotopicMass, product.NeutralMass) && bestEnvelope.charge <= precursorCharge)
                    {
                        matchedFragmentIons.Add(new MatchedFragmentIon(product, bestEnvelope.monoisotopicMass.ToMz(bestEnvelope.charge), bestEnvelope.totalIntensity, bestEnvelope.charge));
                    }
                }

                // look for z=1 fragments without requiring a second isotope peak to be present
                if (commonParameters.AssumeFragmentsAreZ1)
                {
                    // get the closest peak in the spectrum to the theoretical peak assuming z=1
                    int matchedPeakIndex = scan.MassSpectrum.GetClosestPeakIndex(product.NeutralMass.ToMz(1)).Value;

                    double mz = scan.MassSpectrum.XArray[matchedPeakIndex];

                    // is the mass error acceptable and make sure it is not part a deconvoluted isotopic envelope
                    if (commonParameters.ProductMassTolerance.Within(mz, product.NeutralMass.ToMz(1)) 
                        && (deconvolutedPeakMzs == null || !deconvolutedPeakMzs[scan.OneBasedScanNumber].Contains(ClassExtensions.RoundedDouble(mz).Value)))
                    {
                        matchedFragmentIons.Add(new MatchedFragmentIon(product, mz, scan.MassSpectrum.YArray[matchedPeakIndex], 1));
                    }
                }
            }
            if (commonParameters.AddCompIons)//needs to be separate to account for ppm error differences
            {
                double protonMassShift = complementaryIonConversionDictionary[commonParameters.DissociationType].ToMass(1);
                double sumOfCompIonsMz = (precursorMass + protonMassShift).ToMz(1); //FIXME, not valid for all fragmentation (b+y+H = precursor, but c+zdot+2H = precursor)
                foreach (Product product in theoreticalProducts)
                {
                    // unknown fragment mass; this only happens rarely for sequences with unknown amino acids
                    if (double.IsNaN(product.NeutralMass))
                    {
                        continue;
                    }

                    // get the closest peak in the spectrum to the theoretical peak assuming z=1
                    //generate a "comp" product
                    double theoreticalCompMz = sumOfCompIonsMz - product.NeutralMass; //This is NOT the m/z of the product, 
                    //but the m/z of the theoretical complementary that we are looking for in the experimental spectrum.
                    //The complementary to this experimental match will match to the original theoretical product
                    int matchedPeakIndex = scan.MassSpectrum.GetClosestPeakIndex(theoreticalCompMz).Value; //search for the comp ion
                    double mzToCompare = scan.MassSpectrum.XArray[matchedPeakIndex]; //we need the original mz to know the error associated with the comp mz

                    // is the mass error acceptable and has it been counted already?
                    //Need to compare the "noncomplementary" peaks so that the correct mass tolerance is used for Ppm tolerances
                    if (commonParameters.ProductMassTolerance.Within(mzToCompare, theoreticalCompMz))
                    {
                        //the sumOfCompIons - original peak must be converted to mz, because subtracting an mz from an mz creates a mass difference, not an mz (5-3 (m/z) = 2 = 4-2 (mass))
                        matchedFragmentIons.Add(new MatchedFragmentIon(product, (sumOfCompIonsMz - scan.MassSpectrum.XArray[matchedPeakIndex]).ToMz(1),
                            scan.MassSpectrum.YArray[matchedPeakIndex], 1));
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

        private static IsotopicEnvelope BinarySearchForDeconvolutedMass(IsotopicEnvelope[] deconvolutedEnvsForThisScan, Product product)
        {
            if (deconvolutedEnvsForThisScan.Length == 0)
            {
                return null;
            }

            int index = Array.BinarySearch(deconvolutedEnvsForThisScan, new IsotopicEnvelope(null, product.NeutralMass, 0, 0, 0, 0),
                IsotopicEnvelopeMassComparer);

            if (index >= 0)
            {
                return deconvolutedEnvsForThisScan[index];
            }

            index = ~index;

            if (index >= deconvolutedEnvsForThisScan.Length)
            {
                return deconvolutedEnvsForThisScan[index - 1];
            }

            if (index == 0)
            {
                return deconvolutedEnvsForThisScan[index];
            }

            if (product.NeutralMass - deconvolutedEnvsForThisScan[index - 1].monoisotopicMass >
                deconvolutedEnvsForThisScan[index].monoisotopicMass - product.NeutralMass)
            {
                return deconvolutedEnvsForThisScan[index];
            }

            return deconvolutedEnvsForThisScan[index - 1];
        }
    }
}