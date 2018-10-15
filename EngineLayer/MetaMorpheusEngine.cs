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

        protected StoredDeconvolutedMs2Envelopes StoredDeconvolutedMs2Envelopes;

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

        public void DeconvoluteAndStoreMs2(MsDataScan[] ms2Scans)
        {
            StoredDeconvolutedMs2Envelopes = StoredDeconvolutedMs2Envelopes.DeconvoluteAndStoreMs2Scans(ms2Scans, commonParameters);
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

        public static List<MatchedFragmentIon> MatchFragmentIons(MsDataScan scan, List<Product> theoreticalProducts, CommonParameters commonParameters, double precursorMass, 
            int precursorCharge, StoredDeconvolutedMs2Envelopes storedDeconvolutedMs2Envelopes)
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
                if (double.IsNaN(product.NeutralMass))
                {
                    continue;
                }

                // look for higher-charge fragments via deconvolution
                if (commonParameters.DeconvoluteMs2)
                {
                    IsotopicEnvelope bestEnvelope = storedDeconvolutedMs2Envelopes.BinarySearchForDeconvolutedMass(scan.OneBasedScanNumber, product);

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
                        && (storedDeconvolutedMs2Envelopes == null || !storedDeconvolutedMs2Envelopes.HasClaimedMzPeak(scan.OneBasedScanNumber, mz)))
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
    }
}