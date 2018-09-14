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
            // scoring if some fragments get doubled for scoring purposes
            if (maximumMassThatFragmentIonScoreIsDoubled > 0)
            {
                double score = 0;

                foreach (var fragment in matchedFragmentIons)
                {
                    double fragmentScore = 1 + (fragment.Intensity / thisScan.TotalIonCurrent);

                    if (fragment.NeutralTheoreticalProduct.NeutralMass <= maximumMassThatFragmentIonScoreIsDoubled) // TODO: not sure if this is supposed to be neutral mass or mz
                    {
                        score += fragmentScore * 2;
                    }
                    else
                    {
                        score += fragmentScore;
                    }
                }

                return score;
            }

            // normal scoring
            return matchedFragmentIons.Count + (matchedFragmentIons.Sum(v => v.Intensity) / thisScan.TotalIonCurrent);
        }

        public static List<MatchedFragmentIon> MatchFragmentIons(MzSpectrum spectrum, List<Product> theoreticalProducts, CommonParameters commonParameters)
        {
            var matchedFragmentIons = new List<MatchedFragmentIon>();
            var alreadyCountedMzs = new HashSet<double>();

            // if the spectrum has no peaks
            if (spectrum.Size == 0)
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

                // get the closest peak in the spectrum to the theoretical peak assuming z=1
                int matchedPeakIndex = spectrum.GetClosestPeakIndex(product.NeutralMass.ToMz(1)).Value;

                double mz = spectrum.XArray[matchedPeakIndex];
                double intensity = spectrum.YArray[matchedPeakIndex];

                // is the mass error acceptable and has it been counted already?
                if (commonParameters.ProductMassTolerance.Within(mz, product.NeutralMass.ToMz(1)) && !alreadyCountedMzs.Contains(mz))
                {
                    matchedFragmentIons.Add(new MatchedFragmentIon(product, mz, intensity, 1));
                    alreadyCountedMzs.Add(mz);
                }
            }

            return matchedFragmentIons;
        }

        public static MzSpectrum GenerateComplementarySpectrum(MzSpectrum spectrum, double precursorMass, DissociationType dissociationType)
        {
            double protonMassShift = complementaryIonConversionDictionary[dissociationType].ToMass(1);

            double[] newMzSpectrum = new double[spectrum.Size];
            double[] intensity = new double[spectrum.Size];

            for (int i = spectrum.Size - 1; i >= 0; i--)
            {
                int j = spectrum.Size - i - 1;

                double mz = spectrum.XArray[i];
                double compFragmentMass = (precursorMass + protonMassShift) - mz.ToMass(1); //FIXME, not valid for all fragmentation (b+y+H = precursor, but c+zdot+2H = precursor)

                newMzSpectrum[j] = compFragmentMass.ToMz(1);
                intensity[j] = spectrum.YArray[i];
            }

            return new MzSpectrum(newMzSpectrum, intensity, false);
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

        public static List<DissociationType> DetermineDissociationType(List<ProductType> lp)
        {
            List<DissociationType> dissociationTypes = new List<DissociationType>();

            if (lp.Contains(ProductType.b) || lp.Contains(ProductType.y))
            {
                dissociationTypes.Add(DissociationType.HCD);
            }

            if (lp.Contains(ProductType.c) || lp.Contains(ProductType.zPlusOne))
            {
                dissociationTypes.Add(DissociationType.ETD);
            }

            return dissociationTypes;
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