using Chemistry;
using MassSpectrometry;
using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;

namespace EngineLayer
{
    public abstract class MetaMorpheusEngine
    {
        protected readonly CommonParameters CommonParameters;
        protected readonly List<string> NestedIds;

        protected static readonly Dictionary<DissociationType, double> ComplementaryIonConversionDictionary = new Dictionary<DissociationType, double>
        {
            { DissociationType.HCD, Constants.ProtonMass },
            { DissociationType.ETD, 2 * Constants.ProtonMass }
        };

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

        public static void MatchIonsOld(MsDataScan thisScan, Tolerance productMassTolerance, double[] sortedTheoreticalProductMassesForThisPeptide, List<int> matchedIonSeries, List<double> matchedIonMassToChargeRatios, List<double> productMassErrorDa, List<double> productMassErrorPpm, List<double> matchedIonIntensitiesList, double precursorMass, ProductType productType, bool addCompIons)
        {
            var TotalProductsHere = sortedTheoreticalProductMassesForThisPeptide.Length;
            if (TotalProductsHere == 0)
            {
                return;
            }

            int currentTheoreticalIndex = -1;
            double currentTheoreticalMass;
            do
            {
                currentTheoreticalIndex++;
                currentTheoreticalMass = sortedTheoreticalProductMassesForThisPeptide[currentTheoreticalIndex];
            } while (double.IsNaN(currentTheoreticalMass) && currentTheoreticalIndex < sortedTheoreticalProductMassesForThisPeptide.Length - 1);

            if (double.IsNaN(currentTheoreticalMass))
            {
                return;
            }

            double currentTheoreticalMz = currentTheoreticalMass + Constants.ProtonMass;
            int testTheoreticalIndex;
            double testTheoreticalMass;
            double testTheoreticalMz;

            // speed optimizations
            double[] experimental_mzs = thisScan.MassSpectrum.XArray;
            double[] experimental_intensities = thisScan.MassSpectrum.YArray;
            int numExperimentalPeaks = experimental_mzs.Length;

            // Loop over all experimental indices
            for (int experimentalIndex = 0; experimentalIndex < numExperimentalPeaks; experimentalIndex++)
            {
                double currentExperimentalMz = experimental_mzs[experimentalIndex];
                // If found match

                if (productMassTolerance.Within(currentExperimentalMz, currentTheoreticalMz))
                {
                    matchedIonSeries.Add(++currentTheoreticalIndex); //++ because there's no such thing as a y0 ion.
                    matchedIonMassToChargeRatios.Add(currentTheoreticalMz);
                    matchedIonIntensitiesList.Add(experimental_intensities[experimentalIndex]);
                    double currentExperimentalMass = currentExperimentalMz - Constants.ProtonMass;
                    productMassErrorDa.Add(currentExperimentalMass - currentTheoreticalMass);
                    productMassErrorPpm.Add((currentExperimentalMass - currentTheoreticalMass) * 1000000 / currentTheoreticalMass);

                    if (currentTheoreticalIndex == TotalProductsHere)
                    {
                        break;
                    }
                    currentTheoreticalMass = sortedTheoreticalProductMassesForThisPeptide[currentTheoreticalIndex];
                    currentTheoreticalMz = currentTheoreticalMass + Constants.ProtonMass;
                }
                // Else if for sure did not reach the next theoretical yet
                else if (currentExperimentalMz > currentTheoreticalMz)
                {
                    // Move on to next index and never come back!
                    currentTheoreticalIndex++;
                    if (currentTheoreticalIndex == TotalProductsHere)
                    {
                        break;
                    }
                    currentTheoreticalMass = sortedTheoreticalProductMassesForThisPeptide[currentTheoreticalIndex];
                    currentTheoreticalMz = currentTheoreticalMass + Constants.ProtonMass;

                    // Start with the current ones
                    testTheoreticalIndex = currentTheoreticalIndex;
                    testTheoreticalMass = currentTheoreticalMass;
                    testTheoreticalMz = currentTheoreticalMz;
                    // Mark the skipped theoreticals as not found. The last one is not for sure, might be flipped!
                    while (currentExperimentalMz > testTheoreticalMz)
                    {
                        // Store old info for possible reuse
                        currentTheoreticalMz = testTheoreticalMz;
                        currentTheoreticalMass = testTheoreticalMass;
                        currentTheoreticalIndex = testTheoreticalIndex;

                        // Update test stuff!
                        testTheoreticalIndex++;
                        if (testTheoreticalIndex == TotalProductsHere)
                        {
                            break;
                        }
                        testTheoreticalMass = sortedTheoreticalProductMassesForThisPeptide[testTheoreticalIndex];
                        testTheoreticalMz = testTheoreticalMass + Constants.ProtonMass;
                    }
                    experimentalIndex--;
                }
            }
            if (addCompIons)
            {
                double[] complementaryMasses = new double[numExperimentalPeaks];
                double[] complementaryIntensities = new double[numExperimentalPeaks];

                foreach (DissociationType dissociationType in DetermineDissociationType(new List<ProductType> { productType }))
                {
                    if (ComplementaryIonConversionDictionary.TryGetValue(dissociationType, out double protonMassShift))
                    {
                        currentTheoreticalIndex = -1;
                        do
                        {
                            currentTheoreticalIndex++;
                            currentTheoreticalMass = sortedTheoreticalProductMassesForThisPeptide[currentTheoreticalIndex];
                        } while (double.IsNaN(currentTheoreticalMass) && currentTheoreticalIndex < sortedTheoreticalProductMassesForThisPeptide.Length - 1);

                        double massShiftForComplementaryConversion = precursorMass + protonMassShift; //mass shift needed to reobtain the original product ion for calculating tolerance
                        for (int i = numExperimentalPeaks - 1; i >= 0; i--)
                        {
                            complementaryMasses[numExperimentalPeaks - i - 1] = massShiftForComplementaryConversion - experimental_mzs[i];
                            complementaryIntensities[numExperimentalPeaks - i - 1] = experimental_intensities[i];
                        }

                        // Loop over all experimental indices
                        for (int experimentalIndex = 0; experimentalIndex < numExperimentalPeaks; experimentalIndex++)
                        {
                            double currentExperimentalMass = complementaryMasses[experimentalIndex];
                            double originalExperimentalMass = massShiftForComplementaryConversion - currentExperimentalMass;
                            double minBoundary = currentExperimentalMass - originalExperimentalMass + productMassTolerance.GetMinimumValue(originalExperimentalMass);
                            double maxBoundary = currentExperimentalMass - originalExperimentalMass + productMassTolerance.GetMaximumValue(originalExperimentalMass);
                            // If found match
                            if (minBoundary < currentTheoreticalMass && maxBoundary > currentTheoreticalMass)
                            {
                                matchedIonSeries.Add(++currentTheoreticalIndex);
                                matchedIonMassToChargeRatios.Add(currentTheoreticalMass.ToMz(1)); //currentTheoreticalMz is not updated
                                matchedIonIntensitiesList.Add(complementaryIntensities[experimentalIndex]);
                                productMassErrorDa.Add(currentExperimentalMass - currentTheoreticalMass);
                                productMassErrorPpm.Add((currentExperimentalMass - currentTheoreticalMass) * 1000000 / currentTheoreticalMass);

                                if (currentTheoreticalIndex == TotalProductsHere)
                                {
                                    break;
                                }
                                currentTheoreticalMass = sortedTheoreticalProductMassesForThisPeptide[currentTheoreticalIndex];
                            }
                            // Else if for sure passed a theoretical
                            else if (currentExperimentalMass > currentTheoreticalMass)
                            {
                                // Move on to next index and never come back!
                                currentTheoreticalIndex++;
                                if (currentTheoreticalIndex == TotalProductsHere)
                                {
                                    break;
                                }
                                currentTheoreticalMass = sortedTheoreticalProductMassesForThisPeptide[currentTheoreticalIndex];

                                // Start with the current ones
                                testTheoreticalIndex = currentTheoreticalIndex;
                                testTheoreticalMass = currentTheoreticalMass;
                                // Mark the skipped theoreticals as not found. The last one is not for sure, might be flipped!
                                while (currentExperimentalMass > testTheoreticalMass)
                                {
                                    // Store old info for possible reuse
                                    currentTheoreticalMass = testTheoreticalMass;
                                    currentTheoreticalIndex = testTheoreticalIndex;

                                    // Update test stuff!
                                    testTheoreticalIndex++;
                                    if (testTheoreticalIndex == TotalProductsHere)
                                    {
                                        break;
                                    }
                                    testTheoreticalMass = sortedTheoreticalProductMassesForThisPeptide[testTheoreticalIndex];
                                }
                                experimentalIndex--;
                            }
                        }
                    }
                    else
                    {
                        throw new NotImplementedException();
                    }
                }
            }

            //matchedIonSeries assumes 1-n, but bnob1 has no b1, so we need to ++ each series number
            if (productType == ProductType.BnoB1ions)
            {
                for (int i = 0; i < matchedIonSeries.Count; i++)
                {
                    matchedIonSeries[i]++;
                }
            }
        }

        public static double CalculatePeptideScoreOld(MsDataScan thisScan, Tolerance productMassTolerance, double[] sortedTheoreticalProductMassesForThisPeptide, double precursorMass, List<DissociationType> dissociationTypes, bool addCompIons, double maximumMassThatFragmentIonScoreIsDoubled)
        {
            var TotalProductsHere = sortedTheoreticalProductMassesForThisPeptide.Length;
            if (TotalProductsHere == 0)
            {
                return 0;
            }
            int MatchingProductsHere = 0;
            double MatchingIntensityHere = 0;

            int currentTheoreticalIndex = -1;
            double currentTheoreticalMass;
            do
            {
                currentTheoreticalIndex++;
                currentTheoreticalMass = sortedTheoreticalProductMassesForThisPeptide[currentTheoreticalIndex];
            } while (double.IsNaN(currentTheoreticalMass) && currentTheoreticalIndex < sortedTheoreticalProductMassesForThisPeptide.Length - 1);

            if (double.IsNaN(currentTheoreticalMass))
            {
                return 0;
            }

            double currentTheoreticalMz = currentTheoreticalMass + Constants.ProtonMass;
            int testTheoreticalIndex;
            double testTheoreticalMz;

            // speed optimizations
            double[] experimental_mzs = thisScan.MassSpectrum.XArray;
            double[] experimental_intensities = thisScan.MassSpectrum.YArray;
            int numExperimentalPeaks = experimental_mzs.Length;

            // Loop over all experimental indices
            for (int experimentalIndex = 0; experimentalIndex < numExperimentalPeaks; experimentalIndex++)
            {
                double currentExperimentalMz = experimental_mzs[experimentalIndex];
                // If found match
                if (productMassTolerance.Within(currentExperimentalMz, currentTheoreticalMz))
                {
                    MatchingProductsHere++;
                    if (maximumMassThatFragmentIonScoreIsDoubled > currentTheoreticalMz)
                        MatchingProductsHere++;
                    MatchingIntensityHere += experimental_intensities[experimentalIndex];

                    currentTheoreticalIndex++; //prevent multi counting
                    if (currentTheoreticalIndex == TotalProductsHere)
                    {
                        break;
                    }
                    currentTheoreticalMz = sortedTheoreticalProductMassesForThisPeptide[currentTheoreticalIndex] + Constants.ProtonMass;
                }
                // Else if for sure did not reach the next theoretical yet
                else if (currentExperimentalMz > currentTheoreticalMz)
                {
                    // Move on to next index and never come back!
                    currentTheoreticalIndex++;
                    if (currentTheoreticalIndex == TotalProductsHere)
                    {
                        break;
                    }
                    currentTheoreticalMz = sortedTheoreticalProductMassesForThisPeptide[currentTheoreticalIndex] + Constants.ProtonMass;

                    // Start with the current ones
                    testTheoreticalIndex = currentTheoreticalIndex;
                    testTheoreticalMz = currentTheoreticalMz;
                    // Mark the skipped theoreticals as not found. The last one is not for sure, might be flipped!
                    while (currentExperimentalMz > testTheoreticalMz)
                    {
                        // Store old info for possible reuse
                        currentTheoreticalMz = testTheoreticalMz;
                        currentTheoreticalIndex = testTheoreticalIndex;

                        // Update test stuff!
                        testTheoreticalIndex++;
                        if (testTheoreticalIndex == TotalProductsHere)
                        {
                            break;
                        }
                        testTheoreticalMz = sortedTheoreticalProductMassesForThisPeptide[testTheoreticalIndex] + Constants.ProtonMass;
                    }
                    experimentalIndex--;
                }
            }
            if (addCompIons)
            {
                double[] complementaryMasses = new double[numExperimentalPeaks];
                double[] complementaryIntensities = new double[numExperimentalPeaks];

                foreach (DissociationType dissociationType in dissociationTypes)
                {
                    double testTheoreticalMass;
                    if (ComplementaryIonConversionDictionary.TryGetValue(dissociationType, out double protonMassShift))
                    {
                        currentTheoreticalIndex = -1;
                        do
                        {
                            currentTheoreticalIndex++;
                            currentTheoreticalMass = sortedTheoreticalProductMassesForThisPeptide[currentTheoreticalIndex];
                        } while (double.IsNaN(currentTheoreticalMass) && currentTheoreticalIndex < sortedTheoreticalProductMassesForThisPeptide.Length - 1);

                        double massShiftForComplementaryConversion = precursorMass + protonMassShift; //mass shift needed to reobtain the original product ion for calculating tolerance
                        for (int i = numExperimentalPeaks - 1; i >= 0; i--)
                        {
                            complementaryMasses[numExperimentalPeaks - i - 1] = massShiftForComplementaryConversion - experimental_mzs[i];
                            complementaryIntensities[numExperimentalPeaks - i - 1] = experimental_intensities[i];
                        }

                        // Loop over all experimental indices
                        for (int experimentalIndex = 0; experimentalIndex < numExperimentalPeaks; experimentalIndex++)
                        {
                            double currentExperimentalMass = complementaryMasses[experimentalIndex];
                            double originalExperimentalMass = massShiftForComplementaryConversion - currentExperimentalMass;
                            double minBoundary = currentExperimentalMass - originalExperimentalMass + productMassTolerance.GetMinimumValue(originalExperimentalMass);
                            double maxBoundary = currentExperimentalMass - originalExperimentalMass + productMassTolerance.GetMaximumValue(originalExperimentalMass);
                            // If found match
                            if (minBoundary < currentTheoreticalMass && maxBoundary > currentTheoreticalMass)
                            {
                                MatchingProductsHere++;
                                if (maximumMassThatFragmentIonScoreIsDoubled > currentTheoreticalMass)
                                {
                                    MatchingProductsHere++;
                                }
                                MatchingIntensityHere += complementaryIntensities[experimentalIndex];

                                currentTheoreticalIndex++;
                                if (currentTheoreticalIndex == TotalProductsHere)
                                {
                                    break;
                                }
                                currentTheoreticalMass = sortedTheoreticalProductMassesForThisPeptide[currentTheoreticalIndex];
                            }
                            // Else if for sure passed a theoretical
                            else if (currentExperimentalMass > currentTheoreticalMass)
                            {
                                // Move on to next index and never come back!
                                currentTheoreticalIndex++;
                                if (currentTheoreticalIndex == TotalProductsHere)
                                {
                                    break;
                                }
                                currentTheoreticalMass = sortedTheoreticalProductMassesForThisPeptide[currentTheoreticalIndex];

                                // Start with the current ones
                                testTheoreticalIndex = currentTheoreticalIndex;
                                testTheoreticalMass = currentTheoreticalMass;
                                // Mark the skipped theoreticals as not found. The last one is not for sure, might be flipped!
                                while (currentExperimentalMass > testTheoreticalMass)
                                {
                                    // Store old info for possible reuse
                                    currentTheoreticalMass = testTheoreticalMass;
                                    currentTheoreticalIndex = testTheoreticalIndex;

                                    // Update test stuff!
                                    testTheoreticalIndex++;
                                    if (testTheoreticalIndex == TotalProductsHere)
                                    {
                                        break;
                                    }
                                    testTheoreticalMass = sortedTheoreticalProductMassesForThisPeptide[testTheoreticalIndex];
                                }
                                experimentalIndex--;
                            }
                        }
                    }
                    else
                    {
                        throw new NotImplementedException();
                    }
                }
            }
            return (MatchingProductsHere + MatchingIntensityHere / thisScan.TotalIonCurrent);
        }

        public static double CalculatePeptideScore(MsDataScan thisScan, List<MatchedFragmentIon> matchedFragmentIons, double maximumMassThatFragmentIonScoreIsDoubled)
        {
            // scoring if some fragments get doubled for scoring purposes
            if (maximumMassThatFragmentIonScoreIsDoubled > 0)
            {
                double score = 0;

                foreach (var fragment in matchedFragmentIons)
                {
                    double fragmentScore = 1 + (fragment.Intensity / thisScan.TotalIonCurrent);

                    if (fragment.TheoreticalFragmentIon.Mass <= maximumMassThatFragmentIonScoreIsDoubled)
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

        public static List<MatchedFragmentIon> MatchFragmentIons(MzSpectrum spectrum, List<TheoreticalFragmentIon> theoreticalFragmentIons, CommonParameters commonParameters)
        {
            var matchedFragmentIons = new List<MatchedFragmentIon>();
            var alreadyCountedMzs = new HashSet<double>();

            // if the spectrum has no peaks
            if (spectrum.Size == 0)
            {
                return matchedFragmentIons;
            }

            //search for ions in the spectrum
            foreach (var tIon in theoreticalFragmentIons)
            {
                // unknown fragment mass; this only happens rarely for sequences with unknown amino acids
                if (double.IsNaN(tIon.Mass))
                {
                    continue;
                }

                // get the closest peak in the spectrum to the theoretical peak
                int matchedPeakIndex = spectrum.GetClosestPeakIndex(tIon.Mz).Value;

                double mz = spectrum.XArray[matchedPeakIndex];
                double intensity = spectrum.YArray[matchedPeakIndex];

                // is the mass error acceptable and has it been counted already?
                if (commonParameters.ProductMassTolerance.Within(mz.ToMass(tIon.Charge), tIon.Mass) && !alreadyCountedMzs.Contains(mz))
                {
                    matchedFragmentIons.Add(new MatchedFragmentIon(tIon, mz, intensity));
                    alreadyCountedMzs.Add(mz);
                }
            }

            return matchedFragmentIons;
        }

        public static MzSpectrum GenerateComplementarySpectrum(MzSpectrum spectrum, double precursorMass, DissociationType dissociationType)
        {
            double protonMassShift = ComplementaryIonConversionDictionary[dissociationType].ToMass(1);

            double[] newMzSpectrum = new double[spectrum.Size];
            double[] intensity = new double[spectrum.Size];

            for (int i = spectrum.Size - 1; i >= 0; i--)
            {
                int j = spectrum.Size - i - 1;

                double mz = spectrum.XArray[i];
                double compFragmentMass = (precursorMass + protonMassShift) - mz.ToMass(1);

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
            return string.Join(",", NestedIds);
        }

        public static List<DissociationType> DetermineDissociationType(List<ProductType> lp)
        {
            List<DissociationType> dissociationTypes = new List<DissociationType>();

            if (lp.Contains(ProductType.B) || lp.Contains(ProductType.Y) || lp.Contains(ProductType.BnoB1ions))
            {
                dissociationTypes.Add(DissociationType.HCD);
            }
            if (lp.Contains(ProductType.C) || lp.Contains(ProductType.Zdot))
            {
                dissociationTypes.Add(DissociationType.ETD);
            }

            return dissociationTypes;
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