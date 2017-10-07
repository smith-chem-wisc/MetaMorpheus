using Chemistry;
using MassSpectrometry;
using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Diagnostics;

namespace EngineLayer
{
    public abstract class MetaMorpheusEngine
    {
        #region Protected Fields

        protected static readonly Dictionary<DissociationType, double> complementaryIonConversionDictionary = new Dictionary<DissociationType, double>
        {
            { DissociationType.HCD, Constants.protonMass },
            { DissociationType.ETD, 2*Constants.protonMass }
        };

        protected readonly List<string> nestedIds;

        #endregion Protected Fields

        #region Protected Constructors

        protected MetaMorpheusEngine(List<string> nestedIds)
        {
            this.nestedIds = nestedIds;
        }

        #endregion Protected Constructors

        #region Public Events

        public static event EventHandler<SingleEngineEventArgs> StartingSingleEngineHander;

        public static event EventHandler<SingleEngineFinishedEventArgs> FinishedSingleEngineHandler;

        public static event EventHandler<StringEventArgs> OutLabelStatusHandler;

        public static event EventHandler<StringEventArgs> WarnHandler;

        public static event EventHandler<ProgressEventArgs> OutProgressHandler;

        #endregion Public Events

        #region Public Methods

        public static void MatchIons(IMsDataScan<IMzSpectrum<IMzPeak>> thisScan, Tolerance productMassTolerance, double[] sortedTheoreticalProductMassesForThisPeptide, List<double> matchedIonMassesList, List<double> productMassErrorDa, List<double> productMassErrorPpm, double precursorMass, List<DissociationType> dissociationTypes, bool addCompIons)
        {
            var TotalProductsHere = sortedTheoreticalProductMassesForThisPeptide.Length;
            if (TotalProductsHere == 0)
                return;

            int currentTheoreticalIndex = -1;
            double currentTheoreticalMass;
            do
            {
                currentTheoreticalIndex++;
                currentTheoreticalMass = sortedTheoreticalProductMassesForThisPeptide[currentTheoreticalIndex];
            } while (double.IsNaN(currentTheoreticalMass) && currentTheoreticalIndex < sortedTheoreticalProductMassesForThisPeptide.Length - 1);

            if (double.IsNaN(currentTheoreticalMass))
                return;

            double currentTheoreticalMz = currentTheoreticalMass + Constants.protonMass;
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
                    matchedIonMassesList.Add(currentTheoreticalMass);
                    double currentExperimentalMass = currentExperimentalMz - Constants.protonMass;
                    productMassErrorDa.Add(currentExperimentalMass - currentTheoreticalMass);
                    productMassErrorPpm.Add((currentExperimentalMass - currentTheoreticalMass) * 1000000 / currentTheoreticalMass);

                    currentTheoreticalIndex++;
                    if (currentTheoreticalIndex == TotalProductsHere)
                        break;
                    currentTheoreticalMass = sortedTheoreticalProductMassesForThisPeptide[currentTheoreticalIndex];
                    currentTheoreticalMz = currentTheoreticalMass + Constants.protonMass;
                }
                // Else if for sure did not reach the next theoretical yet
                else if (currentExperimentalMz > currentTheoreticalMz)
                {
                    // Move on to next index and never come back!
                    currentTheoreticalIndex++;
                    if (currentTheoreticalIndex == TotalProductsHere)
                        break;
                    currentTheoreticalMass = sortedTheoreticalProductMassesForThisPeptide[currentTheoreticalIndex];
                    currentTheoreticalMz = currentTheoreticalMass + Constants.protonMass;

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
                            break;
                        testTheoreticalMass = sortedTheoreticalProductMassesForThisPeptide[testTheoreticalIndex];
                        testTheoreticalMz = testTheoreticalMass + Constants.protonMass;
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
                    if (complementaryIonConversionDictionary.TryGetValue(dissociationType, out double protonMassShift))
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
                                matchedIonMassesList.Add(currentTheoreticalMass);
                                productMassErrorDa.Add(currentExperimentalMass - currentTheoreticalMass);
                                productMassErrorPpm.Add((currentExperimentalMass - currentTheoreticalMass) * 1000000 / currentTheoreticalMass);

                                currentTheoreticalIndex++;
                                if (currentTheoreticalIndex == TotalProductsHere)
                                    break;
                                currentTheoreticalMass = sortedTheoreticalProductMassesForThisPeptide[currentTheoreticalIndex];
                            }
                            // Else if for sure passed a theoretical
                            else if (currentExperimentalMass > currentTheoreticalMass)
                            {
                                // Move on to next index and never come back!
                                currentTheoreticalIndex++;
                                if (currentTheoreticalIndex == TotalProductsHere)
                                    break;
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
                                        break;
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
        }

        public static double CalculatePeptideScore(IMsDataScan<IMzSpectrum<IMzPeak>> thisScan, Tolerance productMassTolerance, double[] sortedTheoreticalProductMassesForThisPeptide, double precursorMass, List<DissociationType> dissociationTypes, bool addCompIons)
        {
            var TotalProductsHere = sortedTheoreticalProductMassesForThisPeptide.Length;
            if (TotalProductsHere == 0)
                return 0;
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
                return 0;

            double currentTheoreticalMz = currentTheoreticalMass + Constants.protonMass;
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
                    MatchingIntensityHere += experimental_intensities[experimentalIndex];

                    currentTheoreticalIndex++; //prevent multi counting
                    if (currentTheoreticalIndex == TotalProductsHere)
                        break;
                    currentTheoreticalMz = sortedTheoreticalProductMassesForThisPeptide[currentTheoreticalIndex] + Constants.protonMass;
                }
                // Else if for sure did not reach the next theoretical yet
                else if (currentExperimentalMz > currentTheoreticalMz)
                {
                    // Move on to next index and never come back!
                    currentTheoreticalIndex++;
                    if (currentTheoreticalIndex == TotalProductsHere)
                        break;
                    currentTheoreticalMz = sortedTheoreticalProductMassesForThisPeptide[currentTheoreticalIndex] + Constants.protonMass;

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
                            break;
                        testTheoreticalMz = sortedTheoreticalProductMassesForThisPeptide[testTheoreticalIndex] + Constants.protonMass;
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
                    if (complementaryIonConversionDictionary.TryGetValue(dissociationType, out double protonMassShift))
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
                                MatchingIntensityHere += complementaryIntensities[experimentalIndex];

                                currentTheoreticalIndex++;
                                if (currentTheoreticalIndex == TotalProductsHere)
                                    break;
                                currentTheoreticalMass = sortedTheoreticalProductMassesForThisPeptide[currentTheoreticalIndex];
                            }
                            // Else if for sure passed a theoretical
                            else if (currentExperimentalMass > currentTheoreticalMass)
                            {
                                // Move on to next index and never come back!
                                currentTheoreticalIndex++;
                                if (currentTheoreticalIndex == TotalProductsHere)
                                    break;
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
                                        break;
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

        #endregion Public Methods

        #region Protected Methods

        protected static List<DissociationType> DetermineDissociationType(List<ProductType> lp)
        {
            List<DissociationType> dissociationTypes = new List<DissociationType>();

            if (lp.Contains(ProductType.B) || lp.Contains(ProductType.Y))
                dissociationTypes.Add(DissociationType.HCD);

            if (lp.Contains(ProductType.C) || lp.Contains(ProductType.Zdot))
                dissociationTypes.Add(DissociationType.ETD);

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

        #endregion Protected Methods

        #region Private Methods

        private void StartingSingleEngine()
        {
            StartingSingleEngineHander?.Invoke(this, new SingleEngineEventArgs(this));
        }

        private void FinishedSingleEngine(MetaMorpheusEngineResults myResults)
        {
            FinishedSingleEngineHandler?.Invoke(this, new SingleEngineFinishedEventArgs(myResults));
        }

        #endregion Private Methods
    }
}