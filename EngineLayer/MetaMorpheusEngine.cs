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
        #region Protected Fields

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

        public static void MatchIons(IMsDataScan<IMzSpectrum<IMzPeak>> thisScan, Tolerance productMassTolerance, double[] sortedTheoreticalProductMassesForThisPeptide, List<double> matchedIonMassesList, List<double> productMassErrorDa, List<double> productMassErrorPpm, bool addComp, double precursorMass, List<ProductType> lp)
        {
            var TotalProductsHere = sortedTheoreticalProductMassesForThisPeptide.Length;
            if (TotalProductsHere == 0)
                return;
            int MatchingProductsHere = 0;
            double MatchingIntensityHere = 0;

            // speed optimizations
            double[] experimental_mzs = thisScan.MassSpectrum.XArray;
            double[] experimental_intensities = thisScan.MassSpectrum.YArray;

            if(addComp)
                AddComplementaryPeaks(ref experimental_mzs, ref experimental_intensities, precursorMass, lp);

            int numExperimentalPeaks = experimental_mzs.Length;

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
            double testTheoreticalMz;
            double testTheoreticalMass;
            // Loop over all experimental indices
            for (int experimentalIndex = 0; experimentalIndex < numExperimentalPeaks; experimentalIndex++)
            {
                double currentExperimentalMz = experimental_mzs[experimentalIndex];
                // If found match
                if (productMassTolerance.Within(currentExperimentalMz, currentTheoreticalMz))
                {
                    MatchingProductsHere++;
                    MatchingIntensityHere += experimental_intensities[experimentalIndex];

                    matchedIonMassesList.Add(currentTheoreticalMass);
                    productMassErrorDa.Add(currentExperimentalMz - currentTheoreticalMz);
                    productMassErrorPpm.Add((currentExperimentalMz - currentTheoreticalMz) * 1000000 / currentTheoreticalMz);

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
                    testTheoreticalMz = currentTheoreticalMz;
                    testTheoreticalMass = currentTheoreticalMass;
                    // Mark the skipped theoreticals as not found. The last one is not for sure, might be flipped!
                    while (currentExperimentalMz > testTheoreticalMz)
                    {
                        // Store old info for possible reuse
                        currentTheoreticalMass = testTheoreticalMass;
                        currentTheoreticalMz = testTheoreticalMz;
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
            if (addComp)
                MatchIons(thisScan, productMassTolerance, sortedTheoreticalProductMassesForThisPeptide, matchedIonMassesList, productMassErrorDa, productMassErrorPpm, false, precursorMass, lp);
        }

        public static double CalculateClassicScore(IMsDataScan<IMzSpectrum<IMzPeak>> thisScan, Tolerance productMassTolerance, double[] sortedTheoreticalProductMassesForThisPeptide, bool addComp, double precursorMass, List<ProductType> lp)
        {
            var TotalProductsHere = sortedTheoreticalProductMassesForThisPeptide.Length;
            if (TotalProductsHere == 0)
                return 0;
            int MatchingProductsHere = 0;
            double MatchingIntensityHere = 0;
            
            // speed optimizations
            double[] experimental_mzs = thisScan.MassSpectrum.XArray;
            double[] experimental_intensities = thisScan.MassSpectrum.YArray;

            if (addComp)
                AddComplementaryPeaks(ref experimental_mzs, ref experimental_intensities, precursorMass, lp);

            int numExperimentalPeaks = experimental_mzs.Length;

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
            double testTheoreticalMass;
            // Loop over all experimental indices
            for (int experimentalIndex = 0; experimentalIndex < numExperimentalPeaks; experimentalIndex++)
            {
                double currentExperimentalMz = experimental_mzs[experimentalIndex];
                // If found match
                if (productMassTolerance.Within(currentExperimentalMz, currentTheoreticalMz))
                {
                    MatchingProductsHere++;
                    MatchingIntensityHere += experimental_intensities[experimentalIndex];

                    currentTheoreticalIndex++;
                    if (currentTheoreticalIndex == TotalProductsHere)
                        break;
                    currentTheoreticalMass = sortedTheoreticalProductMassesForThisPeptide[currentTheoreticalIndex];
                    currentTheoreticalMz = currentTheoreticalMass + Constants.protonMass;
                }
                // Else if for sure passed a theoretical
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
                    testTheoreticalMz = currentTheoreticalMz;
                    testTheoreticalMass = currentTheoreticalMass;
                    // Mark the skipped theoreticals as not found. The last one is not for sure, might be flipped!
                    while (currentExperimentalMz > testTheoreticalMz)
                    {
                        // Store old info for possible reuse
                        currentTheoreticalMass = testTheoreticalMass;
                        currentTheoreticalMz = testTheoreticalMz;
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
            double addition = 0;
            if (addComp)
            {
                addition = CalculateClassicScore(thisScan, productMassTolerance, sortedTheoreticalProductMassesForThisPeptide, false, precursorMass, lp);
            }

            return (MatchingProductsHere + MatchingIntensityHere / thisScan.TotalIonCurrent)+addition;
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

        #endregion Public Methods

        #region Internal Methods

        internal string GetId()
        {
            return string.Join(",", nestedIds);
        }

        #endregion Internal Methods

        #region Protected Methods

        protected void Warn(string v, List<string> nestedIds)
        {
            WarnHandler?.Invoke(this, new StringEventArgs(v, nestedIds));
        }

        protected void Status(string v, List<string> nestedIds)
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

        private static void AddComplementaryPeaks(ref double[] experimental_mzs, ref double[] experimental_intensities, double precursorMass, List<ProductType> lp)
        {
            List<MzPeak> complementaryPeaks = new List<MzPeak>();

            //If HCD
            if (lp.Contains(ProductType.B) || lp.Contains(ProductType.Y))
            {
                for (int i = 0; i < experimental_mzs.Length; i++)
                {
                    complementaryPeaks.Add(new MzPeak((precursorMass - experimental_mzs[i] + Constants.protonMass * 2), (experimental_intensities[i] / 100)));
                }
            }
            //If ETD
            if (lp.Contains(ProductType.C) || lp.Contains(ProductType.Zdot))
            {
                for (int i = 0; i < experimental_mzs.Length; i++)
                {
                    complementaryPeaks.Add(new MzPeak((precursorMass - experimental_mzs[i] + Constants.protonMass * 3), (experimental_intensities[i] / 100)));
                }
            }

            IEnumerable<MzPeak> sortedPeaksMZ = complementaryPeaks.OrderBy(x => x.Mz);
            experimental_mzs = new double[sortedPeaksMZ.Count()];
            experimental_intensities = new double[sortedPeaksMZ.Count()];
            int index = 0;
            foreach (MzPeak peak in sortedPeaksMZ)
            {
                experimental_mzs[index] = peak.Mz;
                experimental_intensities[index] = peak.Intensity;
                index++;
            }
        }

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