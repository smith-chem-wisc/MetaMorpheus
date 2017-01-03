using Chemistry;
using MassSpectrometry;
using MathNet.Numerics.Statistics;
using Spectra;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace IndexSearchAndAnalyze
{
    public class SearchEngine : MyEngine
    {
        private SearchParams searchParams;

        public SearchEngine(SearchParams searchParams)
        {
            this.searchParams = searchParams;
        }

        protected override MyResults RunSpecific()
        {
            searchParams.OnOutput("In search method!");

            var spectraList = searchParams.myMsDataFile.ToList();
            var totalSpectra = searchParams.myMsDataFile.NumSpectra;

            List<NewPsm>[] newPsms = new List<NewPsm>[searchParams.searchModes.Count];
            for (int i = 0; i < searchParams.searchModes.Count; i++)
                newPsms[i] = new List<NewPsm>(new NewPsm[totalSpectra]);

            int numAllSpectra = 0;
            int numMS2spectra = 0;
            int[] numMS2spectraMatched = new int[searchParams.searchModes.Count];

            var searchModesCount = searchParams.searchModes.Count;

            Parallel.ForEach(Partitioner.Create(0, totalSpectra), fff =>
            {
                for (int i = fff.Item1; i < fff.Item2; i++)
                {
                    var thisScan = spectraList[i];
                    if (thisScan.MsnOrder == 2)
                    {
                        double selectedMZ;
                        int selectedCharge;
                        thisScan.TryGetSelectedIonGuessMonoisotopicMZ(out selectedMZ);
                        thisScan.TryGetSelectedIonGuessChargeStateGuess(out selectedCharge);
                        var scanPrecursorMass = selectedMZ.ToMass(selectedCharge);

                        var fullPeptideScores = CalculatePeptideScores(thisScan, searchParams.peptideIndex, 400, searchParams.keys, searchParams.fragmentIndex, searchParams.fragmentTolerance);

                        CompactPeptide[] bestPeptides = new CompactPeptide[searchModesCount];
                        double[] bestScores = new double[searchModesCount];
                        for (int possibleWinningPeptideIndex = 0; possibleWinningPeptideIndex < fullPeptideScores.Length; possibleWinningPeptideIndex++)
                        {
                            var consideredScore = fullPeptideScores[possibleWinningPeptideIndex];
                            CompactPeptide candidatePeptide = searchParams.peptideIndex[possibleWinningPeptideIndex];
                            for (int j = 0; j < searchModesCount; j++)
                            {
                                // Check if makes sense to add due to peptidescore!
                                var searchMode = searchParams.searchModes[j];
                                double currentBestScore = bestScores[j];
                                if (currentBestScore > 1)
                                {
                                    // Existed! Need to compare with old match
                                    if (Math.Abs(currentBestScore - consideredScore) < 1e-9)
                                    {
                                        // Score is same, need to see if accepts and if prefer the new one
                                        if (searchMode.Accepts(scanPrecursorMass - candidatePeptide.MonoisotopicMass) && FirstIsPreferable(candidatePeptide, bestPeptides[j], scanPrecursorMass))
                                        {
                                            bestPeptides[j] = candidatePeptide;
                                            bestScores[j] = consideredScore;
                                        }
                                    }
                                    else if (currentBestScore < consideredScore)
                                    {
                                        // Score is better, only make sure it is acceptable
                                        if (searchMode.Accepts(scanPrecursorMass - candidatePeptide.MonoisotopicMass))
                                        {
                                            bestPeptides[j] = candidatePeptide;
                                            bestScores[j] = consideredScore;
                                        }
                                    }
                                }
                                else
                                {
                                    // Did not exist! Only make sure that it is acceptable
                                    if (searchMode.Accepts(scanPrecursorMass - candidatePeptide.MonoisotopicMass))
                                    {
                                        bestPeptides[j] = candidatePeptide;
                                        bestScores[j] = consideredScore;
                                    }
                                }
                            }
                        }

                        var psms = new NewPsm[searchModesCount];

                        for (int j = 0; j < searchModesCount; j++)
                        {
                            CompactPeptide theBestPeptide = bestPeptides[j];
                            if (theBestPeptide != null)
                            {
                                newPsms[j][thisScan.OneBasedScanNumber - 1] = new NewPsm(thisScan, searchParams.spectraFileIndex, theBestPeptide, bestScores[j]);
                                //numMS2spectraMatched[j]++;
                            }
                        }
                        //numMS2spectra++;
                    }
                    //numAllSpectra++;
                    //if (numAllSpectra % 100 == 0)
                    //    Console.WriteLine("Spectra: " + numAllSpectra + " / " + totalSpectra);
                }
            });
            return new SearchResults(newPsms, numMS2spectra, numMS2spectraMatched, searchParams);
        }

        private static float[] CalculatePeptideScores(IMsDataScan<IMzSpectrum<MzPeak>> spectrum, List<CompactPeptide> peptides, int maxPeaks, float[] fragmentMassesAscending, List<int>[] fragmentIndex, double fragmentTolerance)
        {
            List<MzPeak> filteredList;
            if (spectrum.MassSpectrum.Count <= maxPeaks)
                filteredList = spectrum.MassSpectrum.ToList();
            else
            {
                var cutoffIntensity = spectrum.MassSpectrum.yArray.Quantile(1.0 - (double)maxPeaks / spectrum.MassSpectrum.Count);
                filteredList = spectrum.MassSpectrum.Where(b => b.Intensity > cutoffIntensity).ToList();
            }
            float[] peptideScores = new float[peptides.Count];
            foreach (var experimentalPeak in filteredList)
            {
                var experimentalPeakInDaltons = experimentalPeak.MZ.ToMass(1);
                float closestPeak = float.NaN;
                var ipos = Array.BinarySearch(fragmentMassesAscending, (float)experimentalPeakInDaltons);
                if (ipos < 0)
                    ipos = ~ipos;

                //Console.WriteLine(" ipos " + ipos);
                if (ipos > 0)
                {
                    var downIpos = ipos - 1;
                    // Try down
                    while (downIpos >= 0)
                    {
                        closestPeak = fragmentMassesAscending[downIpos];
                        // Console.WriteLine("  closestPeak "+ closestPeak);
                        if (Math.Abs(closestPeak - experimentalPeakInDaltons) < fragmentTolerance)
                        {// Console.WriteLine("    ********************************");
                            foreach (var heh in fragmentIndex[downIpos])
                                peptideScores[heh] += (float)(1 + experimentalPeak.Intensity / spectrum.TotalIonCurrent);
                        }
                        else
                            break;
                        downIpos--;
                    }
                }
                if (ipos < fragmentMassesAscending.Length)
                {
                    var upIpos = ipos;
                    // Try here and up
                    while (upIpos < fragmentMassesAscending.Length)
                    {
                        closestPeak = fragmentMassesAscending[upIpos];
                        //Console.WriteLine("  closestPeak " + closestPeak);
                        if (Math.Abs(closestPeak - experimentalPeakInDaltons) < fragmentTolerance)
                        {
                            //Console.WriteLine("    ********************************");
                            foreach (var heh in fragmentIndex[upIpos])
                                peptideScores[heh] += (float)(1 + experimentalPeak.Intensity / spectrum.TotalIonCurrent);
                        }
                        else
                            break;
                        upIpos++;
                    }
                }
            }
            return peptideScores;
        }

        // Want this to return false more!! So less computation is done
        private static bool FirstIsPreferable(CompactPeptide first, CompactPeptide second, double pm)
        {
            if (Math.Abs(first.MonoisotopicMass - pm) < 0.5 && Math.Abs(second.MonoisotopicMass - pm) > 0.5)
                return true;
            if (Math.Abs(first.MonoisotopicMass - pm) > 0.5 && Math.Abs(second.MonoisotopicMass - pm) < 0.5)
                return false;

            if (first.varMod1Type == 0 && second.varMod1Type > 0)
                return true;
            if (first.varMod1Type > 0 && second.varMod1Type == 0)
                return false;
            if (first.varMod2Type == 0 && second.varMod2Type > 0)
                return true;
            if (first.varMod2Type > 0 && second.varMod2Type == 0)
                return false;
            if (first.varMod3Type == 0 && second.varMod3Type > 0)
                return true;
            if (first.varMod3Type > 0 && second.varMod3Type == 0)
                return false;

            return false;
        }
    }
}