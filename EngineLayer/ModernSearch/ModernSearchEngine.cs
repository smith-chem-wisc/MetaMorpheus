﻿using Chemistry;
using MassSpectrometry;
using MzLibUtil;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace EngineLayer.ModernSearch
{
    public class ModernSearchEngine : MetaMorpheusEngine
    {
        #region Protected Fields

        protected readonly List<int>[] fragmentIndex;
        protected readonly Tolerance fragmentTolerance;
        protected readonly float[] keys;
        protected readonly Psm[][] globalPsms;
        protected readonly Ms2ScanWithSpecificMass[] listOfSortedms2Scans;
        protected readonly List<CompactPeptide> peptideIndex;
        protected readonly List<MassDiffAcceptor> searchModes;
        protected readonly bool addCompIons;
        protected readonly List<ProductType> lp;
        protected readonly double scoreCutoff;
        protected readonly int currentPartition;
        protected readonly int totalPartitions;

        #endregion Protected Fields

        #region Private Fields

        private const double tolInDaForPreferringHavingMods = 0.03;

        #endregion Private Fields

        #region Public Constructors

        public ModernSearchEngine(Psm[][] globalPsms, Ms2ScanWithSpecificMass[] listOfSortedms2Scans, List<CompactPeptide> peptideIndex, float[] keys, List<int>[] fragmentIndex, Tolerance fragmentTolerance, List<MassDiffAcceptor> searchModes, List<string> nestedIds, bool addCompIons, List<ProductType> lp, double scoreCutoff, int currentPartition, int totalPartitions) : base(nestedIds)
        {
            this.globalPsms = globalPsms;
            this.listOfSortedms2Scans = listOfSortedms2Scans;
            this.peptideIndex = peptideIndex;
            this.keys = keys;
            this.fragmentIndex = fragmentIndex;
            this.fragmentTolerance = fragmentTolerance;
            this.searchModes = searchModes;
            this.addCompIons = addCompIons;
            this.lp = lp;
            this.scoreCutoff = scoreCutoff;
            this.currentPartition = currentPartition + 1;
            this.totalPartitions = totalPartitions;
        }

        #endregion Public Constructors

        #region Protected Methods

        protected override MetaMorpheusEngineResults RunSpecific()
        {
            Status("In modern search engine..." + currentPartition + "/" + totalPartitions, nestedIds);

            var listOfSortedms2ScansLength = listOfSortedms2Scans.Length;

            var searchModesCount = searchModes.Count;
            var outputObject = new object();
            int scansSeen = 0;
            int old_progress = 0;
            var peptideIndexCount = peptideIndex.Count;
            Parallel.ForEach(Partitioner.Create(0, listOfSortedms2ScansLength), fff =>
            {
                List<CompactPeptide>[] bestPeptides = new List<CompactPeptide>[searchModesCount];
                double[] bestScores = new double[searchModesCount];
                List<int>[] bestNotches = new List<int>[searchModesCount];
                double[] fullPeptideScores = new double[peptideIndexCount];
                for (int i = fff.Item1; i < fff.Item2; i++)
                {
                    var thisScan = listOfSortedms2Scans[i];
                    var thisScanprecursorMass = thisScan.PrecursorMass;
                    Array.Clear(fullPeptideScores, 0, peptideIndexCount);
                    double thePrecursorMass = thisScan.PrecursorMass;
                    CalculatePeptideScores(thisScan.TheScan, fullPeptideScores, thePrecursorMass);

                    Array.Clear(bestPeptides, 0, searchModesCount);
                    Array.Clear(bestScores, 0, searchModesCount);
                    Array.Clear(bestNotches, 0, searchModesCount);

                    for (int possibleWinningPeptideIndex = 0; possibleWinningPeptideIndex < fullPeptideScores.Length; possibleWinningPeptideIndex++)
                    {
                        var consideredScore = fullPeptideScores[possibleWinningPeptideIndex];
                        if (consideredScore > scoreCutoff) //intentionally high. 99.9% of 4-mers are present in a given UniProt database. This saves considerable time
                        {
                            CompactPeptide candidatePeptide = peptideIndex[possibleWinningPeptideIndex];
                            for (int j = 0; j < searchModesCount; j++)
                            {
                                // Check if makes sense to add due to peptidescore!
                                var searchMode = searchModes[j];
                                double currentBestScore = bestScores[j];
                                if (currentBestScore > 1)
                                {
                                    // Existed! Need to compare with old match
                                    if (Math.Abs(currentBestScore - consideredScore) < 1e-9)
                                    {
                                        // Score is same, need to see if accepts and if prefer the new one
                                        int notch = searchMode.Accepts(thisScanprecursorMass, candidatePeptide.MonoisotopicMassIncludingFixedMods);
                                        if (notch >= 0)
                                        {
                                            bestPeptides[j].Add(candidatePeptide);
                                            bestNotches[j].Add(notch);
                                        }
                                    }
                                    else if (currentBestScore < consideredScore)
                                    {
                                        // Score is better, only make sure it is acceptable
                                        int notch = searchMode.Accepts(thisScanprecursorMass, candidatePeptide.MonoisotopicMassIncludingFixedMods);
                                        if (notch >= 0)
                                        {
                                            bestPeptides[j] = new List<CompactPeptide> { candidatePeptide };
                                            bestScores[j] = consideredScore;
                                            bestNotches[j] = new List<int> { notch };
                                        }
                                    }
                                }
                                // Did not exist! Only make sure that it is acceptable
                                else
                                {
                                    int notch = searchMode.Accepts(thisScanprecursorMass, candidatePeptide.MonoisotopicMassIncludingFixedMods);
                                    if (notch >= 0)
                                    {
                                        bestPeptides[j] = new List<CompactPeptide> { candidatePeptide };
                                        bestScores[j] = consideredScore;
                                        bestNotches[j] = new List<int> { notch };
                                    }
                                }
                            }
                        }
                    }
                    for (int j = 0; j < searchModesCount; j++)
                    {
                        if (bestPeptides[j] != null)
                        {
                            globalPsms[j][i] = new Psm(bestPeptides[j][0], bestNotches[j][0], bestScores[j], i, thisScan);
                            for (int k = 1; k < bestPeptides[j].Count; k++)
                            {
                                globalPsms[j][i].AddOrReplace(bestPeptides[j][k], bestScores[j], bestNotches[j][k]);
                            }
                        }
                    }
                }
                lock (outputObject)
                {
                    scansSeen += fff.Item2 - fff.Item1;
                    var new_progress = (int)((double)scansSeen / (listOfSortedms2ScansLength) * 100);
                    if (new_progress > old_progress)
                    {
                        ReportProgress(new ProgressEventArgs(new_progress, "In modern search loop" + currentPartition + "/" + totalPartitions, nestedIds));
                        old_progress = new_progress;
                    }
                }
            });
            return new MetaMorpheusEngineResults(this);
        }

        #endregion Protected Methods

        #region Private Methods

        private void CalculatePeptideScores(IMsDataScan<IMzSpectrum<IMzPeak>> spectrum, double[] peptideScores, double thePrecursorMass)
        {
            if (!addCompIons)
            {
                for (int i = 0; i < spectrum.MassSpectrum.Size; i++)
                {
                    var theAdd = 1 + spectrum.MassSpectrum[i].Intensity / spectrum.TotalIonCurrent;
                    var experimentalPeakInDaltons = spectrum.MassSpectrum[i].Mz - Constants.protonMass;
                    GeneratePeptideScores(theAdd, experimentalPeakInDaltons, peptideScores);
                }
            }
            else
            {
                List<IMzPeak> experimentalPeaks = new List<IMzPeak>();
                for (int i = 0; i < spectrum.MassSpectrum.Size; i++)
                {
                    experimentalPeaks.Add(spectrum.MassSpectrum[i]);
                }

                //If HCD
                if (lp.Contains(ProductType.B) || lp.Contains(ProductType.Y))
                {
                    for (int i = 0; i < spectrum.MassSpectrum.Size; i++)
                    {
                        experimentalPeaks.Add(new MzPeak((thePrecursorMass - spectrum.MassSpectrum[i].Mz + Constants.protonMass * 2), (spectrum.MassSpectrum[i].Intensity / 100)));
                    }
                }
                //If ETD
                if (lp.Contains(ProductType.C) || lp.Contains(ProductType.Zdot))
                {
                    for (int i = 0; i < spectrum.MassSpectrum.Size; i++)
                    {
                        experimentalPeaks.Add(new MzPeak((thePrecursorMass - spectrum.MassSpectrum[i].Mz + Constants.protonMass * 3), (spectrum.MassSpectrum[i].Intensity / 100)));
                    }
                }

                IEnumerable<IMzPeak> sortedPeaksMZ = experimentalPeaks.OrderBy(x => x.Mz);

                foreach (IMzPeak experimentalPeak in sortedPeaksMZ)
                {
                    var theAdd = 1 + experimentalPeak.Intensity / spectrum.TotalIonCurrent;
                    var experimentalPeakInDaltons = experimentalPeak.Mz - Constants.protonMass;
                    GeneratePeptideScores(theAdd, experimentalPeakInDaltons, peptideScores);
                }
            }
        }

        private void GeneratePeptideScores(double theAdd, double experimentalPeakInDaltons, double[] peptideScores)
        {
            float closestPeak;
            var ipos = Array.BinarySearch(keys, (float)experimentalPeakInDaltons);
            if (ipos < 0)
                ipos = ~ipos;

            if (ipos > 0)
            {
                var downIpos = ipos - 1;
                // Try down
                while (downIpos >= 0)
                {
                    closestPeak = keys[downIpos];
                    if (fragmentTolerance.Within(experimentalPeakInDaltons, closestPeak))
                    {
                        foreach (var heh in fragmentIndex[downIpos])
                            peptideScores[heh] += theAdd;
                    }
                    else
                        break;
                    downIpos--;
                }
            }
            if (ipos < keys.Length)
            {
                var upIpos = ipos;
                // Try here and up
                while (upIpos < keys.Length)
                {
                    closestPeak = keys[upIpos];
                    if (fragmentTolerance.Within(experimentalPeakInDaltons, closestPeak))
                    {
                        foreach (var heh in fragmentIndex[upIpos])
                            peptideScores[heh] += theAdd;
                    }
                    else
                        break;
                    upIpos++;
                }
            }
        }

        #endregion Private Methods
    }
}