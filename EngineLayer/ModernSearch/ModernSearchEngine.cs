using Chemistry;
using MassSpectrometry;
using MzLibUtil;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Threading.Tasks;

namespace EngineLayer.ModernSearch
{
    public class ModernSearchEngine : MetaMorpheusEngine
    {

        #region Private Fields

        private const double tolInDaForPreferringHavingMods = 0.03;

        private readonly List<int>[] fragmentIndex;

        private readonly Tolerance fragmentTolerance;

        private readonly float[] keys;

        private readonly Ms2ScanWithSpecificMass[] listOfSortedms2Scans;

        private readonly List<CompactPeptide> peptideIndex;

        private readonly List<MassDiffAcceptor> searchModes;

        #endregion Private Fields

        #region Public Constructors

        public ModernSearchEngine(Ms2ScanWithSpecificMass[] listOfSortedms2Scans, List<CompactPeptide> peptideIndex, float[] keys, List<int>[] fragmentIndex, Tolerance fragmentTolerance, List<MassDiffAcceptor> searchModes, List<string> nestedIds) : base(nestedIds)
        {
            this.listOfSortedms2Scans = listOfSortedms2Scans;
            this.peptideIndex = peptideIndex;
            this.keys = keys;
            this.fragmentIndex = fragmentIndex;
            this.fragmentTolerance = fragmentTolerance;
            this.searchModes = searchModes;
        }

        #endregion Public Constructors

        #region Protected Methods

        protected override MetaMorpheusEngineResults RunSpecific()
        {
            Status("In modern search engine...", nestedIds);

            var listOfSortedms2ScansLength = listOfSortedms2Scans.Length;
            PsmParent[][] newPsms = new PsmParent[searchModes.Count][];
            for (int i = 0; i < searchModes.Count; i++)
                newPsms[i] = new PsmParent[listOfSortedms2Scans.Length];

            var searchModesCount = searchModes.Count;
            var outputObject = new object();
            int scansSeen = 0;
            int old_progress = 0;
            var peptideIndexCount = peptideIndex.Count;
            Parallel.ForEach(Partitioner.Create(0, listOfSortedms2ScansLength), fff =>
            {
                CompactPeptide[] bestPeptides = new CompactPeptide[searchModesCount];
                double[] bestScores = new double[searchModesCount];
                int[] bestNotches = new int[searchModesCount];
                double[] fullPeptideScores = new double[peptideIndexCount];
                for (int i = fff.Item1; i < fff.Item2; i++)
                {
                    var thisScan = listOfSortedms2Scans[i];
                    var thisScanprecursorMass = thisScan.PrecursorMass;
                    Array.Clear(fullPeptideScores, 0, peptideIndexCount);
                    CalculatePeptideScores(thisScan.TheScan, fullPeptideScores);

                    Array.Clear(bestPeptides, 0, searchModesCount);
                    Array.Clear(bestScores, 0, searchModesCount);
                    Array.Clear(bestNotches, 0, searchModesCount);

                    for (int possibleWinningPeptideIndex = 0; possibleWinningPeptideIndex < fullPeptideScores.Length; possibleWinningPeptideIndex++)
                    {
                        var consideredScore = fullPeptideScores[possibleWinningPeptideIndex];
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
                                    int notch = searchMode.Accepts(thisScanprecursorMass, candidatePeptide.MonoisotopicMass);
                                    if (notch >= 0 && FirstIsPreferableWithoutScore(candidatePeptide, bestPeptides[j], thisScanprecursorMass))
                                    {
                                        bestPeptides[j] = candidatePeptide;
                                        bestScores[j] = consideredScore;
                                        bestNotches[j] = notch;
                                    }
                                }
                                else if (currentBestScore < consideredScore)
                                {
                                    // Score is better, only make sure it is acceptable
                                    int notch = searchMode.Accepts(thisScanprecursorMass, candidatePeptide.MonoisotopicMass);
                                    if (notch >= 0)
                                    {
                                        bestPeptides[j] = candidatePeptide;
                                        bestScores[j] = consideredScore;
                                        bestNotches[j] = notch;
                                    }
                                }
                            }
                            // Did not exist! Only make sure that it is acceptable
                            else
                            {
                                int notch = searchMode.Accepts(thisScanprecursorMass, candidatePeptide.MonoisotopicMass);
                                if (notch >= 0)
                                {
                                    bestPeptides[j] = candidatePeptide;
                                    bestScores[j] = consideredScore;
                                    bestNotches[j] = notch;
                                }
                            }
                        }
                    }
                    for (int j = 0; j < searchModesCount; j++)
                    {
                        CompactPeptide theBestPeptide = bestPeptides[j];
                        if (theBestPeptide != null)
                        {
                            newPsms[j][i] = new PsmParent(bestNotches[j], bestScores[j], i, thisScan);
                            newPsms[j][i].Add(theBestPeptide);
                        }
                    }
                }
                lock (outputObject)
                {
                    scansSeen += fff.Item2 - fff.Item1;
                    var new_progress = (int)((double)scansSeen / (listOfSortedms2ScansLength) * 100);
                    if (new_progress > old_progress)
                    {
                        ReportProgress(new ProgressEventArgs(new_progress, "In modern search loop", nestedIds));
                        old_progress = new_progress;
                    }
                }
            });
            return new SearchResults(newPsms, this);
        }

        #endregion Protected Methods

        #region Private Methods

        // Want this to return false more!! So less computation is done. So second is preferable more often.
        private static bool FirstIsPreferableWithoutScore(CompactPeptide first, CompactPeptide second, double pm)
        {
            if (Math.Abs(first.MonoisotopicMass - pm) < tolInDaForPreferringHavingMods && Math.Abs(second.MonoisotopicMass - pm) > tolInDaForPreferringHavingMods)
                return true;
            if (Math.Abs(first.MonoisotopicMass - pm) > tolInDaForPreferringHavingMods && Math.Abs(second.MonoisotopicMass - pm) < tolInDaForPreferringHavingMods)
                return false;

            return false;
        }

        private void CalculatePeptideScores(IMsDataScan<IMzSpectrum<IMzPeak>> spectrum, double[] peptideScores)
        {
            for (int i = 0; i < spectrum.MassSpectrum.Size; i++)
            {
                var theAdd = 1 + spectrum.MassSpectrum[i].Intensity / spectrum.TotalIonCurrent;
                var experimentalPeakInDaltons = spectrum.MassSpectrum[i].Mz - Constants.protonMass;
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
        }

        #endregion Private Methods

    }
}