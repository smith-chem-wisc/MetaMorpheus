using Chemistry;
using MassSpectrometry;
using Spectra;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace InternalLogicEngineLayer
{
    public class ModernSearchEngine : MyEngine
    {
        #region Private Fields

        private readonly List<int>[] fragmentIndex;

        private readonly double fragmentToleranceInDaltons;

        private readonly float[] keys;

        private readonly IMsDataFile<IMzSpectrum<MzPeak>> myMsDataFile;

        private readonly List<CompactPeptide> peptideIndex;

        private readonly List<SearchMode> searchModes;

        private readonly int spectraFileIndex;

        #endregion Private Fields

        #region Public Constructors

        public ModernSearchEngine(IMsDataFile<IMzSpectrum<MzPeak>> myMsDataFile, int spectraFileIndex, List<CompactPeptide> peptideIndex, float[] keys, List<int>[] fragmentIndex, double fragmentToleranceInDaltons, List<SearchMode> searchModes) : base(2)
        {
            this.myMsDataFile = myMsDataFile;
            this.spectraFileIndex = spectraFileIndex;
            this.peptideIndex = peptideIndex;
            this.keys = keys;
            this.fragmentIndex = fragmentIndex;
            this.fragmentToleranceInDaltons = fragmentToleranceInDaltons;
            this.searchModes = searchModes;
        }

        #endregion Public Constructors

        #region Protected Methods

        protected override MyResults RunSpecific()
        {
            var totalSpectra = myMsDataFile.NumSpectra;

            List<ModernSpectrumMatch>[] newPsms = new List<ModernSpectrumMatch>[searchModes.Count];
            for (int i = 0; i < searchModes.Count; i++)
                newPsms[i] = new List<ModernSpectrumMatch>(new ModernSpectrumMatch[totalSpectra]);

            var listOfSortedms2Scans = myMsDataFile.Where(b => b.MsnOrder == 2).Select(b => new LocalMs2Scan(b)).OrderBy(b => b.precursorMass).ToArray();
            var listOfSortedms2ScansLength = listOfSortedms2Scans.Length;
            var searchModesCount = searchModes.Count;
            var outputObject = new object();
            int scansSeen = 0;
            int old_progress = 0;
            var peptideIndexCount = peptideIndex.Count;
            Parallel.ForEach(Partitioner.Create(0, listOfSortedms2ScansLength), fff =>
            {
                CompactPeptide[] bestPeptides = new CompactPeptide[searchModesCount];
                double[] bestScores = new double[searchModesCount];
                double[] fullPeptideScores = new double[peptideIndexCount];
                for (int i = fff.Item1; i < fff.Item2; i++)
                {
                    var thisScan = listOfSortedms2Scans[i];
                    var thisScanprecursorMass = thisScan.precursorMass;
                    Array.Clear(fullPeptideScores, 0, peptideIndexCount);
                    CalculatePeptideScores(thisScan.theScan, fullPeptideScores);

                    Array.Clear(bestPeptides, 0, searchModesCount);
                    Array.Clear(bestScores, 0, searchModesCount);

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
                                    if (searchMode.Accepts(thisScanprecursorMass, candidatePeptide.MonoisotopicMass) && FirstIsPreferableWithoutScore(candidatePeptide, bestPeptides[j], thisScanprecursorMass))
                                    {
                                        bestPeptides[j] = candidatePeptide;
                                        bestScores[j] = consideredScore;
                                    }
                                }
                                else if (currentBestScore < consideredScore)
                                {
                                    // Score is better, only make sure it is acceptable
                                    if (searchMode.Accepts(thisScanprecursorMass, candidatePeptide.MonoisotopicMass))
                                    {
                                        bestPeptides[j] = candidatePeptide;
                                        bestScores[j] = consideredScore;
                                    }
                                }
                            }
                            // Did not exist! Only make sure that it is acceptable
                            else if (searchMode.Accepts(thisScanprecursorMass, candidatePeptide.MonoisotopicMass))
                            {
                                bestPeptides[j] = candidatePeptide;
                                bestScores[j] = consideredScore;
                            }
                        }
                    }
                    for (int j = 0; j < searchModesCount; j++)
                    {
                        CompactPeptide theBestPeptide = bestPeptides[j];
                        if (theBestPeptide != null)
                        {
                            var cool = new ModernSpectrumMatch(thisScan.monoisotopicPrecursorMZ, thisScan.OneBasedScanNumber, thisScan.RetentionTime, thisScan.monoisotopicPrecursorCharge, thisScan.NumPeaks, thisScan.TotalIonCurrent, thisScan.monoisotopicPrecursorIntensity, spectraFileIndex, theBestPeptide, bestScores[j]);
                            newPsms[j][thisScan.OneBasedScanNumber - 1] = cool;
                        }
                    }
                }
                lock (outputObject)
                {
                    scansSeen += fff.Item2 - fff.Item1;
                    var new_progress = (int)((double)scansSeen / (listOfSortedms2ScansLength) * 100);
                    if (new_progress > old_progress)
                    {
                        ReportProgress(new ProgressEventArgs(new_progress, "In modern search loop"));
                        old_progress = new_progress;
                    }
                }
            });
            return new ModernSearchResults(newPsms, this);
        }

        #endregion Protected Methods

        #region Private Methods

        // Want this to return false more!! So less computation is done. So second is preferable more often.
        private static bool FirstIsPreferableWithoutScore(CompactPeptide first, CompactPeptide second, double pm)
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

        private void CalculatePeptideScores(IMsDataScan<IMzSpectrum<MzPeak>> spectrum, double[] peptideScores)
        {
            foreach (var experimentalPeak in spectrum.MassSpectrum)
            {
                var theAdd = 1 + experimentalPeak.Intensity / spectrum.TotalIonCurrent;
                var experimentalPeakInDaltons = experimentalPeak.MZ - Constants.ProtonMass;
                float closestPeak = float.NaN;
                var ipos = Array.BinarySearch(keys, (float)experimentalPeakInDaltons);
                if (ipos < 0)
                    ipos = ~ipos;

                //po.out(" ipos " + ipos);
                if (ipos > 0)
                {
                    var downIpos = ipos - 1;
                    // Try down
                    while (downIpos >= 0)
                    {
                        closestPeak = keys[downIpos];
                        // po.out("  closestPeak "+ closestPeak);
                        if (Math.Abs(closestPeak - experimentalPeakInDaltons) < fragmentToleranceInDaltons)
                        {// po.out("    ********************************");
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
                        //po.out("  closestPeak " + closestPeak);
                        if (Math.Abs(closestPeak - experimentalPeakInDaltons) < fragmentToleranceInDaltons)
                        {
                            //po.out("    ********************************");
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