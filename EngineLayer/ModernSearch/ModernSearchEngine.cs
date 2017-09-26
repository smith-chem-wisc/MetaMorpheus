using Chemistry;
using MassSpectrometry;
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
        protected readonly Psm[] globalPsms;
        protected readonly Ms2ScanWithSpecificMass[] listOfSortedms2Scans;
        protected readonly List<CompactPeptide> peptideIndex;
        protected readonly List<ProductType> lp;
        protected readonly int currentPartition;
        protected readonly CommonParameters CommonParameters;
        protected readonly bool addCompIons;
        protected readonly MassDiffAcceptor massDiffAcceptor;
        protected readonly List<DissociationType> dissociationTypes;
        protected static readonly int fragmentBinsPerDalton = 1000;

        #endregion Protected Fields

        #region Public Constructors

        public ModernSearchEngine(Psm[] globalPsms, Ms2ScanWithSpecificMass[] listOfSortedms2Scans, List<CompactPeptide> peptideIndex, List<int>[] fragmentIndex, List<ProductType> lp, int currentPartition, CommonParameters CommonParameters, bool addCompIons, MassDiffAcceptor massDiffAcceptors, List<string> nestedIds) : base(nestedIds)
        {
            this.globalPsms = globalPsms;
            this.listOfSortedms2Scans = listOfSortedms2Scans;
            this.peptideIndex = peptideIndex;
            this.fragmentIndex = fragmentIndex;
            this.lp = lp;
            this.currentPartition = currentPartition + 1;
            this.CommonParameters = CommonParameters;
            this.addCompIons = addCompIons;
            this.massDiffAcceptor = massDiffAcceptors;
            this.dissociationTypes = DetermineDissociationType(lp);
        }

        #endregion Public Constructors

        #region Protected Methods

        protected override MetaMorpheusEngineResults RunSpecific()
        {
            double progress = 0;
            int oldPercentProgress = 0;
            ReportProgress(new ProgressEventArgs(oldPercentProgress, "Performing modern search... " + currentPartition + "/" + CommonParameters.TotalPartitions, nestedIds));

            int intScoreCutoff = (int)CommonParameters.ScoreCutoff;
            //var roundedPepMasses = peptideIndex.Select(p => (int)Math.Round(p.MonoisotopicMassIncludingFixedMods * 100, 0)).ToList();
            byte byteScoreCutoff = Convert.ToByte(intScoreCutoff);

            Parallel.ForEach(Partitioner.Create(0, listOfSortedms2Scans.Length), new ParallelOptions { MaxDegreeOfParallelism = CommonParameters.MaxThreadsToUse }, range =>
            {
                byte[] scoringTable = new byte[peptideIndex.Count];
                HashSet<int> idsOfPeptidesPossiblyObserved = new HashSet<int>();

                for (int i = range.Item1; i < range.Item2; i++)
                {
                    // empty the scoring table to score the new scan (conserves memory compared to allocating a new array)
                    Array.Clear(scoringTable, 0, scoringTable.Length);
                    idsOfPeptidesPossiblyObserved.Clear();
                    var scan = listOfSortedms2Scans[i];

                    // filter ms2 fragment peaks by intensity
                    int numFragmentsToUse = 0;
                    if (CommonParameters.TopNpeaks != null)
                        numFragmentsToUse = (int)CommonParameters.TopNpeaks;
                    else
                        numFragmentsToUse = scan.NumPeaks;

                    var peaks = scan.TheScan.MassSpectrum.FilterByNumberOfMostIntense(numFragmentsToUse).ToList();
                    double largestIntensity = scan.TheScan.MassSpectrum.YofPeakWithHighestY;

                    // get allowed precursor masses
                    var t = massDiffAcceptor.GetAllowedPrecursorMassIntervals(scan.PrecursorMass);
                    double lowestMassPeptideToLookFor = t.Min(p => p.allowedInterval.Minimum);
                    double highestMassPeptideToLookFor = t.Max(p => p.allowedInterval.Maximum);

                    int numChecksSkipped = 1;
                    int obsPreviousFragmentCeilingMz = 0;

                    // search peaks for matches
                    foreach (IMzPeak peak in peaks)
                    {
                        if (CommonParameters.MinRatio == null || (peak.Intensity / largestIntensity) >= CommonParameters.MinRatio)
                        {
                            // assume charge state 1 to calculate mz tolerance
                            int obsFragmentFloorMz = (int)Math.Floor((CommonParameters.ProductMassTolerance.GetMinimumValue(peak.Mz)) * fragmentBinsPerDalton);
                            if (obsFragmentFloorMz < obsPreviousFragmentCeilingMz)
                                obsFragmentFloorMz = obsPreviousFragmentCeilingMz;
                            int obsFragmentCeilingMz = (int)Math.Ceiling((CommonParameters.ProductMassTolerance.GetMaximumValue(peak.Mz)) * fragmentBinsPerDalton);
                            obsPreviousFragmentCeilingMz = obsFragmentCeilingMz + 1;
                            FirstPassIndexedScoring(obsFragmentFloorMz, obsFragmentCeilingMz, scoringTable, byteScoreCutoff, idsOfPeptidesPossiblyObserved, scan.PrecursorMass, lowestMassPeptideToLookFor, highestMassPeptideToLookFor, ref numChecksSkipped);

                            if (addCompIons)
                            {
                                //okay, we're not actually adding in complementary m/z peaks, we're doing a shortcut and just straight up adding the mass assuming that they're z=1
                                for (int j = 0; j < dissociationTypes.Count; j++)
                                {
                                    int compPrecursor = (int)((scan.PrecursorMass + complementaryIonConversionDictionary[dissociationTypes[j]] + Constants.protonMass) * fragmentBinsPerDalton);
                                    int compFragmentFloorMz = compPrecursor - obsFragmentCeilingMz;
                                    int compFragmentCeilingMz = compPrecursor - obsFragmentFloorMz;
                                    FirstPassIndexedScoring(compFragmentFloorMz, compFragmentCeilingMz, scoringTable, byteScoreCutoff, idsOfPeptidesPossiblyObserved, scan.PrecursorMass, lowestMassPeptideToLookFor, highestMassPeptideToLookFor, ref numChecksSkipped);
                                }
                            }
                        }
                    }
                        

                    // done with initial scoring; refine scores and create PSMs
                    if (idsOfPeptidesPossiblyObserved.Any())
                    {
                        int maxInitialScore = idsOfPeptidesPossiblyObserved.Max(id => scoringTable[id]);

                        foreach (var id in idsOfPeptidesPossiblyObserved.Where(id => scoringTable[id] == maxInitialScore))
                        {
                            var candidatePeptide = peptideIndex[id];
                            double[] fragmentMasses = candidatePeptide.ProductMassesMightHaveDuplicatesAndNaNs(lp).Distinct().Where(p => !Double.IsNaN(p)).OrderBy(p => p).ToArray();
                            double peptideScore = CalculatePeptideScore(scan.TheScan, CommonParameters.ProductMassTolerance, fragmentMasses, scan.PrecursorMass, dissociationTypes, addCompIons);
                            int notch = massDiffAcceptor.Accepts(scan.PrecursorMass, candidatePeptide.MonoisotopicMassIncludingFixedMods);

                            if (globalPsms[i] == null)
                                globalPsms[i] = new Psm(candidatePeptide, notch, peptideScore, i, scan);
                            else
                                globalPsms[i].AddOrReplace(candidatePeptide, peptideScore, notch, CommonParameters.ReportAllAmbiguity);
                        }
                    }

                    // report search progress
                    progress++;
                    var percentProgress = (int)((progress / listOfSortedms2Scans.Length) * 100);

                    if (percentProgress > oldPercentProgress)
                    {
                        oldPercentProgress = percentProgress;
                        ReportProgress(new ProgressEventArgs(percentProgress, "Performing modern search... " + currentPartition + "/" + CommonParameters.TotalPartitions, nestedIds));
                    }
                }
            });

            return new MetaMorpheusEngineResults(this);
        }

        #endregion Protected Methods

        #region Private Methods

        private void FirstPassIndexedScoring(int fragmentFloorMz, int fragmentCeilingMz, byte[] scoringTable, byte byteScoreCutoff, HashSet<int> idsOfPeptidesPossiblyObserved, double scanPrecursorMass, double lowestMassPeptideToLookFor, double highestMassPeptideToLookFor, ref int numChecksSkipped)
        {
            // get all theoretical fragments this experimental fragment could be
            for (int fragmentBin = fragmentFloorMz; fragmentBin <= fragmentCeilingMz; fragmentBin++)
            {
                if (fragmentIndex[fragmentBin] != null)
                {
                    List<int> peptideIdsInThisBin = fragmentIndex[fragmentBin];

                    int m = 0;
                    if (!Double.IsInfinity(lowestMassPeptideToLookFor))
                    {
                        int l = 0;
                        int r = peptideIdsInThisBin.Count - 1;

                        // binary search in the fragment bin for lowest acceptable precursor mass
                        while (l <= r)
                        {
                            m = l + ((r - l) / 2);

                            if (r - l < 2)
                                break;
                            if (peptideIndex[peptideIdsInThisBin[m]].MonoisotopicMassIncludingFixedMods < lowestMassPeptideToLookFor)
                                l = m + 1;
                            else
                                r = m - 1;
                        }
                        if (m > 0)
                            m--;
                    }

                    // add +1 score for each peptide candidate in the scoring table up to the maximum allowed precursor mass
                    if (!Double.IsInfinity(highestMassPeptideToLookFor))
                    {
                        if (numChecksSkipped == byteScoreCutoff)
                        {
                            for (int h = m; h < peptideIdsInThisBin.Count; h++)
                            {
                                int id = peptideIdsInThisBin[h];
                                scoringTable[id]++;

                                // add possible search results to the hashset of id's
                                if (scoringTable[id] == byteScoreCutoff)
                                {
                                    int notch = massDiffAcceptor.Accepts(scanPrecursorMass, peptideIndex[id].MonoisotopicMassIncludingFixedMods);

                                    if (notch >= 0)
                                    {
                                        idsOfPeptidesPossiblyObserved.Add(id);
                                    }
                                }

                                if (peptideIndex[id].MonoisotopicMassIncludingFixedMods > highestMassPeptideToLookFor)
                                    break;
                            }
                        }
                        else
                        {
                            numChecksSkipped++;
                            for (int h = m; h < peptideIdsInThisBin.Count; h++)
                            {
                                int id = peptideIdsInThisBin[h];
                                scoringTable[id]++;

                                if (peptideIndex[id].MonoisotopicMassIncludingFixedMods > highestMassPeptideToLookFor)
                                    break;
                            }
                        }
                    }
                    else
                    {
                        if (numChecksSkipped == byteScoreCutoff)
                        {
                            for (int h = m; h < peptideIdsInThisBin.Count; h++)
                            {
                                int id = peptideIdsInThisBin[h];
                                scoringTable[id]++;

                                // add possible search results to the hashset of id's
                                if (scoringTable[id] == byteScoreCutoff)
                                {
                                    int notch = massDiffAcceptor.Accepts(scanPrecursorMass, peptideIndex[id].MonoisotopicMassIncludingFixedMods);

                                    if (notch >= 0)
                                    {
                                        idsOfPeptidesPossiblyObserved.Add(id);
                                    }
                                }
                            }
                        }
                        else
                        {
                            numChecksSkipped++;
                            for (int h = m; h < peptideIdsInThisBin.Count; h++)
                            {
                                int id = peptideIdsInThisBin[h];
                                scoringTable[id]++;
                            }
                        }
                    }
                }
            }
        }       

        #endregion Private Methods

    }
}