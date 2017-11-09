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

        protected const int fragmentBinsPerDalton = 1000;
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

        #endregion Protected Fields

        #region Public Constructors

        public ModernSearchEngine(Psm[] globalPsms, Ms2ScanWithSpecificMass[] listOfSortedms2Scans, List<CompactPeptide> peptideIndex, List<int>[] fragmentIndex, List<ProductType> lp, int currentPartition, CommonParameters CommonParameters, bool addCompIons, MassDiffAcceptor massDiffAcceptor, List<string> nestedIds) : base(nestedIds)
        {
            this.globalPsms = globalPsms;
            this.listOfSortedms2Scans = listOfSortedms2Scans;
            this.peptideIndex = peptideIndex;
            this.fragmentIndex = fragmentIndex;
            this.lp = lp;
            this.currentPartition = currentPartition + 1;
            this.CommonParameters = CommonParameters;
            this.addCompIons = addCompIons;
            this.massDiffAcceptor = massDiffAcceptor;
            this.dissociationTypes = DetermineDissociationType(lp);
        }

        #endregion Public Constructors

        #region Protected Methods

        protected override MetaMorpheusEngineResults RunSpecific()
        {
            double progress = 0;
            int oldPercentProgress = 0;
            ReportProgress(new ProgressEventArgs(oldPercentProgress, "Performing modern search... " + currentPartition + "/" + CommonParameters.TotalPartitions, nestedIds));

            byte byteScoreCutoff = Convert.ToByte((int)CommonParameters.ScoreCutoff);


            var proteaseString = CommonParameters.DigestionParams.Protease.ToString();
            if (proteaseString.Equals("top-down"))
            {
                Parallel.ForEach(Partitioner.Create(0, listOfSortedms2Scans.Length), new ParallelOptions { MaxDegreeOfParallelism = CommonParameters.MaxThreadsToUsePerFile }, range =>
                {
                    byte[] scoringTable = new byte[peptideIndex.Count];
                    List<int> idsOfProteinsPossiblyObserved = new List<int>();

                    for (int i = range.Item1; i < range.Item2; i++)
                    {
                        Array.Clear(scoringTable, 0, scoringTable.Length);
                        idsOfProteinsPossiblyObserved.Clear();
                        var scan = listOfSortedms2Scans[i];

                        // deconvolute ms2
                        double aggregationTolerancePpm = 5;
                        var decon2 = scan.TheScan.MassSpectrum.Deconvolute(scan.TheScan.ScanWindowRange, 10, 10, 5);
                        List<DeconvolutionFeatureWithMassesAndScans> currentListOfGroups = new List<DeconvolutionFeatureWithMassesAndScans>();

                        foreach (var isotopicEnvelope in decon2)
                        {
                            DeconvolutionFeatureWithMassesAndScans matchingGroup = null;
                            var mass = isotopicEnvelope.monoisotopicMass;
                            foreach (var possibleGroup in currentListOfGroups)
                            {
                                var possibleGroupMass = possibleGroup.Mass;
                                if (Math.Abs(mass - possibleGroupMass) / possibleGroupMass * 1e6 <= aggregationTolerancePpm ||
                                    Math.Abs(mass + 1.002868314 - possibleGroupMass) / possibleGroupMass * 1e6 <= aggregationTolerancePpm ||
                                    Math.Abs(mass + 2.005408917 - possibleGroupMass) / possibleGroupMass * 1e6 <= aggregationTolerancePpm ||
                                    Math.Abs(mass + 3.007841294 - possibleGroupMass) / possibleGroupMass * 1e6 <= aggregationTolerancePpm ||
                                    Math.Abs(mass - 1.002868314 - possibleGroupMass) / possibleGroupMass * 1e6 <= aggregationTolerancePpm ||
                                    Math.Abs(mass - 2.005408917 - possibleGroupMass) / possibleGroupMass * 1e6 <= aggregationTolerancePpm ||
                                    Math.Abs(mass - 3.007841294 - possibleGroupMass) / possibleGroupMass * 1e6 <= aggregationTolerancePpm)
                                {
                                    matchingGroup = possibleGroup;
                                    matchingGroup.AddEnvelope(isotopicEnvelope, scan.TheScan.OneBasedScanNumber, scan.TheScan.RetentionTime);
                                    break;
                                }
                            }

                            if (matchingGroup == null)
                            {
                                var newGroupScans = new DeconvolutionFeatureWithMassesAndScans();
                                newGroupScans.AddEnvelope(isotopicEnvelope, scan.TheScan.OneBasedScanNumber, scan.TheScan.RetentionTime);
                                currentListOfGroups.Add(newGroupScans);
                            }
                        }

                        // get bins to search
                        var allBinsToSearch = new List<int>();
                        foreach (var isotopeCluster in currentListOfGroups)
                        {
                            double fragmentMass = isotopeCluster.Mass;
                            if (fragmentMass > scan.PrecursorMass)
                                continue;

                            // get theoretical fragment bins within mass tolerance
                            int obsFragmentFloorMass = (int)Math.Floor((CommonParameters.ProductMassTolerance.GetMinimumValue(fragmentMass)) * fragmentBinsPerDalton);
                            int obsFragmentCeilingMass = (int)Math.Ceiling((CommonParameters.ProductMassTolerance.GetMaximumValue(fragmentMass)) * fragmentBinsPerDalton);
                            for (int fragmentBin = obsFragmentFloorMass; fragmentBin <= obsFragmentCeilingMass; fragmentBin++)
                            {
                                if (fragmentBin < fragmentIndex.Length)
                                    if (fragmentIndex[fragmentBin] != null)
                                        allBinsToSearch.Add(fragmentBin);
                            }
                        }

                        // get allowed theoretical masses from the known experimental mass
                        // note that this is the OPPOSITE of the classic search (which calculates experimental masses from theoretical values)
                        // this is just PRELIMINARY precursor-mass filtering
                        // additional checks are made later to ensure that the theoretical precursor mass is acceptable
                        var notches = massDiffAcceptor.GetAllowedPrecursorMassIntervals(scan.PrecursorMass);

                        double lowestMassPeptideToLookFor = Double.NegativeInfinity;
                        double highestMassPeptideToLookFor = Double.PositiveInfinity;

                        double largestMassDiff = notches.Max(p => p.allowedInterval.Maximum);
                        double smallestMassDiff = notches.Min(p => p.allowedInterval.Minimum);

                        if (!Double.IsInfinity(largestMassDiff))
                        {
                            double largestOppositeMassDiff = -1 * (notches.Max(p => p.allowedInterval.Maximum) - scan.PrecursorMass);
                            lowestMassPeptideToLookFor = scan.PrecursorMass + largestOppositeMassDiff;
                        }
                        if (!Double.IsNegativeInfinity(smallestMassDiff))
                        {
                            double smallestOppositeMassDiff = -1 * (notches.Min(p => p.allowedInterval.Minimum) - scan.PrecursorMass);
                            highestMassPeptideToLookFor = scan.PrecursorMass + smallestOppositeMassDiff;
                        }

                        // first-pass scoring
                        IndexedScoring(allBinsToSearch, scoringTable, byteScoreCutoff, idsOfProteinsPossiblyObserved, scan.PrecursorMass, lowestMassPeptideToLookFor, highestMassPeptideToLookFor);

                        // create PrSM
                        if (idsOfProteinsPossiblyObserved.Any())
                        {
                            int maxInitialScore = idsOfProteinsPossiblyObserved.Max(id => scoringTable[id]);
                            var possiblyValidIds = idsOfProteinsPossiblyObserved.Where(id => scoringTable[id] == maxInitialScore);

                            foreach (var id in possiblyValidIds)
                            {
                                var peptide = peptideIndex[id];

                                //var productMasses = peptide.ProductMassesMightHaveDuplicatesAndNaNs(lp);
                                //Array.Sort(productMasses);

                                double thePrecursorMass = scan.PrecursorMass;
                                //double score = CalculatePeptideScore(scan.TheScan, CommonParameters.ProductMassTolerance, productMasses, thePrecursorMass, dissociationTypes, addCompIons);
                                double score = scoringTable[id];
                                int notch = massDiffAcceptor.Accepts(scan.PrecursorMass, peptide.MonoisotopicMassIncludingFixedMods);

                                if (score > CommonParameters.ScoreCutoff)
                                {
                                    if (globalPsms[i] == null)
                                        globalPsms[i] = new Psm(peptide, notch, score, i, scan, CommonParameters.ExcelCompatible);
                                    else
                                        globalPsms[i].AddOrReplace(peptide, score, notch, CommonParameters.ReportAllAmbiguity);
                                }
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
            }
            else
            {
                Parallel.ForEach(Partitioner.Create(0, listOfSortedms2Scans.Length), new ParallelOptions { MaxDegreeOfParallelism = CommonParameters.MaxThreadsToUsePerFile }, range =>
                {
                    byte[] scoringTable = new byte[peptideIndex.Count];
                    List<int> idsOfPeptidesPossiblyObserved = new List<int>();

                    for (int i = range.Item1; i < range.Item2; i++)
                    {
                        // empty the scoring table to score the new scan (conserves memory compared to allocating a new array)
                        Array.Clear(scoringTable, 0, scoringTable.Length);
                        idsOfPeptidesPossiblyObserved.Clear();
                        var scan = listOfSortedms2Scans[i];

                        // get fragment bins for this scan
                        List<int> allBinsToSearch = GetBinsToSearch(scan);

                        // get allowed theoretical masses from the known experimental mass
                        // note that this is the OPPOSITE of the classic search (which calculates experimental masses from theoretical values)
                        // this is just PRELIMINARY precursor-mass filtering
                        // additional checks are made later to ensure that the theoretical precursor mass is acceptable
                        var notches = massDiffAcceptor.GetAllowedPrecursorMassIntervals(scan.PrecursorMass);

                        double lowestMassPeptideToLookFor = Double.NegativeInfinity;
                        double highestMassPeptideToLookFor = Double.PositiveInfinity;

                        double largestMassDiff = notches.Max(p => p.allowedInterval.Maximum);
                        double smallestMassDiff = notches.Min(p => p.allowedInterval.Minimum);

                        if (!Double.IsInfinity(largestMassDiff))
                        {
                            double largestOppositeMassDiff = -1 * (notches.Max(p => p.allowedInterval.Maximum) - scan.PrecursorMass);
                            lowestMassPeptideToLookFor = scan.PrecursorMass + largestOppositeMassDiff;
                        }
                        if (!Double.IsNegativeInfinity(smallestMassDiff))
                        {
                            double smallestOppositeMassDiff = -1 * (notches.Min(p => p.allowedInterval.Minimum) - scan.PrecursorMass);
                            highestMassPeptideToLookFor = scan.PrecursorMass + smallestOppositeMassDiff;
                        }

                        // first-pass scoring
                        IndexedScoring(allBinsToSearch, scoringTable, byteScoreCutoff, idsOfPeptidesPossiblyObserved, scan.PrecursorMass, lowestMassPeptideToLookFor, highestMassPeptideToLookFor);

                        // done with indexed scoring; refine scores and create PSMs
                        if (idsOfPeptidesPossiblyObserved.Any())
                        {
                            int maxInitialScore = idsOfPeptidesPossiblyObserved.Max(id => scoringTable[id]);
                            var possiblyValidIds = idsOfPeptidesPossiblyObserved.Where(id => scoringTable[id] == maxInitialScore);

                            foreach (var id in possiblyValidIds)
                            {
                                var peptide = peptideIndex[id];

                                var productMasses = peptide.ProductMassesMightHaveDuplicatesAndNaNs(lp);
                                Array.Sort(productMasses);

                                double thePrecursorMass = scan.PrecursorMass;
                                double score = CalculatePeptideScore(scan.TheScan, CommonParameters.ProductMassTolerance, productMasses, thePrecursorMass, dissociationTypes, addCompIons);
                                int notch = massDiffAcceptor.Accepts(scan.PrecursorMass, peptide.MonoisotopicMassIncludingFixedMods);

                                if (score > CommonParameters.ScoreCutoff)
                                {
                                    if (globalPsms[i] == null)
                                        globalPsms[i] = new Psm(peptide, notch, score, i, scan, CommonParameters.ExcelCompatible);
                                    else
                                        globalPsms[i].AddOrReplace(peptide, score, notch, CommonParameters.ReportAllAmbiguity);
                                }
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
            }

            return new MetaMorpheusEngineResults(this);
        }

        protected List<int> GetBinsToSearch(Ms2ScanWithSpecificMass scan)
        {
            int obsPreviousFragmentCeilingMz = 0;
            List<int> binsToSearch = new List<int>();
            foreach (var peakMz in scan.TheScan.MassSpectrum.XArray)
            {
                // assume charge state 1 to calculate mass tolerance
                double experimentalFragmentMass = ClassExtensions.ToMass(peakMz, 1);

                // get theoretical fragment bins within mass tolerance
                int obsFragmentFloorMass = (int)Math.Floor((CommonParameters.ProductMassTolerance.GetMinimumValue(experimentalFragmentMass)) * fragmentBinsPerDalton);
                if (obsFragmentFloorMass < obsPreviousFragmentCeilingMz)
                    obsFragmentFloorMass = obsPreviousFragmentCeilingMz;
                int obsFragmentCeilingMass = (int)Math.Ceiling((CommonParameters.ProductMassTolerance.GetMaximumValue(experimentalFragmentMass)) * fragmentBinsPerDalton);
                obsPreviousFragmentCeilingMz = obsFragmentCeilingMass + 1;
                for (int fragmentBin = obsFragmentFloorMass; fragmentBin <= obsFragmentCeilingMass; fragmentBin++)
                    if (fragmentIndex[fragmentBin] != null)
                        binsToSearch.Add(fragmentBin);

                if (addCompIons)
                {
                    //okay, we're not actually adding in complementary m/z peaks, we're doing a shortcut and just straight up adding the bins assuming that they're z=1
                    foreach (DissociationType dissociationType in dissociationTypes)
                    {
                        if (complementaryIonConversionDictionary.TryGetValue(dissociationType, out double protonMassShift))
                        {
                            protonMassShift = ClassExtensions.ToMass(protonMassShift, 1);
                            int compFragmentFloorMass = (int)Math.Round(((scan.PrecursorMass + protonMassShift) * fragmentBinsPerDalton)) - obsFragmentCeilingMass;
                            int compFragmentCeilingMass = (int)Math.Round(((scan.PrecursorMass + protonMassShift) * fragmentBinsPerDalton)) - obsFragmentFloorMass;
                            if (compFragmentFloorMass > 0)
                                for (int fragmentBin = compFragmentFloorMass; fragmentBin <= compFragmentCeilingMass; fragmentBin++)
                                    if (fragmentIndex[fragmentBin] != null)
                                        binsToSearch.Add(fragmentBin);
                        }
                        else
                            throw new NotImplementedException();
                    }
                }
            }
            return binsToSearch;
        }

        #endregion Protected Methods

        #region Private Methods

        private int BinarySearchBinForPrecursorIndex(List<int> peptideIdsInThisBin, double peptideMassToLookFor)
        {
            int m = 0;
            int l = 0;
            int r = peptideIdsInThisBin.Count - 1;

            // binary search in the fragment bin for precursor mass
            while (l <= r)
            {
                m = l + ((r - l) / 2);

                if (r - l < 2)
                    break;
                if (peptideIndex[peptideIdsInThisBin[m]].MonoisotopicMassIncludingFixedMods < peptideMassToLookFor)
                    l = m + 1;
                else
                    r = m - 1;
            }
            if (m > 0)
                m--;
            return m;
        }

        private void IndexedScoring(List<int> binsToSearch, byte[] scoringTable, byte byteScoreCutoff, List<int> idsOfPeptidesPossiblyObserved, double scanPrecursorMass, double lowestMassPeptideToLookFor, double highestMassPeptideToLookFor)
        {
            // get all theoretical fragments this experimental fragment could be
            for (int i = 0; i < binsToSearch.Count; i++)
            {
                List<int> peptideIdsInThisBin = fragmentIndex[binsToSearch[i]];

                //get index for minimum monoisotopic allowed
                int lowestPeptideMassIndex = Double.IsInfinity(lowestMassPeptideToLookFor) ? 0 : BinarySearchBinForPrecursorIndex(peptideIdsInThisBin, lowestMassPeptideToLookFor);

                // get index for highest mass allowed
                int highestPeptideMassIndex = peptideIdsInThisBin.Count - 1;

                if (!Double.IsInfinity(highestMassPeptideToLookFor))
                {
                    highestPeptideMassIndex = BinarySearchBinForPrecursorIndex(peptideIdsInThisBin, highestMassPeptideToLookFor);

                    for (int j = highestPeptideMassIndex; j < peptideIdsInThisBin.Count; j++)
                    {
                        int nextId = peptideIdsInThisBin[j];
                        var nextPep = peptideIndex[nextId];
                        if (nextPep.MonoisotopicMassIncludingFixedMods < highestMassPeptideToLookFor)
                            highestPeptideMassIndex = j;
                        else
                            break;
                    }
                }

                // add +1 score for each peptide candidate in the scoring table up to the maximum allowed precursor mass
                for (int j = lowestPeptideMassIndex; j <= highestPeptideMassIndex; j++)
                {
                    int id = peptideIdsInThisBin[j];
                    scoringTable[id]++;

                    // add possible search results to the hashset of id's
                    if (scoringTable[id] == byteScoreCutoff && massDiffAcceptor.Accepts(scanPrecursorMass, peptideIndex[id].MonoisotopicMassIncludingFixedMods) >= 0)
                        idsOfPeptidesPossiblyObserved.Add(id);
                }
            }
        }

        #endregion Private Methods
    }
}