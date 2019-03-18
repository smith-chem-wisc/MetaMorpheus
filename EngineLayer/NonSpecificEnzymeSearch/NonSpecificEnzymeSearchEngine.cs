using Chemistry;
using EngineLayer.FdrAnalysis;
using EngineLayer.ModernSearch;
using Proteomics;
using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace EngineLayer.NonSpecificEnzymeSearch
{
    public class NonSpecificEnzymeSearchEngine : ModernSearchEngine
    {
        private static readonly double WaterMonoisotopicMass = PeriodicTable.GetElement("H").PrincipalIsotope.AtomicMass * 2 + PeriodicTable.GetElement("O").PrincipalIsotope.AtomicMass;

        private readonly List<int>[] PrecursorIndex;
        private readonly int MinimumPeptideLength;
        readonly PeptideSpectralMatch[][] GlobalCategorySpecificPsms;
        readonly CommonParameters ModifiedParametersNoComp;
        readonly List<ProductType> ProductTypesToSearch;

        public NonSpecificEnzymeSearchEngine(PeptideSpectralMatch[][] globalPsms, Ms2ScanWithSpecificMass[] listOfSortedms2Scans, List<PeptideWithSetModifications> peptideIndex, List<int>[] fragmentIndex, List<int>[] precursorIndex, int currentPartition, CommonParameters CommonParameters, MassDiffAcceptor massDiffAcceptor, double maximumMassThatFragmentIonScoreIsDoubled, List<string> nestedIds) : base(null, listOfSortedms2Scans, peptideIndex, fragmentIndex, currentPartition, CommonParameters, massDiffAcceptor, maximumMassThatFragmentIonScoreIsDoubled, nestedIds)
        {
            PrecursorIndex = precursorIndex;
            MinimumPeptideLength = commonParameters.DigestionParams.MinPeptideLength;
            GlobalCategorySpecificPsms = globalPsms;
            ModifiedParametersNoComp = commonParameters.CloneWithNewTerminus(addCompIons: false);
            ProductTypesToSearch = DissociationTypeCollection.ProductsFromDissociationType[commonParameters.DissociationType].Intersect(TerminusSpecificProductTypes.ProductIonTypesFromSpecifiedTerminus[commonParameters.DigestionParams.FragmentationTerminus]).ToList();
        }

        protected override MetaMorpheusEngineResults RunSpecific()
        {
            double progress = 0;
            int oldPercentProgress = 0;
            ReportProgress(new ProgressEventArgs(oldPercentProgress, "Performing nonspecific search... " + CurrentPartition + "/" + commonParameters.TotalPartitions, nestedIds));

            byte byteScoreCutoff = (byte)commonParameters.ScoreCutoff;

            int maxThreadsPerFile = commonParameters.MaxThreadsToUsePerFile;
            int[] threads = Enumerable.Range(0, maxThreadsPerFile).ToArray();
            Parallel.ForEach(threads, (i) =>
            {
                byte[] scoringTable = new byte[PeptideIndex.Count];
                HashSet<int> idsOfPeptidesPossiblyObserved = new HashSet<int>();

                for (; i < ListOfSortedMs2Scans.Length; i += maxThreadsPerFile)
                {
                    // Stop loop if canceled
                    if (GlobalVariables.StopLoops) { return; }

                    // empty the scoring table to score the new scan (conserves memory compared to allocating a new array)
                    Array.Clear(scoringTable, 0, scoringTable.Length);
                    idsOfPeptidesPossiblyObserved.Clear();
                    Ms2ScanWithSpecificMass scan = ListOfSortedMs2Scans[i];

                    //get bins to add points to
                    List<int> allBinsToSearch = GetBinsToSearch(scan);

                    //the entire indexed scoring is done here
                    for (int j = 0; j < allBinsToSearch.Count; j++)
                    {
                        FragmentIndex[allBinsToSearch[j]].ForEach(id => scoringTable[id]++);
                    }

                    //populate ids of possibly observed with those containing allowed precursor masses
                    List<AllowedIntervalWithNotch> validIntervals = MassDiffAcceptor.GetAllowedPrecursorMassIntervalsFromObservedMass(scan.PrecursorMass).ToList(); //get all valid notches
                    foreach (AllowedIntervalWithNotch interval in validIntervals)
                    {
                        int obsPrecursorFloorMz = (int)Math.Floor(interval.AllowedInterval.Minimum * FragmentBinsPerDalton);
                        int obsPrecursorCeilingMz = (int)Math.Ceiling(interval.AllowedInterval.Maximum * FragmentBinsPerDalton);

                        foreach (ProductType pt in ProductTypesToSearch)
                        {
                            int dissociationBinShift = (int)Math.Round((WaterMonoisotopicMass - DissociationTypeCollection.GetMassShiftFromProductType(pt)) * FragmentBinsPerDalton);
                            int lowestBin = obsPrecursorFloorMz - dissociationBinShift;
                            int highestBin = obsPrecursorCeilingMz - dissociationBinShift;
                            for (int bin = lowestBin; bin <= highestBin; bin++)
                            {
                                if (bin < FragmentIndex.Length && FragmentIndex[bin] != null)
                                {
                                    FragmentIndex[bin].ForEach(id => idsOfPeptidesPossiblyObserved.Add(id));
                                }
                            }
                        }

                        for (int bin = obsPrecursorFloorMz; bin <= obsPrecursorCeilingMz; bin++) //no bin shift, since they're precursor masses
                        {
                            if (bin < PrecursorIndex.Length && PrecursorIndex[bin] != null)
                            {
                                PrecursorIndex[bin].ForEach(id => idsOfPeptidesPossiblyObserved.Add(id));
                            }
                        }
                    }

                    // done with initial scoring; refine scores and create PSMs
                    if (idsOfPeptidesPossiblyObserved.Any())
                    {
                        int maxInitialScore = idsOfPeptidesPossiblyObserved.Max(id => scoringTable[id]) + 1;
                        while (maxInitialScore > commonParameters.ScoreCutoff) //go through all until we hit the end
                        {
                            maxInitialScore--;
                            foreach (int id in idsOfPeptidesPossiblyObserved.Where(id => scoringTable[id] == maxInitialScore))
                            {
                                PeptideWithSetModifications peptide = PeptideIndex[id];
                                List<Product> peptideTheorProducts = peptide.Fragment(commonParameters.DissociationType, commonParameters.DigestionParams.FragmentationTerminus).ToList();

                                Tuple<int, PeptideWithSetModifications> notchAndUpdatedPeptide = Accepts(peptideTheorProducts, scan.PrecursorMass, peptide, commonParameters.DigestionParams.FragmentationTerminus, MassDiffAcceptor);
                                int notch = notchAndUpdatedPeptide.Item1;
                                if (notch >= 0)
                                {
                                    peptide = notchAndUpdatedPeptide.Item2;
                                    peptideTheorProducts = peptide.Fragment(commonParameters.DissociationType, FragmentationTerminus.Both).ToList();
                                    List<MatchedFragmentIon> matchedIons = MatchFragmentIons(scan, peptideTheorProducts, ModifiedParametersNoComp);

                                    double thisScore = CalculatePeptideScore(scan.TheScan, matchedIons);
                                    if (thisScore > commonParameters.ScoreCutoff)
                                    {
                                        PeptideSpectralMatch[] localPeptideSpectralMatches = GlobalCategorySpecificPsms[(int)FdrClassifier.GetCleavageSpecificityCategory(peptide.CleavageSpecificityForFdrCategory)];
                                        if (localPeptideSpectralMatches[i] == null)
                                        {
                                            localPeptideSpectralMatches[i] = new PeptideSpectralMatch(peptide, notch, thisScore, i, scan, commonParameters.DigestionParams, matchedIons);
                                        }
                                        else
                                        {
                                            localPeptideSpectralMatches[i].AddOrReplace(peptide, thisScore, notch, commonParameters.ReportAllAmbiguity, matchedIons, 0);
                                        }
                                    }
                                }
                            }
                        }
                    }
                    // report search progress
                    progress++;
                    int percentProgress = (int)((progress / ListOfSortedMs2Scans.Length) * 100);

                    if (percentProgress > oldPercentProgress)
                    {
                        oldPercentProgress = percentProgress;
                        ReportProgress(new ProgressEventArgs(percentProgress, "Performing nonspecific search... " + CurrentPartition + "/" + commonParameters.TotalPartitions, nestedIds));
                    }
                }
            });
            return new MetaMorpheusEngineResults(this);
        }

        private Tuple<int, PeptideWithSetModifications> Accepts(List<Product> fragments, double scanPrecursorMass, PeptideWithSetModifications peptide, FragmentationTerminus fragmentationTerminus, MassDiffAcceptor searchMode)
        {
            //all masses in N and CTerminalMasses are b-ion masses, which are one water away from a full peptide
            int localminPeptideLength = commonParameters.DigestionParams.MinPeptideLength;

            for (int i = localminPeptideLength - 1; i < fragments.Count; i++) //minus one start, because fragment 1 is at index 0
            {
                Product fragment = fragments[i];
                double theoMass = fragment.NeutralMass - DissociationTypeCollection.GetMassShiftFromProductType(fragment.ProductType) + WaterMonoisotopicMass;
                int notch = searchMode.Accepts(scanPrecursorMass, theoMass);
                if (notch >= 0)
                {
                    PeptideWithSetModifications updatedPwsm = null;
                    if (fragmentationTerminus == FragmentationTerminus.N)
                    {
                        int endResidue = peptide.OneBasedStartResidueInProtein + fragment.TerminusFragment.FragmentNumber - 1; //-1 for one based index
                        Dictionary<int, Modification> updatedMods = new Dictionary<int, Modification>();
                        foreach (KeyValuePair<int, Modification> mod in peptide.AllModsOneIsNterminus)
                        {
                            if (mod.Key < endResidue - peptide.OneBasedStartResidueInProtein + 3) //check if we cleaved it off, +1 for N-terminus being mod 1 and first residue being mod 2, +1 again for the -1 on end residue for one based index, +1 (again) for the one-based start residue
                            {
                                updatedMods.Add(mod.Key, mod.Value);
                            }
                        }
                        updatedPwsm = new PeptideWithSetModifications(peptide.Protein, peptide.DigestionParams, peptide.OneBasedStartResidueInProtein, endResidue, CleavageSpecificity.Unknown, "", 0, updatedMods, 0);
                    }
                    else
                    {
                        int startResidue = peptide.OneBasedEndResidueInProtein - fragment.TerminusFragment.FragmentNumber + 1; //plus one for one based index
                        Dictionary<int, Modification> updatedMods = new Dictionary<int, Modification>();  //updateMods
                        int indexShift = startResidue - peptide.OneBasedStartResidueInProtein;
                        foreach (KeyValuePair<int, Modification> mod in peptide.AllModsOneIsNterminus)
                        {
                            if (mod.Key > indexShift + 1) //check if we cleaved it off, +1 for N-terminus being mod 1 and first residue being 2
                            {
                                int key = mod.Key - indexShift;
                                updatedMods.Add(key, mod.Value);
                            }
                        }
                        updatedPwsm = new PeptideWithSetModifications(peptide.Protein, peptide.DigestionParams, startResidue, peptide.OneBasedEndResidueInProtein, CleavageSpecificity.Unknown, "", 0, updatedMods, 0);
                    }
                    return new Tuple<int, PeptideWithSetModifications>(notch, updatedPwsm);
                }
                else if (theoMass > scanPrecursorMass)
                {
                    break;
                }
            }
            //if the theoretical and experimental have the same mass
            if (fragments.Count > localminPeptideLength)
            {
                double totalMass = peptide.MonoisotopicMass;// + Constants.ProtonMass;
                int notch = searchMode.Accepts(scanPrecursorMass, totalMass);
                if (notch >= 0)
                {
                    //need to update so that the cleavage specificity is recorded
                    PeptideWithSetModifications updatedPwsm = new PeptideWithSetModifications(peptide.Protein, peptide.DigestionParams, peptide.OneBasedStartResidueInProtein, peptide.OneBasedEndResidueInProtein, CleavageSpecificity.Unknown, "", 0, peptide.AllModsOneIsNterminus, peptide.NumFixedMods);
                    return new Tuple<int, PeptideWithSetModifications>(notch, updatedPwsm);
                }
            }
            return new Tuple<int, PeptideWithSetModifications>(-1, null);
        }

        public static List<PeptideSpectralMatch> ResolveFdrCategorySpecificPsms(List<PeptideSpectralMatch>[] AllPsms, int numNotches, string taskId, CommonParameters commonParameters)
        {
            //update all psms with peptide info
            AllPsms.ToList()
                .Where(psmArray => psmArray != null).ToList()
                .ForEach(psmArray => psmArray.Where(psm => psm != null).ToList()
                .ForEach(psm => psm.ResolveAllAmbiguities()));

            foreach (List<PeptideSpectralMatch> psmsArray in AllPsms)
            {
                if (psmsArray != null)
                {
                    List<PeptideSpectralMatch> cleanedPsmsArray = psmsArray.Where(b => b != null).OrderByDescending(b => b.Score)
                       .ThenBy(b => b.PeptideMonisotopicMass.HasValue ? Math.Abs(b.ScanPrecursorMass - b.PeptideMonisotopicMass.Value) : double.MaxValue)
                       .GroupBy(b => (b.FullFilePath, b.ScanNumber, b.PeptideMonisotopicMass)).Select(b => b.First()).ToList();

                    new FdrAnalysisEngine(cleanedPsmsArray, numNotches, commonParameters, new List<string> { taskId }).Run();

                    for (int i = 0; i < psmsArray.Count; i++)
                    {
                        if (psmsArray[i] != null)
                        {
                            if (psmsArray[i].FdrInfo == null) //if it was grouped in the cleanedPsmsArray
                            {
                                psmsArray[i] = null;
                            }
                        }
                    }
                }
            }

            int[] ranking = new int[AllPsms.Length]; //high int is good ranking
            List<int> indexesOfInterest = new List<int>();
            for (int i = 0; i < ranking.Length; i++)
            {
                if (AllPsms[i] != null)
                {
                    ranking[i] = AllPsms[i].Where(x => x != null).Count(x => x.FdrInfo.QValue <= 0.01); //set ranking as number of psms above 1% FDR
                    indexesOfInterest.Add(i);
                }
            }

            //get the index of the category with the highest ranking
            int majorCategoryIndex = indexesOfInterest[0];
            for (int i = 1; i < indexesOfInterest.Count; i++)
            {
                int currentCategoryIndex = indexesOfInterest[i];
                if (ranking[currentCategoryIndex] > ranking[majorCategoryIndex])
                {
                    majorCategoryIndex = currentCategoryIndex;
                }
            }

            //update other category q-values
            //There's a chance of weird categories getting a random decoy before a random target, but we don't want to give that target a q value of zero.
            //We can't just take the q of the first decoy, because if the target wasn't random (score = 40), but there are no other targets before the decoy (score = 5), then we're incorrectly dinging the target
            //The current solution is such that if a minor category has a lower q value than it's corresponding score in the major category, then its q-value is changed to what it would be in the major category
            List<PeptideSpectralMatch> majorCategoryPsms = AllPsms[majorCategoryIndex].Where(x => x != null).OrderByDescending(x => x.Score).ToList(); //get sorted major category
            for (int i = 0; i < indexesOfInterest.Count; i++)
            {
                int minorCategoryIndex = indexesOfInterest[i];
                if (minorCategoryIndex != majorCategoryIndex)
                {
                    List<PeptideSpectralMatch> minorCategoryPsms = AllPsms[minorCategoryIndex].Where(x => x != null).OrderByDescending(x => x.Score).ToList(); //get sorted minor category
                    int minorPsmIndex = 0;
                    int majorPsmIndex = 0;
                    while (minorPsmIndex < minorCategoryPsms.Count && majorPsmIndex < majorCategoryPsms.Count) //while in the lists
                    {
                        PeptideSpectralMatch majorPsm = majorCategoryPsms[majorPsmIndex];
                        PeptideSpectralMatch minorPsm = minorCategoryPsms[minorPsmIndex];
                        //major needs to be a lower score than the minor
                        if (majorPsm.Score > minorPsm.Score)
                        {
                            majorPsmIndex++;
                        }
                        else
                        {
                            if (majorPsm.FdrInfo.QValue > minorPsm.FdrInfo.QValue)
                            {
                                minorPsm.FdrInfo.QValue = majorPsm.FdrInfo.QValue;
                            }
                            minorPsmIndex++;
                        }
                    }
                    //wrap up if we hit the end of the major category
                    while (minorPsmIndex < minorCategoryPsms.Count)
                    {
                        PeptideSpectralMatch majorPsm = majorCategoryPsms[majorPsmIndex - 1]; //-1 because it's out of index right now
                        PeptideSpectralMatch minorPsm = minorCategoryPsms[minorPsmIndex];
                        if (majorPsm.FdrInfo.QValue > minorPsm.FdrInfo.QValue)
                        {
                            minorPsm.FdrInfo.QValue = majorPsm.FdrInfo.QValue;
                        }
                        minorPsmIndex++;
                    }
                }
            }

            int numTotalSpectraWithPrecursors = AllPsms[indexesOfInterest[0]].Count;
            List<PeptideSpectralMatch> bestPsmsList = new List<PeptideSpectralMatch>();
            for (int i = 0; i < numTotalSpectraWithPrecursors; i++)
            {
                PeptideSpectralMatch bestPsm = null;
                double lowestQ = double.MaxValue;
                int bestIndex = -1;
                foreach (int index in indexesOfInterest) //foreach category
                {
                    PeptideSpectralMatch currentPsm = AllPsms[index][i];
                    if (currentPsm != null)
                    {
                        double currentQValue = currentPsm.FdrInfo.QValue;
                        if (currentQValue < lowestQ //if the new one is better
                            || (currentQValue == lowestQ && currentPsm.Score > bestPsm.Score))
                        {
                            if (bestIndex != -1)
                            {
                                //remove the old one so we don't use it for fdr later
                                AllPsms[bestIndex][i] = null;
                            }
                            bestPsm = currentPsm;
                            lowestQ = currentQValue;
                            bestIndex = index;
                        }
                        else //remove the old one so we don't use it for fdr later
                        {
                            AllPsms[index][i] = null;
                        }
                    }
                }
                if (bestPsm != null)
                {
                    bestPsmsList.Add(bestPsm);
                }
            }

            //It's probable that psms from some categories were removed by psms from other categories.
            //however, the fdr is still affected by their presence, since it was calculated before their removal.
            foreach (List<PeptideSpectralMatch> psmsArray in AllPsms)
            {
                if (psmsArray != null)
                {
                    List<PeptideSpectralMatch> cleanedPsmsArray = psmsArray.Where(b => b != null).OrderByDescending(b => b.Score)
                       .ThenBy(b => b.PeptideMonisotopicMass.HasValue ? Math.Abs(b.ScanPrecursorMass - b.PeptideMonisotopicMass.Value) : double.MaxValue)
                       .ToList();

                    new FdrAnalysisEngine(cleanedPsmsArray, numNotches, commonParameters, new List<string> { taskId }).Run();
                }
            }

            return bestPsmsList.OrderBy(b => b.FdrInfo.QValue).ThenByDescending(b => b.Score).ToList();
        }
    }
}