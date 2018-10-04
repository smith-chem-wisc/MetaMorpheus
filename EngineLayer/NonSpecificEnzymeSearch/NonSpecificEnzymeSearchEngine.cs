using Chemistry;
using EngineLayer.ModernSearch;
using MassSpectrometry;
using Proteomics;
using Proteomics.AminoAcidPolymer;
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

        private readonly List<int>[] fragmentIndexPrecursor;
        private readonly int MinimumPeptideLength;

        public NonSpecificEnzymeSearchEngine(PeptideSpectralMatch[][] globalPsms, Ms2ScanWithSpecificMass[] listOfSortedms2Scans, List<PeptideWithSetModifications> peptideIndex, List<int>[] fragmentIndex, List<int>[] fragmentIndexPrecursor, int currentPartition, CommonParameters CommonParameters, MassDiffAcceptor massDiffAcceptor, double maximumMassThatFragmentIonScoreIsDoubled, List<string> nestedIds) : base(globalPsms, listOfSortedms2Scans, peptideIndex, fragmentIndex, currentPartition, CommonParameters, massDiffAcceptor, maximumMassThatFragmentIonScoreIsDoubled, nestedIds)
        {
            this.fragmentIndexPrecursor = fragmentIndexPrecursor;
            MinimumPeptideLength = commonParameters.DigestionParams.MinPeptideLength;
        }

        protected override MetaMorpheusEngineResults RunSpecific()
        {
            double progress = 0;
            int oldPercentProgress = 0;
            ReportProgress(new ProgressEventArgs(oldPercentProgress, "Performing nonspecific search... " + CurrentPartition + "/" + commonParameters.TotalPartitions, nestedIds));

            byte byteScoreCutoff = (byte)commonParameters.ScoreCutoff;

            Parallel.ForEach(Partitioner.Create(0, ListOfSortedMs2Scans.Length), new ParallelOptions { MaxDegreeOfParallelism = commonParameters.MaxThreadsToUsePerFile }, range =>
            {
                byte[] scoringTable = new byte[PeptideIndex.Count];
                HashSet<int> idsOfPeptidesPossiblyObserved = new HashSet<int>();

                for (int i = range.Item1; i < range.Item2; i++)
                {
                    // empty the scoring table to score the new scan (conserves memory compared to allocating a new array)
                    Array.Clear(scoringTable, 0, scoringTable.Length);
                    idsOfPeptidesPossiblyObserved.Clear();
                    var scan = ListOfSortedMs2Scans[i];

                    //get bins to add points to
                    List<int> allBinsToSearch = GetBinsToSearch(scan);

                    //the entire indexed scoring is done here
                    for (int j = 0; j < allBinsToSearch.Count; j++)
                    {
                        FragmentIndex[allBinsToSearch[j]].ForEach(id => scoringTable[id]++);
                    }

                    //populate ids of possibly observed with those containing allowed precursor masses
                    List<int> binsToSearch = new List<int>();
                    int obsPrecursorFloorMz = (int)Math.Floor(commonParameters.PrecursorMassTolerance.GetMinimumValue(scan.PrecursorMass) * FragmentBinsPerDalton);
                    int obsPrecursorCeilingMz = (int)Math.Ceiling(commonParameters.PrecursorMassTolerance.GetMaximumValue(scan.PrecursorMass) * FragmentBinsPerDalton);
                    for (int fragmentBin = obsPrecursorFloorMz; fragmentBin <= obsPrecursorCeilingMz; fragmentBin++)
                    {
                        binsToSearch.Add(fragmentBin);
                    }

                    foreach (ProductType pt in DissociationTypeCollection.ProductsFromDissociationType[commonParameters.DissociationType].Intersect(TerminusSpecificProductTypes.ProductIonTypesFromSpecifiedTerminus[commonParameters.DigestionParams.FragmentationTerminus]).ToList())
                    {
                        int binShift = (int)Math.Round((WaterMonoisotopicMass - DissociationTypeCollection.GetMassShiftFromProductType(pt)) * FragmentBinsPerDalton);

                        for (int j = 0; j < binsToSearch.Count; j++)
                        {
                            int bin = binsToSearch[j] - binShift;
                            if (bin < FragmentIndex.Length && FragmentIndex[bin] != null)
                            {
                                FragmentIndex[bin].ForEach(id => idsOfPeptidesPossiblyObserved.Add(id));
                            }
                        }
                    }

                    for (int j = 0; j < binsToSearch.Count; j++)
                    {
                        int bin = binsToSearch[j];
                        if (bin < fragmentIndexPrecursor.Length && fragmentIndexPrecursor[bin] != null)
                        {
                            fragmentIndexPrecursor[bin].ForEach(id => idsOfPeptidesPossiblyObserved.Add(id));
                        }
                    }

                    // done with initial scoring; refine scores and create PSMs
                    if (idsOfPeptidesPossiblyObserved.Any())
                    {
                        int maxInitialScore = idsOfPeptidesPossiblyObserved.Max(id => scoringTable[id]) + 1;
                        while (maxInitialScore > commonParameters.ScoreCutoff) //go through all until we hit the end
                        {
                            maxInitialScore--;
                            foreach (var id in idsOfPeptidesPossiblyObserved.Where(id => scoringTable[id] == maxInitialScore))
                            {
                                PeptideWithSetModifications peptide = PeptideIndex[id];
                                List<Product> peptideTheorProducts = peptide.Fragment(commonParameters.DissociationType, commonParameters.DigestionParams.FragmentationTerminus).ToList();

                                List<MatchedFragmentIon> matchedIons = MatchFragmentIons(scan.TheScan.MassSpectrum, peptideTheorProducts, commonParameters, scan.PrecursorMass);

                                double thisScore = CalculatePeptideScore(scan.TheScan, matchedIons, MaxMassThatFragmentIonScoreIsDoubled);
                                if (thisScore > commonParameters.ScoreCutoff)
                                {
                                    Tuple<int, PeptideWithSetModifications> notchAndUpdatedPeptide = Accepts(peptideTheorProducts, scan.PrecursorMass, peptide, commonParameters.DigestionParams.FragmentationTerminus, MassDiffAcceptor);
                                    if (notchAndUpdatedPeptide.Item1 >= 0)
                                    {
                                        PeptideSpectralMatch[] localPeptideSpectralMatches = PeptideSpectralMatches[(int)FdrClassifier.GetCleavageSpecificityCategory(notchAndUpdatedPeptide.Item2.CleavageSpecificityForFdrCategory)];
                                        if (localPeptideSpectralMatches[i] == null)
                                        {
                                            localPeptideSpectralMatches[i] = new PeptideSpectralMatch(notchAndUpdatedPeptide.Item2, notchAndUpdatedPeptide.Item1, thisScore, i, scan, commonParameters.DigestionParams, matchedIons);
                                        }
                                        else
                                        {
                                            localPeptideSpectralMatches[i].AddOrReplace(notchAndUpdatedPeptide.Item2, thisScore, notchAndUpdatedPeptide.Item1, commonParameters.ReportAllAmbiguity, matchedIons);
                                        }                        
                                    }
                                }
                            }
                        }
                    }
                    // report search progress
                    progress++;
                    var percentProgress = (int)((progress / ListOfSortedMs2Scans.Length) * 100);

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

            for (int i = localminPeptideLength; i < fragments.Count; i++)
            {
                Product fragment = fragments[i];
                double theoMass = fragment.NeutralMass - DissociationTypeCollection.GetMassShiftFromProductType(fragment.ProductType) + WaterMonoisotopicMass;
                int notch = searchMode.Accepts(scanPrecursorMass, theoMass);
                if (notch >= 0)
                {
                    PeptideWithSetModifications updatedPwsm = null;
                    if(fragmentationTerminus == FragmentationTerminus.N)
                    {
                        int endResidue = peptide.OneBasedStartResidueInProtein + fragment.TerminusFragment.FragmentNumber -1; //-1 for one based index
                        Dictionary<int, Modification> updatedMods = new Dictionary<int, Modification>();
                        foreach(var mod in peptide.AllModsOneIsNterminus)
                        {
                            if(mod.Key<endResidue-peptide.OneBasedStartResidueInProtein) //check if we cleaved it off
                            {
                                updatedMods.Add(mod.Key, mod.Value);
                            }
                        }
                        updatedPwsm = new PeptideWithSetModifications(peptide.Protein, peptide.DigestionParams, peptide.OneBasedStartResidueInProtein, endResidue, CleavageSpecificity.Unknown, "", 0, updatedMods, 0);
                    }
                    else
                    {
                        int startResidue = peptide.OneBasedEndResidueInProtein - fragment.TerminusFragment.FragmentNumber+1; //plus one for one based index
                        Dictionary<int, Modification> updatedMods = new Dictionary<int, Modification>();  //updateMods
                        int indexShift = startResidue - peptide.OneBasedStartResidueInProtein;
                        foreach (var mod in peptide.AllModsOneIsNterminus)
                        {
                            if (mod.Key > indexShift) //check if we cleaved it off
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
                    var updatedPwsm = new PeptideWithSetModifications(peptide.Protein, peptide.DigestionParams, peptide.OneBasedStartResidueInProtein, peptide.OneBasedEndResidueInProtein, CleavageSpecificity.Unknown, "", 0, peptide.AllModsOneIsNterminus, peptide.NumFixedMods);
                    return new Tuple<int, PeptideWithSetModifications>(notch, updatedPwsm);
                }
            }
            return new Tuple<int, PeptideWithSetModifications>(-1, null);
        }
    }
}