using Chemistry;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace EngineLayer.NonSpecificEnzymeSearch
{
    public class NonSpecificEnzymeSearchEngine : ModernSearch.ModernSearchEngine
    {
        #region Private Fields

        private static readonly double nitrogenAtomMonoisotopicMass = PeriodicTable.GetElement("N").PrincipalIsotope.AtomicMass;
        private static readonly double oxygenAtomMonoisotopicMass = PeriodicTable.GetElement("O").PrincipalIsotope.AtomicMass;
        private static readonly double hydrogenAtomMonoisotopicMass = PeriodicTable.GetElement("H").PrincipalIsotope.AtomicMass;
        private static readonly double waterMonoisotopicMass = PeriodicTable.GetElement("H").PrincipalIsotope.AtomicMass * 2 + PeriodicTable.GetElement("O").PrincipalIsotope.AtomicMass;
        private static readonly int bBinShift = (int)Math.Round((waterMonoisotopicMass) * fragmentBinsPerDalton);
        private static readonly int cBinShift = (int)Math.Round((nitrogenAtomMonoisotopicMass + 3 * hydrogenAtomMonoisotopicMass) * fragmentBinsPerDalton);
        private static readonly int zdotBinShift = (int)Math.Round((oxygenAtomMonoisotopicMass - nitrogenAtomMonoisotopicMass) * fragmentBinsPerDalton);
        private readonly List<int>[] fragmentIndexPrecursor;

        #endregion Private Fields

        #region Public Constructors

        public NonSpecificEnzymeSearchEngine(Psm[] globalPsms, Ms2ScanWithSpecificMass[] listOfSortedms2Scans, List<CompactPeptide> peptideIndex, List<int>[] fragmentIndex, List<int>[] fragmentIndexPrecursor, List<ProductType> lp, int currentPartition, ICommonParameters CommonParameters, bool addCompIons, MassDiffAcceptor massDiffAcceptors, double maximumMassThatFragmentIonScoreIsDoubled, List<string> nestedIds) : base(globalPsms, listOfSortedms2Scans, peptideIndex, fragmentIndex, lp, currentPartition, CommonParameters, addCompIons, massDiffAcceptors, maximumMassThatFragmentIonScoreIsDoubled, nestedIds)
        {
            this.fragmentIndexPrecursor = fragmentIndexPrecursor;
        }

        #endregion Public Constructors

        #region Protected Methods

        protected override MetaMorpheusEngineResults RunSpecific()
        {
            double progress = 0;
            int oldPercentProgress = 0;
            ReportProgress(new ProgressEventArgs(oldPercentProgress, "Performing nonspecific search... " + currentPartition + "/" + CommonParameters.TotalPartitions, nestedIds));
            TerminusType terminusType = ProductTypeMethod.IdentifyTerminusType(lp);

            byte byteScoreCutoff = (byte)CommonParameters.ScoreCutoff;

            Parallel.ForEach(Partitioner.Create(0, listOfSortedms2Scans.Length), new ParallelOptions { MaxDegreeOfParallelism = CommonParameters.MaxThreadsToUsePerFile }, range =>
            {
                byte[] scoringTable = new byte[peptideIndex.Count];
                HashSet<int> idsOfPeptidesPossiblyObserved = new HashSet<int>();

                for (int i = range.Item1; i < range.Item2; i++)
                {
                    // empty the scoring table to score the new scan (conserves memory compared to allocating a new array)
                    Array.Clear(scoringTable, 0, scoringTable.Length);
                    idsOfPeptidesPossiblyObserved.Clear();
                    var scan = listOfSortedms2Scans[i];

                    //get bins to add points to
                    List<int> allBinsToSearch = GetBinsToSearch(scan);

                    for (int j = 0; j < allBinsToSearch.Count; j++)
                        fragmentIndex[allBinsToSearch[j]].ForEach(id => scoringTable[id]++);

                    //populate ids of possibly observed with those containing allowed precursor masses
                    List<int> binsToSearch = new List<int>();
                    int obsPrecursorFloorMz = (int)Math.Floor(CommonParameters.PrecursorMassTolerance.GetMinimumValue(scan.PrecursorMass) * fragmentBinsPerDalton);
                    int obsPrecursorCeilingMz = (int)Math.Ceiling(CommonParameters.PrecursorMassTolerance.GetMaximumValue(scan.PrecursorMass) * fragmentBinsPerDalton);
                    for (int fragmentBin = obsPrecursorFloorMz; fragmentBin <= obsPrecursorCeilingMz; fragmentBin++)
                        binsToSearch.Add(fragmentBin);

                    foreach (ProductType pt in lp)
                    {
                        int binShift;
                        switch (pt)
                        {
                            case ProductType.B:
                                binShift = bBinShift;
                                break;

                            case ProductType.Y:
                                binShift = 0;
                                break;

                            case ProductType.C:
                                binShift = cBinShift;
                                break;

                            case ProductType.Zdot:
                                binShift = zdotBinShift;
                                break;

                            default:
                                throw new NotImplementedException();
                        }
                        for (int j = 0; j < binsToSearch.Count; j++)
                        {
                            int bin = binsToSearch[j] - binShift;
                            if (bin < fragmentIndex.Length && fragmentIndex[bin] != null)
                                fragmentIndex[bin].ForEach(id => idsOfPeptidesPossiblyObserved.Add(id));
                        }
                    }

                    for (int j = 0; j < binsToSearch.Count; j++)
                    {
                        int bin = binsToSearch[j];
                        if (bin < fragmentIndexPrecursor.Length && fragmentIndexPrecursor[bin] != null)
                            fragmentIndexPrecursor[bin].ForEach(id => idsOfPeptidesPossiblyObserved.Add(id));
                    }

                    // done with initial scoring; refine scores and create PSMs
                    if (idsOfPeptidesPossiblyObserved.Any())
                    {
                        int maxInitialScore = idsOfPeptidesPossiblyObserved.Max(id => scoringTable[id]) + 1;
                        while (maxInitialScore > CommonParameters.ScoreCutoff)
                        {
                            maxInitialScore--;
                            foreach (var id in idsOfPeptidesPossiblyObserved.Where(id => scoringTable[id] == maxInitialScore))
                            {
                                var candidatePeptide = peptideIndex[id];
                                double[] fragmentMasses = candidatePeptide.ProductMassesMightHaveDuplicatesAndNaNs(lp).Distinct().Where(p => !Double.IsNaN(p)).OrderBy(p => p).ToArray();
                                double peptideScore = CalculatePeptideScore(scan.TheScan, CommonParameters.ProductMassTolerance, fragmentMasses, scan.PrecursorMass, dissociationTypes, addCompIons, 0);
                                Tuple<int, double> notchAndPrecursor = Accepts(scan.PrecursorMass, candidatePeptide, terminusType, massDiffAcceptor);
                                if (notchAndPrecursor.Item1 >= 0)
                                {
                                    CompactPeptideWithModifiedMass cp = new CompactPeptideWithModifiedMass(candidatePeptide, notchAndPrecursor.Item2);

                                    if (globalPsms[i] == null)
                                        globalPsms[i] = new Psm(cp, notchAndPrecursor.Item1, peptideScore, i, scan);
                                    else
                                        globalPsms[i].AddOrReplace(cp, peptideScore, notchAndPrecursor.Item1, CommonParameters.ReportAllAmbiguity);
                                }
                            }
                            if (globalPsms[i] != null)
                                break;
                        }
                    }

                    // report search progress
                    progress++;
                    var percentProgress = (int)((progress / listOfSortedms2Scans.Length) * 100);

                    if (percentProgress > oldPercentProgress)
                    {
                        oldPercentProgress = percentProgress;
                        ReportProgress(new ProgressEventArgs(percentProgress, "Performing nonspecific search... " + currentPartition + "/" + CommonParameters.TotalPartitions, nestedIds));
                    }
                }
            });
            return new MetaMorpheusEngineResults(this);
        }

        #endregion Protected Methods

        #region Private Methods

        private Tuple<int, double> Accepts(double scanPrecursorMass, CompactPeptide peptide, TerminusType terminusType, MassDiffAcceptor searchMode)
        {
            //all masses in N and CTerminalMasses are b-ion masses, which are one water away from a full peptide
            int localminPeptideLength = CommonParameters.DigestionParams.MinPeptideLength ?? 0;
            if (terminusType == TerminusType.N)
            {
                for (int i = localminPeptideLength; i < peptide.NTerminalMasses.Count(); i++)
                {
                    double theoMass = peptide.NTerminalMasses[i] + waterMonoisotopicMass;
                    int notch = searchMode.Accepts(scanPrecursorMass, theoMass);
                    if (notch >= 0)
                    {
                        return new Tuple<int, double>(notch, theoMass);
                    }
                    else if (theoMass > scanPrecursorMass)
                    {
                        break;
                    }
                }
                //if the theoretical and experimental have the same mass
                if (peptide.NTerminalMasses.Count() > localminPeptideLength)
                {
                    double totalMass = peptide.MonoisotopicMassIncludingFixedMods;// + Constants.protonMass;
                    int notch = searchMode.Accepts(scanPrecursorMass, totalMass);
                    if (notch >= 0)
                    {
                        return new Tuple<int, double>(notch, totalMass);
                    }
                }
            }
            else//if (terminusType==TerminusType.C)
            {
                for (int i = localminPeptideLength; i < peptide.CTerminalMasses.Count(); i++)
                {
                    double theoMass = peptide.CTerminalMasses[i] + waterMonoisotopicMass;
                    int notch = searchMode.Accepts(scanPrecursorMass, theoMass);
                    if (notch >= 0)
                    {
                        return new Tuple<int, double>(notch, theoMass);
                    }
                    else if (theoMass > scanPrecursorMass)
                    {
                        break;
                    }
                }
                //if the theoretical and experimental have the same mass
                if (peptide.CTerminalMasses.Count() > localminPeptideLength)
                {
                    double totalMass = peptide.MonoisotopicMassIncludingFixedMods;// + Constants.protonMass;
                    int notch = searchMode.Accepts(scanPrecursorMass, totalMass);
                    if (notch >= 0)
                    {
                        return new Tuple<int, double>(notch, totalMass);
                    }
                }
            }
            return new Tuple<int, double>(-1, -1);
        }

        #endregion Private Methods
    }
}