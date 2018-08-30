using Chemistry;
using EngineLayer.ModernSearch;
using MassSpectrometry;
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
        private static readonly double nitrogenAtomMonoisotopicMass = PeriodicTable.GetElement("N").PrincipalIsotope.AtomicMass;
        private static readonly double oxygenAtomMonoisotopicMass = PeriodicTable.GetElement("O").PrincipalIsotope.AtomicMass;
        private static readonly double hydrogenAtomMonoisotopicMass = PeriodicTable.GetElement("H").PrincipalIsotope.AtomicMass;
        private static readonly double waterMonoisotopicMass = PeriodicTable.GetElement("H").PrincipalIsotope.AtomicMass * 2 + PeriodicTable.GetElement("O").PrincipalIsotope.AtomicMass;
        private static readonly int bBinShift = (int)Math.Round((waterMonoisotopicMass) * FragmentBinsPerDalton);
        private static readonly int cBinShift = (int)Math.Round((nitrogenAtomMonoisotopicMass + 3 * hydrogenAtomMonoisotopicMass) * FragmentBinsPerDalton);
        private static readonly int zdotBinShift = (int)Math.Round((oxygenAtomMonoisotopicMass - nitrogenAtomMonoisotopicMass) * FragmentBinsPerDalton);
        private readonly List<int>[] fragmentIndexPrecursor;

        public NonSpecificEnzymeSearchEngine(PeptideSpectralMatch[] globalPsms, Ms2ScanWithSpecificMass[] listOfSortedms2Scans, List<PeptideWithSetModifications> peptideIndex, List<int>[] fragmentIndex, List<int>[] fragmentIndexPrecursor, List<ProductType> lp, int currentPartition, CommonParameters CommonParameters, MassDiffAcceptor massDiffAcceptor, double maximumMassThatFragmentIonScoreIsDoubled, List<string> nestedIds) : base(globalPsms, listOfSortedms2Scans, peptideIndex, fragmentIndex, lp, currentPartition, CommonParameters, massDiffAcceptor, maximumMassThatFragmentIonScoreIsDoubled, nestedIds)
        {
            this.fragmentIndexPrecursor = fragmentIndexPrecursor;
        }

        protected override MetaMorpheusEngineResults RunSpecific()
        {
            double progress = 0;
            int oldPercentProgress = 0;
            ReportProgress(new ProgressEventArgs(oldPercentProgress, "Performing nonspecific search... " + CurrentPartition + "/" + commonParameters.TotalPartitions, nestedIds));

            byte byteScoreCutoff = (byte)commonParameters.ScoreCutoff;

            Parallel.ForEach(Partitioner.Create(0, ListOfSortedms2Scans.Length), new ParallelOptions { MaxDegreeOfParallelism = commonParameters.MaxThreadsToUsePerFile }, range =>
            {
                byte[] scoringTable = new byte[PeptideIndex.Count];
                HashSet<int> idsOfPeptidesPossiblyObserved = new HashSet<int>();

                for (int i = range.Item1; i < range.Item2; i++)
                {
                    // empty the scoring table to score the new scan (conserves memory compared to allocating a new array)
                    Array.Clear(scoringTable, 0, scoringTable.Length);
                    idsOfPeptidesPossiblyObserved.Clear();
                    var scan = ListOfSortedms2Scans[i];

                    //get bins to add points to
                    List<int> allBinsToSearch = GetBinsToSearch(scan, commonParameters, FragmentIndex);

                    for (int j = 0; j < allBinsToSearch.Count; j++)
                        base.FragmentIndex[allBinsToSearch[j]].ForEach(id => scoringTable[id]++);

                    //populate ids of possibly observed with those containing allowed precursor masses
                    List<int> binsToSearch = new List<int>();
                    int obsPrecursorFloorMz = (int)Math.Floor(commonParameters.PrecursorMassTolerance.GetMinimumValue(scan.PrecursorMass) * FragmentBinsPerDalton);
                    int obsPrecursorCeilingMz = (int)Math.Ceiling(commonParameters.PrecursorMassTolerance.GetMaximumValue(scan.PrecursorMass) * FragmentBinsPerDalton);
                    for (int fragmentBin = obsPrecursorFloorMz; fragmentBin <= obsPrecursorCeilingMz; fragmentBin++)
                        binsToSearch.Add(fragmentBin);

                    foreach (ProductType pt in ProductTypes)
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
                            if (bin < base.FragmentIndex.Length && base.FragmentIndex[bin] != null)
                                base.FragmentIndex[bin].ForEach(id => idsOfPeptidesPossiblyObserved.Add(id));
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
                        while (maxInitialScore > commonParameters.ScoreCutoff)
                        {
                            maxInitialScore--;
                            foreach (var id in idsOfPeptidesPossiblyObserved.Where(id => scoringTable[id] == maxInitialScore))
                            {
                                PeptideWithSetModifications peptide = PeptideIndex[id];

                                List<Product> peptideTheorProducts = peptide.Fragment(commonParameters.DissociationType, FragmentationTerminus.Both).ToList();
                                List<MatchedFragmentIon> matchedIons = MatchFragmentIons(scan.TheScan.MassSpectrum, peptideTheorProducts, commonParameters);

                                if (commonParameters.AddCompIons)
                                {
                                    MzSpectrum complementarySpectrum = GenerateComplementarySpectrum(scan.TheScan.MassSpectrum, scan.PrecursorMass, commonParameters.DissociationType);
                                    matchedIons.AddRange(MatchFragmentIons(complementarySpectrum, peptideTheorProducts, commonParameters));
                                }

                                double thisScore = CalculatePeptideScore(scan.TheScan, matchedIons, MaxMassThatFragmentIonScoreIsDoubled);

                                Tuple<int, double> notchAndPrecursor = Accepts(scan.PrecursorMass, peptide, commonParameters.FragmentationTerminus, MassDiffAcceptor);
                                if (notchAndPrecursor.Item1 >= 0)
                                {
                                    if (PeptideSpectralMatches[i] == null)
                                        PeptideSpectralMatches[i] = new PeptideSpectralMatch(peptide, notchAndPrecursor.Item1, thisScore, i, scan, commonParameters.DigestionParams, matchedIons);
                                    else
                                        PeptideSpectralMatches[i].AddOrReplace(peptide, thisScore, notchAndPrecursor.Item1, commonParameters.ReportAllAmbiguity, matchedIons);
                                }
                            }
                            if (PeptideSpectralMatches[i] != null)
                                break;
                        }
                    }

                    // report search progress
                    progress++;
                    var percentProgress = (int)((progress / ListOfSortedms2Scans.Length) * 100);

                    if (percentProgress > oldPercentProgress)
                    {
                        oldPercentProgress = percentProgress;
                        ReportProgress(new ProgressEventArgs(percentProgress, "Performing nonspecific search... " + CurrentPartition + "/" + commonParameters.TotalPartitions, nestedIds));
                    }
                }
            });
            return new MetaMorpheusEngineResults(this);
        }

        private Tuple<int, double> Accepts(double scanPrecursorMass, PeptideWithSetModifications peptide, FragmentationTerminus fragmentationTerminus, MassDiffAcceptor searchMode)
        {
            //all masses in N and CTerminalMasses are b-ion masses, which are one water away from a full peptide
            int localminPeptideLength = commonParameters.DigestionParams.MinPeptideLength;
            if (fragmentationTerminus == FragmentationTerminus.N)
            {
                double[] NterminalMasses = peptide.Fragment(commonParameters.DissociationType, FragmentationTerminus.N).Select(m => m.NeutralMass).ToArray();

                for (int i = localminPeptideLength; i < NterminalMasses.Count(); i++)
                {
                    double theoMass = NterminalMasses[i] + waterMonoisotopicMass;
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
                if (NterminalMasses.Count() > localminPeptideLength)
                {
                    double totalMass = peptide.MonoisotopicMass;// + Constants.ProtonMass;
                    int notch = searchMode.Accepts(scanPrecursorMass, totalMass);
                    if (notch >= 0)
                    {
                        return new Tuple<int, double>(notch, totalMass);
                    }
                }
            }
            else//if (terminusType==TerminusType.C)
            {
                double[] CterminalMasses = peptide.Fragment(commonParameters.DissociationType, FragmentationTerminus.C).Select(m => m.NeutralMass).ToArray();

                for (int i = localminPeptideLength; i < CterminalMasses.Count(); i++)
                {
                    double theoMass = CterminalMasses[i] + waterMonoisotopicMass;
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
                if (CterminalMasses.Count() > localminPeptideLength)
                {
                    double totalMass = peptide.MonoisotopicMass;// + Constants.ProtonMass;
                    int notch = searchMode.Accepts(scanPrecursorMass, totalMass);
                    if (notch >= 0)
                    {
                        return new Tuple<int, double>(notch, totalMass);
                    }
                }
            }
            return new Tuple<int, double>(-1, -1);
        }
    }
}