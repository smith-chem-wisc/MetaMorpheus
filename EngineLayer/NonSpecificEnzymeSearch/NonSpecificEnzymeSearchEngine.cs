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

        #endregion Private Fields

        #region Public Constructors

        public NonSpecificEnzymeSearchEngine(Psm[] globalPsms, Ms2ScanWithSpecificMass[] listOfSortedms2Scans, List<CompactPeptide> peptideIndex, List<int>[] fragmentIndex, List<ProductType> lp, int currentPartition, CommonParameters CommonParameters, bool addCompIons, MassDiffAcceptor massDiffAcceptors, List<string> nestedIds) : base(globalPsms, listOfSortedms2Scans, peptideIndex, fragmentIndex, lp, currentPartition, CommonParameters, addCompIons, massDiffAcceptors, nestedIds)
        {
        }

        #endregion Public Constructors

        #region Protected Methods

        protected override MetaMorpheusEngineResults RunSpecific()
        {
            Status("In nonspecific search engine..." + currentPartition + "/" + CommonParameters.TotalPartitions, nestedIds);
            TerminusType terminusType = ProductTypeMethod.IdentifyTerminusType(lp);
            var listOfSortedms2ScansLength = listOfSortedms2Scans.Length;

            var outputObject = new object();
            int scansSeen = 0;
            int old_progress = 0;
            var peptideIndexCount = peptideIndex.Count;
            Parallel.ForEach(Partitioner.Create(0, listOfSortedms2ScansLength), fff =>
            {
                List<CompactPeptideBase> bestPeptides;
                double bestScores;
                List<int> bestNotches;
                double[] fullPeptideScores = new double[peptideIndexCount];
                for (int i = fff.Item1; i < fff.Item2; i++)
                {
                    var thisScan = listOfSortedms2Scans[i];
                    var thisScanprecursorMass = thisScan.PrecursorMass;
                    Array.Clear(fullPeptideScores, 0, peptideIndexCount);
                    double thePrecursorMass = thisScan.PrecursorMass;
                    //CalculatePeptideScores(thisScan.TheScan, fullPeptideScores, thePrecursorMass);

                    bestPeptides = null;
                    bestScores = 0;
                    bestNotches = null;

                    for (int possibleWinningPeptideIndex = 0; possibleWinningPeptideIndex < fullPeptideScores.Length; possibleWinningPeptideIndex++)
                    {
                        var consideredScore = fullPeptideScores[possibleWinningPeptideIndex];
                        if (consideredScore > CommonParameters.ScoreCutoff) //intentionally high. 99.9% of 4-mers are present in a given UniProt database. This saves considerable time
                        {
                            CompactPeptide candidatePeptide = peptideIndex[possibleWinningPeptideIndex];
                            // Check if makes sense to add due to peptidescore!
                            var searchMode = massDiffAcceptors;
                            double currentBestScore = bestScores;
                            if (currentBestScore > 1)
                            {
                                // Existed! Need to compare with old match
                                if ((Math.Abs(currentBestScore - consideredScore) < 1e-9) && (CommonParameters.ReportAllAmbiguity || bestPeptides != null))
                                {
                                    // Score is same, need to see if accepts and if prefer the new one
                                    double precursorMass = Accepts(thisScanprecursorMass, candidatePeptide, terminusType, searchMode);
                                    if (precursorMass > 1)
                                    {
                                        CompactPeptideWithModifiedMass cp = new CompactPeptideWithModifiedMass(candidatePeptide, precursorMass);
                                        cp.SwapMonoisotopicMassWithModifiedMass();
                                        if (bestPeptides == null) //have to check, because current best score may have been found in a previous partition
                                        {
                                            bestPeptides = new List<CompactPeptideBase> { cp };
                                            bestScores = consideredScore;
                                            bestNotches = new List<int> { 0 };
                                        }
                                        else if (!bestPeptides.Contains(cp) && (globalPsms[i] == null || !globalPsms[i].CompactPeptidesContainsKey(cp)))
                                        {
                                            bestPeptides.Add(cp);
                                            bestNotches.Add(0);
                                        }
                                    }
                                }
                                else if (currentBestScore < consideredScore)
                                {
                                    // Score is better, only make sure it is acceptable
                                    double precursorMass = Accepts(thisScanprecursorMass, candidatePeptide, terminusType, searchMode);
                                    if (precursorMass > 1)
                                    {
                                        CompactPeptideWithModifiedMass cp = new CompactPeptideWithModifiedMass(candidatePeptide, precursorMass);
                                        cp.SwapMonoisotopicMassWithModifiedMass();
                                        bestPeptides = new List<CompactPeptideBase> { cp };
                                        bestScores = consideredScore;
                                        bestNotches = new List<int> { 0 };
                                        currentBestScore = consideredScore;
                                    }
                                }
                            }
                            // Did not exist! Only make sure that it is acceptable
                            else
                            {
                                double precursorMass = Accepts(thisScanprecursorMass, candidatePeptide, terminusType, searchMode);
                                if (precursorMass > 1)
                                {
                                    CompactPeptideWithModifiedMass cp = new CompactPeptideWithModifiedMass(candidatePeptide, precursorMass);
                                    cp.SwapMonoisotopicMassWithModifiedMass();
                                    bestPeptides = new List<CompactPeptideBase> { cp };
                                    bestScores = consideredScore;
                                    bestNotches = new List<int> { 0 };
                                    currentBestScore = consideredScore;
                                }
                            }
                        }
                    }
                    if (bestPeptides != null)
                    {
                        foreach (CompactPeptideBase cpb in bestPeptides)
                            (cpb as CompactPeptideWithModifiedMass).SwapMonoisotopicMassWithModifiedMass();

                        int startIndex = 0;

                        if (globalPsms[i] == null)
                        {
                            globalPsms[i] = new Psm(bestPeptides[0], bestNotches[0], bestScores, i, thisScan, CommonParameters.ExcelCompatible);
                            startIndex = 1;
                        }

                        for (int k = startIndex; k < bestPeptides.Count; k++)
                        {
                            globalPsms[i].AddOrReplace(bestPeptides[k], bestScores, bestNotches[k], CommonParameters.ReportAllAmbiguity);
                        }
                    }
                }
                lock (outputObject)
                {
                    scansSeen += fff.Item2 - fff.Item1;
                    var new_progress = (int)((double)scansSeen / (listOfSortedms2ScansLength) * 100);
                    if (new_progress > old_progress)
                    {
                        ReportProgress(new ProgressEventArgs(new_progress, "In nonspecific search loop " + currentPartition + "/" + CommonParameters.TotalPartitions, nestedIds));
                        old_progress = new_progress;
                    }
                }
            });
            return new MetaMorpheusEngineResults(this);
        }

        #endregion Protected Methods

        #region Private Methods

        private double Accepts(double scanPrecursorMass, CompactPeptide peptide, TerminusType terminusType, MassDiffAcceptor searchMode)
        {
            //all masses in N and CTerminalMasses are b-ion masses, which are one water away from a full peptide
            int localminPeptideLength = CommonParameters.DigestionParams.MinPeptideLength ?? 0;
            if (terminusType == TerminusType.N)
            {
                for (int i = localminPeptideLength; i < peptide.NTerminalMasses.Count(); i++)
                {
                    double theoMass = peptide.NTerminalMasses[i] + waterMonoisotopicMass;
                    if (searchMode.Accepts(scanPrecursorMass, theoMass) >= 0)
                    {
                        return theoMass;
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
                    if (searchMode.Accepts(scanPrecursorMass, totalMass) >= 0)
                    {
                        return totalMass;
                    }
                }
            }
            else//if (terminusType==TerminusType.C)
            {
                for (int i = localminPeptideLength; i < peptide.CTerminalMasses.Count(); i++)
                {
                    double theoMass = peptide.CTerminalMasses[i] + waterMonoisotopicMass;
                    if (searchMode.Accepts(scanPrecursorMass, theoMass) >= 0)
                    {
                        return theoMass;
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
                    if (searchMode.Accepts(scanPrecursorMass, totalMass) >= 0)
                    {
                        return totalMass;
                    }
                }
            }
            return 0;
        }

        #endregion Private Methods
    }
}