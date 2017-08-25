﻿using Chemistry;
using MassSpectrometry;
using MzLibUtil;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace EngineLayer.NonSpecificEnzymeSearch
{
    public class NonSpecificEnzymeEngine : ModernSearch.ModernSearchEngine
    {
        #region Private Fields

        private static readonly double nitrogenAtomMonoisotopicMass = PeriodicTable.GetElement("N").PrincipalIsotope.AtomicMass;
        private static readonly double oxygenAtomMonoisotopicMass = PeriodicTable.GetElement("O").PrincipalIsotope.AtomicMass;
        private static readonly double hydrogenAtomMonoisotopicMass = PeriodicTable.GetElement("H").PrincipalIsotope.AtomicMass;
        private static readonly double waterMonoisotopicMass = PeriodicTable.GetElement("H").PrincipalIsotope.AtomicMass * 2 + PeriodicTable.GetElement("O").PrincipalIsotope.AtomicMass;
        private TerminusType terminusType;

        #endregion Private Fields

        #region Public Constructors

        public NonSpecificEnzymeEngine(Psm[][] globalPsms, Ms2ScanWithSpecificMass[] listOfSortedms2Scans, List<CompactPeptide> peptideIndex, float[] keys, List<int>[] fragmentIndex, List<ProductType> lp, int currentPartition, CommonParameters CommonParameters, bool addCompIons, List<MassDiffAcceptor> massDiffAcceptors, TerminusType terminusType, List<string> nestedIds) : base(globalPsms, listOfSortedms2Scans, peptideIndex, keys, fragmentIndex, lp, currentPartition, CommonParameters, addCompIons, massDiffAcceptors, nestedIds)
        {
            this.terminusType = terminusType;
        }

        #endregion Public Constructors

        #region Protected Methods

        protected override MetaMorpheusEngineResults RunSpecific()
        {
            Status("In nonspecific search engine..." + currentPartition + "/" + CommonParameters.TotalPartitions, nestedIds);
            bool classicAntigens = false;
            double precursorToleranceDouble = 5;//default 5ppm
            int openSearchIndex = 0;
            if (massDiffAcceptors.Count() > 1)
            {
                if (massDiffAcceptors[0].ToString().Contains("ppmAroundZero"))
                {
                    string name = massDiffAcceptors[0].ToString();
                    int index = name.IndexOf("ppmAroundZero");
                    precursorToleranceDouble = Convert.ToDouble(name.Substring(0, index));
                    classicAntigens = true;
                    openSearchIndex = 1;
                }
                else if (massDiffAcceptors[1].ToString().Contains("ppmAroundZero"))
                {
                    string name = massDiffAcceptors[1].ToString();
                    int index = name.IndexOf("ppmAroundZero");
                    precursorToleranceDouble = Convert.ToDouble(name.Substring(0, index));
                    classicAntigens = true;
                }
            }
            PpmTolerance precursorTolerance = new PpmTolerance(precursorToleranceDouble);
            var listOfSortedms2ScansLength = listOfSortedms2Scans.Length;

            var searchModesCount = massDiffAcceptors.Count;
            var outputObject = new object();
            int scansSeen = 0;
            int old_progress = 0;
            var peptideIndexCount = peptideIndex.Count;
            Parallel.ForEach(Partitioner.Create(0, listOfSortedms2ScansLength), fff =>
            {
                List<CompactPeptideBase>[] bestPeptides = new List<CompactPeptideBase>[searchModesCount];
                double[] bestScores = new double[searchModesCount];
                List<int>[] bestNotches = new List<int>[searchModesCount];
                double[] fullPeptideScores = new double[peptideIndexCount];
                for (int i = fff.Item1; i < fff.Item2; i++)
                {
                    var thisScan = listOfSortedms2Scans[i];
                    var thisScanprecursorMass = thisScan.PrecursorMass;
                    Array.Clear(fullPeptideScores, 0, peptideIndexCount);
                    double thePrecursorMass = thisScan.PrecursorMass;
                    CalculatePeptideScores(thisScan.TheScan, fullPeptideScores, thePrecursorMass, precursorTolerance);

                    Array.Clear(bestPeptides, 0, searchModesCount);
                    Array.Clear(bestScores, 0, searchModesCount);
                    Array.Clear(bestNotches, 0, searchModesCount);

                    if (classicAntigens)
                    {
                        double currentBestScore = bestScores[openSearchIndex];
                        if (globalPsms[openSearchIndex][i] != null)
                            currentBestScore = globalPsms[openSearchIndex][i].Score;
                        for (int possibleWinningPeptideIndex = 0; possibleWinningPeptideIndex < fullPeptideScores.Length; possibleWinningPeptideIndex++)
                        {
                            var consideredScore = fullPeptideScores[possibleWinningPeptideIndex];
                            if (consideredScore > CommonParameters.ScoreCutoff) //intentionally high. 99.9% of 4-mers are present in a given UniProt database. This saves considerable time
                            {
                                CompactPeptide candidatePeptide = peptideIndex[possibleWinningPeptideIndex];
                                // Check if makes sense to add due to peptidescore!
                                if (currentBestScore > 1)
                                {
                                    // Existed! Need to compare with old match
                                    if (Math.Abs(currentBestScore - consideredScore) < 1e-9)
                                    {
                                        // Score is same, need to see if accepts and if prefer the new one
                                        double precursorMass = Accepts(thisScanprecursorMass, candidatePeptide, precursorTolerance, terminusType);
                                        if (precursorMass > 1)
                                        {
                                            CompactPeptideWithModifiedMass cp = new CompactPeptideWithModifiedMass(candidatePeptide, precursorMass);
                                            cp.SwapMonoisotopicMassWithModifiedMass();
                                            if (bestPeptides[openSearchIndex] == null)
                                            {
                                                bestPeptides[openSearchIndex] = new List<CompactPeptideBase> { cp };
                                                bestScores[openSearchIndex] = consideredScore;
                                                bestNotches[openSearchIndex] = new List<int> { 0 };
                                            }
                                            else if (!bestPeptides[openSearchIndex].Contains(cp) && (globalPsms[openSearchIndex][i] == null || !globalPsms[openSearchIndex][i].CompactPeptidesContainsKey(cp)))
                                            {
                                                bestPeptides[openSearchIndex].Add(cp);
                                                bestNotches[openSearchIndex].Add(0);
                                            }
                                        }
                                    }
                                    else if (currentBestScore < consideredScore)
                                    {
                                        // Score is better, only make sure it is acceptable
                                        double precursorMass = Accepts(thisScanprecursorMass, candidatePeptide, precursorTolerance, terminusType);
                                        if (precursorMass > 1)
                                        {
                                            CompactPeptideWithModifiedMass cp = new CompactPeptideWithModifiedMass(candidatePeptide, precursorMass);
                                            cp.SwapMonoisotopicMassWithModifiedMass();
                                            bestPeptides[openSearchIndex] = new List<CompactPeptideBase> { cp };
                                            bestScores[openSearchIndex] = consideredScore;
                                            bestNotches[openSearchIndex] = new List<int> { 0 };
                                            currentBestScore = consideredScore;
                                        }
                                    }
                                }
                                // Did not exist! Only make sure that it is acceptable
                                else
                                {
                                    double precursorMass = Accepts(thisScanprecursorMass, candidatePeptide, precursorTolerance, terminusType);
                                    if (precursorMass > 1)
                                    {
                                        CompactPeptideWithModifiedMass cp = new CompactPeptideWithModifiedMass(candidatePeptide, precursorMass);
                                        cp.SwapMonoisotopicMassWithModifiedMass();
                                        bestPeptides[openSearchIndex] = new List<CompactPeptideBase> { cp };
                                        bestScores[openSearchIndex] = consideredScore;
                                        bestNotches[openSearchIndex] = new List<int> { 0 };
                                        currentBestScore = consideredScore;
                                    }
                                }
                            }
                        }
                        if (bestPeptides[openSearchIndex] != null)
                            foreach (CompactPeptideBase cpb in bestPeptides[openSearchIndex])
                                (cpb as CompactPeptideWithModifiedMass).SwapMonoisotopicMassWithModifiedMass();
                    }
                    else //(assumes open search)
                    {
                        double currentBestScore = bestScores[openSearchIndex];
                        var searchMode = massDiffAcceptors[openSearchIndex];
                        for (int possibleWinningPeptideIndex = 0; possibleWinningPeptideIndex < fullPeptideScores.Length; possibleWinningPeptideIndex++)
                        {
                            var consideredScore = fullPeptideScores[possibleWinningPeptideIndex];
                            if (consideredScore > CommonParameters.ScoreCutoff) //intentionally high. 99.9% of 4-mers are present in a given UniProt database. This saves considerable time
                            {
                                CompactPeptide candidatePeptide = peptideIndex[possibleWinningPeptideIndex];
                                // Check if makes sense to add due to peptidescore!

                                if (currentBestScore > 1)
                                {
                                    // Existed! Need to compare with old match
                                    if (Math.Abs(currentBestScore - consideredScore) < 1e-9)
                                    {
                                        // Score is same, need to see if accepts and if prefer the new one
                                        int notch = searchMode.Accepts(thisScanprecursorMass, candidatePeptide.MonoisotopicMassIncludingFixedMods);
                                        if (notch >= 0)
                                        {
                                            bestPeptides[openSearchIndex].Add(candidatePeptide);
                                            bestNotches[openSearchIndex].Add(notch);
                                        }
                                    }
                                    else if (currentBestScore < consideredScore)
                                    {
                                        // Score is better, only make sure it is acceptable
                                        int notch = searchMode.Accepts(thisScanprecursorMass, candidatePeptide.MonoisotopicMassIncludingFixedMods);
                                        if (notch >= 0)
                                        {
                                            bestPeptides[openSearchIndex] = new List<CompactPeptideBase> { candidatePeptide };
                                            bestScores[openSearchIndex] = consideredScore;
                                            bestNotches[openSearchIndex] = new List<int> { notch };
                                            currentBestScore = consideredScore;
                                        }
                                    }
                                }
                                // Did not exist! Only make sure that it is acceptable
                                else
                                {
                                    int notch = searchMode.Accepts(thisScanprecursorMass, candidatePeptide.MonoisotopicMassIncludingFixedMods);
                                    if (notch >= 0)
                                    {
                                        bestPeptides[openSearchIndex] = new List<CompactPeptideBase> { candidatePeptide };
                                        bestScores[openSearchIndex] = consideredScore;
                                        bestNotches[openSearchIndex] = new List<int> { notch };
                                        currentBestScore = consideredScore;
                                    }
                                }
                            }
                        }
                    }
                    for (int j = 0; j < searchModesCount; j++)
                    {
                        if (bestPeptides[j] != null)
                        {
                            int startIndex = 0;

                            if (globalPsms[j][i] == null)
                            {
                                globalPsms[j][i] = new Psm(bestPeptides[j][0], bestNotches[j][0], bestScores[j], i, thisScan);
                                startIndex = 1;
                            }

                            for (int k = startIndex; k < bestPeptides[j].Count; k++)
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
                        ReportProgress(new ProgressEventArgs(new_progress, "In nonspecific search loop " + currentPartition + "/" + CommonParameters.TotalPartitions, nestedIds));
                        old_progress = new_progress;
                    }
                }
            });
            return new MetaMorpheusEngineResults(this);
        }

        #endregion Protected Methods

        #region Private Methods

        private double Accepts(double scanPrecursorMass, CompactPeptide peptide, PpmTolerance precursorTolerance, TerminusType terminusType)
        {
            //all masses in N and CTerminalMasses are b-ion masses, which are one water away from a full peptide
            int localminPeptideLength = CommonParameters.MinPeptideLength ?? 0;
            if (terminusType == TerminusType.N)
            {
                for (int i = localminPeptideLength; i < peptide.NTerminalMasses.Count(); i++)
                {
                    double theoMass = peptide.NTerminalMasses[i] + waterMonoisotopicMass;
                    if (Math.Abs((scanPrecursorMass - theoMass) / (theoMass) * 1e6) < precursorTolerance.Value)
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
                    if (Math.Abs((scanPrecursorMass - totalMass) / (totalMass) * 1e6) < precursorTolerance.Value)
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
                    if (Math.Abs((scanPrecursorMass - theoMass) / (theoMass) * 1e6) < precursorTolerance.Value)
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
                    if (Math.Abs((scanPrecursorMass - totalMass) / (totalMass) * 1e6) < precursorTolerance.Value)
                    {
                        return totalMass;
                    }
                }
            }
            return 0;
        }

        private void CalculatePeptideScores(IMsDataScan<IMzSpectrum<IMzPeak>> spectrum, double[] peptideScores, double thePrecursorMass, PpmTolerance precursorTolerance)
        {
            for (int i = 0; i < spectrum.MassSpectrum.Size; i++)
            {
                var theAdd = 1 + spectrum.MassSpectrum[i].Intensity / spectrum.TotalIonCurrent;
                var experimentalPeakInDaltons = spectrum.MassSpectrum[i].Mz - Constants.protonMass;
                GeneratePeptideScores(theAdd, experimentalPeakInDaltons, peptideScores, CommonParameters.ProductMassTolerance);
            }
            if (addCompIons)
            {
                List<IMzPeak> experimentalPeaks = new List<IMzPeak>();
                //If HCD
                if (lp.Contains(ProductType.B) | lp.Contains(ProductType.Y))
                {
                    for (int i = 0; i < spectrum.MassSpectrum.Size; i++)
                    {
                        experimentalPeaks.Add(new MzPeak((thePrecursorMass - spectrum.MassSpectrum[i].Mz + Constants.protonMass), (spectrum.MassSpectrum[i].Intensity)));
                    }
                }
                //If ETD
                if (lp.Contains(ProductType.C) | lp.Contains(ProductType.Zdot))
                {
                    for (int i = 0; i < spectrum.MassSpectrum.Size; i++)
                    {
                        experimentalPeaks.Add(new MzPeak((thePrecursorMass - spectrum.MassSpectrum[i].Mz + Constants.protonMass * 2), (spectrum.MassSpectrum[i].Intensity)));
                    }
                }

                IEnumerable<IMzPeak> sortedPeaksMZ = experimentalPeaks.OrderBy(x => x.Mz);
                //propogation of error from precursor mass and complementary product mass
                AbsoluteTolerance expandedFragmentTolerance = new AbsoluteTolerance(Math.Sqrt(Math.Pow(CommonParameters.ProductMassTolerance.Value, 2) + Math.Pow(thePrecursorMass / 1000000 * precursorTolerance.Value, 2)));
                foreach (IMzPeak experimentalPeak in sortedPeaksMZ)
                {
                    var theAdd = 1 + experimentalPeak.Intensity / spectrum.TotalIonCurrent;
                    GeneratePeptideScores(theAdd, experimentalPeak.Mz, peptideScores, expandedFragmentTolerance);
                }
            }
        }

        private void GeneratePeptideScores(double theAdd, double experimentalPeakInDaltons, double[] peptideScores, Tolerance fragmentTolerance)
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