using Chemistry;
using MassSpectrometry;
using MzLibUtil;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Threading.Tasks;
using System.Linq;

namespace EngineLayer.NonSpecificEnzymeSearch
{
    class NonSpecificEnzymeEngine : ModernSearch.ModernSearchEngine
    {
        private static readonly double nitrogenAtomMonoisotopicMass = PeriodicTable.GetElement("N").PrincipalIsotope.AtomicMass;
        private static readonly double oxygenAtomMonoisotopicMass = PeriodicTable.GetElement("O").PrincipalIsotope.AtomicMass;
        private static readonly double hydrogenAtomMonoisotopicMass = PeriodicTable.GetElement("H").PrincipalIsotope.AtomicMass;
        private static readonly double waterMonoisotopicMass = PeriodicTable.GetElement("H").PrincipalIsotope.AtomicMass * 2 + PeriodicTable.GetElement("O").PrincipalIsotope.AtomicMass;

        public NonSpecificEnzymeEngine(Ms2ScanWithSpecificMass[] listOfSortedms2Scans, List<CompactPeptide> peptideIndex, float[] keys, List<int>[] fragmentIndex, Tolerance fragmentTolerance, List<MassDiffAcceptor> searchModes, List<string> nestedIds, bool addCompIons, List<ProductType> lp, Protease protease, int? minPeptideLength) : base(listOfSortedms2Scans, peptideIndex, keys, fragmentIndex, fragmentTolerance, searchModes, nestedIds, addCompIons, lp, protease, minPeptideLength) { }

        protected override MetaMorpheusEngineResults RunSpecific()
        {
            Status("In modern search engine...", nestedIds);
            bool classicAntigens = false;
            double precursorToleranceDouble = 5;//default 5ppm
            int openSearchIndex = 0;
            if (searchModes.Count() > 1)
            {
                if (searchModes[0].ToString().Contains("ppmAroundZero"))
                {
                    string name = searchModes[0].ToString();
                    int index = name.IndexOf("ppmAroundZero");
                    precursorToleranceDouble = Convert.ToDouble(name.Substring(0, index));
                    classicAntigens = true;
                    openSearchIndex = 1;
                }
                else if (searchModes[1].ToString().Contains("ppmAroundZero"))
                {
                    string name = searchModes[1].ToString();
                    int index = name.IndexOf("ppmAroundZero");
                    precursorToleranceDouble = Convert.ToDouble(name.Substring(0, index));
                    classicAntigens = true;
                }
            }
            PpmTolerance precursorTolerance = new PpmTolerance(precursorToleranceDouble);
            var listOfSortedms2ScansLength = listOfSortedms2Scans.Length;
            Psm[][] newPsms = new Psm[searchModes.Count][];
            for (int i = 0; i < searchModes.Count; i++)
                newPsms[i] = new Psm[listOfSortedms2Scans.Length];

            var searchModesCount = searchModes.Count;
            var outputObject = new object();
            int scansSeen = 0;
            int old_progress = 0;
            var peptideIndexCount = peptideIndex.Count;
            Parallel.ForEach(Partitioner.Create(0, listOfSortedms2ScansLength), fff =>
            {
                List<CompactPeptide>[] bestPeptides = new List<CompactPeptide>[searchModesCount];
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
                        var searchMode = searchModes[openSearchIndex];
                        double currentBestScore = bestScores[openSearchIndex];
                        for (int possibleWinningPeptideIndex = 0; possibleWinningPeptideIndex < fullPeptideScores.Length; possibleWinningPeptideIndex++)
                        {
                            var consideredScore = fullPeptideScores[possibleWinningPeptideIndex];
                            if (consideredScore > 4)
                            {
                                CompactPeptide candidatePeptide = peptideIndex[possibleWinningPeptideIndex];
                                // Check if makes sense to add due to peptidescore!

                                //if (thisScan.TheScan.OneBasedScanNumber == 20467)
                                //{ if (possibleWinningPeptideIndex == 11)
                                //    { }
                                //}
                                if (currentBestScore > 1)
                                {
                                    // Existed! Need to compare with old match
                                    if (Math.Abs(currentBestScore - consideredScore) < 1e-9)
                                    {
                                        // Score is same, need to see if accepts and if prefer the new one
                                        double precursorMass = Accepts(thisScanprecursorMass, candidatePeptide, precursorTolerance, protease);
                                        if (precursorMass > 1)
                                        {
                                            //CompactPeptide cp = new CompactPeptide(candidatePeptide.NTerminalMasses, candidatePeptide.CTerminalMasses, precursorMass, addCompIons);
                                            bestPeptides[openSearchIndex].Add(candidatePeptide);
                                            bestNotches[openSearchIndex].Add(0);
                                        }
                                    }
                                    else if (currentBestScore < consideredScore)
                                    {
                                        // Score is better, only make sure it is acceptable
                                        double precursorMass = Accepts(thisScanprecursorMass, candidatePeptide, precursorTolerance, protease);
                                        if (precursorMass > 1)
                                        {
                                            // CompactPeptide cp = new CompactPeptide(candidatePeptide.NTerminalMasses, candidatePeptide.CTerminalMasses, precursorMass, addCompIons);
                                            bestPeptides[openSearchIndex] = new List<CompactPeptide> { candidatePeptide };
                                            bestScores[openSearchIndex] = consideredScore;
                                            bestNotches[openSearchIndex] = new List<int> { 0 };
                                            currentBestScore = consideredScore;
                                        }
                                    }
                                }
                                // Did not exist! Only make sure that it is acceptable
                                else
                                {
                                    double precursorMass = Accepts(thisScanprecursorMass, candidatePeptide, precursorTolerance, protease);
                                    if (precursorMass > 1)
                                    {
                                        //CompactPeptide cp = new CompactPeptide(candidatePeptide.NTerminalMasses, candidatePeptide.CTerminalMasses, precursorMass, addCompIons);
                                        bestPeptides[openSearchIndex] = new List<CompactPeptide> { candidatePeptide };
                                        bestScores[openSearchIndex] = consideredScore;
                                        bestNotches[openSearchIndex] = new List<int> { 0 };
                                        currentBestScore = consideredScore;
                                    }
                                }
                            }
                        }
                    }
                    else //if looking for fusion peptides
                    {
                        double currentBestScore = bestScores[openSearchIndex];
                        var searchMode = searchModes[openSearchIndex];
                        for (int possibleWinningPeptideIndex = 0; possibleWinningPeptideIndex < fullPeptideScores.Length; possibleWinningPeptideIndex++)
                        {
                            var consideredScore = fullPeptideScores[possibleWinningPeptideIndex];
                            if (consideredScore > 4)
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
                                            bestPeptides[openSearchIndex] = new List<CompactPeptide> { candidatePeptide };
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
                                        bestPeptides[openSearchIndex] = new List<CompactPeptide> { candidatePeptide };
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
                            newPsms[j][i] = new Psm(bestPeptides[j][0], bestNotches[j][0], bestScores[j], i, thisScan);
                            for (int k = 1; k < bestPeptides[j].Count; k++)
                            {
                                newPsms[j][i].Add(bestPeptides[j][k], bestNotches[j][k]);
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
                        ReportProgress(new ProgressEventArgs(new_progress, "In modern search loop", nestedIds));
                        old_progress = new_progress;
                    }
                }
            });
            return new SearchResults(newPsms, this);
        }

        private double Accepts(double scanPrecursorMass, CompactPeptide peptide, PpmTolerance precursorTolerance, Protease protease)
        {
            //all masses in N and CTerminalMasses are b-ion masses, which are one water away from a full peptide
            int localminPeptideLength = minPeptideLength ?? 0;
            if (protease.Name.Equals("singleN"))
            {
                for (int i = localminPeptideLength; i < peptide.NTerminalMasses.Count(); i++)
                {
                    double theoMass = peptide.NTerminalMasses[i] + waterMonoisotopicMass;
                    if (Math.Abs((scanPrecursorMass - theoMass) / (theoMass) * 1e6) < precursorTolerance.Value)
                    {
                        return theoMass;
                    }
                }

            }
            else//if singleC
            {
                for (int i = localminPeptideLength; i < peptide.CTerminalMasses.Count(); i++)
                {
                    double theoMass = peptide.CTerminalMasses[i] + waterMonoisotopicMass;
                    if (Math.Abs((scanPrecursorMass - theoMass) / (theoMass) * 1e6) < precursorTolerance.Value)
                    {
                        return theoMass;
                    }
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
            return 0;
        }

        private void CalculatePeptideScores(IMsDataScan<IMzSpectrum<IMzPeak>> spectrum, double[] peptideScores, double thePrecursorMass, PpmTolerance precursorTolerance)
        {
            for (int i = 0; i < spectrum.MassSpectrum.Size; i++)
            {
                var theAdd = 1 + spectrum.MassSpectrum[i].Intensity / spectrum.TotalIonCurrent;
                var experimentalPeakInDaltons = spectrum.MassSpectrum[i].Mz - Constants.protonMass;
                GeneratePeptideScores(theAdd, experimentalPeakInDaltons, peptideScores, fragmentTolerance);
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
                AbsoluteTolerance expandedFragmentTolerance = new AbsoluteTolerance(Math.Sqrt(Math.Pow(fragmentTolerance.Value,2) + Math.Pow(thePrecursorMass / 1000000 * precursorTolerance.Value,2)));
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
    }
}
