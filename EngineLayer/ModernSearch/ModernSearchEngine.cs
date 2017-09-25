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
        protected readonly MassDiffAcceptor massDiffAcceptors;
        protected readonly List<DissociationType> dissociationTypes;

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
            this.massDiffAcceptors = massDiffAcceptors;
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

            int fragmentBinsPerDalton = 1000;
            
            Parallel.ForEach(Partitioner.Create(0, listOfSortedms2Scans.Length), new ParallelOptions { MaxDegreeOfParallelism = CommonParameters.MaxThreadsToUse }, range =>
            {
                byte[] scoringTable = new byte[peptideIndex.Count];
                HashSet<int> idsOfHighestScoringPeptides = new HashSet<int>();
                int highestScore;

                for (int i = range.Item1; i < range.Item2; i++)
                {
                    // empty the scoring table to score the new scan (conserves memory compared to allocating a new array)
                    Array.Clear(scoringTable, 0, scoringTable.Length);
                    idsOfHighestScoringPeptides.Clear();
                    var scan = listOfSortedms2Scans[i];
                    highestScore = intScoreCutoff;

                    // filter ms2 fragment peaks by intensity
                    int numFragmentsToUse = 0;
                    if (CommonParameters.TopNpeaks != null)
                        numFragmentsToUse = (int)CommonParameters.TopNpeaks;
                    else
                        numFragmentsToUse = scan.NumPeaks;

                    var peaks = scan.TheScan.MassSpectrum.FilterByNumberOfMostIntense(numFragmentsToUse).ToList();
                    double largestIntensity = scan.TheScan.MassSpectrum.YofPeakWithHighestY;

                    // get allowed precursor masses
                    var t = massDiffAcceptors.GetAllowedPrecursorMassIntervals(scan.PrecursorMass);
                    double lowestMassPeptideToLookFor = t.Min(p => p.allowedInterval.Minimum);
                    double highestMassPeptideToLookFor = t.Max(p => p.allowedInterval.Maximum);
                    
                    // search peaks for matches
                    for (int k = 0; k < peaks.Count; k++)
                    {
                        if (CommonParameters.MinRatio == null || (peaks[k].Intensity / largestIntensity) >= CommonParameters.MinRatio)
                        {
                            // assume charge state 1 to calculate mz tolerance
                            var mzTolerance = (CommonParameters.ProductMassTolerance.Value / 1e6) * peaks[k].Mz;
                            int fragmentFloorMz = (int)Math.Floor((peaks[k].Mz - mzTolerance) * fragmentBinsPerDalton);
                            int fragmentCeilingMz = (int)Math.Ceiling((peaks[k].Mz + mzTolerance) * fragmentBinsPerDalton);

                            // get all theoretical fragments this experimental fragment could be
                            for (int fragmentBin = fragmentFloorMz; fragmentBin <= fragmentCeilingMz; fragmentBin++)
                            {
                                if (fragmentIndex[fragmentBin].Count != 0)
                                {
                                    List<int> peptideIdsInThisBin = fragmentIndex[fragmentBin];
                                    int m = 0;
                                    if (lowestMassPeptideToLookFor !=0)
                                    {
                                        int y = 0;
                                        int r = peptideIdsInThisBin.Count - 1;

                                        // binary search in the fragment bin for lowest acceptable precursor mass
                                        while (y <= r)
                                        {
                                            m = y + ((r - y) / 2);

                                            if (r - y < 2)
                                                break;
                                            if (peptideIndex[peptideIdsInThisBin[m]].MonoisotopicMassIncludingFixedMods < lowestMassPeptideToLookFor)
                                                y = m + 1;
                                            else
                                                r = m - 1;
                                        }
                                        if (m > 0)
                                            m--;
                                    }

                                    // add +1 score for each peptide candidate in the scoring table up to the maximum allowed precursor mass
                                    if (highestMassPeptideToLookFor > 0)
                                    {
                                        for (int h = m; h < peptideIdsInThisBin.Count; h++)
                                        {
                                            int pepId = peptideIdsInThisBin[h];
                                            scoringTable[pepId]++;
                                            if (peptideIndex[pepId].MonoisotopicMassIncludingFixedMods > highestMassPeptideToLookFor)
                                                break;
                                        }
                                    }
                                    else
                                    {
                                        for (int h = m; h < peptideIdsInThisBin.Count; h++)
                                        {
                                            scoringTable[peptideIdsInThisBin[h]]++;
                                        }
                                    }
                                }
                            }
                        }
                    }

                    //filter theoreticals into only the highests scoring ones
                    for (int j = 0; j < scoringTable.Length; j++)
                    {
                        if(scoringTable[j]<highestScore)
                        {
                        }
                        else if(scoringTable[j]==highestScore)
                        {
                            idsOfHighestScoringPeptides.Add(j);
                        }
                        else
                        {
                            idsOfHighestScoringPeptides.Clear();
                            idsOfHighestScoringPeptides.Add(j);
                            highestScore = scoringTable[j];
                        }
                    }

                    foreach(int id in idsOfHighestScoringPeptides)
                    { 
                        var candidatePeptide = peptideIndex[id];
                        int notch = massDiffAcceptors.Accepts(scan.PrecursorMass, candidatePeptide.MonoisotopicMassIncludingFixedMods);

                        if (notch >= 0)
                        {
                            double[] fragmentMasses = candidatePeptide.ProductMassesMightHaveDuplicatesAndNaNs(lp).Distinct().Where(p => !Double.IsNaN(p)).OrderBy(p => p).ToArray();
                            double peptideScore = CalculatePeptideScore(scan.TheScan, CommonParameters.ProductMassTolerance, fragmentMasses, scan.PrecursorMass, dissociationTypes, addCompIons);

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

        /*
        protected void CalculatePeptideScores(IMsDataScan<IMzSpectrum<IMzPeak>> spectrum, double[] peptideScores, double thePrecursorMass)
        {
            //create previous variables to determine if peaks can be sequestered
            double previousTheAdd = 1 + spectrum.MassSpectrum.YArray[0] / spectrum.TotalIonCurrent;
            double previousExperimentalPeakInDaltons = spectrum.MassSpectrum.XArray[0] - Constants.protonMass;
            double previousMinRange = CommonParameters.ProductMassTolerance.GetMinimumValue(previousExperimentalPeakInDaltons);
            double previousMaxRange;
            //search observed peaks
            for (int i = 1; i < spectrum.MassSpectrum.Size; i++)
            {
                double experimentalPeakInDaltons = spectrum.MassSpectrum.XArray[i] - Constants.protonMass;
                if (CommonParameters.ProductMassTolerance.Within(previousExperimentalPeakInDaltons, experimentalPeakInDaltons))
                {
                    previousTheAdd += spectrum.MassSpectrum.YArray[i] / spectrum.TotalIonCurrent; //open to debate, currently sum intensities of all peaks within tolerance like it was low res
                }
                else
                {
                    previousMaxRange = CommonParameters.ProductMassTolerance.GetMaximumValue(previousExperimentalPeakInDaltons);
                    FindPeakMatches(previousTheAdd, previousMinRange, previousMaxRange, peptideScores);
                    previousTheAdd = 1 + spectrum.MassSpectrum.YArray[i] / spectrum.TotalIonCurrent;
                    previousMinRange = CommonParameters.ProductMassTolerance.GetMinimumValue(experimentalPeakInDaltons);
                }
                previousExperimentalPeakInDaltons = experimentalPeakInDaltons;
            }
            previousMaxRange = CommonParameters.ProductMassTolerance.GetMaximumValue(previousExperimentalPeakInDaltons);
            FindPeakMatches(previousTheAdd, previousMinRange, previousMaxRange, peptideScores);

            //generate experimental complementary ions if specified
            if (addCompIons)
            {
                //okay, we're not actually adding in complementary m/z peaks, we're doing a shortcut and just straight up adding the mass assuming that they're z=1
                int numCompIons = spectrum.MassSpectrum.Size;
                (double mass, double intensity)[] complementaryIons = new(double mass, double intensity)[numCompIons];

                foreach (DissociationType dissociationType in dissociationTypes)
                {
                    if (complementaryIonConversionDictionary.TryGetValue(dissociationType, out double protonMassShift))
                    {
                        double massShiftForComplementaryConversion = thePrecursorMass + protonMassShift; //mass shift needed to reobtain the original product ion for calculating tolerance
                        for (int i = numCompIons - 1; i >= 0; i--)
                            complementaryIons[numCompIons - i - 1] = (massShiftForComplementaryConversion - spectrum.MassSpectrum.XArray[i], spectrum.MassSpectrum.YArray[i]);

                        //propogation of error from precursor mass and complementary product mass
                        //IMPLEMENT AbsoluteTolerance expandedFragmentTolerance = new AbsoluteTolerance(Math.Sqrt(Math.Pow(CommonParameters.ProductMassTolerance.Value, 2) + Math.Pow(thePrecursorMass / 1000000 * precursorTolerance.Value, 2)));
                        previousTheAdd = 1 + complementaryIons[0].intensity / spectrum.TotalIonCurrent;
                        //we already subtracted that proton, so don't add it again (unit test should break if you do!)
                        previousExperimentalPeakInDaltons = complementaryIons[0].mass;
                        //need to use original tolerance since it's mass based.
                        double previousOriginalMassInDaltons = massShiftForComplementaryConversion - previousExperimentalPeakInDaltons;
                        previousMinRange = previousExperimentalPeakInDaltons - previousOriginalMassInDaltons + CommonParameters.ProductMassTolerance.GetMinimumValue(previousOriginalMassInDaltons);
                        for (int i = 1; i < complementaryIons.Length; i++)
                        {
                            //we already subtracted that proton when making comp ions, so don't add it again (unit test should break if you do!)
                            double experimentalPeakInDaltons = complementaryIons[i].mass;
                            double originalMassInDaltons = massShiftForComplementaryConversion - experimentalPeakInDaltons;
                            if (CommonParameters.ProductMassTolerance.Within(previousOriginalMassInDaltons, originalMassInDaltons))
                            {
                                previousTheAdd += complementaryIons[i].intensity / spectrum.TotalIonCurrent; //open to debate, currently sum intensities of all peaks within tolerance like it was low res. Classic search takes first intensity.
                            }
                            else
                            {
                                previousMaxRange = previousExperimentalPeakInDaltons - previousOriginalMassInDaltons + CommonParameters.ProductMassTolerance.GetMaximumValue(previousOriginalMassInDaltons);
                                FindPeakMatches(previousTheAdd, previousMinRange, previousMaxRange, peptideScores);
                                previousTheAdd = 1 + complementaryIons[i].intensity / spectrum.TotalIonCurrent;
                                previousMinRange = experimentalPeakInDaltons - originalMassInDaltons + CommonParameters.ProductMassTolerance.GetMinimumValue(originalMassInDaltons);
                            }
                            previousExperimentalPeakInDaltons = experimentalPeakInDaltons;
                            previousOriginalMassInDaltons = massShiftForComplementaryConversion - previousExperimentalPeakInDaltons;
                        }
                        previousMaxRange = previousExperimentalPeakInDaltons - previousOriginalMassInDaltons + CommonParameters.ProductMassTolerance.GetMaximumValue(previousOriginalMassInDaltons);
                        FindPeakMatches(previousTheAdd, previousMinRange, previousMaxRange, peptideScores);
                    }
                    else
                    {
                        throw new NotImplementedException();
                    }
                }
            }
        }
        */

        /*
        protected void FindPeakMatches(double theAdd, double min, double max, double[] peptideScores)
        {
            float closestPeak;
            int ipos = Array.BinarySearch(keys, (float)min);
            if (ipos < 0)
                ipos = ~ipos;

            while (ipos < keys.Length)
            {
                closestPeak = keys[ipos];
                if (closestPeak < max)
                {
                    foreach (int heh in fragmentIndex[ipos])
                        peptideScores[heh] += theAdd;
                }
                else
                    break;
                ipos++;
            }
        }
        */

        #endregion Protected Methods
    }
}