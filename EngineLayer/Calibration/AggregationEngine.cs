using MassSpectrometry;
using System;
using System.Collections.Generic;
using System.Linq;
using MzLibUtil;
using Chemistry;

namespace EngineLayer.Calibration
{
    public class AggregationEngine : MetaMorpheusEngine
    {
        private readonly double MaxRetentionTimeDifferenceAllowedInMinutes;
        private readonly double MinCosineScoreAllowed;
        private readonly Tolerance PrecursorTolerance;
        private readonly Tolerance ProductTolerance;
        private MsDataFile OriginalFile;
        private const int NumberOfStrikesBeforeOut = 2;

        public Tolerance SuggestedPrecursorTolerance { get; private set; }
        public Tolerance SuggestedProductTolerance { get; private set; }
        public MsDataFile AggregatedDataFile { get; private set; }
        public string OriginalFilePath { get; private set; }

        public AggregationEngine(MsDataFile originalFile, string originalFilePath, CommonParameters commonParameters, List<string> nestedIds, double maxRetentionTimeDifferenceAllowedInMinutes, double minCosineScoreAllowed) : base(commonParameters, nestedIds)
        {
            OriginalFile = originalFile;
            MaxRetentionTimeDifferenceAllowedInMinutes = maxRetentionTimeDifferenceAllowedInMinutes;
            MinCosineScoreAllowed = minCosineScoreAllowed;
            PrecursorTolerance = commonParameters.PrecursorMassTolerance;
            ProductTolerance = commonParameters.ProductMassTolerance;
        }

        protected override MetaMorpheusEngineResults RunSpecific()
        {
            //MS1: this engine matches adjacent MS1 peaks that fall within a specified tolerance
            //If an MS1 peak does not have an adjacent match, then it is removed           
            //average the m/zs found across scans and update each value with the average

            //MS2: this engine finds multiple MS2 scans that share a precursor as determined by precursor mass (not m/z) and dot product of the ms2 scans
            //ms2 scans are deconvoluted (i.e. the same ms2 scan is present multiple times for coisolation) prior to this process.
            //Matches are used to create a consensus spectrum, where unshared peaks are removed. Shared peaks are averaged for intensity and m/z.

            int oldPercentProgress = 0; //progress keeper

            Status("Averaging MS1 spectra");
            MsDataScan[] allScans = OriginalFile.GetAllScansList().ToArray();
            //we are ONLY averaging m/z and NOT intensity. Intensity averaging could convolute downstream quantification
            MsDataScan[] ms1scans = OriginalFile.GetMS1Scans().ToArray();
            //we have a set of peaks in the ms1 scan, and we'll cycle through until:
            //-we have two consecutive ms1 scans that do not contain a peak, 
            //-we reach the end of the file
            double[][] ms1mzs = new double[ms1scans.Length][]; //quickly have mzs on hand            
            List<double>[][] allMs1PeaksFound = new List<double>[ms1scans.Length][];
            //^this is a tricky index. Each ms1 scan has an index of List<double>[], 
            //where each peak of the ms1 scan has a List<double> that contains all the grouped mzs for that peak.

            for (int i = 0; i < ms1scans.Length; i++) //populate arrays
            {
                MzSpectrum spectrum = ms1scans[i].MassSpectrum;
                double[] mzArray = spectrum.XArray;
                ms1mzs[i] = mzArray;
                allMs1PeaksFound[i] = new List<double>[mzArray.Length];
            }

            //go through every scan
            for (int seedScanIndex = 0; seedScanIndex < ms1scans.Length; seedScanIndex++) //foreach ms1 scan
            {
                //check if canceled
                if (GlobalVariables.StopLoops)
                {
                    return new MetaMorpheusEngineResults(this);
                }

                double[] seedMzs = ms1mzs[seedScanIndex]; //grab current ms1scan
                int[] numberOfStrikesForEachPeak = new int[seedMzs.Length]; //keep track of the number of missed ms1 scans for each peak 

                //see what peaks have been found already and populate or mark them
                List<double>[] seedPeaksFound = allMs1PeaksFound[seedScanIndex];
                for (int seedPeakIndex = 0; seedPeakIndex < seedMzs.Length; seedPeakIndex++)
                {
                    if (seedPeaksFound[seedPeakIndex] == null) //this is a new peak for us that was not previously seen
                    {
                        seedPeaksFound[seedPeakIndex] = new List<double> { seedMzs[seedPeakIndex] };
                    }
                    else //we saw this already, just use that existing list
                    {
                        numberOfStrikesForEachPeak[seedPeakIndex] = NumberOfStrikesBeforeOut;
                    }
                }

                //start cycling through other ms1 scans
                for (int branchScanIndex = seedScanIndex + 1; branchScanIndex < ms1mzs.Length; branchScanIndex++)
                {
                    bool done = true; //this is used to determine if all of the seed peaks for the current scan have been completed
                    double[] branchMzs = ms1mzs[branchScanIndex];

                    int seedPeakIndex = 0;
                    int branchPeakIndex = 0;
                    //foreach peak in the seed spectrum, find matches in the branch spectrum
                    for (; seedPeakIndex < seedMzs.Length; seedPeakIndex++)
                    {
                        if (numberOfStrikesForEachPeak[seedPeakIndex] != NumberOfStrikesBeforeOut) //if this seed is not dead to us
                        {
                            double seedMz = seedMzs[seedPeakIndex];
                            //see if it has a buddy! Find the index of the peak (could use a binary search for this)
                            for (; branchPeakIndex < branchMzs.Length; branchPeakIndex++)
                            {
                                if (branchMzs[branchPeakIndex] > seedMz) //we found it, or went one too far
                                {
                                    break;
                                }
                            }

                            //backup the branchIndex if needed (if we went one too far)
                            if (branchPeakIndex == branchMzs.Length
                                || (branchPeakIndex != 0 && branchMzs[branchPeakIndex] - seedMz > seedMz - branchMzs[branchPeakIndex - 1]))
                            {
                                branchPeakIndex--;
                            }

                            double branchMz = branchMzs[branchPeakIndex];
                            if (ProductTolerance.Within(seedMz, branchMz) && SimilarPeak(seedPeaksFound[seedPeakIndex], branchMz)) //if a match
                            {
                                done = false; //some seeds are still going
                                numberOfStrikesForEachPeak[seedPeakIndex] = 0; //reset
                                seedPeaksFound[seedPeakIndex].Add(branchMz);
                                allMs1PeaksFound[branchScanIndex][branchPeakIndex] = seedPeaksFound[seedPeakIndex]; //update index with the reference
                                seedMzs[seedPeakIndex] = (branchMz + seedMz) / 2; //update seed Mz so that the Mz we look for is closest to the most recent Mz (instead of the first). This helps if there's drift over time. Average prevents a single mistaken peak from throwing us off too terribly
                            }
                            else
                            {
                                if (++numberOfStrikesForEachPeak[seedPeakIndex] != NumberOfStrikesBeforeOut)
                                {
                                    done = false;
                                }
                                else //hey, this is dead to us now. Let's do some post processing. (ie average)
                                {
                                    //average the grouped mzs
                                    AggregateMzs(seedPeaksFound[seedPeakIndex]);
                                }
                            }
                        }
                    }

                    if (done) //if none of the peaks are still viable. 
                    {
                        break; //let's move on to the next seed
                    }
                }

                //Aggregate peaks that made it to the end of the file
                for (int peakIndex = 0; peakIndex < numberOfStrikesForEachPeak.Length; peakIndex++)
                {
                    if (numberOfStrikesForEachPeak[peakIndex] != NumberOfStrikesBeforeOut) //if we haven't finished this already
                    {
                        AggregateMzs(seedPeaksFound[peakIndex]);
                    }
                }


                //All grouping is done for this seed, now update the scan with the aggregated mzs
                List<double> mzsToAdd = new List<double>();
                List<double> intensitiesToAdd = new List<double>();
                MsDataScan originalScan = ms1scans[seedScanIndex];

                double[] unchangedIntensities = originalScan.MassSpectrum.YArray; //get the old intensities
                for (int peakIndex = 0; peakIndex < seedPeaksFound.Length; peakIndex++)
                {
                    List<double> mz = seedPeaksFound[peakIndex];
                    if (mz.Count == 1) // if we didn't remove the peak
                    {
                        mzsToAdd.Add(mz[0]); //[0] because there's only one value now that we've averaged them
                        intensitiesToAdd.Add(unchangedIntensities[peakIndex]);
                    }
                }

                double[] mzArray = mzsToAdd.ToArray();
                double[] intensityArray = intensitiesToAdd.ToArray();
                MzSpectrum syntheticSpectrum = new MzSpectrum(mzArray, intensityArray, false);
                ms1scans[seedScanIndex] = CloneDataScanWithUpdatedFields(originalScan, syntheticSpectrum);

                int percentProgress = (int)((1d * seedScanIndex / ms1mzs.Length) * 100);

                if (percentProgress > oldPercentProgress)
                {
                    oldPercentProgress = percentProgress;
                    ReportProgress(new ProgressEventArgs(percentProgress, "Averaging MS1 spectra... ", NestedIds));
                }
            }


            //update the datafile with the new scans
            int ms1Index = 0;
            MsDataScan currentSyntheticMS1 = ms1scans[ms1Index];
            for (int i = 0; i < allScans.Length; i++)
            {
                if (currentSyntheticMS1.OneBasedScanNumber == allScans[i].OneBasedScanNumber)
                {
                    allScans[i] = currentSyntheticMS1;
                    ms1Index++;
                    if (ms1Index == ms1scans.Length)
                    {
                        break;
                    }
                    else
                    {
                        currentSyntheticMS1 = ms1scans[ms1Index];
                    }
                }
            }

            //estimate optimal retention time tolerance
            //we previously recorded the difference between scan indexes for elutions, now use those
            //double totalRunTime = ms1scans.Last().RetentionTime;
            //double estimatedTimeBetweenMs1Scans = totalRunTime / ms1scans.Length;
            //double averageElution = elutionProfileWidthsInScans.Average() * estimatedTimeBetweenMs1Scans;
            //double innerQuartile = Statistics.InterquartileRange(elutionProfileWidthsInScans) * estimatedTimeBetweenMs1Scans;
            //currently not doing anything with this info


            Status("Getting ms2 scans...");
            Ms2ScanWithSpecificMass[] MS2Scans = GetMs2Scans(OriginalFile, OriginalFilePath, CommonParameters).ToArray();
            {/*
                #region Identifying MS2 groups

                Status("Identifying MS2 groups");
                //need to group together which scans to compare
                //index ms2scans by precursor masses
                //Could change this to the max MS1 range observed
                int maxMass = 30000;
                int binsPerDalton = 1000;
                List<int>[] massIndex = new List<int>[maxMass * binsPerDalton]; //scan number possesing a precursor mass of the index
                for (int i = 0; i < MS2Scans.Length; i++)
                {
                    int mass = (int)Math.Round(MS2Scans[i].PrecursorMonoisotopicPeakMz * binsPerDalton);
                    if (massIndex[mass] == null)
                    {
                        massIndex[mass] = new List<int> { i };
                    }
                    else
                    {
                        massIndex[mass].Add(i);
                    }
                }
                //somewhat tricky. We want to group a scan with all the other scans, but there's a chance that the tolerance falls outside for some but not for others.
                //we can just expand the range for tolerance whenever another is added.

                //want to group all scans to compare and later subgroup those where the comparison is above a threshold
                bool[] seen = new bool[MS2Scans.Length]; //don't double group stuff
                List<List<Ms2ScanWithSpecificMass>> groups = new List<List<Ms2ScanWithSpecificMass>>();
                for (int i = 0; i < MS2Scans.Length; i++)
                {
                    //get indexes to compare
                    if (!seen[i])
                    {
                        if (GlobalVariables.StopLoops)
                        {
                            return new MetaMorpheusEngineResults(this);
                        }
                        seen[i] = true; //we've seen it, so don't use it again
                        var scan = MS2Scans[i]; //get the scan
                        int obsFragmentFloorMass = (int)Math.Floor((PrecursorTolerance.GetMinimumValue(scan.PrecursorMonoisotopicPeakMz)) * binsPerDalton);
                        int obsFragmentCeilingMass = (int)Math.Ceiling((PrecursorTolerance.GetMaximumValue(scan.PrecursorMonoisotopicPeakMz)) * binsPerDalton);
                        int minBinObs = -1; //save outer bounds so we can expand tolerances if needed
                        int maxBinObs = -1; //save outer bounds so we can expand tolerances if needed
                        List<Ms2ScanWithSpecificMass> groupToAdd = new List<Ms2ScanWithSpecificMass>(); //current group
                                                                                                        //foreach mz in range, expand range if necessary
                        for (int bin = obsFragmentFloorMass; bin <= obsFragmentCeilingMass; bin++) //go through bins and add to the group
                        {
                            List<int> scans = massIndex[bin];
                            if (scans != null) //FIXME: It's still possible to group things twice since large mz ppm can hit a smaller mz where the small ppm wouldn't
                            {
                                if (minBinObs == -1)
                                {
                                    minBinObs = bin;
                                }
                                maxBinObs = bin;
                                foreach (int scanIndex in scans)
                                {
                                    seen[scanIndex] = true;
                                    groupToAdd.Add(MS2Scans[scanIndex]);
                                }
                            }
                        }

                        //get lower bound if it's expanded
                        int decreasingValues = obsFragmentFloorMass - 1;
                        int minimumValue = (int)Math.Floor(PrecursorTolerance.GetMinimumValue(minBinObs));
                        while (decreasingValues > minimumValue)
                        {
                            List<int> scans = massIndex[decreasingValues];
                            if (scans != null)
                            {
                                minimumValue = (int)Math.Floor(PrecursorTolerance.GetMinimumValue(decreasingValues));

                                foreach (int scanIndex in scans)
                                {
                                    seen[scanIndex] = true;
                                    groupToAdd.Add(MS2Scans[scanIndex]);
                                }
                            }
                            decreasingValues--;
                        }

                        //get upper bound if it's expanded
                        int increasingValues = obsFragmentCeilingMass + 1;
                        int maximumValue = (int)Math.Ceiling(PrecursorTolerance.GetMaximumValue(maxBinObs));
                        while (increasingValues < maximumValue)
                        {
                            List<int> scans = massIndex[increasingValues];
                            if (scans != null)
                            {
                                maximumValue = (int)Math.Ceiling(PrecursorTolerance.GetMaximumValue(increasingValues));

                                foreach (int scanIndex in scans)
                                {
                                    seen[scanIndex] = true;
                                    groupToAdd.Add(MS2Scans[scanIndex]);
                                }
                            }
                            increasingValues++;
                        }
                        groups.Add(groupToAdd);
                    }
                }

                //Now that we've separated groups by mass, let's try to separate based on retention time
                List<List<Ms2ScanWithSpecificMass>> retGroups = new List<List<Ms2ScanWithSpecificMass>>();
                foreach (List<Ms2ScanWithSpecificMass> group in groups) //go over all the previously made groups
                {
                    if (GlobalVariables.StopLoops)
                    {
                        return new MetaMorpheusEngineResults(this);
                    }
                    List<List<Ms2ScanWithSpecificMass>> subGroups = new List<List<Ms2ScanWithSpecificMass>>(); //local subgroups go here, needed so you don't regroup previous classifications
                    foreach (Ms2ScanWithSpecificMass scan in group) //iterate through each scan in the previous group
                    {
                        bool foundSpot = false;
                        foreach (List<Ms2ScanWithSpecificMass> subGroup in subGroups) //see if the scan fits somewhere, if not make a new subgroup
                        {
                            Ms2ScanWithSpecificMass scanInSubGroup = subGroup.Last(); //only need to check the lawst
                            {
                                if (scan.RetentionTime > scanInSubGroup.RetentionTime - MaxRetentionTimeDifferenceAllowedInMinutes
                                    && scan.RetentionTime < scanInSubGroup.RetentionTime + MaxRetentionTimeDifferenceAllowedInMinutes) //if a match, add it
                                {
                                    subGroup.Add(scan);
                                    foundSpot = true;
                                    break;
                                }
                            }
                            if (foundSpot)
                            {
                                break;
                            }
                        }
                        if (!foundSpot) //if not a match, create a new category
                        {
                            subGroups.Add(new List<Ms2ScanWithSpecificMass> { scan });
                        }
                    }
                    retGroups.AddRange(subGroups);
                }
                groups = retGroups; //save reassignments

                //Now let's separate based on cosine score
                List<List<Ms2ScanWithSpecificMass>> scoredGroups = new List<List<Ms2ScanWithSpecificMass>>();
                foreach (List<Ms2ScanWithSpecificMass> group in groups) //go over all the previously made groups
                {
                    if (GlobalVariables.StopLoops)
                    {
                        return new MetaMorpheusEngineResults(this);
                    }
                    List<List<Ms2ScanWithSpecificMass>> subGroups = new List<List<Ms2ScanWithSpecificMass>>(); //local subgroups go here, needed so you don't regroup previous classifications
                    foreach (Ms2ScanWithSpecificMass scan in group) //iterate through each scan in the previous group
                    {
                        bool foundSpot = false;
                        foreach (List<Ms2ScanWithSpecificMass> subGroup in subGroups) //see if the scan fits somewhere, if not make a new subgroup
                        {
                            foreach (Ms2ScanWithSpecificMass scanInSubGroup in subGroup) //iterate through each member of each previously found group
                            {
                                if (MinCosineScoreAllowed <= CosineScore(scan.TheScan.MassSpectrum, scanInSubGroup.TheScan.MassSpectrum, ProductTolerance)) //if a match, add it
                                {
                                    subGroup.Add(scan);
                                    foundSpot = true;
                                    break;
                                }
                            }
                            if (foundSpot)
                            {
                                break;
                            }
                        }
                        if (!foundSpot) //if not a match, create a new category
                        {
                            subGroups.Add(new List<Ms2ScanWithSpecificMass> { scan });
                        }
                    }
                    scoredGroups.AddRange(subGroups);
                }
                groups = scoredGroups; //save

                //order by the middle scan number
                groups = groups.OrderBy(x => x[x.Count / 2].OneBasedScanNumber).ToList();

                #endregion Identifying MS2 groups

                #region Averaging MS2 spectra

                Status("Averaging MS2 spectra");
                //Each MS2 spectra is going to be assigned a new scan number that's placed at the middle of the group. 
                //There's the possibility that more spectra are generated because of chimeras, or fewer spectra because of aggregation (regardless, fewer comparisons afterward)
                int assignedScanNumber = 0; //this gets ++ immediately, so there shouldn't be any 0 scan numbers
                int ms1Index = 0;
                int mostRecentPrecursorNumber = 0;

                List<MsDataScan> syntheticSpectra = new List<MsDataScan>(); //include MS1 as we go through this
                foreach (List<Ms2ScanWithSpecificMass> group in groups)
                {
                    if (GlobalVariables.StopLoops)
                    {
                        return new MetaMorpheusEngineResults(this);
                    }
                    assignedScanNumber++;
                    Ms2ScanWithSpecificMass representativeScan = group[group.Count / 2];

                    while (ms1Index < ms1scans.Length
                        && representativeScan.OneBasedScanNumber > ms1scans[ms1Index].OneBasedScanNumber) //if there should be a precursor at this point
                    {
                        mostRecentPrecursorNumber = assignedScanNumber; //track so that ms2 precursor scan numbers can be updated
                        syntheticSpectra.Add(CloneDataScanWithUpdatedFields(ms1scans[ms1Index], scanNumber: assignedScanNumber)); //update scan number
                        ms1Index++; //update to the next precursor
                        assignedScanNumber++; //update the current scan
                    }

                    if (group.Count == 1) //if nothing to aggregate, just save it after updating scan number and precursor info
                    {
                        syntheticSpectra.Add(CloneDataScanWithUpdatedFields(
                            representativeScan.TheScan,
                            scanNumber: assignedScanNumber, //update scan number
                            precursorScanNumber: mostRecentPrecursorNumber,
                            precursorMonoisotopicMz: representativeScan.PrecursorMonoisotopicPeakMz,
                            precursorCharge: representativeScan.PrecursorCharge,
                            precursorMass: representativeScan.PrecursorMass
                            ));
                    }
                    else
                    {
                        List<double> syntheticMZs = new List<double>();
                        List<double> syntheticIntensities = new List<double>();

                        //get values
                        double[][] mzsForGroupedScans = group.Select(x => x.TheScan.MassSpectrum.XArray).ToArray();
                        double[][] intensitiesForGroupedScans = group.Select(x => x.TheScan.MassSpectrum.YArray).ToArray();

                        List<double> allMzs = new List<double>();
                        List<double> groupedMzMaxRanges = new List<double>();

                        //Sort array of all mzs sorted low to high
                        allMzs.Sort();

                        //group mzs
                        int seedIndex = 0;
                        for (int i = 0; i < allMzs.Count; i++)
                        {
                            seedIndex = i; //the actual seed
                            double maxValue = ProductTolerance.GetMaximumValue(allMzs[i]);
                            i++;
                            for (; i < allMzs.Count; i++)
                            {
                                double currentMz = allMzs[i];
                                if (currentMz > maxValue) //stop if out of range
                                {
                                    break;
                                }
                                else //update range
                                {
                                    maxValue = ProductTolerance.GetMaximumValue(currentMz);
                                }
                            }
                            //record
                            groupedMzMaxRanges.Add(maxValue);
                        }

                        //we're done grouping peaks, see if there's enough info (half of all grouped Ms2s must contain this peak group)
                        int[] peakPositionArray = new int[group.Count]; //save position so we don't have to iterate through each time

                        foreach (double maxValue in groupedMzMaxRanges) //foreach peak group
                        {
                            int numScansNeeded = group.Count / 2;
                            List<(double mz, double intensity)> groupedPeaks = new List<(double mz, double intensity)>();
                            for (int i = 0; i < group.Count; i++) //foreach scan we're looking at here
                            {
                                double[] currentScanMzs = mzsForGroupedScans[i]; //get mzs
                                double[] currentScanIntensities = intensitiesForGroupedScans[i]; //get intensities
                                int j = peakPositionArray[i]; //previous position we left off at
                                for (; j < currentScanMzs.Length; j++)
                                {
                                    if (currentScanMzs[j] < maxValue) //it's in range, so save it!
                                    {
                                        groupedPeaks.Add((currentScanMzs[j], currentScanIntensities[j]));
                                    }
                                    else //not a match
                                    {
                                        break;
                                    }
                                }
                                if (peakPositionArray[i] != j) //check if we found anything
                                {
                                    peakPositionArray[i] = j; //update position
                                    numScansNeeded--; //we had a hit, so we need one fewer scan now!
                                }
                            }
                            if (numScansNeeded <= 0) //if we saw this peak group enough times to believe it
                            {
                                (double mz, double intensity) syntheticPeak = AverageMzsAndIntensities(groupedPeaks);
                                syntheticMZs.Add(syntheticPeak.mz);
                                syntheticIntensities.Add(syntheticPeak.intensity);
                            }
                        }

                        MzSpectrum syntheticSpectrum = new MzSpectrum(syntheticMZs.ToArray(), syntheticIntensities.ToArray(), false);
                        syntheticSpectra.Add(CloneDataScanWithUpdatedFields(
                            representativeScan.TheScan,
                            syntheticSpectrum,
                            scanNumber: assignedScanNumber,
                            precursorScanNumber: mostRecentPrecursorNumber,
                            precursorMonoisotopicMz: representativeScan.PrecursorMonoisotopicPeakMz,
                            precursorCharge: representativeScan.PrecursorCharge,
                            precursorMass: representativeScan.PrecursorMass
                            ));
                    }
                }
                //wrap up any precursor scans that didn't make it in
                for (; ms1Index < ms1scans.Length; ms1Index++)
                {
                    assignedScanNumber++; //update the current scan
                    syntheticSpectra.Add(CloneDataScanWithUpdatedFields(ms1scans[ms1Index], scanNumber: assignedScanNumber)); //update scan number
                }

                #endregion Averaging MS2 spectra            
                */
            }
            ////////////// TEST WRITE
            List<MsDataScan> syntheticSpectra = new List<MsDataScan>(); //include MS1 as we go through this
            {
                int assignedScanNumber = 0; //this gets ++ immediately, so there shouldn't be any 0 scan numbers
                ms1Index = 0;
                int mostRecentPrecursorNumber = 0;
                //Print diagnostics
                List<string> diagnosticLines = new List<string>();
                foreach (Ms2ScanWithSpecificMass ms2 in MS2Scans)
                {
                    assignedScanNumber++;
                    while (ms1Index < ms1scans.Length
                        && ms2.OneBasedScanNumber > ms1scans[ms1Index].OneBasedScanNumber) //if there should be a precursor at this point
                    {
                        mostRecentPrecursorNumber = assignedScanNumber; //track so that ms2 precursor scan numbers can be updated
                        syntheticSpectra.Add(CloneDataScanWithUpdatedFields(ms1scans[ms1Index], scanNumber: assignedScanNumber)); //update scan number
                        diagnosticLines.Add(assignedScanNumber.ToString() +
                            '\t' + "0" +//ms2.PrecursorMonoisotopicPeakMz.ToString() +
                            '\t' + "0" +//ms2.PrecursorCharge.ToString() +
                            '\t' + ms1scans[ms1Index].OneBasedScanNumber);
                        ms1Index++; //update to the next precursor
                        assignedScanNumber++; //update the current scan
                    }
                    syntheticSpectra.Add(CloneDataScanWithUpdatedFields(
        ms2.TheScan,
        scanNumber: assignedScanNumber,
        precursorScanNumber: mostRecentPrecursorNumber,
        precursorMonoisotopicMz: ms2.PrecursorMonoisotopicPeakMz,
        precursorCharge: ms2.PrecursorCharge,
        precursorMass: ms2.PrecursorMass
        ));
                    diagnosticLines.Add(assignedScanNumber.ToString() +
        '\t' + ms2.PrecursorMonoisotopicPeakMz.ToString() +
        '\t' + ms2.PrecursorCharge.ToString() +
        '\t' + ms2.OneBasedScanNumber);
                }
                for (; ms1Index < ms1scans.Length; ms1Index++)
                {
                    assignedScanNumber++; //update the current scan
                    syntheticSpectra.Add(CloneDataScanWithUpdatedFields(ms1scans[ms1Index], scanNumber: assignedScanNumber)); //update scan number
                    diagnosticLines.Add(assignedScanNumber.ToString() +
        '\t' + "0" +//ms2.PrecursorMonoisotopicPeakMz.ToString() +
        '\t' + "0" +//ms2.PrecursorCharge.ToString() +
        '\t' + ms1scans[ms1Index].OneBasedScanNumber);
                }
            }
            ////////////////////
            AggregatedDataFile = new MsDataFile(syntheticSpectra.ToArray(), OriginalFile.SourceFile);
            return new MetaMorpheusEngineResults(this);
        }


        public double CosineScore(MzSpectrum scan1, MzSpectrum scan2, Tolerance tolerance)
        {
            double[] mz1 = scan1.XArray;
            double[] intensity1 = scan1.YArray;
            double[] mz2 = scan2.XArray;
            double[] intensity2 = scan2.YArray;
            return CosineScore(mz1, mz2, intensity1, intensity2, tolerance);
        }

        public double CosineScore(MsDataScan scan1, MsDataScan scan2, Tolerance tolerance)
        {
            double[] mz1 = scan1.MassSpectrum.XArray;
            double[] intensity1 = scan1.MassSpectrum.YArray;
            double[] mz2 = scan2.MassSpectrum.XArray;
            double[] intensity2 = scan2.MassSpectrum.YArray;
            return CosineScore(mz1, mz2, intensity1, intensity2, tolerance);
        }

        public double CosineScore(double[] mz1, double[] mz2, double[] intensity1, double[] intensity2, Tolerance tolerance)
        {
            //convert spectra to vectors
            List<double> vector1 = new List<double>();
            List<double> vector2 = new List<double>();
            int i = 0;
            int j = 0;

            while (i != mz1.Length && j != mz2.Length)
            {
                double one = mz1[i];
                double two = mz2[j];
                if (tolerance.Within(one, two))
                {
                    vector1.Add(intensity1[i]);
                    vector2.Add(intensity2[j]);
                    i++;
                    j++;
                }
                else if (one > two)
                {
                    vector1.Add(0);
                    vector2.Add(intensity2[j]);
                    j++;
                }
                else //two>one
                {
                    vector1.Add(intensity1[i]);
                    vector2.Add(0);
                    i++;
                }
            }
            //wrap up leftover peaks
            for (; i < mz1.Length; i++)
            {
                vector1.Add(intensity1[i]);
                vector2.Add(0);
            }
            for (; j < mz2.Length; j++)
            {
                vector1.Add(0);
                vector2.Add(intensity2[j]);
            }

            //cosine score of vectors
            //numerator
            double numerator = 0;
            for (i = 0; i < vector1.Count; i++)
            {
                numerator += vector1[i] * vector2[i];
            }

            //denominator
            double denominator = Math.Sqrt(vector1.Sum(x => x * x)) * Math.Sqrt(vector2.Sum(x => x * x));

            //calculate cosine score
            return Math.Round(numerator / denominator * 1000) / 1000;
        }

        public bool SimilarPeak(List<double> currentMzs, double putativeMz)
        {
            //determine if the putativeMz is the same as the previous ones.
            if (currentMzs.Count < 3)
            {
                return true;
            }
            else
            {
                //Compute the Average      
                double avg = currentMzs.Average();
                //Perform the Sum of (value-avg)_2_2      
                double sumOfSquares = currentMzs.Sum(d => Math.Pow(d - avg, 2));
                //Put it all together      
                double stdev = Math.Sqrt((sumOfSquares) / (currentMzs.Count - 1));
                return (Math.Abs(avg - putativeMz) < stdev * 1.96);
            }
        }

        public void AggregateMzs(List<double> referenceListOfMzs) //use for MS1 scans
        {
            //Currently NOT using intensity for weighting. Oddly, intensity weighting yielded poorer results than the median
            if (referenceListOfMzs.Count != 1) //if it's worth aggregating
            {
                referenceListOfMzs.Sort();
                double aggregateMz = referenceListOfMzs[referenceListOfMzs.Count / 2]; //get the median
                referenceListOfMzs.Clear(); //need to clear to keep the reference
                referenceListOfMzs.Add(aggregateMz); //add the median as the only peak
            }
            else //this peak was only found once... isn't that a little odd? like maybe it's noise?
            {
                referenceListOfMzs.Clear(); //Let's wipe it, but keep the entry to know that it's blank now.
            }
        }

        public (double mz, double intensity) AverageMzsAndIntensities(List<(double mz, double intensity)> peaks) //use for MS2 scans
        {
            double normalizedMzSum = peaks.Sum(x => x.mz * x.intensity); //weight by intensity
            double intensitySum = peaks.Sum(x => x.intensity);
            double averageMz = normalizedMzSum / intensitySum;
            double averageIntensity = intensitySum / peaks.Count;
            return (averageMz, averageIntensity);
        }

        public MsDataScan CloneDataScanWithUpdatedFields(MsDataScan oldScan,
            MzSpectrum updatedSpectrum = null,
            int scanNumber = -1,
            int? precursorScanNumber = -1,
            double? precursorMonoisotopicMz = -1,
            int? precursorCharge = -1,
            double? precursorMass = -1
            )
        {
            if (scanNumber == -1)
            {
                scanNumber = oldScan.OneBasedScanNumber;
            }
            if (precursorScanNumber == -1)
            {
                precursorScanNumber = oldScan.OneBasedPrecursorScanNumber;
            }
            if (updatedSpectrum == null)
            {
                updatedSpectrum = oldScan.MassSpectrum;
            }
            if (precursorMonoisotopicMz.Value == -1)
            {
                precursorMonoisotopicMz = oldScan.SelectedIonMonoisotopicGuessMz;
            }
            if (precursorCharge.Value == -1)
            {
                precursorCharge = oldScan.SelectedIonChargeStateGuess;
            }
            double? selectedIonMz = precursorMass == -1 ? oldScan.SelectedIonMZ : precursorMass.Value.ToMz(precursorCharge.Value);

            var scan = new MsDataScan(
                updatedSpectrum, //changed?
                scanNumber, //changed?
                oldScan.MsnOrder,
                oldScan.IsCentroid,
                oldScan.Polarity,
                oldScan.RetentionTime,
                oldScan.ScanWindowRange,
                oldScan.ScanFilter,
                oldScan.MzAnalyzer,
                oldScan.TotalIonCurrent,
                oldScan.InjectionTime,
                oldScan.NoiseData,
                oldScan.NativeId,
                selectedIonMz, //changed?
                precursorCharge, //changed?
                oldScan.SelectedIonIntensity,
                selectedIonMz, //changed?
                oldScan.IsolationWidth,
                oldScan.DissociationType,
                precursorScanNumber, //changed?
                precursorMonoisotopicMz //changed?
                );

            return scan;
        }
    }
}