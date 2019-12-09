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
        private readonly Tolerance ProductTolerance;
        private readonly MsDataFile OriginalFile;
        private const int NumberOfStrikesBeforeOut = 2;

        public MsDataFile AggregatedDataFile { get; private set; }

        public AggregationEngine(MsDataFile originalFile, CommonParameters commonParameters, List<string> nestedIds) : base(commonParameters, nestedIds)
        {
            OriginalFile = originalFile;            
            ProductTolerance = commonParameters.ProductMassTolerance;
        }

        protected override MetaMorpheusEngineResults RunSpecific()
        {
            //MS1: this engine matches adjacent MS1 peaks that fall within a specified tolerance
            //If an MS1 peak does not have an adjacent match, then it is removed           
            //average the m/zs found across scans and update each value with the average

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

            AggregatedDataFile = new MsDataFile(allScans, OriginalFile.SourceFile);
            return new MetaMorpheusEngineResults(this);
        }

        public static bool SimilarPeak(List<double> currentMzs, double putativeMz)
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

        public static void AggregateMzs(List<double> referenceListOfMzs) //use for MS1 scans
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

        public static MsDataScan CloneDataScanWithUpdatedFields(MsDataScan oldScan,
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
                precursorMonoisotopicMz, //changed?
                oldScan.HcdEnergy
                );

            return scan;
        }
    }
}