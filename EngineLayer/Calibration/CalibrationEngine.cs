using MassSpectrometry;
using System;
using System.Collections.Generic;
using System.Linq;

namespace EngineLayer.Calibration
{
    public class CalibrationEngine : MetaMorpheusEngine
    {
        private const int NumberOfScansUsedForSmoothingOnEachSide = 100;
        private readonly MsDataFile MyMsDataFile;
        private readonly DataPointAquisitionResults Datapoints;
        public MsDataFile CalibratedDataFile { get; private set; }

        public CalibrationEngine(MsDataFile myMSDataFile, DataPointAquisitionResults datapoints, CommonParameters commonParameters, List<string> nestedIds) : base(commonParameters, nestedIds)
        {
            MyMsDataFile = myMSDataFile;
            Datapoints = datapoints;
        }

        protected override MetaMorpheusEngineResults RunSpecific()
        {
            Status("Calibrating spectra");
            List<LabeledDataPoint> ms1Points = Datapoints.Ms1List;
            List<LabeledDataPoint> ms2Points = Datapoints.Ms2List;
            List<MsDataScan> originalScans = MyMsDataFile.GetAllScansList();

            List<MsDataScan> ms1Scans = new List<MsDataScan>();
            List<MsDataScan> ms2Scans = new List<MsDataScan>();
            //separate scans by msnOrder, because the mass accuracy varies between the two
            foreach (MsDataScan scan in originalScans)
            {
                if (scan.MsnOrder == 1)
                {
                    ms1Scans.Add(scan);
                }
                else
                {
                    ms2Scans.Add(scan);
                }
            }

            //create a way to go from the scan number to the scan order. This can be a single array, since there shouldn't be any overlap of scan numbers between ms1 and ms2 scans
            int[] scanNumberToScanPlacement = new int[originalScans.Max(x => x.OneBasedScanNumber) + 1];
            for (int i = 0; i < ms1Scans.Count; i++)
            {
                scanNumberToScanPlacement[ms1Scans[i].OneBasedScanNumber] = i;
            }
            for (int i = 0; i < ms2Scans.Count; i++)
            {
                scanNumberToScanPlacement[ms2Scans[i].OneBasedScanNumber] = i;
            }

            //Populate the weighted average relative error for each scan, where index of returned array is the placement
            double[] ms1RelativeErrors = PopulateErrors(ms1Points, scanNumberToScanPlacement, ms1Scans.Count);
            double[] ms2RelativeErrors = new double[ms2Scans.Count];
            if(!(commonParameters.DissociationType == DissociationType.LowCID))
            {
                ms2RelativeErrors = PopulateErrors(ms2Points, scanNumberToScanPlacement, ms2Scans.Count);
            }
                
            //generate new scans
            MsDataScan[] calibratedScans = new MsDataScan[originalScans.Count];

            //hard copy original scans
            for (int i = 0; i < originalScans.Count; i++)
            {
                calibratedScans[i] = originalScans[i];
            }

            //apply a smoothing function, so that outlier scans aren't wildly shifted
            double[] ms1SmoothedErrors = SmoothErrors(ms1RelativeErrors);
            double[] ms2SmoothedErrors = new double[ms2RelativeErrors.Length];
            if (!(commonParameters.DissociationType == DissociationType.LowCID))
            {
                ms2SmoothedErrors = SmoothErrors(ms2RelativeErrors);
            }
            

            //calibrate the data
            int ms1Index = 0;
            int ms2Index = 0;
            double mostRecentMS1SmoothedError = ms1SmoothedErrors.FirstOrDefault(); //this is needed to update the precursor mass error of MS2 scans
            for (int scanIndex = 0; scanIndex < calibratedScans.Length; scanIndex++) //foreach scan
            {
                MsDataScan originalScan = originalScans[scanIndex]; //get original scan
                if (originalScan.MsnOrder == 1) //if ms1
                {
                    mostRecentMS1SmoothedError = ms1SmoothedErrors[ms1Index]; //update the mass error
                    calibratedScans[scanIndex] = CalibrateScan(originalScan, mostRecentMS1SmoothedError);
                    ms1Index++;
                }
                else if (!(commonParameters.DissociationType == DissociationType.LowCID))
                {
                    calibratedScans[scanIndex] = CalibrateScan(originalScan, ms2SmoothedErrors[ms2Index], mostRecentMS1SmoothedError);
                    ms2Index++;
                }
                else
                {
                    //no change for low res
                    calibratedScans[scanIndex] = originalScan;
                }
            }

            CalibratedDataFile = new MsDataFile(calibratedScans, MyMsDataFile.SourceFile);
            return new MetaMorpheusEngineResults(this);
        }

        private static double[] PopulateErrors(List<LabeledDataPoint> datapoints, int[] scanNumberToScanPlacement, int arrayLength)
        {
            //return an array of weighted average relative errors for each scan
            double[] averageRelativeErrors = new double[arrayLength];

            int currentScanNumber = datapoints.First().ScanNumber; //get the first scan number. This should be an ordered list
            List<(double massError, double logIntensity)> localRelativeErrors = new List<(double massError, double logIntensity)>(); //create a list to store the mass error

            foreach (LabeledDataPoint datapoint in datapoints)
            {
                if (datapoint.ScanNumber == currentScanNumber)
                {
                    localRelativeErrors.Add((datapoint.RelativeMzError, datapoint.LogIntensity));
                }
                else
                {
                    //wrap up the old set
                    averageRelativeErrors[scanNumberToScanPlacement[currentScanNumber]] = CalculateAverageRelativeErrors(localRelativeErrors);

                    //update
                    currentScanNumber = datapoint.ScanNumber;
                    localRelativeErrors = new List<(double massError, double logIntensity)> { (datapoint.RelativeMzError, datapoint.LogIntensity) };
                }
            }
            //finish loop
            averageRelativeErrors[scanNumberToScanPlacement[currentScanNumber]] = CalculateAverageRelativeErrors(localRelativeErrors);

            return averageRelativeErrors;
        }

        private static double CalculateAverageRelativeErrors(List<(double massError, double logIntensity)> localRelativeErrors)
        {
            //double logIntensityToSubtract = localRelativeErrors.Min(x => x.logIntensity) - 1; // normalize each log intensity so that the minimum log intensity is 1. 
            //Convert from log to actual intensity to more heavily weight intensities.
            double weightedSumOfErrors = localRelativeErrors.Sum(x => x.massError * Math.Pow(10, x.logIntensity));
            double sumOfIntensities = localRelativeErrors.Sum(x => Math.Pow(10, x.logIntensity));
            return weightedSumOfErrors / sumOfIntensities;
        }

        private static double[] SmoothErrors(double[] relativeErrors)
        {
            //impute missing values
            //not all scans are guarenteed to contain data points. We can infer these data point with nearby points.
            for (int index = 0; index < relativeErrors.Length; index++)
            {
                if (relativeErrors[index] == 0) //if there were no points, then the value should be perfectly zero (the double default)
                {
                    int startingBlankIndex = index;
                    //increase the index until we find the next scan containing a data point
                    while (index < relativeErrors.Length && relativeErrors[index] == 0)
                    {
                        index++;
                    }
                    double nextError = index == relativeErrors.Length ? relativeErrors[startingBlankIndex - 1] : relativeErrors[index]; //can't go all the way through without any data points, the original function checks for enough data points (where enough is more than zero)
                    double previousError = startingBlankIndex > 0 ? relativeErrors[startingBlankIndex - 1] : nextError;
                    int numberOfConsecutiveScansWithoutDataPoints = index - startingBlankIndex;
                    for (int tempIndex = 0; tempIndex < numberOfConsecutiveScansWithoutDataPoints; tempIndex++)
                    {
                        relativeErrors[startingBlankIndex + tempIndex] = ((tempIndex + 1) * nextError + (numberOfConsecutiveScansWithoutDataPoints - tempIndex - 1) * previousError) / numberOfConsecutiveScansWithoutDataPoints;
                    }
                }
            }

            //get correction factor for each scan, where half of the smooth comes from both directions
            double[] smoothedErrors = new double[relativeErrors.Length];
            int leftIndex = 0; //starting scan index used for numbers less than the current scan
            int rightIndex = 1; //
            double smoothedCorrectionFactor = 0; //this variable is the sum of all nearby errors, to be later divided by the number of summed errors for a smoothed average error
            //for scan #1 (index 0)
            //no left scans, because we're at the beginning of the file. Just populate the first "numberOfScansUsedForSmoothingOnEachSide" scan errors on the right side

            while (rightIndex < NumberOfScansUsedForSmoothingOnEachSide && rightIndex < relativeErrors.Length)
            {
                smoothedCorrectionFactor += relativeErrors[rightIndex];
                rightIndex++;
            }

            //for each scan
            for (int i = 0; i < relativeErrors.Length; i++)
            {
                //for left index, remove an error if 
                if (i > NumberOfScansUsedForSmoothingOnEachSide)
                {
                    smoothedCorrectionFactor -= relativeErrors[leftIndex];
                    leftIndex++;
                }

                if (rightIndex < relativeErrors.Length) //need to check so we don't run off the end of the file
                {
                    smoothedCorrectionFactor += relativeErrors[rightIndex];
                    rightIndex++;
                }

                smoothedErrors[i] = smoothedCorrectionFactor / (rightIndex - leftIndex);
            }
            return smoothedErrors;
        }

        private static MsDataScan CalibrateScan(MsDataScan oldScan, double smoothedRelativeError, double? precursorSmoothedRelativeError = null)
        {
            double correctionFactor = 1 - smoothedRelativeError; //create the multiplier. Positive mass errors mean that the experimental mass was greater than the theoretical, so we want to shift the experimental DOWN
            double[] originalMzs = oldScan.MassSpectrum.XArray;
            double[] calibratedMzs = new double[originalMzs.Length];
            //calibrate the mzs. Because peaks are in mz and we are making a mz shift, we don't need to deconvolute
            for (int i = 0; i < originalMzs.Length; i++)
            {
                calibratedMzs[i] = originalMzs[i] * correctionFactor;
            }

            //update precursor values (if applicable)
            double? selectedIonMz = null;
            double? isolationMz = null;
            double? selectedIonMonoisotopicGuessMz = null;
            if (precursorSmoothedRelativeError != null)
            {
                correctionFactor = 1 - precursorSmoothedRelativeError.Value;
                if (oldScan.SelectedIonMZ.HasValue)
                {
                    selectedIonMz = oldScan.SelectedIonMZ * correctionFactor;
                }
                if (oldScan.IsolationMz.HasValue)
                {
                    isolationMz = oldScan.IsolationMz * correctionFactor;
                }
                if (oldScan.SelectedIonMZ.HasValue)
                {
                    selectedIonMonoisotopicGuessMz = oldScan.SelectedIonMonoisotopicGuessMz * correctionFactor;
                }
            }

            //create new calibrated spectrum
            MzSpectrum calibratedSpectrum = new MzSpectrum(calibratedMzs, oldScan.MassSpectrum.YArray, false);
            return new MsDataScan(
                calibratedSpectrum, //changed
                oldScan.OneBasedScanNumber,
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
                oldScan.SelectedIonChargeStateGuess,
                oldScan.SelectedIonIntensity,
                isolationMz, //changed?
                oldScan.IsolationWidth,
                oldScan.DissociationType,
                oldScan.OneBasedPrecursorScanNumber,
                selectedIonMonoisotopicGuessMz //changed?
                );
        }
    }
}