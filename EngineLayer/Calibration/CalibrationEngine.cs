using MassSpectrometry;
using System;
using System.Collections.Generic;
using System.Linq;

namespace EngineLayer.Calibration
{
    public class CalibrationEngine : MetaMorpheusEngine
    {
        private const int ScansUsedForSmoothingOnEachSide = 50;
        private readonly MsDataFile MyMsDataFile;
        private readonly DataPointAquisitionResults Datapoints;
        public MsDataFile CalibratedDataFile { get; private set; }

        public CalibrationEngine(MsDataFile myMSDataFile, DataPointAquisitionResults datapoints, CommonParameters commonParameters, List<string> nestedIds) : base(commonParameters, nestedIds)
        {
            MyMsDataFile = myMSDataFile;
            Datapoints = datapoints;
            int numMs1Scans = MyMsDataFile.GetAllScansList().Where(x => x.MsnOrder == 1).Count();
        }

        protected override MetaMorpheusEngineResults RunSpecific()
        {
            Status("Calibrating spectra");
            List<LabeledDataPoint> ms1Points = Datapoints.Ms1List;
            List<LabeledDataPoint> ms2Points = Datapoints.Ms2List;
            List<MsDataScan> originalScans = MyMsDataFile.GetAllScansList();

            List<MsDataScan> ms1Scans = new List<MsDataScan>();
            List<MsDataScan> ms2Scans = new List<MsDataScan>();
            //separate scans by msnOrder
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

            //index for scanNumber to scan placement and vice versa
            int[] ms1PlacementToScanNumber = ms1Scans.Select(x => x.OneBasedScanNumber).ToArray();
            int[] ms2PlacementToScanNumber = ms2Scans.Select(x => x.OneBasedScanNumber).ToArray();

            int[] scanNumberToScanPlacement = new int[originalScans.Max(x => x.OneBasedScanNumber) + 1];
            for (int i = 0; i < ms1PlacementToScanNumber.Length; i++)
            {
                scanNumberToScanPlacement[ms1PlacementToScanNumber[i]] = i;
            }
            for (int i = 0; i < ms2PlacementToScanNumber.Length; i++)
            {
                scanNumberToScanPlacement[ms2PlacementToScanNumber[i]] = i;
            }

            //Populate the average relative error for each scan, where index of returned array is the placement
            double[] ms1RelativeErrors = PopulateErrors(ms1Points, scanNumberToScanPlacement, ms1PlacementToScanNumber.Length);
            double[] ms2RelativeErrors = PopulateErrors(ms2Points, scanNumberToScanPlacement, ms2PlacementToScanNumber.Length);

            //current scans
            MsDataScan[] calibratedScans = new MsDataScan[originalScans.Count];

            //saved scans
            MsDataScan[] ms1Saves = new MsDataScan[originalScans.Count];
            MsDataScan[] ms2Saves = new MsDataScan[originalScans.Count];

            //hard copy original scans
            for (int i=0; i<originalScans.Count; i++)
            {
                MsDataScan originalScan = originalScans[i];
                calibratedScans[i] = originalScans[i];
                ms1Saves[i] = originalScan;
                ms2Saves[i] = originalScan;
            }

            //record the best performance
            double bestMS1Performance = ms1RelativeErrors.Average(x => Math.Abs(x));
            double bestMS2Performance = ms2RelativeErrors.Average(x => Math.Abs(x));
            int[] bestParams = new int[4];



            for (int i = 0; i < 5; i++)
            {
                double[] ms1SmoothedErrors = SmoothErrors(ms1RelativeErrors, ScansUsedForSmoothingOnEachSide);
                double[] ms2SmoothedErrors = SmoothErrors(ms2RelativeErrors, ScansUsedForSmoothingOnEachSide);
                int previousScanNumber = ms1Points.FirstOrDefault().ScanNumber;
                double previousSmoothedError = ms1SmoothedErrors[scanNumberToScanPlacement[previousScanNumber]];

                int ms1Index = 0;
                int ms2Index = 0;
                double mostRecentMS1SmoothedError = ms1SmoothedErrors.FirstOrDefault();
                for (int scanIndex = 0; scanIndex < calibratedScans.Length; scanIndex++)
                {
                    MsDataScan originalScan = calibratedScans[scanIndex];
                    if (originalScan.MsnOrder == 1)
                    {
                        mostRecentMS1SmoothedError = ms1SmoothedErrors[ms1Index];
                        ms1RelativeErrors[ms1Index] -= mostRecentMS1SmoothedError;
                        ms1Index++;
                        calibratedScans[scanIndex] = CalibrateScan(originalScan, mostRecentMS1SmoothedError);
                    }
                    else //if ms2
                    {
                        double smoothedRelativeError = ms2SmoothedErrors[ms2Index];
                        ms2RelativeErrors[ms2Index] -= smoothedRelativeError;
                        ms2Index++;
                        calibratedScans[scanIndex] = CalibrateScan(originalScan, smoothedRelativeError, mostRecentMS1SmoothedError);
                    }
                }
                double currentMS1Performance = ms1RelativeErrors.Average(x => Math.Abs(x));
                double currentMS2Performance = ms2RelativeErrors.Average(x => Math.Abs(x));
                if (currentMS1Performance < bestMS1Performance) //low is better!
                {
                    bestMS1Performance = currentMS1Performance;
                    for (int k = 0; k < calibratedScans.Length; k++)
                    {
                        ms1Saves[k] = calibratedScans[k];
                    }
                }
                if (currentMS2Performance < bestMS2Performance) //low is better!
                {
                    bestMS2Performance = currentMS2Performance;
                    for (int k = 0; k < calibratedScans.Length; k++)
                    {
                        ms2Saves[k] = calibratedScans[k];
                    }
                }
            }
            

            //get best of ms1 and ms2 saves
            for (int scanIndex = 0; scanIndex < ms2Saves.Length; scanIndex++)
            {
                MsDataScan genericScan = ms2Saves[scanIndex];
                if (genericScan.MsnOrder == 1)
                {
                    calibratedScans[scanIndex] = ms1Saves[scanIndex];
                }
                else //if ms2
                {
                    calibratedScans[scanIndex] = genericScan;
                }
            }
            CalibratedDataFile = new MsDataFile(calibratedScans, MyMsDataFile.SourceFile);
            return new MetaMorpheusEngineResults(this);
        }

        private double[] PopulateErrors(List<LabeledDataPoint> datapoints, int[] scanNumberToScanPlacement, int arrayLength)
        {
            double[] averageRelativeErrors = new double[arrayLength];
            int currentScanNumber = datapoints.First().ScanNumber;
            List<double> localRelativeErrors = new List<double>();
            foreach (LabeledDataPoint datapoint in datapoints)
            {
                if (datapoint.ScanNumber == currentScanNumber)
                {
                    localRelativeErrors.Add(datapoint.RelativeMzError);
                }
                else
                {
                    //wrap up
                    averageRelativeErrors[scanNumberToScanPlacement[currentScanNumber]] = localRelativeErrors.Average();
                    //update
                    currentScanNumber = datapoint.ScanNumber;
                    localRelativeErrors = new List<double> { datapoint.RelativeMzError };
                }
            }
            //finish loop
            averageRelativeErrors[scanNumberToScanPlacement[currentScanNumber]] = localRelativeErrors.Average();
            return averageRelativeErrors;
        }

        private double[] SmoothErrors(double[] relativeErrors, int numberOfScansUsedForSmoothingOnEachSide)
        {
            //get correction factor for each scan, where half of the smooth comes from both directions
            double[] smoothedErrors = new double[relativeErrors.Length];
            int leftIndex = 0;
            int rightIndex = 1;
            int emptyIndexesCount = 0;
            double smoothedCorrectionFactor = 0;
            int leftValuesUsed = 0;
            int rightValuesUsed = 0;
            //for scan #1 (index 0)
            //no left scans, because we're at the beginning of the file
            for (; rightIndex < (numberOfScansUsedForSmoothingOnEachSide + emptyIndexesCount) && rightIndex < relativeErrors.Length; rightIndex++)
            {
                double addedFactor = relativeErrors[rightIndex];
                if (addedFactor == 0)
                {
                    emptyIndexesCount++;
                }
                else
                {
                    smoothedCorrectionFactor += addedFactor;
                    rightValuesUsed++;
                }
            }

            //for all other scans
            for (int i = 0; i < relativeErrors.Length; i++)
            {
                //make previous right value a left value. (the current value is classified as a left value)
                double currentError = relativeErrors[i];
                if (currentError != 0)
                {
                    leftValuesUsed++;
                    rightValuesUsed--;
                }
                else
                {
                    emptyIndexesCount--;
                }
                //for left index, remove an error if applicable
                while (leftValuesUsed > numberOfScansUsedForSmoothingOnEachSide)
                {
                    double errorToRemove = relativeErrors[leftIndex];
                    if (errorToRemove != 0)
                    {
                        smoothedCorrectionFactor -= errorToRemove;
                        leftValuesUsed--;
                    }
                    leftIndex++;
                }
                for (; rightIndex < (i + numberOfScansUsedForSmoothingOnEachSide + emptyIndexesCount) && rightIndex < relativeErrors.Length; rightIndex++)
                {
                    double addedFactor = relativeErrors[rightIndex];
                    if (addedFactor == 0)
                    {
                        emptyIndexesCount++;
                    }
                    else
                    {
                        smoothedCorrectionFactor += addedFactor;
                        rightValuesUsed++;
                    }
                }
                smoothedErrors[i] = smoothedCorrectionFactor / (rightValuesUsed + leftValuesUsed);
            }
            return smoothedErrors;
        }

        private MsDataScan CalibrateScan(MsDataScan oldScan, double smoothedRelativeError, double? precursorSmoothedRelativeError = null)
        {
            double correctionFactor = 1 - smoothedRelativeError; //wasn't in ppm before
            double[] originalMzs = oldScan.MassSpectrum.XArray;
            double[] calibratedMzs = new double[originalMzs.Length];
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
            MzSpectrum calibratedSpectrum = new MzSpectrum(calibratedMzs, oldScan.MassSpectrum.YArray, false);
            var scan = new MsDataScan(
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

            return scan;
        }
    }
}