using Chemistry;
using InternalLogic;
using MassSpectrometry;
using MathNet.Numerics.Statistics;
using Proteomics;
using Spectra;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text.RegularExpressions;

namespace InternalLogicCalibration
{
    public class CalibrationEngine : MyEngine
    {
        #region isotopologue parameters

        // THIS PARAMETER IS FRAGILE!!!
        // TUNED TO CORRESPOND TO SPECTROMETER OUTPUT
        // BETTER SPECTROMETERS WOULD HAVE BETTER (LOWER) RESOLUIONS
        // Parameter for isotopolouge distribution searching
        public double fineResolution = 0.1;

        #endregion isotopologue parameters

        public HashSet<int> MS2spectraToWatch;
        public HashSet<int> MS1spectraToWatch;
        public DoubleRange mzRange;

        public IMsDataFile<IMzSpectrum<MzPeak>> myMsDataFile;
        public List<InternalLogic.NewPsmWithFDR> identifications;


        public bool calibrateSpectra = true;
        internal int randomSeed;
        public string paramString = "";
        public int minMS2 = 2;
        public int minMS1 = 3;
        public double toleranceInMZforMS2Search;
        internal double toleranceInMZforMS1Search = 0.01;
        public HashSet<int> matchesToExclude;

        public CalibrationEngine(IMsDataFile<IMzSpectrum<MzPeak>> myMsDataFile, int randomSeed, double toleranceInMZforMS2Search)
        {
            this.myMsDataFile = myMsDataFile;
            MS1spectraToWatch = new HashSet<int>();
            MS2spectraToWatch = new HashSet<int>();
            this.randomSeed = randomSeed;
            this.toleranceInMZforMS2Search = toleranceInMZforMS2Search;
        }
        protected override MyResults RunSpecific()
        {
            output("Welcome to my software lock mass implementation");
            output("Calibrating " + Path.GetFileName(myMsDataFile.FilePath));

            List<int> trainingPointCounts = new List<int>();
            List<LabeledDataPoint> pointList;
            for (int calibrationRound = 1; ; calibrationRound++)
            {
                output("Calibration round " + calibrationRound);

                output("Getting Training Points");

                pointList = GetDataPoints(myMsDataFile, identifications, matchesToExclude);

                if (calibrationRound >= 2 && pointList.Count <= trainingPointCounts[calibrationRound - 2])
                    break;

                trainingPointCounts.Add(pointList.Count);

                var pointList1 = pointList.Where((b) => b.inputs[0] == 1).ToList();
                if (pointList1.Count == 0)
                {
                    output("Not enough MS1 training points, identification quality is poor");
                    return null;
                }
                WriteDataToFiles(pointList1, "pointList1" + myMsDataFile.Name + calibrationRound);
                output("pointList1.Count() = " + pointList1.Count());
                var pointList2 = pointList.Where((b) => b.inputs[0] == 2).ToList();
                if (pointList2.Count == 0)
                {
                    output("Not enough MS2 training points, identification quality is poor");
                    return null;
                }
                WriteDataToFiles(pointList2, "pointList2" + myMsDataFile.Name + calibrationRound);
                output("pointList2.Count() = " + pointList2.Count());

                CalibrationFunction identityPredictor = new IdentityCalibrationFunction();
                output("Uncalibrated MSE, " + identityPredictor.getMSE(pointList1) + "," + identityPredictor.getMSE(pointList2) + "," + identityPredictor.getMSE(pointList));

                CalibrationFunction combinedCalibration = Calibrate(pointList);

                combinedCalibration.writePredictedLables(pointList1, "pointList1predictedLabels" + myMsDataFile.Name + "CalibrationRound" + calibrationRound);
                combinedCalibration.writePredictedLables(pointList2, "pointList2predictedLabels" + myMsDataFile.Name + "CalibrationRound" + calibrationRound);
            }

            output("Post-processing");

            //postProcessing(p);

            output("Finished running my software lock mass implementation");

            return new CalibrationResults(this);
        }

        private CalibrationFunction Calibrate(List<LabeledDataPoint> trainingPoints)
        {
            var rnd = new Random(randomSeed);
            var shuffledTrainingPoints = trainingPoints.OrderBy(item => rnd.Next()).ToArray();

            var trainList = shuffledTrainingPoints.Take(trainingPoints.Count * 3 / 4).ToList();
            var testList = shuffledTrainingPoints.Skip(trainingPoints.Count * 3 / 4).ToList();

            var trainList1 = trainList.Where((b) => b.inputs[0] == 1).ToList();
            WriteDataToFiles(trainList1, "train1" + myMsDataFile.Name);
            output("trainList1.Count() = " + trainList1.Count());
            var trainList2 = trainList.Where((b) => b.inputs[0] == 2).ToList();
            WriteDataToFiles(trainList2, "train2" + myMsDataFile.Name);
            output("trainList2.Count() = " + trainList2.Count());
            var testList1 = testList.Where((b) => b.inputs[0] == 1).ToList();
            WriteDataToFiles(testList1, "test1" + myMsDataFile.Name);
            output("testList1.Count() = " + testList1.Count());
            var testList2 = testList.Where((b) => b.inputs[0] == 2).ToList();
            WriteDataToFiles(testList2, "test2" + myMsDataFile.Name);
            output("testList2.Count() = " + testList2.Count());

            CalibrationFunction bestMS1predictor = new IdentityCalibrationFunction();
            CalibrationFunction bestMS2predictor = new IdentityCalibrationFunction();
            CalibrationFunction combinedCalibration = new SeparateCalibrationFunction(bestMS1predictor, bestMS2predictor);
            double bestMS1MSE = bestMS1predictor.getMSE(testList1);
            double bestMS2MSE = bestMS2predictor.getMSE(testList2);
            double combinedMSE = combinedCalibration.getMSE(testList);
            output("Uncalibrated MSE, " + bestMS1MSE + "," + bestMS2MSE + "," + combinedMSE);

            CalibrationFunction ms1regressor = new ConstantCalibrationFunction();
            CalibrationFunction ms2regressor = new ConstantCalibrationFunction();
            ms1regressor.Train(trainList1);
            ms2regressor.Train(trainList2);
            combinedCalibration = new SeparateCalibrationFunction(ms1regressor, ms2regressor);
            combinedCalibration.writePredictedLables(trainList1, "trainList1Constant" + myMsDataFile.Name);
            combinedCalibration.writePredictedLables(trainList2, "trainList2Constant" + myMsDataFile.Name);
            combinedCalibration.writePredictedLables(testList1, "testList1Constant" + myMsDataFile.Name);
            combinedCalibration.writePredictedLables(testList2, "testList2Constant" + myMsDataFile.Name);
            double MS1mse = ms1regressor.getMSE(testList1);
            double MS2mse = ms2regressor.getMSE(testList2);
            combinedMSE = combinedCalibration.getMSE(testList);
            //OnOutput("Constant calibration MSE, " + MS1mse + "," + MS2mse + "," + combinedMSE));
            if (MS1mse < bestMS1MSE)
            {
                bestMS1MSE = MS1mse;
                bestMS1predictor = ms1regressor;
            }
            if (MS2mse < bestMS2MSE)
            {
                bestMS2MSE = MS2mse;
                bestMS2predictor = ms2regressor;
            }

            //ms1regressor = new ByHandCalibrationFunction(OnOutput, trainList1);
            //ms2regressor = new ByHandCalibrationFunction(OnOutput, trainList2);
            //combinedCalibration = new SeparateCalibrationFunction(ms1regressor, ms2regressor);
            //combinedCalibration.writeNewLabels(trainList1, "trainList1byHand" + myMsDataFile.Name);
            //combinedCalibration.writeNewLabels(trainList2, "trainList2byHand" + myMsDataFile.Name);
            //combinedCalibration.writeNewLabels(testList1, "testList1byHand" + myMsDataFile.Name);
            //combinedCalibration.writeNewLabels(testList2, "testList2byHand" + myMsDataFile.Name);
            //MS1mse = ms1regressor.getMSE(testList1);
            //MS2mse = ms2regressor.getMSE(testList2);
            //combinedMSE = combinedCalibration.getMSE(testList);
            //OnOutput(new OutputHandlerEventArgs("By hand calibration MSE, " + MS1mse + "," + MS2mse + "," + combinedMSE));
            //if (MS1mse < bestMS1MSE)
            //{
            //    bestMS1MSE = MS1mse;
            //    bestMS1predictor = ms1regressor;
            //}
            //if (MS2mse < bestMS2MSE)
            //{
            //    bestMS2MSE = MS2mse;
            //    bestMS2predictor = ms2regressor;
            //}

            List<TransformFunction> transforms = new List<TransformFunction>();

            transforms.Add(new TransformFunction(b => new double[1] { b[1] }, 1, "TFFFF"));
            transforms.Add(new TransformFunction(b => new double[1] { b[2] }, 1, "FTFFF"));
            transforms.Add(new TransformFunction(b => new double[1] { Math.Log(b[3]) }, 1, "FFTFF"));
            transforms.Add(new TransformFunction(b => new double[1] { Math.Log(b[4]) }, 1, "FFFTF"));
            transforms.Add(new TransformFunction(b => new double[1] { Math.Log(b[5]) }, 1, "FFFFT"));

            transforms.Add(new TransformFunction(b => new double[2] { b[1], b[2] }, 2, "TTFFF"));
            transforms.Add(new TransformFunction(b => new double[2] { b[1], Math.Log(b[3]) }, 2, "TFTFF"));
            transforms.Add(new TransformFunction(b => new double[2] { b[1], Math.Log(b[4]) }, 2, "TFFTF"));
            transforms.Add(new TransformFunction(b => new double[2] { b[1], Math.Log(b[5]) }, 2, "TFFFT"));
            transforms.Add(new TransformFunction(b => new double[2] { b[2], Math.Log(b[3]) }, 2, "FTTFF"));
            transforms.Add(new TransformFunction(b => new double[2] { b[2], Math.Log(b[4]) }, 2, "FTFTF"));
            transforms.Add(new TransformFunction(b => new double[2] { b[2], Math.Log(b[5]) }, 2, "FTFFT"));
            transforms.Add(new TransformFunction(b => new double[2] { Math.Log(b[3]), Math.Log(b[4]) }, 2, "FFTTF"));
            transforms.Add(new TransformFunction(b => new double[2] { Math.Log(b[3]), Math.Log(b[5]) }, 2, "FFTFT"));
            transforms.Add(new TransformFunction(b => new double[2] { Math.Log(b[4]), Math.Log(b[5]) }, 2, "FFFTT"));

            transforms.Add(new TransformFunction(b => new double[3] { b[1], b[2], Math.Log(b[3]) }, 3, "TTTFF"));
            transforms.Add(new TransformFunction(b => new double[3] { b[1], b[2], Math.Log(b[4]) }, 3, "TTFTF"));
            transforms.Add(new TransformFunction(b => new double[3] { b[1], b[2], Math.Log(b[5]) }, 3, "TTFFT"));
            transforms.Add(new TransformFunction(b => new double[3] { b[1], Math.Log(b[3]), Math.Log(b[4]) }, 3, "TFTTF"));
            transforms.Add(new TransformFunction(b => new double[3] { b[1], Math.Log(b[3]), Math.Log(b[5]) }, 3, "TFTFT"));
            transforms.Add(new TransformFunction(b => new double[3] { b[1], Math.Log(b[4]), Math.Log(b[5]) }, 3, "TFFTT"));
            transforms.Add(new TransformFunction(b => new double[3] { b[2], Math.Log(b[3]), Math.Log(b[4]) }, 3, "FTTTF"));
            transforms.Add(new TransformFunction(b => new double[3] { b[2], Math.Log(b[3]), Math.Log(b[5]) }, 3, "FTTFT"));
            transforms.Add(new TransformFunction(b => new double[3] { b[2], Math.Log(b[4]), Math.Log(b[5]) }, 3, "FTFTT"));
            transforms.Add(new TransformFunction(b => new double[3] { Math.Log(b[3]), Math.Log(b[4]), Math.Log(b[5]) }, 3, "FFTTT"));

            transforms.Add(new TransformFunction(b => new double[4] { b[1], b[2], Math.Log(b[3]), Math.Log(b[4]) }, 4, "TTTTF"));
            transforms.Add(new TransformFunction(b => new double[4] { b[1], b[2], Math.Log(b[3]), Math.Log(b[5]) }, 4, "TTTFT"));
            transforms.Add(new TransformFunction(b => new double[4] { b[1], b[2], Math.Log(b[4]), Math.Log(b[5]) }, 4, "TTFTT"));
            transforms.Add(new TransformFunction(b => new double[4] { b[1], Math.Log(b[3]), Math.Log(b[4]), Math.Log(b[5]) }, 4, "TFTTT"));
            transforms.Add(new TransformFunction(b => new double[4] { b[2], Math.Log(b[3]), Math.Log(b[4]), Math.Log(b[5]) }, 4, "FTTTT"));

            transforms.Add(new TransformFunction(b => new double[5] { b[1], b[2], Math.Log(b[3]), Math.Log(b[4]), Math.Log(b[5]) }, 5, "TTTTT"));

            try
            {
                foreach (var transform in transforms)
                {
                    ms1regressor = new LinearCalibrationFunctionMathNet(transform);
                    ms2regressor = new LinearCalibrationFunctionMathNet(transform);
                    ms1regressor.Train(trainList1);
                    ms2regressor.Train(trainList2);
                    MS1mse = ms1regressor.getMSE(testList1);
                    MS2mse = ms2regressor.getMSE(testList2);
                    //OnOutput(ms1regressor.name + transform.name + " MSE, " + MS1mse + "," + MS2mse));
                    if (MS1mse < bestMS1MSE)
                    {
                        bestMS1MSE = MS1mse;
                        bestMS1predictor = ms1regressor;
                        ms1regressor.writePredictedLables(trainList1, "train1" + ms1regressor.name + transform.name + myMsDataFile.Name);
                        ms1regressor.writePredictedLables(testList1, "test1" + ms1regressor.name + transform.name + myMsDataFile.Name);
                    }
                    if (MS2mse < bestMS2MSE)
                    {
                        bestMS2MSE = MS2mse;
                        bestMS2predictor = ms2regressor;
                        ms2regressor.writePredictedLables(trainList2, "train2" + ms2regressor.name + transform.name + myMsDataFile.Name);
                        ms2regressor.writePredictedLables(testList2, "test2" + ms2regressor.name + transform.name + myMsDataFile.Name);
                    }
                }
                foreach (var transform in transforms)
                {
                    ms1regressor = new QuadraticCalibrationFunctionMathNet(transform);
                    ms2regressor = new QuadraticCalibrationFunctionMathNet(transform);
                    ms1regressor.Train(trainList1);
                    ms2regressor.Train(trainList2);
                    MS1mse = ms1regressor.getMSE(testList1);
                    MS2mse = ms2regressor.getMSE(testList2);
                    //OnOutput(ms1regressor.name + transform.name + " MSE, " + MS1mse + "," + MS2mse));
                    if (MS1mse < bestMS1MSE)
                    {
                        bestMS1MSE = MS1mse;
                        bestMS1predictor = ms1regressor;
                        ms1regressor.writePredictedLables(trainList1, "train1" + ms1regressor.name + transform.name + myMsDataFile.Name);
                        ms1regressor.writePredictedLables(testList1, "test1" + ms1regressor.name + transform.name + myMsDataFile.Name);
                    }
                    if (MS2mse < bestMS2MSE)
                    {
                        bestMS2MSE = MS2mse;
                        bestMS2predictor = ms2regressor;
                        ms2regressor.writePredictedLables(trainList2, "train2" + ms2regressor.name + transform.name + myMsDataFile.Name);
                        ms2regressor.writePredictedLables(testList2, "test2" + ms2regressor.name + transform.name + myMsDataFile.Name);
                    }
                }
            }
            catch (ArgumentException e)
            {
                output("Could not calibrate: " + e.Message);
            }

            CalibrationFunction bestCf = new SeparateCalibrationFunction(bestMS1predictor, bestMS2predictor);

            status("Calibrating Spectra");

            CalibrateSpectra(bestCf);

            return bestCf;
        }

        private void CalibrateSpectra(CalibrationFunction bestCf)
        {
            foreach (var a in myMsDataFile)
            {
                if (a.MsnOrder == 2)
                {
                    if (MS2spectraToWatch.Contains(a.OneBasedScanNumber))
                    {
                        output("Calibrating scan number " + a.OneBasedScanNumber);
                        output(" before calibration:");
                        output(" " + string.Join(",", a.MassSpectrum.newSpectrumExtract(mzRange).xArray));
                    }

                    int oneBasedScanNumber;
                    a.TryGetPrecursorOneBasedScanNumber(out oneBasedScanNumber);
                    var precursorScan = myMsDataFile.GetOneBasedScan(oneBasedScanNumber);

                    double precursorMZ;
                    a.TryGetSelectedIonGuessMZ(out precursorMZ);
                    double precursorIntensity;
                    a.TryGetSelectedIonGuessIntensity(out precursorIntensity);
                    double newSelectedMZ = precursorMZ - bestCf.Predict(new double[6] { 1, precursorMZ, precursorScan.RetentionTime, precursorIntensity, precursorScan.TotalIonCurrent, precursorScan.InjectionTime });

                    double monoisotopicMZ;
                    a.TryGetSelectedIonGuessMonoisotopicMZ(out monoisotopicMZ);
                    double monoisotopicIntensity;
                    a.TryGetSelectedIonGuessMonoisotopicIntensity(out monoisotopicIntensity);

                    if (double.IsNaN(monoisotopicIntensity))
                        monoisotopicIntensity = precursorScan.MassSpectrum.GetClosestPeak(monoisotopicMZ).Intensity;

                    double newMonoisotopicMZ = monoisotopicMZ - bestCf.Predict(new double[6] { 1, monoisotopicMZ, precursorScan.RetentionTime, monoisotopicIntensity, precursorScan.TotalIonCurrent, precursorScan.InjectionTime });

                    int SelectedIonGuessChargeStateGuess;
                    a.TryGetSelectedIonGuessChargeStateGuess(out SelectedIonGuessChargeStateGuess);
                    double IsolationMZ;
                    a.TryGetIsolationMZ(out IsolationMZ);

                    Func<MzPeak, double> theFunc = x => x.MZ - bestCf.Predict(new double[9] { 2, x.MZ, a.RetentionTime, x.Intensity, a.TotalIonCurrent, a.InjectionTime, SelectedIonGuessChargeStateGuess, IsolationMZ, (x.MZ - a.ScanWindowRange.Minimum) / (a.ScanWindowRange.Maximum - a.ScanWindowRange.Minimum) });
                    a.tranformByApplyingFunctionsToSpectraAndReplacingPrecursorMZs(theFunc, newSelectedMZ, newMonoisotopicMZ);

                    if (MS2spectraToWatch.Contains(a.OneBasedScanNumber))
                    {
                        output(" after calibration:");
                        output(" precursorMZ:" + precursorMZ);
                        output(" monoisotopicMZ:" + monoisotopicMZ);
                        output(" newSelectedMZ:" + newSelectedMZ);
                        output(" newMonoisotopicMZ:" + newMonoisotopicMZ);
                        output(" " + string.Join(",", a.MassSpectrum.newSpectrumExtract(mzRange).xArray));
                    }
                }
                else
                {
                    if (MS1spectraToWatch.Contains(a.OneBasedScanNumber))
                    {
                        output("Calibrating scan number " + a.OneBasedScanNumber);
                        output(" before calibration:");
                        output(" " + string.Join(",", a.MassSpectrum.newSpectrumExtract(mzRange).xArray));
                    }
                    Func<MzPeak, double> theFUnc = x => x.MZ - bestCf.Predict(new double[6] { 1, x.MZ, a.RetentionTime, x.Intensity, a.TotalIonCurrent, a.InjectionTime });
                    a.tranformByApplyingFunctionsToSpectraAndReplacingPrecursorMZs(theFUnc, double.NaN, double.NaN); if (MS1spectraToWatch.Contains(a.OneBasedScanNumber))
                    {
                        output(" after calibration:");
                        output(string.Join(",", a.MassSpectrum.newSpectrumExtract(mzRange).xArray));
                    }
                }
            }
        }

        public void WriteDataToFiles(IEnumerable<LabeledDataPoint> trainingPoints, string prefix)
        {
            if (trainingPoints.Count() == 0)
                return;
            var fullFileName = Path.Combine(@"DataPoints", prefix + ".dat");
            Directory.CreateDirectory(Path.GetDirectoryName(fullFileName));

            using (StreamWriter file = new StreamWriter(fullFileName))
            {
                if (trainingPoints.First().inputs.Count() == 9)
                    file.WriteLine("MS, MZ, RetentionTime, Intensity,TotalIonCurrent, InjectionTime, SelectedIonGuessChargeStateGuess, IsolationMZ, relativeMZ, label");
                else
                    file.WriteLine("MS, MZ, RetentionTime, Intensity,TotalIonCurrent, InjectionTime, label");
                foreach (LabeledDataPoint d in trainingPoints)
                    file.WriteLine(string.Join(",", d.inputs) + "," + d.output);
            }
        }

        public override void ValidateParams()
        {
            throw new NotImplementedException();
        }
        private const int numFragmentsNeeded = 10;

        public List<LabeledDataPoint> GetDataPoints(IMsDataFile<IMzSpectrum<MzPeak>> myMsDataFile, List<InternalLogic.NewPsmWithFDR> identifications, HashSet<int> matchesToExclude)
        {
            output("Extracting data points:");
            // The final training point list
            List<LabeledDataPoint> trainingPointsToReturn = new List<LabeledDataPoint>();

            // Set of peaks, identified by m/z and retention time. If a peak is in here, it means it has been a part of an accepted identification, and should be rejected
            HashSet<Tuple<double, double>> peaksAddedFromMS1HashSet = new HashSet<Tuple<double, double>>();

            int numIdentifications = identifications.Count;
            // Loop over identifications
            for (int matchIndex = 0; matchIndex < numIdentifications; matchIndex++)
            {
                var identification = identifications[matchIndex];
                if (identification.QValue > 0.01)
                    break;

                // Progress
                if (numIdentifications < 100 || matchIndex % (numIdentifications / 100) == 0)
                    ReportProgress(new InternalLogic.ProgressEventArgs(100 * matchIndex / numIdentifications, "Looking at identifications..."));

                // Skip decoys, they are for sure not there!
                if (identification.isDecoy)
                    continue;

                // Skip exclusions! These are for testing
                if (matchesToExclude.Contains(matchIndex))
                    continue;

                // Each identification has an MS2 spectrum attached to it.
                int ms2spectrumIndex = identification.thisPSM.newPsm.scanNumber;

                // Get the peptide, don't forget to add the modifications!!!!
                var SequenceWithChemicalFormulas = identification.thisPSM.SequenceWithChemicalFormulas;
                int peptideCharge = identification.thisPSM.newPsm.scanPrecursorCharge;

                #region watch

                if (MS2spectraToWatch.Contains(ms2spectrumIndex))
                {
                    output("ms2spectrumIndex: " + ms2spectrumIndex);
                    output(" peptide: " + SequenceWithChemicalFormulas);
                }

                #endregion watch

                int numFragmentsIdentified = -1;
                List<LabeledDataPoint> candidateTrainingPointsForPeptide = new List<LabeledDataPoint>();
                Proteomics.Peptide coolPeptide = null;
                try
                {
                    coolPeptide = new Peptide(SequenceWithChemicalFormulas);
                }
                catch
                {
                    continue;
                }

                candidateTrainingPointsForPeptide = SearchMS2Spectrum(myMsDataFile.GetOneBasedScan(ms2spectrumIndex), coolPeptide, peptideCharge, out numFragmentsIdentified);

                //SoftwareLockMassRunner.WriteDataToFiles(candidateTrainingPointsForPeptide, ms2spectrumIndex.ToString());

                //RTBoutput(new OutputHandlerEventArgs(numFragmentsIdentified.ToString()));

                // If MS2 has low evidence for peptide, skip and go to next one
                if (numFragmentsIdentified < numFragmentsNeeded)
                    continue;

                // Calculate isotopic distribution of the full peptide

                IsotopicDistribution dist = new IsotopicDistribution(coolPeptide.GetChemicalFormula(), fineResolution, 0.001);

                double[] masses = new double[dist.Masses.Count];
                double[] intensities = new double[dist.Intensities.Count];
                for (int i = 0; i < dist.Masses.Count; i++)
                {
                    masses[i] = dist.Masses[i];
                    intensities[i] = dist.Intensities[i];
                }
                Array.Sort(intensities, masses, Comparer<double>.Create((x, y) => y.CompareTo(x)));

                if (MS2spectraToWatch.Contains(ms2spectrumIndex))
                {
                    output(" Isotopologue distribution: ");
                    output(" masses = " + string.Join(", ", masses) + "...");
                    output(" intensities = " + string.Join(", ", intensities) + "...");
                }

                int lowestMS1ind = SearchMS1Spectra(myMsDataFile, masses, intensities, candidateTrainingPointsForPeptide, ms2spectrumIndex, -1, peaksAddedFromMS1HashSet, peptideCharge);

                int highestMS1ind = SearchMS1Spectra(myMsDataFile, masses, intensities, candidateTrainingPointsForPeptide, ms2spectrumIndex, 1, peaksAddedFromMS1HashSet, peptideCharge);

                if (MS2spectraToWatch.Contains(ms2spectrumIndex))
                {
                    output(" ms1range: " + lowestMS1ind + " to " + highestMS1ind);
                }
                trainingPointsToReturn.AddRange(candidateTrainingPointsForPeptide);
            }
            output("Number of training points: " + trainingPointsToReturn.Count());
            return trainingPointsToReturn;
        }

        private int GetOneBasedSpectrumIndexFromID(string spectrumID, IMsDataFile<IMzSpectrum<MzPeak>> myMsDataFile)
        {
            var lastInt = Convert.ToInt32(Regex.Match(spectrumID, @"\d+$").Value);
            if (lastInt > 0 && lastInt <= myMsDataFile.Count() && myMsDataFile.GetOneBasedScan(lastInt).id.Equals(spectrumID))
                return lastInt;
            if (lastInt - 1 > 0 && lastInt - 1 <= myMsDataFile.Count() && myMsDataFile.GetOneBasedScan(lastInt - 1).id.Equals(spectrumID))
                return lastInt - 1;
            if (lastInt + 1 > 0 && lastInt + 1 <= myMsDataFile.Count() && myMsDataFile.GetOneBasedScan(lastInt + 1).id.Equals(spectrumID))
                return lastInt + 1;

            int oneBasedScanNumber = 0;
            do
            {
                oneBasedScanNumber++;
            } while (!spectrumID.Equals(myMsDataFile.GetOneBasedScan(oneBasedScanNumber).id));
            return oneBasedScanNumber;
        }

        private List<LabeledDataPoint> SearchMS2Spectrum(IMsDataScan<IMzSpectrum<MzPeak>> ms2DataScan, Proteomics.Peptide peptide, int peptideCharge, out int candidateFragmentsIdentified)
        {
            List<LabeledDataPoint> myCandidatePoints = new List<LabeledDataPoint>();

            // Key: mz value, Value: error
            Dictionary<double, double> addedPeaks = new Dictionary<double, double>();

            int SelectedIonGuessChargeStateGuess;
            ms2DataScan.TryGetSelectedIonGuessChargeStateGuess(out SelectedIonGuessChargeStateGuess);
            double IsolationMZ;
            ms2DataScan.TryGetIsolationMZ(out IsolationMZ);

            int ms2spectrumIndex = ms2DataScan.OneBasedScanNumber;

            var countForThisMS2 = 0;
            var countForThisMS2a = 0;
            var numFragmentsIdentified = 0;

            var scanWindowRange = ms2DataScan.ScanWindowRange;

            Fragment[] fragmentList = peptide.Fragment(FragmentTypes.b | FragmentTypes.y, true).ToArray();

            if (MS2spectraToWatch.Contains(ms2spectrumIndex))
            {
                output(" Considering individual fragments:");
            }

            foreach (IHasChemicalFormula fragment in fragmentList)
            {
                bool fragmentIdentified = false;
                bool computedIsotopologues = false;
                double[] masses = new double[0];
                double[] intensities = new double[0];
                // First look for monoisotopic masses, do not compute distribution spectrum!
                if (MS2spectraToWatch.Contains(ms2spectrumIndex))
                {
                    output("  Considering fragment " + (fragment as Fragment).Sequence + " with formula " + fragment.ThisChemicalFormula.Formula);
                    //if ((fragment as Fragment).Modifications.Count() > 0)
                    //RTBoutput("  Modifications: " + string.Join(", ", (fragment as Fragment).Modifications));
                }

                #region loop to determine if need to compute isotopologue distribution

                for (int chargeToLookAt = 1; chargeToLookAt <= peptideCharge; chargeToLookAt++)
                {
                    var monoisotopicMZ = fragment.MonoisotopicMass.ToMassToChargeRatio(chargeToLookAt);
                    if (monoisotopicMZ > scanWindowRange.Maximum)
                        continue;
                    if (monoisotopicMZ < scanWindowRange.Minimum)
                        break;
                    var closestPeakMZ = ms2DataScan.MassSpectrum.GetClosestPeakXvalue(monoisotopicMZ);
                    if (Math.Abs(closestPeakMZ - monoisotopicMZ) < toleranceInMZforMS2Search)
                    {
                        if (!computedIsotopologues)
                        {
                            if (MS2spectraToWatch.Contains(ms2spectrumIndex))
                            {
                                output("    Computing isotopologues because absolute error " + Math.Abs(closestPeakMZ - monoisotopicMZ) + " is smaller than tolerance " + toleranceInMZforMS2Search);
                                output("    Charge was = " + chargeToLookAt + "  closestPeakMZ = " + closestPeakMZ + " while monoisotopicMZ = " + monoisotopicMZ);
                            }

                            IsotopicDistribution dist = new IsotopicDistribution(fragment.ThisChemicalFormula, fineResolution, 0.001);

                            masses = new double[dist.Masses.Count];
                            intensities = new double[dist.Intensities.Count];
                            for (int i = 0; i < dist.Masses.Count; i++)
                            {
                                masses[i] = dist.Masses[i];
                                intensities[i] = dist.Intensities[i];
                            }
                            Array.Sort(intensities, masses, Comparer<double>.Create((x, y) => y.CompareTo(x)));
                            computedIsotopologues = true;
                            if (MS2spectraToWatch.Contains(ms2spectrumIndex))
                            {
                                output("    Isotopologue distribution: ");
                                output("    masses = " + string.Join(", ", masses) + "...");
                                output("    intensities = " + string.Join(", ", intensities) + "...");
                            }

                            break;
                        }
                    }
                }

                #endregion loop to determine if need to compute isotopologue distribution

                if (computedIsotopologues)
                {
                    #region actually add training points

                    if (MS2spectraToWatch.Contains(ms2spectrumIndex))
                    {
                        output("   Considering individual charges, to get training points:");
                    }
                    bool startingToAdd = false;
                    for (int chargeToLookAt = 1; chargeToLookAt <= peptideCharge; chargeToLookAt++)
                    {
                        if (MS2spectraToWatch.Contains(ms2spectrumIndex))
                        {
                            output("    Considering charge " + chargeToLookAt);
                        }
                        if (masses.First().ToMassToChargeRatio(chargeToLookAt) > scanWindowRange.Maximum)
                        {
                            if (MS2spectraToWatch.Contains(ms2spectrumIndex))
                            {
                                output("    Out of range: too high");
                            }
                            continue;
                        }
                        if (masses.Last().ToMassToChargeRatio(chargeToLookAt) < scanWindowRange.Minimum)
                        {
                            if (MS2spectraToWatch.Contains(ms2spectrumIndex))
                            {
                                output("    Out of range: too low");
                            }
                            break;
                        }
                        List<TrainingPoint> trainingPointsToAverage = new List<TrainingPoint>();
                        foreach (double a in masses)
                        {
                            double theMZ = a.ToMassToChargeRatio(chargeToLookAt);
                            var npwr = ms2DataScan.MassSpectrum.NumPeaksWithinRange(theMZ - toleranceInMZforMS2Search, theMZ + toleranceInMZforMS2Search);
                            if (npwr == 0)
                            {
                                if (MS2spectraToWatch.Contains(ms2spectrumIndex))
                                    output("     Breaking because extracted.Count = " + npwr);
                                break;
                            }
                            if (npwr > 1)
                            {
                                if (MS2spectraToWatch.Contains(ms2spectrumIndex))
                                    output("     Not looking for " + theMZ + " because extracted.Count = " + npwr);
                                continue;
                            }
                            var closestPeak = ms2DataScan.MassSpectrum.GetClosestPeak(theMZ);
                            var closestPeakMZ = closestPeak.MZ;
                            if (MS2spectraToWatch.Contains(ms2spectrumIndex))
                            {
                                output("     Found       " + closestPeakMZ + "   Error is    " + (closestPeakMZ - theMZ));
                            }
                            if (!addedPeaks.ContainsKey(closestPeakMZ))
                            {
                                addedPeaks.Add(closestPeakMZ, Math.Abs(closestPeakMZ - theMZ));
                                trainingPointsToAverage.Add(new TrainingPoint(new DataPoint(closestPeakMZ, double.NaN, 0, closestPeak.Intensity, double.NaN, double.NaN), closestPeakMZ - theMZ));
                            }
                            else
                            {
                                if (MS2spectraToWatch.Contains(ms2spectrumIndex))
                                {
                                    output("     Not using because already added peak");
                                }
                            }
                        }
                        // If started adding and suddnely stopped, go to next one, no need to look at higher charges
                        if (trainingPointsToAverage.Count == 0 && startingToAdd == true)
                            break;
                        if (trainingPointsToAverage.Count < Math.Min(minMS2, intensities.Count()))
                        {
                            if (MS2spectraToWatch.Contains(ms2spectrumIndex))
                            {
                                output("    Not adding, since not enought isotopes were seen");
                            }
                        }
                        else
                        {
                            startingToAdd = true;
                            if (!fragmentIdentified)
                            {
                                fragmentIdentified = true;
                                numFragmentsIdentified += 1;
                            }

                            countForThisMS2 += trainingPointsToAverage.Count;
                            countForThisMS2a += 1;

                            double addedMZ = trainingPointsToAverage.Select(b => b.dp.mz).Average();
                            double relativeMZ = (addedMZ - ms2DataScan.ScanWindowRange.Minimum) / (ms2DataScan.ScanWindowRange.Maximum - ms2DataScan.ScanWindowRange.Minimum);
                            double[] inputs = new double[9] { 2, addedMZ, ms2DataScan.RetentionTime, trainingPointsToAverage.Select(b => b.dp.intensity).Average(), ms2DataScan.TotalIonCurrent, ms2DataScan.InjectionTime, SelectedIonGuessChargeStateGuess, IsolationMZ, relativeMZ };
                            var a = new LabeledDataPoint(inputs, trainingPointsToAverage.Select(b => b.l).Median());

                            if (MS2spectraToWatch.Contains(ms2spectrumIndex))
                            {
                                output("    Adding aggregate of " + trainingPointsToAverage.Count + " points FROM MS2 SPECTRUM");
                                output("    a.dmz " + a.inputs[1]);
                                output("    a.drt " + a.inputs[2]);
                                output("    a.l     " + a.output);
                            }
                            myCandidatePoints.Add(a);
                        }
                    }

                    #endregion actually add training points
                }
            }

            //RTBoutput(new OutputHandlerEventArgs("ind = " + ms2spectrumIndex + " count = " + countForThisMS2 + " count2 = " + countForThisMS2a + " fragments = " + numFragmentsIdentified));
            if (MS2spectraToWatch.Contains(ms2spectrumIndex))
            {
                output(" countForThisMS2 = " + countForThisMS2);
                output(" countForThisMS2a = " + countForThisMS2a);
                output(" numFragmentsIdentified = " + numFragmentsIdentified);
            }
            candidateFragmentsIdentified = numFragmentsIdentified;
            return myCandidatePoints;
        }

        private int SearchMS1Spectra(IMsDataFile<IMzSpectrum<MzPeak>> myMsDataFile, double[] originalMasses, double[] originalIntensities, List<LabeledDataPoint> myCandidatePoints, int ms2spectrumIndex, int direction, HashSet<Tuple<double, double>> peaksAddedHashSet, int peptideCharge)
        {
            int goodIndex = -1;
            List<int> scores = new List<int>();
            var theIndex = -1;
            if (direction == 1)
                theIndex = ms2spectrumIndex;
            else
                theIndex = ms2spectrumIndex - 1;

            bool addedAscan = true;

            int highestKnownChargeForThisPeptide = peptideCharge;
            while (theIndex >= 1 && theIndex <= myMsDataFile.NumSpectra && addedAscan == true)
            {
                int countForThisScan = 0;
                if (myMsDataFile.GetOneBasedScan(theIndex).MsnOrder > 1)
                {
                    theIndex += direction;
                    continue;
                }
                addedAscan = false;
                if (MS2spectraToWatch.Contains(ms2spectrumIndex) || MS1spectraToWatch.Contains(theIndex))
                {
                    output(" Looking in MS1 spectrum " + theIndex + " because of MS2 spectrum " + ms2spectrumIndex);
                }
                List<LabeledDataPoint> myCandidatePointsForThisMS1scan = new List<LabeledDataPoint>();
                var fullMS1scan = myMsDataFile.GetOneBasedScan(theIndex);
                double ms1RetentionTime = fullMS1scan.RetentionTime;
                var scanWindowRange = fullMS1scan.ScanWindowRange;
                var fullMS1spectrum = fullMS1scan.MassSpectrum;
                if (fullMS1spectrum.Count == 0)
                    break;

                bool startingToAddCharges = false;
                int chargeToLookAt = 1;
                do
                {
                    if (MS2spectraToWatch.Contains(ms2spectrumIndex) || MS1spectraToWatch.Contains(theIndex))
                    {
                        output("  Looking at charge " + chargeToLookAt);
                    }
                    if (originalMasses[0].ToMassToChargeRatio(chargeToLookAt) > scanWindowRange.Maximum)
                    {
                        chargeToLookAt++;
                        continue;
                    }
                    if (originalMasses[0].ToMassToChargeRatio(chargeToLookAt) < scanWindowRange.Minimum)
                        break;
                    List<TrainingPoint> trainingPointsToAverage = new List<TrainingPoint>();
                    foreach (double a in originalMasses)
                    {
                        double theMZ = a.ToMassToChargeRatio(chargeToLookAt);

                        var npwr = fullMS1spectrum.NumPeaksWithinRange(theMZ - toleranceInMZforMS1Search, theMZ + toleranceInMZforMS1Search);
                        if (npwr == 0)
                        {
                            if (MS1spectraToWatch.Contains(theIndex))
                                output("      Breaking because extracted.Count = " + npwr);
                            break;
                        }
                        if (npwr > 1)
                        {
                            if (MS1spectraToWatch.Contains(theIndex))
                                output("      Not looking for " + theMZ + " because extracted.Count = " + npwr);
                            continue;
                        }

                        var closestPeak = fullMS1spectrum.GetClosestPeak(theMZ);
                        var closestPeakMZ = closestPeak.MZ;

                        //if (closestPeak.Intensity / fullMS1scan.TotalIonCurrent < 2e-4)
                        //{
                        //    if (MS2spectraToWatch.Contains(ms2spectrumIndex) || MS1spectraToWatch.Contains(theIndex))
                        //        RTBoutput("      Breaking because intensity fraction is " + closestPeak.Intensity / fullMS1scan.TotalIonCurrent);
                        //    break;
                        //}
                        if ((MS2spectraToWatch.Contains(ms2spectrumIndex) || MS1spectraToWatch.Contains(theIndex)) && mzRange.Contains(theMZ))
                        {
                            output("      Looking for " + theMZ + " found " + closestPeakMZ + " error is " + (closestPeakMZ - theMZ));
                        }

                        var theTuple = Tuple.Create(closestPeakMZ, ms1RetentionTime);
                        if (!peaksAddedHashSet.Contains(theTuple))
                        {
                            peaksAddedHashSet.Add(theTuple);
                            highestKnownChargeForThisPeptide = Math.Max(highestKnownChargeForThisPeptide, chargeToLookAt);
                            //if ((MS2spectraToWatch.Contains(ms2spectrumIndex) || MS1spectraToWatch.Contains(theIndex)) && mzRange.Contains(theMZ))
                            if ((MS2spectraToWatch.Contains(ms2spectrumIndex) || MS1spectraToWatch.Contains(theIndex)))
                            {
                                output("      Found " + closestPeakMZ + ", was looking for " + theMZ + ", e=" + (closestPeakMZ - theMZ) + ", if=" + closestPeak.Intensity / fullMS1scan.TotalIonCurrent);
                            }
                            trainingPointsToAverage.Add(new TrainingPoint(new DataPoint(closestPeakMZ, double.NaN, 1, closestPeak.Intensity, double.NaN, double.NaN), closestPeakMZ - theMZ));
                        }
                        else
                            break;
                    }
                    // If started adding and suddnely stopped, go to next one, no need to look at higher charges
                    if (trainingPointsToAverage.Count == 0 && startingToAddCharges == true)
                    {
                        if (MS2spectraToWatch.Contains(ms2spectrumIndex) || MS1spectraToWatch.Contains(theIndex))
                        {
                            output("    Started adding and suddnely stopped, no need to look at higher charges");
                        }
                        break;
                    }
                    if ((trainingPointsToAverage.Count == 0 || (trainingPointsToAverage.Count == 1 && originalIntensities[0] < 0.65)) && (peptideCharge <= chargeToLookAt))
                    {
                        if (MS2spectraToWatch.Contains(ms2spectrumIndex) || MS1spectraToWatch.Contains(theIndex))
                        {
                            output("    Did not find (or found without isotopes) charge " + chargeToLookAt + ", no need to look at higher charges");
                        }
                        break;
                    }
                    if (trainingPointsToAverage.Count == 1 && originalIntensities[0] < 0.65)
                    {
                        if ((MS2spectraToWatch.Contains(ms2spectrumIndex) || MS1spectraToWatch.Contains(theIndex)))
                        {
                            output("    Not adding, since originalIntensities[0] is " + originalIntensities[0] + " which is too low");
                        }
                    }
                    else if (trainingPointsToAverage.Count < Math.Min(minMS1, originalIntensities.Count()))
                    {
                        if ((MS2spectraToWatch.Contains(ms2spectrumIndex) || MS1spectraToWatch.Contains(theIndex)))
                        {
                            output("    Not adding, since count " + trainingPointsToAverage.Count + " is too low");
                        }
                    }
                    else
                    {
                        addedAscan = true;
                        startingToAddCharges = true;
                        countForThisScan += 1;
                        double[] inputs = new double[6] { 1, trainingPointsToAverage.Select(b => b.dp.mz).Average(), fullMS1scan.RetentionTime, trainingPointsToAverage.Select(b => b.dp.intensity).Average(), fullMS1scan.TotalIonCurrent, fullMS1scan.InjectionTime };
                        var a = new LabeledDataPoint(inputs, trainingPointsToAverage.Select(b => b.l).Median());
                        //if (a.output > 0)
                        //    RTBoutput(theIndex + "," + ms2spectrumIndex);
                        if (MS2spectraToWatch.Contains(ms2spectrumIndex) || MS1spectraToWatch.Contains(theIndex))
                        {
                            output("    Adding aggregate of " + trainingPointsToAverage.Count + " points FROM MS1 SPECTRUM");
                            output("    a.dmz " + a.inputs[1]);
                            output("    a.drt " + a.inputs[2]);
                            output("    a.l     " + a.output);
                        }
                        myCandidatePointsForThisMS1scan.Add(a);
                    }
                    chargeToLookAt++;
                } while (chargeToLookAt <= highestKnownChargeForThisPeptide + 1);

                if (myCandidatePointsForThisMS1scan.Count > 0)
                    goodIndex = theIndex;
                //    SoftwareLockMassRunner.WriteDataToFiles(myCandidatePointsForThisMS1scan, theIndex.ToString());
                myCandidatePoints.AddRange(myCandidatePointsForThisMS1scan);

                scores.Add(countForThisScan);
                theIndex += direction;
            }
            return goodIndex;
        }

    }
}