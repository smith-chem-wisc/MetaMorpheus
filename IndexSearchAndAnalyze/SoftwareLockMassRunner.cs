using Spectra;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace IndexSearchAndAnalyze
{
    public static class SoftwareLockMassRunner
    {
        public static void Run(SoftwareLockMassParams p)
        {
            p.po.output("Welcome to my software lock mass implementation");
            p.po.output("Calibrating " + Path.GetFileName(p.myMsDataFile.FilePath));

            List<int> trainingPointCounts = new List<int>();
            List<LabeledDataPoint> pointList;
            for (int calibrationRound = 1; ; calibrationRound++)
            {
                p.po.output("Calibration round " + calibrationRound);

                p.po.output("Getting Training Points");

                pointList = TrainingPointsExtractor.GetDataPoints(p.myMsDataFile, p.identifications, p, p.matchesToExclude);

                if (calibrationRound >= 2 && pointList.Count <= trainingPointCounts[calibrationRound - 2])
                    break;

                trainingPointCounts.Add(pointList.Count);

                var pointList1 = pointList.Where((b) => b.inputs[0] == 1).ToList();
                if (pointList1.Count == 0)
                {
                    p.po.output("Not enough MS1 training points, identification quality is poor");
                    return;
                }
                WriteDataToFiles(pointList1, "pointList1" + p.myMsDataFile.Name + calibrationRound);
                p.po.output("pointList1.Count() = " + pointList1.Count());
                var pointList2 = pointList.Where((b) => b.inputs[0] == 2).ToList();
                if (pointList2.Count == 0)
                {
                    p.po.output("Not enough MS2 training points, identification quality is poor");
                    return;
                }
                WriteDataToFiles(pointList2, "pointList2" + p.myMsDataFile.Name + calibrationRound);
                p.po.output("pointList2.Count() = " + pointList2.Count());

                CalibrationFunction identityPredictor = new IdentityCalibrationFunction();
                p.po.output("Uncalibrated MSE, " + identityPredictor.getMSE(pointList1) + "," + identityPredictor.getMSE(pointList2) + "," + identityPredictor.getMSE(pointList));

                CalibrationFunction combinedCalibration = Calibrate(pointList, p);

                combinedCalibration.writePredictedLables(pointList1, "pointList1predictedLabels" + p.myMsDataFile.Name + "CalibrationRound" + calibrationRound);
                combinedCalibration.writePredictedLables(pointList2, "pointList2predictedLabels" + p.myMsDataFile.Name + "CalibrationRound" + calibrationRound);
            }

            p.po.output("Post-processing");

            p.postProcessing(p);

            p.po.output("Finished running my software lock mass implementation");
        }

        private static CalibrationFunction Calibrate(List<LabeledDataPoint> trainingPoints, SoftwareLockMassParams p)
        {
            var rnd = new Random(p.randomSeed);
            var shuffledTrainingPoints = trainingPoints.OrderBy(item => rnd.Next()).ToArray();

            var trainList = shuffledTrainingPoints.Take(trainingPoints.Count * 3 / 4).ToList();
            var testList = shuffledTrainingPoints.Skip(trainingPoints.Count * 3 / 4).ToList();

            var trainList1 = trainList.Where((b) => b.inputs[0] == 1).ToList();
            WriteDataToFiles(trainList1, "train1" + p.myMsDataFile.Name);
            p.po.output("trainList1.Count() = " + trainList1.Count());
            var trainList2 = trainList.Where((b) => b.inputs[0] == 2).ToList();
            WriteDataToFiles(trainList2, "train2" + p.myMsDataFile.Name);
            p.po.output("trainList2.Count() = " + trainList2.Count());
            var testList1 = testList.Where((b) => b.inputs[0] == 1).ToList();
            WriteDataToFiles(testList1, "test1" + p.myMsDataFile.Name);
            p.po.output("testList1.Count() = " + testList1.Count());
            var testList2 = testList.Where((b) => b.inputs[0] == 2).ToList();
            WriteDataToFiles(testList2, "test2" + p.myMsDataFile.Name);
            p.po.output("testList2.Count() = " + testList2.Count());

            CalibrationFunction bestMS1predictor = new IdentityCalibrationFunction();
            CalibrationFunction bestMS2predictor = new IdentityCalibrationFunction();
            CalibrationFunction combinedCalibration = new SeparateCalibrationFunction(bestMS1predictor, bestMS2predictor);
            double bestMS1MSE = bestMS1predictor.getMSE(testList1);
            double bestMS2MSE = bestMS2predictor.getMSE(testList2);
            double combinedMSE = combinedCalibration.getMSE(testList);
            p.po.output("Uncalibrated MSE, " + bestMS1MSE + "," + bestMS2MSE + "," + combinedMSE);

            CalibrationFunction ms1regressor = new ConstantCalibrationFunction();
            CalibrationFunction ms2regressor = new ConstantCalibrationFunction();
            ms1regressor.Train(trainList1);
            ms2regressor.Train(trainList2);
            combinedCalibration = new SeparateCalibrationFunction(ms1regressor, ms2regressor);
            combinedCalibration.writePredictedLables(trainList1, "trainList1Constant" + p.myMsDataFile.Name);
            combinedCalibration.writePredictedLables(trainList2, "trainList2Constant" + p.myMsDataFile.Name);
            combinedCalibration.writePredictedLables(testList1, "testList1Constant" + p.myMsDataFile.Name);
            combinedCalibration.writePredictedLables(testList2, "testList2Constant" + p.myMsDataFile.Name);
            double MS1mse = ms1regressor.getMSE(testList1);
            double MS2mse = ms2regressor.getMSE(testList2);
            combinedMSE = combinedCalibration.getMSE(testList);
            //p.OnOutput("Constant calibration MSE, " + MS1mse + "," + MS2mse + "," + combinedMSE));
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

            //ms1regressor = new ByHandCalibrationFunction(p.OnOutput, trainList1);
            //ms2regressor = new ByHandCalibrationFunction(p.OnOutput, trainList2);
            //combinedCalibration = new SeparateCalibrationFunction(ms1regressor, ms2regressor);
            //combinedCalibration.writeNewLabels(trainList1, "trainList1byHand" + p.myMsDataFile.Name);
            //combinedCalibration.writeNewLabels(trainList2, "trainList2byHand" + p.myMsDataFile.Name);
            //combinedCalibration.writeNewLabels(testList1, "testList1byHand" + p.myMsDataFile.Name);
            //combinedCalibration.writeNewLabels(testList2, "testList2byHand" + p.myMsDataFile.Name);
            //MS1mse = ms1regressor.getMSE(testList1);
            //MS2mse = ms2regressor.getMSE(testList2);
            //combinedMSE = combinedCalibration.getMSE(testList);
            //p.OnOutput(new OutputHandlerEventArgs("By hand calibration MSE, " + MS1mse + "," + MS2mse + "," + combinedMSE));
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
                    //p.OnOutput(ms1regressor.name + transform.name + " MSE, " + MS1mse + "," + MS2mse));
                    if (MS1mse < bestMS1MSE)
                    {
                        bestMS1MSE = MS1mse;
                        bestMS1predictor = ms1regressor;
                        ms1regressor.writePredictedLables(trainList1, "train1" + ms1regressor.name + transform.name + p.myMsDataFile.Name);
                        ms1regressor.writePredictedLables(testList1, "test1" + ms1regressor.name + transform.name + p.myMsDataFile.Name);
                    }
                    if (MS2mse < bestMS2MSE)
                    {
                        bestMS2MSE = MS2mse;
                        bestMS2predictor = ms2regressor;
                        ms2regressor.writePredictedLables(trainList2, "train2" + ms2regressor.name + transform.name + p.myMsDataFile.Name);
                        ms2regressor.writePredictedLables(testList2, "test2" + ms2regressor.name + transform.name + p.myMsDataFile.Name);
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
                    //p.OnOutput(ms1regressor.name + transform.name + " MSE, " + MS1mse + "," + MS2mse));
                    if (MS1mse < bestMS1MSE)
                    {
                        bestMS1MSE = MS1mse;
                        bestMS1predictor = ms1regressor;
                        ms1regressor.writePredictedLables(trainList1, "train1" + ms1regressor.name + transform.name + p.myMsDataFile.Name);
                        ms1regressor.writePredictedLables(testList1, "test1" + ms1regressor.name + transform.name + p.myMsDataFile.Name);
                    }
                    if (MS2mse < bestMS2MSE)
                    {
                        bestMS2MSE = MS2mse;
                        bestMS2predictor = ms2regressor;
                        ms2regressor.writePredictedLables(trainList2, "train2" + ms2regressor.name + transform.name + p.myMsDataFile.Name);
                        ms2regressor.writePredictedLables(testList2, "test2" + ms2regressor.name + transform.name + p.myMsDataFile.Name);
                    }
                }
            }
            catch (ArgumentException e)
            {
                p.po.output("Could not calibrate: " + e.Message);
            }

            CalibrationFunction bestCf = new SeparateCalibrationFunction(bestMS1predictor, bestMS2predictor);

            p.po.status("Calibrating Spectra");

            CalibrateSpectra(p, bestCf);

            return bestCf;
        }

        private static void CalibrateSpectra(SoftwareLockMassParams p, CalibrationFunction bestCf)
        {
            foreach (var a in p.myMsDataFile)
            {
                if (a.MsnOrder == 2)
                {
                    if (p.MS2spectraToWatch.Contains(a.OneBasedScanNumber))
                    {
                        p.po.output("Calibrating scan number " + a.OneBasedScanNumber);
                        p.po.output(" before calibration:");
                        p.po.output(" " + string.Join(",", a.MassSpectrum.newSpectrumExtract(p.mzRange).xArray));
                    }

                    int oneBasedScanNumber;
                    a.TryGetPrecursorOneBasedScanNumber(out oneBasedScanNumber);
                    var precursorScan = p.myMsDataFile.GetOneBasedScan(oneBasedScanNumber);

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

                    if (p.MS2spectraToWatch.Contains(a.OneBasedScanNumber))
                    {
                        p.po.output(" after calibration:");
                        p.po.output(" precursorMZ:" + precursorMZ);
                        p.po.output(" monoisotopicMZ:" + monoisotopicMZ);
                        p.po.output(" newSelectedMZ:" + newSelectedMZ);
                        p.po.output(" newMonoisotopicMZ:" + newMonoisotopicMZ);
                        p.po.output(" " + string.Join(",", a.MassSpectrum.newSpectrumExtract(p.mzRange).xArray));
                    }
                }
                else
                {
                    if (p.MS1spectraToWatch.Contains(a.OneBasedScanNumber))
                    {
                        p.po.output("Calibrating scan number " + a.OneBasedScanNumber);
                        p.po.output(" before calibration:");
                        p.po.output(" " + string.Join(",", a.MassSpectrum.newSpectrumExtract(p.mzRange).xArray));
                    }
                    Func<MzPeak, double> theFUnc = x => x.MZ - bestCf.Predict(new double[6] { 1, x.MZ, a.RetentionTime, x.Intensity, a.TotalIonCurrent, a.InjectionTime });
                    a.tranformByApplyingFunctionsToSpectraAndReplacingPrecursorMZs(theFUnc, double.NaN, double.NaN); if (p.MS1spectraToWatch.Contains(a.OneBasedScanNumber))
                    {
                        p.po.output(" after calibration:");
                        p.po.output(string.Join(",", a.MassSpectrum.newSpectrumExtract(p.mzRange).xArray));
                    }
                }
            }
        }

        public static void WriteDataToFiles(IEnumerable<LabeledDataPoint> trainingPoints, string prefix)
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
    }
}