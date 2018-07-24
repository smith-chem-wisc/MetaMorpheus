using MassSpectrometry;
using System;
using System.Collections.Generic;

namespace EngineLayer.Calibration
{
    public class CalibrationEngine : MetaMorpheusEngine
    {
        private const double FractionOfFileUsedForSmoothing = 0.05;
        private readonly MsDataFile MyMsDataFile;
        private readonly DataPointAquisitionResults Datapoints;
        private int ScansUsedForSmoothing;

        public CalibrationEngine(MsDataFile myMSDataFile, DataPointAquisitionResults datapoints, CommonParameters commonParameters, List<string> nestedIds) : base(commonParameters, nestedIds)
        {
            MyMsDataFile = myMSDataFile;
            Datapoints = datapoints;
            ScansUsedForSmoothing = (int)Math.Round(FractionOfFileUsedForSmoothing * MyMsDataFile.NumSpectra);
        }

        protected override MetaMorpheusEngineResults RunSpecific()
        {
            Status("Generating MS1 calibration function");
            //var ms1Model = GetRandomForestModel(myMs1DataPoints, ms1fracForTraining);
            //var ms1Model = GetGradientBoostModel(myMs1DataPoints, ms1fracForTraining);

            Status("Generating MS2 calibration function");
            //var ms2Model = GetRandomForestModel(myMs2DataPoints, ms2fracForTraining);
            //var ms2Model = GetGradientBoostModel(myMs2DataPoints, ms2fracForTraining);

            Status("Calibrating spectra");

            CalibrateSpectra(Datapoints.Ms1List, 1);
            CalibrateSpectra(Datapoints.Ms1List, 2);

            return new MetaMorpheusEngineResults(this);
        }
    }
}