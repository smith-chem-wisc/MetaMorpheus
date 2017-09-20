using MassSpectrometry;
using MathNet.Numerics.Statistics;
using SharpLearning.Common.Interfaces;
using System.Linq;
using System.Text;
using System;

namespace EngineLayer.Calibration
{
    public class CalibrationResults : MetaMorpheusEngineResults
    {
        //public readonly List<Tuple<double, double>> ms1meanSds;
        //public readonly List<Tuple<double, double>> ms2meanSds;

        #region Public Fields

        //private readonly List<int> numMs1MassChargeCombinationsConsideredList;
        //private readonly List<int> numMs1MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaksList;
        //private readonly List<int> numMs2MassChargeCombinationsConsideredList;
        //private readonly List<int> numMs2MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaksList;
        //private readonly List<int> countList;
        //private readonly List<IPredictorModel<double>> ms1calibrationFunctions;
        //private readonly List<IPredictorModel<double>> ms2calibrationFunctions;
        public readonly Tuple<double, double> ms1Info;
        public readonly Tuple<double, double> ms2Info;

        #endregion Public Fields

        #region Public Constructors
        
        public CalibrationResults(Tuple<double, double> ms1Info, Tuple<double, double> ms2Info, CalibrationEngine s) : base(s)
        {
            this.ms1Info = ms1Info;
            this.ms2Info = ms2Info;
        }

        #endregion Public Constructors

        //public CalibrationResults(IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> myMSDataFile, CalibrationEngine s) : base(s)
        //{
        //    this.MyMSDataFile = myMSDataFile;
        //    ms1calibrationFunctions = new List<IPredictorModel<double>>();
        //    ms2calibrationFunctions = new List<IPredictorModel<double>>();
        //    ms1meanSds = new List<Tuple<double, double>>();
        //    ms2meanSds = new List<Tuple<double, double>>();
        //    numMs1MassChargeCombinationsConsideredList = new List<int>();
        //    numMs1MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaksList = new List<int>();
        //    numMs2MassChargeCombinationsConsideredList = new List<int>();
        //    numMs2MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaksList = new List<int>();
        //    countList = new List<int>();
        //}

        #region Public Properties

        public IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> MyMSDataFile { get; private set; }

        #endregion Public Properties

        #region Public Methods

        public override string ToString()
        {
            var sb = new StringBuilder();
            sb.AppendLine(base.ToString());
            return sb.ToString();
        }

        #endregion Public Methods

        //#region Internal Methods

        //internal void Add(DataPointAquisitionResults res)
        //{
        //    numMs1MassChargeCombinationsConsideredList.Add(res.numMs1MassChargeCombinationsConsidered);
        //    numMs1MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaksList.Add(res.numMs1MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks);
        //    numMs2MassChargeCombinationsConsideredList.Add(res.numMs2MassChargeCombinationsConsidered);
        //    numMs2MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaksList.Add(res.numMs2MassChargeCombinationsThatAreIgnoredBecauseOfTooManyPeaks);
        //    countList.Add(res.Count);
        //    ms1meanSds.Add(res.Ms1List.Select(b => b.label).MeanStandardDeviation());
        //    ms2meanSds.Add(res.Ms2List.Select(b => b.label).MeanStandardDeviation());
        //}

        //internal void Add(IPredictorModel<double> ms1calib, IPredictorModel<double> ms2calib)
        //{
        //    ms1calibrationFunctions.Add(ms1calib);
        //    ms2calibrationFunctions.Add(ms2calib);
        //}

        //#endregion Internal Methods
    }
}