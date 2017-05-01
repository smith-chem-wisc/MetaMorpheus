using Chemistry;
using MassSpectrometry;
using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;

namespace EngineLayer
{
    public abstract class MetaMorpheusEngine
    {

        #region Public Events

        public static event EventHandler<SingleEngineEventArgs> StartingSingleEngineHander;

        public static event EventHandler<SingleEngineFinishedEventArgs> FinishedSingleEngineHandler;

        public static event EventHandler<StringEventArgs> OutLabelStatusHandler;

        public static event EventHandler<StringEventArgs> WarnHandler;

        public static event EventHandler<ProgressEventArgs> OutProgressHandler;

        #endregion Public Events

        #region Public Methods

        public static IEnumerable<Ms2ScanWithSpecificMass> GetMs2Scans(IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> myMSDataFile)
        {
            foreach (var ms2scan in myMSDataFile.OfType<IMsDataScanWithPrecursor<IMzSpectrum<IMzPeak>>>())
            {
                var precursorSpectrum = myMSDataFile.GetOneBasedScan(ms2scan.OneBasedPrecursorScanNumber);

                ms2scan.RefineSelectedMzAndIntensity(precursorSpectrum.MassSpectrum);

                if (ms2scan.SelectedIonMonoisotopicGuessMz.HasValue)
                    ms2scan.ComputeMonoisotopicPeakIntensity(precursorSpectrum.MassSpectrum);

                // By now ms2scan may or may not have SelectedIonMonoisotopicGuessMz

                int maxAssumedChargeState = 10;
                var massTolerance = new Tolerance("10 PPM");
                int intensityRatio = 4;
                var isolatedStuff = ms2scan.GetIsolatedMassesAndCharges(precursorSpectrum.MassSpectrum, maxAssumedChargeState, massTolerance, intensityRatio).ToList();

                if (ms2scan.SelectedIonChargeStateGuess.HasValue)
                {
                    var PrecursorCharge = ms2scan.SelectedIonChargeStateGuess.Value;
                    if (ms2scan.SelectedIonMonoisotopicGuessMz.HasValue)
                    {
                        var PrecursorMZ = ms2scan.SelectedIonMonoisotopicGuessMz.Value;
                        if (!isolatedStuff.Any(b => massTolerance.Within(PrecursorMZ.ToMass(PrecursorCharge), b.Item1.First().ToMass(b.Item2))))
                        {
                            isolatedStuff.Add(new Tuple<List<double>, int>(new List<double> { PrecursorMZ }, PrecursorCharge));
                        }
                    }
                    else
                    {
                        var PrecursorMZ = ms2scan.SelectedIonMZ;
                        if (!isolatedStuff.Any(b => massTolerance.Within(PrecursorMZ.ToMass(PrecursorCharge), b.Item1.First().ToMass(b.Item2))))
                        {
                            isolatedStuff.Add(new Tuple<List<double>, int>(new List<double> { PrecursorMZ }, PrecursorCharge));
                        }
                    }
                }
                foreach (var heh in isolatedStuff)
                    yield return new Ms2ScanWithSpecificMass(ms2scan, heh.Item1.First(), heh.Item2);
            }
        }

        public MetaMorpheusEngineResults Run()
        {
            StartingSingleEngine();
            var stopWatch = new Stopwatch();
            stopWatch.Start();
            var myResults = RunSpecific();
            stopWatch.Stop();
            myResults.Time = stopWatch.Elapsed;
            FinishedSingleEngine(myResults);
            return myResults;
        }

        #endregion Public Methods

        #region Protected Methods

        protected void Warn(string v, List<string> nestedIds)
        {
            WarnHandler?.Invoke(this, new StringEventArgs(v, nestedIds));
        }

        protected void Status(string v, List<string> nestedIds)
        {
            OutLabelStatusHandler?.Invoke(this, new StringEventArgs(v, nestedIds));
        }

        protected void ReportProgress(ProgressEventArgs v)
        {
            OutProgressHandler?.Invoke(this, v);
        }

        protected abstract MetaMorpheusEngineResults RunSpecific();

        #endregion Protected Methods

        #region Private Methods

        private void StartingSingleEngine()
        {
            StartingSingleEngineHander?.Invoke(this, new SingleEngineEventArgs(this));
        }

        private void FinishedSingleEngine(MetaMorpheusEngineResults myResults)
        {
            FinishedSingleEngineHandler?.Invoke(this, new SingleEngineFinishedEventArgs(myResults));
        }

        #endregion Private Methods

    }
}