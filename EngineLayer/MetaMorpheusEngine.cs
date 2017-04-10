using MassSpectrometry;
using System;
using System.Collections.Generic;
using System.Diagnostics;

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

        public static IEnumerable<LocalMS2Scan> GetMs2Scans(IMsDataFile<IMsDataScan<IMzSpectrum<IMzPeak>>> myMSDataFile)
        {
            foreach (var heh in myMSDataFile)
            {
                var ms2scan = heh as IMsDataScanWithPrecursor<IMzSpectrum<IMzPeak>>;
                if (ms2scan != null)
                {
                    var precursorSpectrum = myMSDataFile.GetOneBasedScan(ms2scan.OneBasedPrecursorScanNumber);
                    if (!ms2scan.SelectedIonGuessChargeStateGuess.HasValue)
                        ms2scan.RecomputeChargeState(precursorSpectrum.MassSpectrum, 0.01, 10);

                    if (!ms2scan.SelectedIonGuessIntensity.HasValue && !ms2scan.SelectedIonGuessMZ.HasValue)
                        ms2scan.RecomputeSelectedPeak(precursorSpectrum.MassSpectrum);

                    if (!ms2scan.SelectedIonGuessIntensity.HasValue)
                        ms2scan.ComputeSelectedPeakIntensity(precursorSpectrum.MassSpectrum);

                    if (!ms2scan.SelectedIonGuessMonoisotopicIntensity.HasValue && !ms2scan.SelectedIonGuessMonoisotopicMZ.HasValue)
                        ms2scan.RecomputeMonoisotopicPeak(precursorSpectrum.MassSpectrum, 0.01, 0.3);

                    if (!ms2scan.SelectedIonGuessMonoisotopicIntensity.HasValue)
                        ms2scan.ComputeMonoisotopicPeakIntensity(precursorSpectrum.MassSpectrum);

                    yield return new LocalMS2Scan(ms2scan);
                }
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