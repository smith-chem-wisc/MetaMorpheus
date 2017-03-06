using MassSpectrometry;
using MzLibUtil;
using Proteomics;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Net;
using System.Reflection;

namespace EngineLayer
{
    public abstract class MyEngine
    {

        #region Private Fields

        private static readonly string elementsLocation = Path.Combine(@"Data", @"elements.dat");
        private static readonly string unimodLocation = Path.Combine(@"Data", @"unimod.xml");
        private static readonly string uniprotLocation = Path.Combine(@"Mods", @"ptmlist.txt");

        #endregion Private Fields

        #region Public Constructors

        static MyEngine()
        {
            try
            {
                UsefulProteomicsDatabases.Loaders.LoadElements(elementsLocation);
                UnimodDeserialized = UsefulProteomicsDatabases.Loaders.LoadUnimod(unimodLocation).ToList();
                UniprotDeseralized = UsefulProteomicsDatabases.Loaders.LoadUniprot(uniprotLocation).ToList();
            }
            catch (WebException)
            {
            }

            MetaMorpheusVersion = Assembly.GetExecutingAssembly().GetName().Version.ToString();
            SearchModesKnown = LoadSearchModesFromFile().ToList();
        }

        #endregion Public Constructors

        #region Public Events

        public static event EventHandler<SingleEngineEventArgs> StartingSingleEngineHander;

        public static event EventHandler<SingleEngineFinishedEventArgs> FinishedSingleEngineHandler;

        public static event EventHandler<StringEventArgs> OutLabelStatusHandler;

        public static event EventHandler<StringEventArgs> WarnHandler;

        public static event EventHandler<ProgressEventArgs> OutProgressHandler;

        #endregion Public Events

        #region Public Properties

        public static string MetaMorpheusVersion { get; private set; }
        public static IEnumerable<Modification> UnimodDeserialized { get; private set; }
        public static IEnumerable<Modification> UniprotDeseralized { get; private set; }

        public static List<SearchMode> SearchModesKnown { get; private set; }

        #endregion Public Properties

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

        public MyResults Run()
        {
            startingSingleEngine();
            var stopWatch = new Stopwatch();
            stopWatch.Start();
            var myResults = RunSpecific();
            stopWatch.Stop();
            myResults.Time = stopWatch.Elapsed;
            finishedSingleEngine(myResults);
            return myResults;
        }

        #endregion Public Methods

        #region Protected Methods

        protected void Warn(string v)
        {
            WarnHandler?.Invoke(this, new StringEventArgs(v));
        }

        protected void Status(string v)
        {
            OutLabelStatusHandler?.Invoke(this, new StringEventArgs(v));
        }

        protected void ReportProgress(ProgressEventArgs v)
        {
            OutProgressHandler?.Invoke(this, v);
        }

        protected abstract MyResults RunSpecific();

        #endregion Protected Methods

        #region Private Methods

        private static IEnumerable<SearchMode> LoadSearchModesFromFile()
        {
            yield return new SinglePpmAroundZeroSearchMode(10);
            yield return new SinglePpmAroundZeroSearchMode(5);
            yield return new SinglePpmAroundZeroSearchMode(20);
            yield return new SingleAbsoluteAroundZeroSearchMode(0.05);
            yield return new DotSearchMode("3mm", new double[] { 0, 1.003, 2.006, 3.009 }, new Tolerance(ToleranceUnit.PPM, 5));
            yield return new IntervalSearchMode(new List<DoubleRange>() { new DoubleRange(-2.1, 2.1) });
            yield return new OpenSearchMode();
            yield return new IntervalSearchMode(new List<DoubleRange> { new DoubleRange(-0.005, 0.005), new DoubleRange(21.981943 - 0.005, 21.981943 + 0.005) });
            yield return new IntervalSearchMode(new List<DoubleRange> { new DoubleRange(-187, double.PositiveInfinity) });
            foreach (var sm in GetResidueInclusionExclusionSearchModes(new DoubleRange(-187, double.PositiveInfinity), 0.0075))
                yield return sm;
        }

        /// <summary>
        /// Ideally v is less than 0.00168565165, so no overlaps happen
        /// </summary>
        /// <param name="doubleRange"></param>
        /// <param name="v"></param>
        /// <returns></returns>
        private static IEnumerable<SearchMode> GetResidueInclusionExclusionSearchModes(DoubleRange doubleRange, double v)
        {
            List<double> massesToExclude = new List<double>();
            for (char c = 'A'; c <= 'Z'; c++)
            {
                Residue residue;
                if (Residue.TryGetResidue(c, out residue))
                {
                    massesToExclude.Add(residue.MonoisotopicMass);
                    massesToExclude.Add(-residue.MonoisotopicMass);
                    for (char cc = 'A'; cc <= 'Z'; cc++)
                    {
                        Residue residueCC;
                        if (Residue.TryGetResidue(cc, out residueCC))
                        {
                            massesToExclude.Add(residue.MonoisotopicMass + residueCC.MonoisotopicMass);
                            massesToExclude.Add(residue.MonoisotopicMass - residueCC.MonoisotopicMass);
                            massesToExclude.Add(-residue.MonoisotopicMass + residueCC.MonoisotopicMass);
                            massesToExclude.Add(-residue.MonoisotopicMass - residueCC.MonoisotopicMass);
                        }
                    }
                }
            }
            List<double> filteredMasses = massesToExclude.GroupBy(b => Math.Round(b, 6)).Select(b => b.FirstOrDefault()).OrderBy(b => b).ToList();

            yield return new DotSearchMode("Only AAs", filteredMasses, new Tolerance(ToleranceUnit.Absolute, v));

            List<DoubleRange> doubleRanges = new List<DoubleRange>();

            var prevGoodMin = double.NegativeInfinity;

            for (int i = 0; i < filteredMasses.Count; i++)
            {
                if (Math.Round(filteredMasses[i], 6) == 0)
                    continue;
                doubleRanges.Add(new DoubleRange(prevGoodMin, filteredMasses[i] - v));
                prevGoodMin = filteredMasses[i] + v;
            }
            doubleRanges.Add(new DoubleRange(prevGoodMin, double.PositiveInfinity));

            doubleRanges = doubleRanges.Where(b => b.Minimum <= doubleRange.Maximum && b.Maximum >= doubleRange.Minimum).Select(b => new DoubleRange(Math.Max(doubleRange.Minimum, b.Minimum), Math.Min(doubleRange.Maximum, b.Maximum))).ToList();

            yield return new IntervalSearchMode("Exclude AAs", doubleRanges);
        }

        private void startingSingleEngine()
        {
            StartingSingleEngineHander?.Invoke(this, new SingleEngineEventArgs(this));
        }

        private void finishedSingleEngine(MyResults myResults)
        {
            FinishedSingleEngineHandler?.Invoke(this, new SingleEngineFinishedEventArgs(myResults));
        }

        #endregion Private Methods

    }
}