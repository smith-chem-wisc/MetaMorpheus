using Chemistry;
using MassSpectrometry;
using MzLibUtil;
using Proteomics.Fragmentation;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Threading.Tasks;

namespace EngineLayer
{
    public abstract class MetaMorpheusEngine
    {
        protected static readonly Dictionary<DissociationType, double> complementaryIonConversionDictionary = new Dictionary<DissociationType, double>
        {
            { DissociationType.HCD, Constants.ProtonMass },
            { DissociationType.ETD, 2 * Constants.ProtonMass }, //presence of zplusone (zdot) makes this two instead of one
            { DissociationType.CID, Constants.ProtonMass },
            //TODO: refactor such that complementary ions are generated specifically for their complementary pair.
            //TODO: create a method to auto-determine the conversion
        };

        protected readonly CommonParameters CommonParameters;

        protected readonly List<string> NestedIds;

        protected MetaMorpheusEngine(CommonParameters commonParameters, List<string> nestedIds)
        {
            CommonParameters = commonParameters;
            NestedIds = nestedIds;
        }

        public static event EventHandler<SingleEngineEventArgs> StartingSingleEngineHander;

        public static event EventHandler<SingleEngineFinishedEventArgs> FinishedSingleEngineHandler;

        public static event EventHandler<StringEventArgs> OutLabelStatusHandler;

        public static event EventHandler<StringEventArgs> WarnHandler;

        public static event EventHandler<ProgressEventArgs> OutProgressHandler;

        public static double CalculatePeptideScore(MsDataScan thisScan, List<MatchedFragmentIon> matchedFragmentIons)
        {
            double score = 0;

            if (thisScan.MassSpectrum.XcorrProcessed)
            {
                foreach (var fragment in matchedFragmentIons)
                {
                    switch (fragment.NeutralTheoreticalProduct.ProductType)
                    {
                        case ProductType.aDegree:
                        case ProductType.aStar:
                        case ProductType.bDegree:
                        case ProductType.bStar:
                        case ProductType.yDegree:
                        case ProductType.yStar:
                            score += 0.01 * fragment.Intensity;
                            break;

                        default:
                            score += 1 * fragment.Intensity;
                            break;
                    }
                }
            }
            else
            {
                foreach (var fragment in matchedFragmentIons)
                {
                    double fragmentScore = 1 + (fragment.Intensity / thisScan.TotalIonCurrent);
                    score += fragmentScore;
                }
            }

            return score;
        }

        public static List<MatchedFragmentIon> MatchFragmentIons(Ms2ScanWithSpecificMass scan, List<Product> theoreticalProducts, CommonParameters commonParameters)
        {
            var matchedFragmentIons = new List<MatchedFragmentIon>();

            // if the spectrum has no peaks
            if (!scan.ExperimentalFragments.Any())
            {
                return matchedFragmentIons;
            }

            if (scan.TheScan.MassSpectrum.XcorrProcessed)
            {
                foreach (Product product in theoreticalProducts)
                {
                    // unknown fragment mass; this only happens rarely for sequences with unknown amino acids
                    if (double.IsNaN(product.NeutralMass))
                    {
                        continue;
                    }

                    double theoreticalFragmentMz = Math.Round(product.NeutralMass.ToMz(1) / 1.0005079, 0) * 1.0005079;
                    var closestMzIndex = scan.TheScan.MassSpectrum.GetClosestPeakIndex(theoreticalFragmentMz).Value;

                    if (commonParameters.ProductMassTolerance.Within(scan.TheScan.MassSpectrum.XArray[closestMzIndex], theoreticalFragmentMz))
                    {
                        matchedFragmentIons.Add(new MatchedFragmentIon(product, theoreticalFragmentMz, scan.TheScan.MassSpectrum.YArray[closestMzIndex], 1));
                    }
                }

                return matchedFragmentIons;
            }

            // search for ions in the spectrum
            foreach (Product product in theoreticalProducts)
            {
                // unknown fragment mass; this only happens rarely for sequences with unknown amino acids
                if (double.IsNaN(product.NeutralMass))
                {
                    continue;
                }

                // get the closest peak in the spectrum to the theoretical peak
                var closestExperimentalMass = scan.GetClosestExperimentalFragmentMass(product.NeutralMass);

                // is the mass error acceptable?
                if (commonParameters.ProductMassTolerance.Within(closestExperimentalMass.monoisotopicMass, product.NeutralMass) && closestExperimentalMass.charge <= scan.PrecursorCharge)
                {
                    matchedFragmentIons.Add(new MatchedFragmentIon(product, closestExperimentalMass.monoisotopicMass.ToMz(closestExperimentalMass.charge),
                        closestExperimentalMass.peaks.First().intensity, closestExperimentalMass.charge));
                }
            }
            if (commonParameters.AddCompIons)
            {
                double protonMassShift = complementaryIonConversionDictionary[commonParameters.DissociationType].ToMass(1);

                foreach (Product product in theoreticalProducts)
                {
                    // unknown fragment mass or diagnostic ion or precursor; skip those
                    if (double.IsNaN(product.NeutralMass) || product.ProductType == ProductType.D || product.ProductType == ProductType.M)
                    {
                        continue;
                    }

                    double compIonMass = scan.PrecursorMass + protonMassShift - product.NeutralMass;

                    // get the closest peak in the spectrum to the theoretical peak
                    var closestExperimentalMass = scan.GetClosestExperimentalFragmentMass(compIonMass);

                    // is the mass error acceptable?
                    if (commonParameters.ProductMassTolerance.Within(closestExperimentalMass.monoisotopicMass, compIonMass) && closestExperimentalMass.charge <= scan.PrecursorCharge)
                    {
                        matchedFragmentIons.Add(new MatchedFragmentIon(product, closestExperimentalMass.monoisotopicMass.ToMz(closestExperimentalMass.charge),
                            closestExperimentalMass.totalIntensity, closestExperimentalMass.charge));
                    }
                }
            }

            return matchedFragmentIons;
        }


        public static IEnumerable<Ms2ScanWithSpecificMass> GetMs2Scans(MsDataFile myMSDataFile, string fullFilePath, CommonParameters commonParameters)
        {
            var msNScans = myMSDataFile.GetAllScansList().Where(x => x.MsnOrder > 1).ToArray();
            var ms2Scans = msNScans.Where(p => p.MsnOrder == 2).ToArray();
            var ms3Scans = msNScans.Where(p => p.MsnOrder == 3).ToArray();
            List<Ms2ScanWithSpecificMass>[] scansWithPrecursors = new List<Ms2ScanWithSpecificMass>[ms2Scans.Length];

            if (!ms2Scans.Any())
            {
                return new List<Ms2ScanWithSpecificMass>();
            }

            Parallel.ForEach(Partitioner.Create(0, ms2Scans.Length), new ParallelOptions { MaxDegreeOfParallelism = commonParameters.MaxThreadsToUsePerFile },
                (partitionRange, loopState) =>
                {
                    for (int i = partitionRange.Item1; i < partitionRange.Item2; i++)
                    {
                        if (GlobalVariables.StopLoops) { break; }

                        MsDataScan ms2scan = ms2Scans[i];

                        List<(double, int)> precursors = new List<(double, int)>();
                        if (ms2scan.OneBasedPrecursorScanNumber.HasValue)
                        {
                            var precursorSpectrum = myMSDataFile.GetOneBasedScan(ms2scan.OneBasedPrecursorScanNumber.Value);

                            try
                            {
                                ms2scan.RefineSelectedMzAndIntensity(precursorSpectrum.MassSpectrum);
                            }
                            catch (MzLibException ex)
                            {
                                Warn("Could not get precursor ion for MS2 scan #" + ms2scan.OneBasedScanNumber + "; " + ex.Message);
                                continue;
                            }

                            if (ms2scan.SelectedIonMonoisotopicGuessMz.HasValue)
                            {
                                ms2scan.ComputeMonoisotopicPeakIntensity(precursorSpectrum.MassSpectrum);
                            }

                            if (commonParameters.DoPrecursorDeconvolution)
                            {
                                foreach (var envelope in ms2scan.GetIsolatedMassesAndCharges(
                                    precursorSpectrum.MassSpectrum, 1,
                                    commonParameters.DeconvolutionMaxAssumedChargeState,
                                    commonParameters.DeconvolutionMassTolerance.Value,
                                    commonParameters.DeconvolutionIntensityRatio))
                                {
                                    var monoPeakMz = envelope.monoisotopicMass.ToMz(envelope.charge);
                                    precursors.Add((monoPeakMz, envelope.charge));
                                }
                            }
                        }

                        if (commonParameters.UseProvidedPrecursorInfo && ms2scan.SelectedIonChargeStateGuess.HasValue)
                        {
                            var precursorCharge = ms2scan.SelectedIonChargeStateGuess.Value;
                            if (ms2scan.SelectedIonMonoisotopicGuessMz.HasValue)
                            {
                                double precursorMZ = ms2scan.SelectedIonMonoisotopicGuessMz.Value;
                                if (!precursors.Any(b =>
                                    commonParameters.DeconvolutionMassTolerance.Within(
                                        precursorMZ.ToMass(precursorCharge), b.Item1.ToMass(b.Item2))))
                                {
                                    precursors.Add((precursorMZ, precursorCharge));
                                }
                            }
                            else
                            {
                                double precursorMZ = ms2scan.SelectedIonMZ.Value;
                                if (!precursors.Any(b =>
                                    commonParameters.DeconvolutionMassTolerance.Within(
                                        precursorMZ.ToMass(precursorCharge), b.Item1.ToMass(b.Item2))))
                                {
                                    precursors.Add((precursorMZ, precursorCharge));
                                }
                            }
                        }

                        scansWithPrecursors[i] = new List<Ms2ScanWithSpecificMass>();
                        IsotopicEnvelope[] neutralExperimentalFragments = null;

                        if (commonParameters.DissociationType != DissociationType.LowCID)
                        {
                            neutralExperimentalFragments = Ms2ScanWithSpecificMass.GetNeutralExperimentalFragments(ms2scan, commonParameters);
                        }

                        // get child scans
                        List<MsDataScan> ms2ChildScans = new List<MsDataScan>();
                        List<MsDataScan> ms3ChildScans = new List<MsDataScan>();
                        if (commonParameters.ChildScanDissociationType != DissociationType.Unknown)
                        {
                            ms3ChildScans = ms3Scans.Where(p => p.OneBasedPrecursorScanNumber == ms2scan.OneBasedScanNumber).ToList();

                            ms2ChildScans = ms2Scans.Where(p => p.OneBasedPrecursorScanNumber == ms2scan.OneBasedPrecursorScanNumber
                                && p.OneBasedScanNumber > ms2scan.OneBasedScanNumber
                                && Math.Abs(p.IsolationMz.Value - ms2scan.IsolationMz.Value) < 0.01).ToList();
                        }

                        foreach (var precursor in precursors)
                        {
                            var scan = new Ms2ScanWithSpecificMass(ms2scan, precursor.Item1,
                                precursor.Item2, fullFilePath, commonParameters, neutralExperimentalFragments);

                            foreach (var ms2ChildScan in ms2ChildScans)
                            {
                                IsotopicEnvelope[] childNeutralExperimentalFragments = null;

                                if (commonParameters.ChildScanDissociationType != DissociationType.LowCID)
                                {
                                    childNeutralExperimentalFragments = Ms2ScanWithSpecificMass.GetNeutralExperimentalFragments(ms2ChildScan, commonParameters);
                                }

                                scan.ChildScans.Add(new Ms2ScanWithSpecificMass(ms2ChildScan, precursor.Item1,
                                    precursor.Item2, fullFilePath, commonParameters, childNeutralExperimentalFragments));
                            }

                            foreach (var ms3ChildScan in ms3ChildScans)
                            {
                                int precursorCharge = 1;
                                var precursorSpectrum = ms2scan;

                                try
                                {
                                    ms3ChildScan.RefineSelectedMzAndIntensity(precursorSpectrum.MassSpectrum);
                                }
                                catch (MzLibException ex)
                                {
                                    Warn("Could not get precursor ion for MS3 scan #" + ms3ChildScan.OneBasedScanNumber + "; " + ex.Message);
                                    continue;
                                }

                                if (ms3ChildScan.SelectedIonMonoisotopicGuessMz.HasValue)
                                {
                                    ms3ChildScan.ComputeMonoisotopicPeakIntensity(precursorSpectrum.MassSpectrum);
                                }

                                if (ms3ChildScan.SelectedIonChargeStateGuess.HasValue)
                                {
                                    precursorCharge = ms3ChildScan.SelectedIonChargeStateGuess.Value;
                                }
                                if (!ms3ChildScan.SelectedIonMonoisotopicGuessMz.HasValue)
                                {
                                    Warn("Could not get precursor ion m/z for MS3 scan #" + ms3ChildScan.OneBasedScanNumber);
                                    continue;
                                }

                                IsotopicEnvelope[] childNeutralExperimentalFragments = null;

                                if (commonParameters.ChildScanDissociationType != DissociationType.LowCID)
                                {
                                    childNeutralExperimentalFragments = Ms2ScanWithSpecificMass.GetNeutralExperimentalFragments(ms3ChildScan, commonParameters);
                                }

                                scan.ChildScans.Add(new Ms2ScanWithSpecificMass(ms3ChildScan, ms3ChildScan.SelectedIonMonoisotopicGuessMz.Value,
                                    ms3ChildScan.SelectedIonChargeStateGuess.Value, fullFilePath, commonParameters, childNeutralExperimentalFragments));
                            }

                            scansWithPrecursors[i].Add(scan);
                        }
                    }
                });

            var childScanNumbers = new HashSet<int>(scansWithPrecursors.SelectMany(p => p.SelectMany(v => v.ChildScans.Select(x => x.OneBasedScanNumber))));
            var parentScans = scansWithPrecursors.Where(p => p.Any() && !childScanNumbers.Contains(p.First().OneBasedScanNumber)).SelectMany(v => v);

            // XCorr pre-processing for low-res data. this is here because the parent/child scans may have different
            // resolutions, so this pre-processing must take place after the parent/child scans have been determined
            foreach (var parentScan in parentScans)
            {
                if (commonParameters.DissociationType == DissociationType.LowCID && !parentScan.TheScan.MassSpectrum.XcorrProcessed)
                {
                    parentScan.TheScan.MassSpectrum.XCorrPrePreprocessing(0, 1969, parentScan.TheScan.IsolationMz.Value);
                }

                foreach (var childScan in parentScan.ChildScans)
                {
                    if (commonParameters.ChildScanDissociationType == DissociationType.LowCID && !childScan.TheScan.MassSpectrum.XcorrProcessed)
                    {
                        childScan.TheScan.MassSpectrum.XCorrPrePreprocessing(0, 1969, childScan.TheScan.IsolationMz.Value);
                    }
                }
            }

            return parentScans;
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

        public string GetId()
        {
            return string.Join(",", NestedIds);
        }

        protected static void Warn(string v)
        {
            WarnHandler?.Invoke(null, new StringEventArgs(v, null));
        }

        protected void Status(string v)
        {
            OutLabelStatusHandler?.Invoke(this, new StringEventArgs(v, NestedIds));
        }

        protected void ReportProgress(ProgressEventArgs v)
        {
            OutProgressHandler?.Invoke(this, v);
        }

        protected abstract MetaMorpheusEngineResults RunSpecific();

        private void StartingSingleEngine()
        {
            StartingSingleEngineHander?.Invoke(this, new SingleEngineEventArgs(this));
        }

        private void FinishedSingleEngine(MetaMorpheusEngineResults myResults)
        {
            FinishedSingleEngineHandler?.Invoke(this, new SingleEngineFinishedEventArgs(myResults));
        }
    }
}