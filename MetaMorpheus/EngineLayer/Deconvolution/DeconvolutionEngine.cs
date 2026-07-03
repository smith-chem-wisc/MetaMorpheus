using EngineLayer.DIA;
using EngineLayer.Util;
using MassSpectrometry;
using MzLibUtil;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace EngineLayer.Deconvolution;

public class DeconvolutionEngineResults : MetaMorpheusEngineResults
{
    public DeconvolutionEngineResults(DeconvolutionEngine engine, IEnumerable<Ms2ScanWithSpecificMass> ms2Scans) : base(engine)
    {
        Ms2Scans = ms2Scans;
    }
    public IEnumerable<Ms2ScanWithSpecificMass> Ms2Scans { get; }
}


public class DeconvolutionEngine : MetaMorpheusEngine
{
    private readonly MsDataFile myMSDataFile;

    public DeconvolutionEngine(MsDataFile toDeconvolute, CommonParameters commonParameters, List<(string FileName, CommonParameters Parameters)> fileSpecificParameters, List<string> nestedIds) 
        : base(commonParameters, fileSpecificParameters, nestedIds)
    {
        myMSDataFile = toDeconvolute;
    }

    protected override MetaMorpheusEngineResults RunSpecific()
    {
        IEnumerable<Ms2ScanWithSpecificMass>? scans = GetMs2Scans(myMSDataFile, myMSDataFile.FilePath, CommonParameters);

        if (scans == null)
            throw new ArgumentNullException(nameof(scans), "Scans could not be deconvoluted.");

        return new DeconvolutionEngineResults(this, scans);
    }

    public static IEnumerable<Ms2ScanWithSpecificMass> GetMs2Scans(MsDataFile myMSDataFile, string fullFilePath, CommonParameters commonParameters)
    {
        if (commonParameters.DIAparameters != null)
        {
            return commonParameters.DIAparameters.AanalysisType switch
            {
                DIAanalysisType.DIA => new DIAEngine(myMSDataFile, commonParameters).GetPseudoMs2Scans(),
                DIAanalysisType.ISD => new ISDEngine(myMSDataFile, commonParameters).GetPseudoMs2Scans(),
                _ => throw new NotImplementedException("DIA analysis type not implemented."),
            };
        }

        var scansWithPrecursors = GetGroupedMs2Scans(myMSDataFile, fullFilePath, commonParameters);

        if (scansWithPrecursors.Length == 0)
        {
            return new List<Ms2ScanWithSpecificMass>();
        }

        var childScanNumbers = new HashSet<int>(scansWithPrecursors.SelectMany(p => p.SelectMany(v => v.ChildScans.Select(x => x.OneBasedScanNumber))));
        var parentScans = scansWithPrecursors.Where(p => p.Any() && !childScanNumbers.Contains(p.First().OneBasedScanNumber))
            .SelectMany(v => v)
            .OrderBy(p => p.OneBasedScanNumber)
            .ToArray();

        // XCorr pre-processing for low-res data. this is here because the parent/child scans may have different
        // resolutions, so this pre-processing must take place after the parent/child scans have been determined
        if (commonParameters.DissociationType == DissociationType.LowCID || commonParameters.MS2ChildScanDissociationType == DissociationType.LowCID || commonParameters.MS3ChildScanDissociationType == DissociationType.LowCID)
        {
            Parallel.ForEach(Partitioner.Create(0, parentScans.Length), new ParallelOptions { MaxDegreeOfParallelism = commonParameters.MaxThreadsToUsePerFile },
                (partitionRange, loopState) =>
                {
                    for (int i = partitionRange.Item1; i < partitionRange.Item2; i++)
                    {
                        if (GlobalVariables.StopLoops) { break; }

                        var parentScan = parentScans[i];

                        if (commonParameters.DissociationType == DissociationType.LowCID && !parentScan.TheScan.MassSpectrum.XcorrProcessed)
                        {
                            lock (parentScan.TheScan)
                            {
                                if (!parentScan.TheScan.MassSpectrum.XcorrProcessed)
                                {
                                    parentScan.TheScan.MassSpectrum.XCorrPrePreprocessing(0, 1969, parentScan.TheScan.IsolationMz.Value);
                                }
                            }
                        }

                        foreach (var childScan in parentScan.ChildScans)
                        {
                            if (((childScan.TheScan.MsnOrder == 2 && commonParameters.MS2ChildScanDissociationType == DissociationType.LowCID)
                                || (childScan.TheScan.MsnOrder == 3 && commonParameters.MS3ChildScanDissociationType == DissociationType.LowCID))
                            && !childScan.TheScan.MassSpectrum.XcorrProcessed)
                            {
                                lock (childScan.TheScan)
                                {
                                    if (!childScan.TheScan.MassSpectrum.XcorrProcessed)
                                    {
                                        childScan.TheScan.MassSpectrum.XCorrPrePreprocessing(0, 1969, childScan.TheScan.IsolationMz.Value);
                                    }
                                }
                            }
                        }
                    }
                });
            }

        return parentScans;
    }

    public static List<Ms2ScanWithSpecificMass> GetMs2ScansWrapByScanNum(MsDataFile myMSDataFile, string fullFilePath,
        CommonParameters commonParameters, out List<List<(double, int, double)>> precursors)
    {
        precursors = new List<List<(double, int, double)>>();

        if (commonParameters.DIAparameters != null)
        {
            var diaParentScans = GetMs2Scans(myMSDataFile, fullFilePath, commonParameters).OrderBy(p => p.OneBasedScanNumber).ToList();
            foreach (var scan in diaParentScans)
            {
                precursors.Add(new List<(double, int, double)>
                {
                    (scan.PrecursorMass, scan.PrecursorCharge, scan.PrecursorMonoisotopicPeakMz)
                });
            }

            return diaParentScans;
        }

        var scansWithPrecursors = GetGroupedMs2Scans(myMSDataFile, fullFilePath, commonParameters);
        var parentScansByGroup = GetParentScansByGroup(scansWithPrecursors);
        var parentScans = new List<Ms2ScanWithSpecificMass>();

        foreach (var group in parentScansByGroup)
        {
            parentScans.Add(group[0]);
            precursors.Add(group.Select(p => (p.PrecursorMass, p.PrecursorCharge, p.PrecursorMonoisotopicPeakMz)).ToList());
        }

        return parentScans;
    }

    public static List<Ms2ScanWithSpecificMass>[] GetGroupedMs2Scans(MsDataFile myMSDataFile, string fullFilePath, CommonParameters commonParameters)
    {
        var msNScans = myMSDataFile.GetAllScansList().Where(x => x.MsnOrder > 1).ToArray();
        var ms2Scans = msNScans.Where(p => p.MsnOrder == 2).ToArray();
        var ms3Scans = msNScans.Where(p => p.MsnOrder == 3).ToArray();
        List<Ms2ScanWithSpecificMass>[] scansWithPrecursors = new List<Ms2ScanWithSpecificMass>[ms2Scans.Length];

        if (!ms2Scans.Any())
        {
            return scansWithPrecursors;
        }

        Parallel.ForEach(Partitioner.Create(0, ms2Scans.Length), new ParallelOptions { MaxDegreeOfParallelism = commonParameters.MaxThreadsToUsePerFile },
            (partitionRange, loopState) =>
            {
                var precursorSet = new PrecursorSet(commonParameters.DeconvolutionMassTolerance);

                for (int i = partitionRange.Item1; i < partitionRange.Item2; i++)
                {
                    if (GlobalVariables.StopLoops) { break; }

                    precursorSet.Clear();
                    MsDataScan ms2scan = ms2Scans[i];

                    if (ms2scan.OneBasedPrecursorScanNumber.HasValue)
                    {
                        MsDataScan precursorSpectrum = myMSDataFile.GetOneBasedScan(ms2scan.OneBasedPrecursorScanNumber.Value);

                        if (precursorSpectrum is null)
                            goto PrecursorFromScanHeader;

                        try
                        {
                            ms2scan.RefineSelectedMzAndIntensity(precursorSpectrum.MassSpectrum);
                        }
                        catch (MzLibException ex)
                        {
                            Warn("Could not get precursor ion for MS2 scan #" + ms2scan.OneBasedScanNumber + "; " + ex.Message, null);
                            continue;
                        }

                        if (ms2scan.SelectedIonMonoisotopicGuessMz.HasValue)
                        {
                            ms2scan.ComputeMonoisotopicPeakIntensity(precursorSpectrum.MassSpectrum);
                        }

                        if (commonParameters.DoPrecursorDeconvolution)
                        {
                            foreach (IsotopicEnvelope envelope in ms2scan.GetIsolatedMassesAndCharges(
                                precursorSpectrum.MassSpectrum, commonParameters.PrecursorDeconvolutionParameters))
                            {
                                double? intensity = null;
                                if (commonParameters.UseMostAbundantPrecursorIntensity)
                                    intensity = envelope.Peaks.Max(p => p.intensity);

                                var fractionalIntensity = envelope.TotalIntensity /
                                      precursorSpectrum.MassSpectrum.YArray
                                      [
                                          precursorSpectrum.MassSpectrum.GetClosestPeakIndex(ms2scan.IsolationRange.Minimum)..precursorSpectrum.MassSpectrum.GetClosestPeakIndex(ms2scan.IsolationRange.Maximum)
                                      ].Sum();

                                // Method-agnostic envelope-quality score from mzLib (idempotent: caches on the
                                // envelope, so re-asking the same envelope is cheap).
                                double genericScore = envelope.GetOrComputeGenericScore(
                                    commonParameters.PrecursorDeconvolutionParameters);

                                precursorSet.Add(new Precursor(envelope, intensity, fractionalIntensity)
                                {
                                    DeconvolutionScore = genericScore
                                });
                            }
                        }
                    }

                    // If using precursor info from scan header and scan header has charge state.
                    // MsAlign uses this conditional to construct its precursors. 
                PrecursorFromScanHeader:
                    if (commonParameters.UseProvidedPrecursorInfo && ms2scan.SelectedIonChargeStateGuess.HasValue && ms2scan.SelectedIonChargeStateGuess != 0)
                    {
                        int precursorCharge = ms2scan.SelectedIonChargeStateGuess.Value;
                        double precursorIntensity = ms2scan.SelectedIonMonoisotopicGuessIntensity ?? ms2scan.SelectedIonIntensity ?? 1.0;
                        double precursorMz = ms2scan.SelectedIonMonoisotopicGuessMz ?? ms2scan.SelectedIonMZ.Value;

                        precursorSet.Add(new Precursor(precursorMz, precursorCharge, precursorIntensity, 1, null));
                    }

                    scansWithPrecursors[i] = new List<Ms2ScanWithSpecificMass>();
                    IsotopicEnvelope[] neutralExperimentalFragments = null;

                    if (commonParameters.DissociationType != DissociationType.LowCID)
                    {
                        neutralExperimentalFragments = Ms2ScanWithSpecificMass.GetNeutralExperimentalFragments(ms2scan, commonParameters);
                    }

                    // get child scans
                    List<MsDataScan> ms2ChildScans = null;
                    List<MsDataScan> ms3ChildScans = null;
                    if (commonParameters.MS2ChildScanDissociationType != DissociationType.Unknown
                        || commonParameters.MS3ChildScanDissociationType != DissociationType.Unknown)
                    {
                        ms3ChildScans = ms3Scans.Where(p => p.OneBasedPrecursorScanNumber == ms2scan.OneBasedScanNumber).ToList();

                        ms2ChildScans = ms2Scans.Where(p => p.OneBasedPrecursorScanNumber == ms2scan.OneBasedScanNumber ||
                        (p.OneBasedPrecursorScanNumber == ms2scan.OneBasedPrecursorScanNumber
                            && p.OneBasedScanNumber > ms2scan.OneBasedScanNumber
                            && Math.Abs(p.IsolationMz.Value - ms2scan.IsolationMz.Value) < 0.01)).ToList();
                    }

                    foreach (var precursor in precursorSet)
                    {
                        // assign precursor for this MS2 scan
                        var scan = new Ms2ScanWithSpecificMass(ms2scan, precursor.MonoisotopicPeakMz,
                            precursor.Charge, fullFilePath, commonParameters, neutralExperimentalFragments,
                            precursor.Intensity, precursor.EnvelopePeakCount, precursor.FractionalIntensity,
                            precursor.DeconvolutionScore);

                        // assign precursors for MS2 child scans
                        if (ms2ChildScans != null)
                        {
                            foreach (var ms2ChildScan in ms2ChildScans)
                            {
                                IsotopicEnvelope[] childNeutralExperimentalFragments = null;

                                if (commonParameters.MS2ChildScanDissociationType != DissociationType.LowCID)
                                {
                                    childNeutralExperimentalFragments = Ms2ScanWithSpecificMass.GetNeutralExperimentalFragments(ms2ChildScan, commonParameters);
                                }
                                var theChildScan = new Ms2ScanWithSpecificMass(ms2ChildScan, precursor.MonoisotopicPeakMz,
                                    precursor.Charge, fullFilePath, commonParameters, childNeutralExperimentalFragments,
                                    precursor.Intensity, precursor.EnvelopePeakCount, precursor.FractionalIntensity);
                                scan.ChildScans.Add(theChildScan);
                            }
                        }

                        // assign precursors for MS3 child scans
                        if (ms3ChildScans != null)
                        {
                            foreach (var ms3ChildScan in ms3ChildScans)
                            {
                                int precursorCharge = 1;
                                double precursorMz = 0;

                                if (ms3ChildScan.SelectedIonMonoisotopicGuessMz.HasValue)
                                {
                                    precursorMz = ms3ChildScan.SelectedIonMonoisotopicGuessMz.Value;
                                }
                                else if (ms3ChildScan.SelectedIonMZ.HasValue)
                                {
                                    precursorMz = ms3ChildScan.SelectedIonMZ.Value;
                                }

                                if (ms3ChildScan.SelectedIonChargeStateGuess.HasValue)
                                {
                                    precursorCharge = ms3ChildScan.SelectedIonChargeStateGuess.Value;
                                }

                                IsotopicEnvelope[] childNeutralExperimentalFragments = null;

                                if (commonParameters.MS3ChildScanDissociationType != DissociationType.LowCID)
                                {
                                    childNeutralExperimentalFragments = Ms2ScanWithSpecificMass.GetNeutralExperimentalFragments(ms3ChildScan, commonParameters);
                                }
                                var theChildScan = new Ms2ScanWithSpecificMass(ms3ChildScan, precursorMz,
                                    precursorCharge, fullFilePath, commonParameters, childNeutralExperimentalFragments);

                                scan.ChildScans.Add(theChildScan);
                            }
                        }

                        scansWithPrecursors[i].Add(scan);
                    }
                }
            });

        return scansWithPrecursors;
    }

    private static List<List<Ms2ScanWithSpecificMass>> GetParentScansByGroup(List<Ms2ScanWithSpecificMass>[] scansWithPrecursors)
    {
        if (scansWithPrecursors.Length == 0)
        {
            return new List<List<Ms2ScanWithSpecificMass>>();
        }

        var childScanNumbers = new HashSet<int>(scansWithPrecursors.SelectMany(p => p.SelectMany(v => v.ChildScans.Select(x => x.OneBasedScanNumber))));
        return scansWithPrecursors
            .Where(p => p.Any() && !childScanNumbers.Contains(p.First().OneBasedScanNumber))
            .OrderBy(p => p.First().OneBasedScanNumber)
            .Select(p => p.ToList())
            .ToList();
    }
}
