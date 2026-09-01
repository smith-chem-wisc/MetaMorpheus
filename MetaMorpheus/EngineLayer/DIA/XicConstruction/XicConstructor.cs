using MassSpectrometry;
using MzLibUtil;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Threading.Tasks;

namespace EngineLayer.DIA
{
    /// <summary>
    /// XicConstructor defines the process of constructing extracted ion chromatograms (XICs) for a given set of scans. It could operate on spectral peaks
    /// or deconvoluted masses depending on the implementation. Any XicConstructor should have a method that returns all XICs that can be found in the set of scans.
    /// </summary>
    public abstract class XicConstructor
    {
        public Tolerance PeakFindingTolerance { get; set; } 
        public int MaxMissedScansAllowed { get; set; } 
        public double MaxPeakHalfWidth { get; set; }
        public int MinNumberOfPeaks { get; set; }
        public XicSpline? XicSplineEngine { get; set; }

        public XicConstructor(Tolerance peakFindingTolerance, int maxMissedScansAllowed, double maxPeakHalfWidth, int minNumberOfPeaks, XicSpline? xicSpline = null)
        {
            this.PeakFindingTolerance = peakFindingTolerance;
            MaxMissedScansAllowed = maxMissedScansAllowed;
            MaxPeakHalfWidth = maxPeakHalfWidth;
            MinNumberOfPeaks = minNumberOfPeaks;
            XicSplineEngine = xicSpline;
        }

        public abstract List<ExtractedIonChromatogram> GetAllXics(MsDataScan[] scans, MzRange isolationRange = null);

        public List<ExtractedIonChromatogram> GetAllXicsWithXicSpline(MsDataScan[] scans, MzRange isolationRange = null, int? numberOfThreads = null)
        {
            var allXics = GetAllXics(scans, isolationRange);
            if (XicSplineEngine != null)
            {
                Parallel.ForEach(Partitioner.Create(0, allXics.Count), new ParallelOptions { MaxDegreeOfParallelism = numberOfThreads.HasValue ? numberOfThreads.Value : Environment.ProcessorCount - 2 },
                (partitionRange, loopState) =>
                {
                    for (int i = partitionRange.Item1; i < partitionRange.Item2; i++)
                    {
                        var xic = allXics[i];
                        XicSplineEngine.SetXicSplineXYData(xic);
                    }
                });
            }
            return allXics;
        }
    }
}
