//using MassSpectrometry;
//using System;
//using System.Collections.Generic;
//using System.Linq;
//using System.Text;
//using MzLibUtil;
//using FlashLFQ;

//namespace EngineLayer.DIA.XicConstruction
//{
    //public class HighestIsotopePeakXicConstructor : XicConstructor
    //{
    //    public DeconvolutionParameters DeconParameters { get; set; }
    //    public Dictionary<IIndexedPeak, ExtractedIonChromatogram> PeakXicDictionary { get; set; } 

    //    public HighestIsotopePeakXicConstructor(Tolerance peakFindingTolerance, int maxMissedScansAllowed, double maxPeakHalfWidth, int minNumberOfPeaks, DeconvolutionParameters deconParameters, XicSpline? xicSpline = null)
    //        : base(peakFindingTolerance, maxMissedScansAllowed, maxPeakHalfWidth, minNumberOfPeaks, xicSpline)
    //    {
    //        DeconParameters = deconParameters;
    //    }

        //public override List<ExtractedIonChromatogram> GetAllXics(MsDataScan[] scans, MzRange isolationRange = null)
        //{
        //    var mzPeakIndexingEngine = PeakIndexingEngine.InitializeIndexingEngine(scans);
        //    if (PeakXicDictionary == null)
        //    {
        //        var allMzXics = mzPeakIndexingEngine.GetAllXics(PeakFindingTolerance, MaxMissedScansAllowed, MaxPeakHalfWidth, MinNumberOfPeaks, out var matchedPeaks);
        //        PeakXicDictionary = matchedPeaks;
        //    }
        //    var foundedXics = new HashSet<ExtractedIonChromatogram>();
        //    for (int i = 0; i < scans.Length; i++)
        //    {
        //        var envelopes = Deconvoluter.Deconvolute(scans[i], DeconParameters, isolationRange);
        //        foreach (var envelope in envelopes)
        //        {
        //            var highestPeak = envelope.Peaks.MaxBy(p => p.intensity);
        //            var indexedPeak = mzPeakIndexingEngine.GetIndexedPeak(highestPeak.mz, i, PeakFindingTolerance);
        //            if (indexedPeak != null && PeakXicDictionary.ContainsKey(indexedPeak))
        //            {
        //                var foundXic = PeakXicDictionary[indexedPeak];
        //                if (foundXic != null)
        //                {
        //                    foundXic.AveragedMassOrMz = envelope.MonoisotopicMass;
        //                    foundedXics.Add(foundXic);
        //                }
        //            }
        //        }
        //    }
        //    return foundedXics.ToList();
        //}
//    }
//}
