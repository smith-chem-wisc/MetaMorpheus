using Chemistry;
using MassSpectrometry;
using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Linq;

namespace EngineLayer
{
    public class Ms2ScanWithSpecificMass
    {
        public Ms2ScanWithSpecificMass(MsDataScan mzLibScan, double precursorMonoisotopicPeakMz, int precursorCharge, string fullFilePath, CommonParameters commonParam,
            IsotopicEnvelope[] neutralExperimentalFragments = null, double? precursorIntensity = null, int? envelopePeakCount = null, double? precursorFractionalIntensity = null)
        {
            PrecursorMonoisotopicPeakMz = precursorMonoisotopicPeakMz;
            PrecursorCharge = precursorCharge;
            PrecursorMass = PrecursorMonoisotopicPeakMz.ToMass(precursorCharge);
            PrecursorIntensity = precursorIntensity ?? 1;
            PrecursorEnvelopePeakCount = envelopePeakCount ?? 1;
            PrecursorFractionalIntensity = precursorFractionalIntensity ?? -1;
            FullFilePath = fullFilePath;
            ChildScans = new List<Ms2ScanWithSpecificMass>();
            NativeId = mzLibScan.NativeId;

            TheScan = mzLibScan;

            if (commonParam.DissociationType != DissociationType.LowCID)
            {
                ExperimentalFragments = neutralExperimentalFragments ?? GetNeutralExperimentalFragments(mzLibScan, commonParam);
            }
            if (ExperimentalFragments != null && ExperimentalFragments.Any())
            {
                DeconvolutedMonoisotopicMasses = ExperimentalFragments.Select(p => p.MonoisotopicMass).ToArray();
            }
            else
            {
                DeconvolutedMonoisotopicMasses = new double[0];
            }
        }

        public MsDataScan TheScan { get; }
        public double PrecursorMonoisotopicPeakMz { get; }
        public double PrecursorMass { get; }
        public int PrecursorCharge { get; }
        public double PrecursorIntensity { get; }
        public int PrecursorEnvelopePeakCount { get; }
        public double PrecursorFractionalIntensity { get; }
        public string FullFilePath { get; }
        public IsotopicEnvelope[] ExperimentalFragments { get; private set; }
        public List<Ms2ScanWithSpecificMass> ChildScans { get; set; } // MS2/MS3 scans that are children of this MS2 scan
        private double[] DeconvolutedMonoisotopicMasses;
        public string NativeId { get; }
        public int OneBasedScanNumber => TheScan.OneBasedScanNumber;
        public int? OneBasedPrecursorScanNumber => TheScan.OneBasedPrecursorScanNumber;
        public double RetentionTime => TheScan.RetentionTime;
        public int NumPeaks => TheScan.MassSpectrum.Size;
        public double TotalIonCurrent => TheScan.TotalIonCurrent;
        /// <summary>
        /// An array containing the intensities of the reporter ions for isobaric mass tags. 
        /// If multiplex quantification wasn't performed, this will be null
        /// </summary>
        public double[]? IsobaricMassTagReporterIonIntensities { get; private set; }

        public static IsotopicEnvelope[] GetNeutralExperimentalFragments(MsDataScan scan, CommonParameters commonParam)
        {
            var neutralExperimentalFragmentMasses =
                Deconvoluter.Deconvolute(scan, commonParam.ProductDeconvolutionParameters, scan.MassSpectrum.Range).ToList();

            if (!commonParam.AssumeOrphanPeaksAreZ1Fragments || scan.MassSpectrum is NeutralMassSpectrum)
                return neutralExperimentalFragmentMasses.OrderBy(p => p.MonoisotopicMass).ToArray();

            HashSet<double> alreadyClaimedMzs = new HashSet<double>(neutralExperimentalFragmentMasses
                .SelectMany(p => p.Peaks.Select(v => v.mz.RoundedDouble()!.Value)));

            int charge = scan.Polarity == Polarity.Positive ? 1 : -1;
            for (int i = 0; i < scan.MassSpectrum.XArray.Length; i++)
            {
                double mz = scan.MassSpectrum.XArray[i];
                double intensity = scan.MassSpectrum.YArray[i];

                if (!alreadyClaimedMzs.Contains(mz.RoundedDouble()!.Value))
                {
                    neutralExperimentalFragmentMasses.Add(new IsotopicEnvelope(
                        new List<(double mz, double intensity)> { (mz, intensity) },
                        mz.ToMass(charge), charge, intensity, 0));
                }
            }

            return neutralExperimentalFragmentMasses.OrderBy(p => p.MonoisotopicMass).ToArray();
        }


        private bool ChildScansHaveHigherHcdEnergy(out Ms2ScanWithSpecificMass mostIntenseChildScan)
        {
            mostIntenseChildScan = null;
            if (ChildScans.IsNullOrEmpty()) return false;
            double parentScanHcdEnergy = Double.TryParse(TheScan.HcdEnergy, out var thisHcdEnergy) ? thisHcdEnergy : 0;
            double childScanMaxHcdEnergy = ChildScans.Max(s => Double.TryParse(s.TheScan.HcdEnergy, out var hcdEnergy) ? hcdEnergy : 0);

            if (childScanMaxHcdEnergy > parentScanHcdEnergy)
            {
                mostIntenseChildScan = ChildScans
                return true;
            }

            return false;
        }

        /// <summary>
        /// Writes the reporter ion intensities into the IsobaricMassTagReporterIonIntensities property
        /// If the scan has
        /// </summary>
        /// <param name="massTag"></param>
        public void GetIsobaricMassTagReporterIonIntensities(IsobaricMassTag massTag)
        {
            if (TheScan.ScanWindowRange.Minimum > massTag.ReporterIonMzs[^1] // If the scan window is above the highest reporter ion mass, check if MS3 was used for the reporter ions
                || (ChildScans.Any() 
                && ChildScans.Max(s => Double.TryParse(s.TheScan.HcdEnergy, out var hcdEnergy) ? hcdEnergy : 0) > 
                (Double.TryParse(TheScan.HcdEnergy, out var thisHcdEnergy) ? thisHcdEnergy : 0))) // If any child scan has higher HCD energy than this scan, assume reporter ions were measured in MS3
            {
                if (ChildScans.IsNullOrEmpty()) return;
                if (ChildScans.All(s => s.TheScan.ScanWindowRange.Minimum > massTag.ReporterIonMzs[^1])) return;
                // If there are multiple child scans, we'll just pick the most intense one (In general, there should only be one MS3 scan for reporter ions)
                var ms3Scan = ChildScans.Where(s => s.TheScan.MsnOrder == 3)
                    .MaxBy(s => s.TotalIonCurrent);
                IsobaricMassTagReporterIonIntensities = massTag.GetReporterIonIntensities(ms3Scan.TheScan.MassSpectrum);
                return;
            }
            IsobaricMassTagReporterIonIntensities = massTag.GetReporterIonIntensities(TheScan.MassSpectrum);
        }
       
        public IsotopicEnvelope GetClosestExperimentalIsotopicEnvelope(double theoreticalNeutralMass)
        {
            if (DeconvolutedMonoisotopicMasses.Length == 0)
            {
                return null;
            }
            return ExperimentalFragments[GetClosestFragmentMass(theoreticalNeutralMass)];
        }

        public int GetClosestFragmentMass(double mass)
        {
            int index = Array.BinarySearch(DeconvolutedMonoisotopicMasses, mass);
            if (index >= 0)
            {
                return index;
            }
            index = ~index;

            if (index == DeconvolutedMonoisotopicMasses.Length)
            {
                return index - 1;
            }
            if (index == 0 || mass - DeconvolutedMonoisotopicMasses[index - 1] > DeconvolutedMonoisotopicMasses[index] - mass)
            {
                return index;
            }

            return index - 1;
        }

        //look for IsotopicEnvelopes which are in the range of acceptable mass 
        public IsotopicEnvelope[] GetClosestExperimentalIsotopicEnvelopeList(double minimumMass, double maxMass)
        {

            if (DeconvolutedMonoisotopicMasses.Length == 0)
            {
                return null;
            }

            //if no mass is in this range, then return null
            if (DeconvolutedMonoisotopicMasses[0] > maxMass || DeconvolutedMonoisotopicMasses.Last() < minimumMass)
            {
                return null;
            }

            int startIndex = GetClosestFragmentMass(minimumMass);
            int endIndex = GetClosestFragmentMass(maxMass);

            //the index we get from GetClosestFragmentMass is the closest mass, while the acceptable mass we need is between minimumMass and maxMass
            //so the startIndex mass is supposed to be larger than minimumMass, if not , then startIndex increases by 1;
            //the endIndex mass is supposed to be smaller than maxMass, if not , then endIndex decreases by 1;
            if (DeconvolutedMonoisotopicMasses[startIndex]<minimumMass)
            {
                startIndex = startIndex+1;
            }
            if(DeconvolutedMonoisotopicMasses[endIndex] > maxMass)
            {
                endIndex = endIndex - 1;
            }
            int length = endIndex - startIndex + 1;

            if (length < 1)
            {
                return null;
            }
            IsotopicEnvelope[] isotopicEnvelopes = ExperimentalFragments.Skip(startIndex).Take(length).ToArray();
            return isotopicEnvelopes;
        }
    }
}