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

        /// <summary>
        /// Writes the reporter ion intensities into the IsobaricMassTagReporterIonIntensities property
        /// If the scan has
        /// </summary>
        /// <param name="massTag"></param>
        public void SetIsobaricMassTagReporterIonIntensities(IsobaricMassTag massTag)
        {
            if (UseChildScansForIsobaricQuant(out var mostIntenseChildScan))
            {
                IsobaricMassTagReporterIonIntensities = massTag.GetReporterIonIntensities(mostIntenseChildScan.TheScan.MassSpectrum);
                return;
            }
            IsobaricMassTagReporterIonIntensities = massTag.GetReporterIonIntensities(TheScan.MassSpectrum);
        }

        /// <summary>
        /// Helper method to determine if child scans contain isobaric mass tag (e.g., TMT) reporter ions
        /// In most cases, if MS3 scans exist, there will only be one MS3 scan per MS2 scan, and the ChildScans list will contain that one MS3 scan
        /// If there are methods that generate multiple child scans, then this method will return the most intense child scan
        /// However, if we ever have multiple child scans with isobaric mass tag reporter ions, we should revisit this logic
        /// </summary>
        /// <param name="mostIntenseChildScan">The child scan that has the highest ion current </param>
        private bool UseChildScansForIsobaricQuant(out Ms2ScanWithSpecificMass mostIntenseChildScan)
        {
            mostIntenseChildScan = null;
            if (ChildScans.IsNullOrEmpty()) return false;
            mostIntenseChildScan = ChildScans.MaxBy(s => s.TotalIonCurrent);
            if(mostIntenseChildScan != null)
            {
                return true;
            }
            return false;
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