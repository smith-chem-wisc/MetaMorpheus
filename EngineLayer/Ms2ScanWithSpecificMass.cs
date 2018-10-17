using System;
using System.Collections.Generic;
using System.Linq;
using Chemistry;
using MassSpectrometry;

namespace EngineLayer
{
    public class Ms2ScanWithSpecificMass : IScan
    {
        public Ms2ScanWithSpecificMass(MsDataScan mzLibScan, double precursorMonoisotopicPeakMz, int precursorCharge, string fullFilePath, CommonParameters commonParam, IsotopicEnvelope[] neutralExperimentalFragments = null)
        {
            PrecursorMonoisotopicPeakMz = precursorMonoisotopicPeakMz;
            PrecursorCharge = precursorCharge;
            PrecursorMass = PrecursorMonoisotopicPeakMz.ToMass(precursorCharge);
            FullFilePath = fullFilePath;

            TheScan = mzLibScan;

            ExperimentalFragments = neutralExperimentalFragments ?? GetNeutralExperimentalFragments(mzLibScan, commonParam);

            if (ExperimentalFragments.Any())
            {
                DeconvolutedMonoisotopicMasses = ExperimentalFragments.Select(p => p.monoisotopicMass).ToArray();
            }
        }

        public MsDataScan TheScan { get; }
        public double PrecursorMonoisotopicPeakMz { get; }
        public double PrecursorMass { get; }
        public int PrecursorCharge { get; }
        public string FullFilePath { get; }
        public IsotopicEnvelope[] ExperimentalFragments { get; private set; }
        private double[] DeconvolutedMonoisotopicMasses;

        public int OneBasedScanNumber => TheScan.OneBasedScanNumber;

        public int? OneBasedPrecursorScanNumber => TheScan.OneBasedPrecursorScanNumber;

        public double RetentionTime => TheScan.RetentionTime;

        public int NumPeaks => TheScan.MassSpectrum.Size;

        public double TotalIonCurrent => TheScan.TotalIonCurrent;

        public static IsotopicEnvelope[] GetNeutralExperimentalFragments(MsDataScan scan, CommonParameters commonParam)
        {
            int minZ = 1;
            int maxZ = 10;

            var neutralExperimentalFragmentMasses = scan.MassSpectrum.Deconvolute(scan.MassSpectrum.Range, 
                minZ, maxZ, commonParam.DeconvolutionMassTolerance.Value, commonParam.DeconvolutionIntensityRatio).ToList();
            
            if (commonParam.AssumeOrphanPeaksAreZ1Fragments)
            {
                HashSet<double> alreadyClaimedMzs = new HashSet<double>(neutralExperimentalFragmentMasses
                    .SelectMany(p => p.peaks.Select(v => ClassExtensions.RoundedDouble(v.mz).Value)));

                for (int i = 0; i < scan.MassSpectrum.XArray.Length; i++)
                {
                    double mz = scan.MassSpectrum.XArray[i];
                    double intensity = scan.MassSpectrum.YArray[i];

                    if (!alreadyClaimedMzs.Contains(ClassExtensions.RoundedDouble(mz).Value))
                    {
                        neutralExperimentalFragmentMasses.Add(new IsotopicEnvelope(
                            new List<(double mz, double intensity)> { (mz, intensity) },
                            mz.ToMass(1), 1, intensity, 0, 0));
                    }
                }
            }

            return neutralExperimentalFragmentMasses.OrderBy(p => p.monoisotopicMass).ToArray();
        }

        public IsotopicEnvelope GetClosestExperimentalFragmentMass(double theoreticalNeutralMass)
        {
            if (DeconvolutedMonoisotopicMasses.Length == 0)
            {
                return null;
            }
            return ExperimentalFragments[GetClosestFragmentMass(theoreticalNeutralMass).Value];
        }

        private int? GetClosestFragmentMass(double mass)
        {
            if (DeconvolutedMonoisotopicMasses.Length == 0)
            {
                return null;
            }
            int index = Array.BinarySearch(DeconvolutedMonoisotopicMasses, mass);
            if (index >= 0)
            {
                return index;
            }
            index = ~index;

            if (index >= DeconvolutedMonoisotopicMasses.Length)
            {
                return index - 1;
            }
            if (index == 0)
            {
                return index;
            }

            if (mass - DeconvolutedMonoisotopicMasses[index - 1] > DeconvolutedMonoisotopicMasses[index] - mass)
            {
                return index;
            }
            return index - 1;
        }
    }
}