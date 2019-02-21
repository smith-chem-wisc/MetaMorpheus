using System;
using System.Collections.Generic;
using System.Linq;
using Chemistry;
using MassSpectrometry;
using Proteomics.Fragmentation;

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
            TheExperimentalFragments = ExperimentalFragments;
            if (ExperimentalFragments.Any())
            {
                DeconvolutedMonoisotopicMasses = ExperimentalFragments.Select(p => p.monoisotopicMass).ToArray();
                TheDeconvolutedMonoisotopicMasses = DeconvolutedMonoisotopicMasses;
            }
        }

        public MsDataScan TheScan { get; }
        public double PrecursorMonoisotopicPeakMz { get; }
        public double PrecursorMass { get; }
        public int PrecursorCharge { get; }
        public string FullFilePath { get; }
        public IsotopicEnvelope[] ExperimentalFragments { get; private set; }
        public double[] DeconvolutedMonoisotopicMasses { get; private set; }

        public IsotopicEnvelope[] TheExperimentalFragments { get; private set; }
        public double[] TheDeconvolutedMonoisotopicMasses { get; private set; }

        public List<Ms2ScanWithSpecificMass> childMs2ScanWithSpecificMass { get; set; }
        //Get all experimentalFragments from MS2 scan and its children scan.
        public void GetAllNeutralExperimentalFragments(MsDataScan scan, CommonParameters commonParam)
        {
                     
            if (childMs2ScanWithSpecificMass.Count > 0)
            {
                foreach (var aMsScan in childMs2ScanWithSpecificMass)
                {
                    aMsScan.ExperimentalFragments = GetNeutralExperimentalFragments(aMsScan.TheScan, commonParam);
                    ExperimentalFragments = ExperimentalFragments.Concat(aMsScan.ExperimentalFragments).OrderBy(p => p.monoisotopicMass).ToArray();
                }

                if (ExperimentalFragments.Any())
                {
                    DeconvolutedMonoisotopicMasses = ExperimentalFragments.Select(p => p.monoisotopicMass).ToArray();
                }
            }
        }

        //TO DO: Optimize the filter. (Whether a product ion can be found in this scan.)
        public bool AllowProductType(Product product)
        {
            if (product.ProductType == ProductType.M || !TheScan.DissociationType.HasValue || TheScan.DissociationType.Value == DissociationType.Unknown) //!TheScan.DissociationType.HasValue or DissociationType.Unknown is used in XLtest
            {
                return true;
            }
            if(TheScan.DissociationType.Value == DissociationType.CID || TheScan.DissociationType.Value == DissociationType.HCD || TheScan.DissociationType.Value == DissociationType.EThcD)
            {
                if (product.ProductType == ProductType.b || product.ProductType == ProductType.y)
                {
                    return true;
                }
            }
            if (TheScan.DissociationType.Value == DissociationType.ETD || TheScan.DissociationType.Value == DissociationType.EThcD)
            {
                if (product.ProductType == ProductType.c || product.ProductType == ProductType.zDot)
                {
                    return true;
                }
            }
            return false;
        }

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

        public IsotopicEnvelope GetClosestExperimentalFragmentMass(IsotopicEnvelope[] ExperimentalFragments, double[] DeconvolutedMonoisotopicMasses, double theoreticalNeutralMass)
        {
            if (DeconvolutedMonoisotopicMasses.Length == 0)
            {
                return null;
            }
            return ExperimentalFragments[GetClosestFragmentMass(DeconvolutedMonoisotopicMasses, theoreticalNeutralMass).Value];
        }

        private int? GetClosestFragmentMass(double[] DeconvolutedMonoisotopicMasses, double mass)
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