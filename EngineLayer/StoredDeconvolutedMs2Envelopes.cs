using System;
using System.Collections.Generic;
using System.Linq;
using Chemistry;
using MassSpectrometry;
using Proteomics.Fragmentation;

namespace EngineLayer
{
    public class StoredDeconvolutedMs2Envelopes
    {
        private Dictionary<int, IsotopicEnvelope[]> ScanToEnvelopes;
        private Dictionary<int, double[]> ScanToDeconvolutedMonoisotopicMasses;
        private Dictionary<int, HashSet<double>> ScanToDeconvolutedMzs;

        public StoredDeconvolutedMs2Envelopes(Dictionary<int, IsotopicEnvelope[]> scanToEnvelopes, Dictionary<int, double[]> scanToDeconvolutedMonoisotopicMasses,
            Dictionary<int, HashSet<double>> scanToDeconvolutedMzs)
        {
            ScanToEnvelopes = scanToEnvelopes;
            ScanToDeconvolutedMonoisotopicMasses = scanToDeconvolutedMonoisotopicMasses;
            ScanToDeconvolutedMzs = scanToDeconvolutedMzs;
        }

        public IsotopicEnvelope BinarySearchForDeconvolutedMass(int oneBasedScan, Product product)
        {
            var deconvolutedEnvsForThisScan = ScanToEnvelopes[oneBasedScan];
            var deconvolutedMonoisotopicMasses = ScanToDeconvolutedMonoisotopicMasses[oneBasedScan];

            if (deconvolutedEnvsForThisScan.Length == 0)
            {
                return null;
            }

            int index = Array.BinarySearch(deconvolutedMonoisotopicMasses, product.NeutralMass);

            if (index >= 0)
            {
                return deconvolutedEnvsForThisScan[index];
            }

            index = ~index;

            if (index >= deconvolutedEnvsForThisScan.Length)
            {
                return deconvolutedEnvsForThisScan[index - 1];
            }

            if (index == 0)
            {
                return deconvolutedEnvsForThisScan[index];
            }

            if (product.NeutralMass - deconvolutedEnvsForThisScan[index - 1].monoisotopicMass >
                deconvolutedEnvsForThisScan[index].monoisotopicMass - product.NeutralMass)
            {
                return deconvolutedEnvsForThisScan[index];
            }

            return deconvolutedEnvsForThisScan[index - 1];
        }

        public bool HasClaimedMzPeak(int oneBasedScan, double mz)
        {
            return ScanToDeconvolutedMzs[oneBasedScan].Contains(ClassExtensions.RoundedDouble(mz).Value);
        }

        public IEnumerable<IsotopicEnvelope> IsotopicEnvelopesInScan(int oneBasedScan)
        {
            return ScanToEnvelopes[oneBasedScan];
        }

        public static StoredDeconvolutedMs2Envelopes DeconvoluteAndStoreMs2Scans(MsDataScan[] ms2Scans, CommonParameters commonParameters)
        {
            var deconvolutedMs2IsotopicEnvelopes = new Dictionary<int, IsotopicEnvelope[]>();
            var deconvolutedPeakMzs = new Dictionary<int, HashSet<double>>();
            var deconvolutedMonoisotopicMasses = new Dictionary<int, double[]>();

            double ms2DeconvolutionPpmTolerance = 5.0;
            int minZ = 1;
            int maxZ = 10;

            foreach (MsDataScan scan in ms2Scans)
            {
                if (deconvolutedMs2IsotopicEnvelopes.ContainsKey(scan.OneBasedScanNumber))
                {
                    continue;
                }

                // deconvolute the scan
                var isotopicEnvelopes = scan.MassSpectrum.Deconvolute(scan.MassSpectrum.Range, minZ, maxZ,
                    ms2DeconvolutionPpmTolerance, commonParameters.DeconvolutionIntensityRatio).OrderBy(p => p.monoisotopicMass).ToArray();

                // store the scan's deconvoluted envelopes
                deconvolutedMs2IsotopicEnvelopes.Add(scan.OneBasedScanNumber, isotopicEnvelopes);

                // store the scan's deconvoluted envelope peaks
                HashSet<double> deconvolutedMzForThisScan = new HashSet<double>();

                foreach (IsotopicEnvelope isotopicEnvelope in isotopicEnvelopes)
                {
                    foreach (var peak in isotopicEnvelope.peaks)
                    {
                        deconvolutedMzForThisScan.Add(ClassExtensions.RoundedDouble(peak.mz).Value);
                    }
                }

                deconvolutedPeakMzs.Add(scan.OneBasedScanNumber, deconvolutedMzForThisScan);
                deconvolutedMonoisotopicMasses.Add(scan.OneBasedScanNumber, isotopicEnvelopes.Select(p => p.monoisotopicMass).OrderBy(p => p).ToArray());
            }

            return new StoredDeconvolutedMs2Envelopes(deconvolutedMs2IsotopicEnvelopes, deconvolutedMonoisotopicMasses, deconvolutedPeakMzs);
        }
    }
}
