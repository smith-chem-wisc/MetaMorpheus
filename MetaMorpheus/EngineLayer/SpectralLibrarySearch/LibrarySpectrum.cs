using Proteomics.Fragmentation;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System;
using MassSpectrometry.MzSpectra;

namespace EngineLayer
{
    /// <summary>
    /// This class describes one element of a spectral library.
    /// </summary>
    public class LibrarySpectrum
    {

        public string Sequence { get; set; }
        public double RetentionTime { get; set; }
        public double PrecursorMz { get; set; }
        public int ChargeState { get; set; }
        public List<MatchedFragmentIon> MatchedFragmentIons { get; set; }
        public bool IsDecoy { get; set; }

        public double[] XArray { get; private set; }
        public double[] YArray { get; private set; }
        public string Name
        {
            get { return Sequence + "/" + ChargeState; }
        }

        public LibrarySpectrum(string sequence, double precursorMz, int chargeState, List<MatchedFragmentIon> peaks, double rt, bool isDecoy = false)
        {
            Sequence = sequence;
            PrecursorMz = precursorMz;
            MatchedFragmentIons = peaks;
            ChargeState = chargeState;
            IsDecoy = isDecoy;
            RetentionTime = rt;
            XArray = peaks.Select(p => p.Mz).ToArray();
            YArray = peaks.Select(p => p.Intensity).ToArray();
            Array.Sort(XArray, YArray);
        }

        /// <summary>
        /// This function enables the spectrum angle to be computed between an individual experimental spectrum and the loaded library spectrum within MetaDraw
        /// </summary>
        /// <param name="librarySpectrum"></param>
        /// <returns></returns>
        public string CalculateSpectralAngleOnTheFly(List<MatchedFragmentIon> spectrumMatchFragments)
        {

            SpectralSimilarity spectraComparison = new SpectralSimilarity(
                spectrumMatchFragments.Select(f => f.Mz).ToArray(),
                spectrumMatchFragments.Select(f => f.Intensity).ToArray(),
                MatchedFragmentIons.Select(f => f.Mz).ToArray(),
                MatchedFragmentIons.Select(f => f.Intensity).ToArray(),
                SpectralSimilarity.SpectrumNormalizationScheme.mostAbundantPeak,
                toleranceInPpm: 20,
                allPeaks: true);
            double? spectralContrastAngle = spectraComparison.SpectralContrastAngle();

            return spectralContrastAngle == null
                ? "N/A"
                : ((double)spectralContrastAngle).ToString("F4");
        }

        public override string ToString()
        {
            StringBuilder spectrum = new StringBuilder();
            spectrum.Append("Name: " + Name);
            spectrum.Append("\nMW: " + PrecursorMz);
            spectrum.Append("\nComment: ");
            spectrum.Append("Parent=" + PrecursorMz);
            spectrum.Append(" RT=" + RetentionTime);
            spectrum.Append("\nNum peaks: " + MatchedFragmentIons.Count);

            double maxIntensity = MatchedFragmentIons.Select(b => b.Intensity).Max();

            foreach (MatchedFragmentIon matchedIon in MatchedFragmentIons)
            {
                double intensityFraction = matchedIon.Intensity / maxIntensity;

                string neutralLoss = null;
                if (matchedIon.NeutralTheoreticalProduct.NeutralLoss != 0)
                {
                    neutralLoss = "-" + matchedIon.NeutralTheoreticalProduct.NeutralLoss;
                }

                spectrum.Append("\n"+matchedIon.Mz + "\t" + intensityFraction + "\t" + "\"" +
                    matchedIon.NeutralTheoreticalProduct.ProductType.ToString() +
                    matchedIon.NeutralTheoreticalProduct.FragmentNumber.ToString() + "^" +
                    matchedIon.Charge + neutralLoss + "/" + 0 + "ppm" + "\"");
            }

            return spectrum.ToString();
        }
    }
}