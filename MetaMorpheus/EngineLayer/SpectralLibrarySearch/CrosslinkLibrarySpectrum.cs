using Omics.Fragmentation;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Text.RegularExpressions;
using EngineLayer.CrosslinkSearch;
using System.Threading.Tasks;

namespace EngineLayer
{
    public class CrosslinkLibrarySpectrum : LibrarySpectrum
    {
        public CrosslinkLibrarySpectrum BetaPeptideSpectrum { get; }
        public string AlphaPeptideSequence { get; private set; }
        public string BetaPeptideSequence { get; private set; }
        public string UniqueSequence { get; private set; }
        public bool IsBetaPeptide { get; }
        public static Regex CrosslinkRegex = new Regex(@"\(\d+\)");
        public override string Name => UniqueSequence + "/" + ChargeState; 

        public CrosslinkLibrarySpectrum(
            string uniqueSequence,
            double precursorMz,
            int precursorCharge,
            List<MatchedFragmentIon> peaks,
            double rt,
            CrosslinkSpectralMatch betaPeptide,
            bool isDecoy = false) : this(
                uniqueSequence,
                precursorMz, 
                precursorCharge,
                peaks, 
                rt, 
                new CrosslinkLibrarySpectrum(uniqueSequence, precursorMz, precursorCharge, betaPeptide.MatchedFragmentIons, rt), 
                isDecoy) 
        { }

        public CrosslinkLibrarySpectrum(
            string uniqueSequence,
            double precursorMz,
            int precursorCharge,
            List<MatchedFragmentIon> peaks,
            double rt,
            List<MatchedFragmentIon> betaPeaks,
            bool isDecoy = false) : this(
                uniqueSequence,
                precursorMz,
                precursorCharge,
                peaks,
                rt,
                new CrosslinkLibrarySpectrum(uniqueSequence, precursorMz, precursorCharge, betaPeaks, rt),
                isDecoy)   
        { }

        public CrosslinkLibrarySpectrum(
            string uniqueSequence,
            double precursorMz,
            int precursorCharge,
            List<MatchedFragmentIon> peaks,
            double rt,
            CrosslinkLibrarySpectrum betaSpectrum = null,
            bool isDecoy = false) : base(uniqueSequence, precursorMz, precursorCharge, peaks, rt, isDecoy)
        {
            UniqueSequence = uniqueSequence;
            if (betaSpectrum == null)
            {
                IsBetaPeptide = true;
            }
            else
            {
                BetaPeptideSpectrum = betaSpectrum;
            }
            SetAlphaBetaSequence();
        }

        private void SetAlphaBetaSequence()
        {
            string[] uniqueSequenceSplit = CrosslinkRegex.Split(UniqueSequence);
            if (uniqueSequenceSplit.Length >= 2)
            {
                AlphaPeptideSequence = uniqueSequenceSplit[0];
                BetaPeptideSequence = uniqueSequenceSplit[1];
            }
            else
            {
                AlphaPeptideSequence = null;
                BetaPeptideSequence = null;
            }
        }

        public override string ToString()
        {
            StringBuilder spectrum = new();
            spectrum.AppendLine("Name: " + Name);
            spectrum.AppendLine("MW: " + PrecursorMz);
            spectrum.Append("Comment: ");
            spectrum.Append("Parent=" + PrecursorMz);
            spectrum.AppendLine(" RT=" + RetentionTime);
            spectrum.Append("Num alpha peaks: " + MatchedFragmentIons.Count);
            spectrum.AppendLine(", Num beta peaks: " + BetaPeptideSpectrum.MatchedFragmentIons.Count);

            double maxIntensity = Math.Max(MatchedFragmentIons.Max(b => b.Intensity),
                BetaPeptideSpectrum.MatchedFragmentIons.Max(s => s.Intensity));

            foreach (MatchedFragmentIon matchedIon in MatchedFragmentIons)
            {
                double intensityFraction = matchedIon.Intensity / maxIntensity;

                string neutralLoss = null;
                if (matchedIon.NeutralTheoreticalProduct.NeutralLoss != 0)
                {
                    neutralLoss = "-" + matchedIon.NeutralTheoreticalProduct.NeutralLoss;
                }

                spectrum.AppendLine(matchedIon.Mz + "\t" + intensityFraction + "\t" + "\"" +
                    matchedIon.NeutralTheoreticalProduct.ProductType.ToString() +
                    matchedIon.NeutralTheoreticalProduct.FragmentNumber.ToString() + "^" +
                    matchedIon.Charge + neutralLoss + "/" + 0 + "ppm" + "\"");
            }

            foreach (MatchedFragmentIon matchedIon in BetaPeptideSpectrum.MatchedFragmentIons)
            {
                double intensityFraction = matchedIon.Intensity / maxIntensity;

                string neutralLoss = null;
                if (matchedIon.NeutralTheoreticalProduct.NeutralLoss != 0)
                {
                    neutralLoss = "-" + matchedIon.NeutralTheoreticalProduct.NeutralLoss;
                }

                spectrum.AppendLine(matchedIon.Mz + "\t" + intensityFraction + "\t" + "\"" +
                    matchedIon.NeutralTheoreticalProduct.ProductType.ToString() +
                    matchedIon.NeutralTheoreticalProduct.FragmentNumber.ToString() + "^" +
                    matchedIon.Charge + neutralLoss + "/" + 0 + "ppm" + "\"" +
                    "\t" + "BetaPeptideIon");
            }

            return spectrum.ToString().Trim();
        }
    }
}
