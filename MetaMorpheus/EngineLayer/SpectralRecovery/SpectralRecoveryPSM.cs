using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text;
using Chemistry;
using Easy.Common.Extensions;
using FlashLFQ;
using MassSpectrometry;
using MassSpectrometry.MzSpectra;
using Proteomics;
using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using ThermoFisher.CommonCore.Data.Business;
using Peptide = Proteomics.AminoAcidPolymer.Peptide;

namespace EngineLayer.SpectralRecovery
{
    public class SpectralRecoveryPSM : PeptideSpectralMatch
    {
        public ChromatographicPeak AcceptorPeak { get; set; }
        public PeptideSpectralMatch OriginalSpectralMatch { get; private set; }
        /// <summary>
        /// In spectral recovery, we see situations where a mass was selected for fragmentation
        /// but could not be deconvoluted. As a result, base.ScanPrecursorMass == 0. However, the
        /// RecoveredMs2ScanWithSpecificMass contains information on the precursor peak that was closest
        /// to the theoretical precursor m/z. ClosestPrecursorPeak stores that info
        /// </summary>
        public MzPeak ClosestPrecursorPeak { get; }

        public bool DeconvolutablePrecursor { get; }
        public PeptideWithSetModifications DonorPeptide { get; }
        public double? RetentionTimeShift { get; private set; }
        public double? PrecursorIsotpicEnvelopeAngle { get; private set; }

        public SpectralRecoveryPSM(
            ChromatographicPeak acceptorPeak,
            PeptideWithSetModifications peptide,
            double score,
            RecoveredMs2ScanWithSpecificMass scan,
            CommonParameters commonParameters,
            List<MatchedFragmentIon> matchedFragmentIons,
            int notch = -1,
            int scanIndex = -1,
            double xcorr = 0) :
            base(peptide, notch, score, scanIndex, scan, commonParameters, matchedFragmentIons, xcorr)
        {
            AcceptorPeak = acceptorPeak;
            DonorPeptide = peptide;
            ClosestPrecursorPeak = scan.PeakClosestToDonor;
            RetentionTimeShift = SetRetentionTimeShift(acceptorPeak);
            PrecursorIsotpicEnvelopeAngle = GetIsotopeCorrelation(peptide, acceptorPeak);
            DeconvolutablePrecursor = true;

            //In cases where the precursor wasn't deconvoluted, the precursor mz and mass are set to 0 and 0 - z*1.008, respectively.
            // This updates the precursor mz and mass using the ClosestPrecursorPeak mz
            if (Math.Abs(ScanPrecursorMonoisotopicPeakMz) < 0.1)
            {
                DeconvolutablePrecursor = false;

                ScanPrecursorMonoisotopicPeakMz = ClosestPrecursorPeak.Mz;
                ScanPrecursorMass = ScanPrecursorMonoisotopicPeakMz.ToMass(AcceptorPeak.Apex.ChargeState);
            }
        }

        // <summary>
        /// Calculates the difference between the actual and expected retention times in minutes.
        /// Calculates the "Z-Score" by dividing this difference by the standard deviation of the alignment between
        /// peaks in the donor and acceptor files
        /// </summary>
        public static double? SetRetentionTimeShift(ChromatographicPeak acceptorPeak)
        {
            return (acceptorPeak?.RtHypothesis == null || acceptorPeak.Apex == null)
                ? null
                : acceptorPeak.Apex.IndexedPeak.RetentionTime - (double)acceptorPeak.RtHypothesis;
        }

        public void FindOriginalPsm(List<PeptideSpectralMatch> originalSearchPsms)
        {
            OriginalSpectralMatch = originalSearchPsms.FirstOrDefault(
                p => p.FullFilePath.Equals(FullFilePath) &&
                     p.MsDataScan.OneBasedScanNumber == MsDataScan.OneBasedScanNumber &&
                     (FullSequence.Equals(p.FullSequence) ||
                      p.FullSequence.Split('|').Any(seq => seq.Equals(FullSequence))));
        }

        public static string TabSeparatedHeader => string.Join('\t', SpectralRecoveryDataDictionary(null, null).Keys);

        public override string ToString()
        {
            return ToString(new Dictionary<string, int>());
        }

        public override string ToString(IReadOnlyDictionary<string, int> modsToWritePruned)
        {
            return string.Join('\t', SpectralRecoveryDataDictionary(this, modsToWritePruned).Values);
        }

        public static Dictionary<string, string> SpectralRecoveryDataDictionary(SpectralRecoveryPSM srPsm,
            IReadOnlyDictionary<string, int> modsToWritePruned)
        {
            // Get information from base PSM class
            Dictionary<string, string> psmDictionary = new();

            PsmTsvWriter.AddBasicMatchData(psmDictionary, srPsm);
            AddSpectralRecoveryData(psmDictionary, srPsm);
            PsmTsvWriter.AddPeptideSequenceData(psmDictionary, srPsm, modsToWritePruned);
            PsmTsvWriter.AddMatchedIonsData(psmDictionary, srPsm?.MatchedFragmentIons);
            PsmTsvWriter.AddMatchScoreData(psmDictionary, srPsm);

            return psmDictionary;
        }

        /// <summary>
        /// Populate fields specific to RecoveredPSMs
        /// </summary>
        /// <param name="psmDictionary"></param>
        /// <param name="srPsm"></param>
        /// <returns></returns>
        private static void AddSpectralRecoveryData(Dictionary<string, string> psmDictionary,
            SpectralRecoveryPSM srPsm)
        {
            // Chromatographic Peak Info
            psmDictionary[PsmTsvHeader_SpectralRecovery.PeakApexRt] =
                srPsm?.AcceptorPeak?.Apex == null
                    ? " "
                    : srPsm.AcceptorPeak.Apex.IndexedPeak.RetentionTime.ToString(CultureInfo.InvariantCulture);
            psmDictionary[PsmTsvHeader_SpectralRecovery.PeakShift] =
                srPsm?.RetentionTimeShift == null
                    ? " "
                    : srPsm.RetentionTimeShift.NullableToString(CultureInfo.InvariantCulture);

            // Scan Isolation Window info
            psmDictionary[PsmTsvHeader_SpectralRecovery.PrecursorDeconvolutedBool] =
                srPsm == null ? " "
                : srPsm.DeconvolutablePrecursor ? "Y" : "N";
            psmDictionary[PsmTsvHeader_SpectralRecovery.PrecursorIsotopicEnvelopeAngle] =
                srPsm == null
                    ? " "
                    : srPsm.PrecursorIsotpicEnvelopeAngle.NullableToString(CultureInfo.InvariantCulture);
            psmDictionary[PsmTsvHeader_SpectralRecovery.IsolationWindowCenter] =
                srPsm == null
                    ? " "
                    : srPsm.MsDataScan.IsolationMz.NullableToString(CultureInfo.InvariantCulture);
            psmDictionary[PsmTsvHeader_SpectralRecovery.PrecursorOffset] =
                srPsm == null
                    ? " "
                    : CalculatePrecursorOffset(srPsm).NullableToString(CultureInfo.InvariantCulture);
            psmDictionary[PsmTsvHeader_SpectralRecovery.IsolationWindowWidth] =
                srPsm == null
                    ? " "
                    : srPsm.MsDataScan.IsolationWidth.NullableToString(CultureInfo.InvariantCulture);

            // Original Psm Info
            psmDictionary[PsmTsvHeader_SpectralRecovery.OriginalPsmQ] =
                srPsm?.OriginalSpectralMatch == null
                    ? " "
                    : srPsm.OriginalSpectralMatch.FdrInfo.QValue.ToString(CultureInfo.InvariantCulture);
            psmDictionary[PsmTsvHeader_SpectralRecovery.OriginalPsmPEP] =
                srPsm?.OriginalSpectralMatch == null
                    ? " "
                    : srPsm.OriginalSpectralMatch.FdrInfo.PEP.ToString(CultureInfo.InvariantCulture);
            psmDictionary[PsmTsvHeader_SpectralRecovery.OriginalPsmPEP_QValue] =
                srPsm?.OriginalSpectralMatch == null
                    ? " "
                    : srPsm.OriginalSpectralMatch.FdrInfo.PEP_QValue.ToString(CultureInfo.InvariantCulture);
        }

        private static double? CalculatePrecursorOffset(SpectralRecoveryPSM srPsm)
        {
            if (srPsm?.MsDataScan.IsolationMz == null) return null;
            double precursorMz = Math.Abs(srPsm.ScanPrecursorMonoisotopicPeakMz) > 0.1
                ? srPsm.ScanPrecursorMonoisotopicPeakMz
                : srPsm.ClosestPrecursorPeak.Mz;
            return precursorMz - srPsm.MsDataScan.IsolationMz;
        }

        private static double? AveragineMass { get; set; }

        public static double? GetIsotopeCorrelation(PeptideWithSetModifications selectedPeptide, ChromatographicPeak peak)
        {
            if (peak == null) return null;
            double? isotopeAngle = null;
            double fineResolution = 0.01;
            double minimumProbability = 0.005;
            double massCorrection = 0; // Accounts for mass differences caused by mods with unknown chemical formulas.
            
            ChemicalFormula peptideFormula;
            try
            {
                peptideFormula = selectedPeptide.FullChemicalFormula;
            }
            // If the full chemical formula field is null, it returns a null reference exception. i.e.,
            // you can't check to see if the field is null w/o getting an exception. This should be fixed at some point
            catch (NullReferenceException e)
            {
                peptideFormula = null;
            }

            if (peptideFormula == null || peptideFormula.AtomCount > 0)
            {
                // calculate averagine (used for isotopic distributions for unknown modifications)
                double averageC = 4.9384;
                double averageH = 7.7583;
                double averageO = 1.4773;
                double averageN = 1.3577;
                double averageS = 0.0417;

                // If averagine mass hasn't been calculated yet, it's stored
                // as private field _averagine mass
                double averagineMass = AveragineMass ??
                                       PeriodicTable.GetElement("C").AverageMass * averageC +
                                       PeriodicTable.GetElement("H").AverageMass * averageH +
                                       PeriodicTable.GetElement("O").AverageMass * averageO +
                                       PeriodicTable.GetElement("N").AverageMass * averageN +
                                       PeriodicTable.GetElement("S").AverageMass * averageS;
                AveragineMass ??= averagineMass;

                if (!String.IsNullOrEmpty(selectedPeptide.BaseSequence))
                {
                    Peptide baseSequence = new Peptide(selectedPeptide.BaseSequence);
                    peptideFormula = baseSequence.GetChemicalFormula();
                    var modifications = selectedPeptide.AllModsOneIsNterminus.Values;
                    foreach (Modification mod in modifications)
                    {
                        if (mod.ChemicalFormula != null)
                        {
                            peptideFormula.Add(mod.ChemicalFormula);
                        }
                        else if(mod.MonoisotopicMass != null)
                        {
                            massCorrection += (double)mod.MonoisotopicMass;
                        }
                    }

                    // add averagine for any unknown mass difference (i.e., modification w/o Chemical Formula)
                    massCorrection = Math.Max(peptideFormula.MonoisotopicMass - baseSequence.MonoisotopicMass, massCorrection);
                    if (Math.Abs(massCorrection) >= 20) // 20 Da difference is pulled directly from FlashLfq
                    {
                        double averagines = massCorrection / averagineMass;
                        ChemicalFormula unknownMod = new ChemicalFormula();

                        unknownMod.Add("C", (int)Math.Round(averagines * averageC, 0));
                        unknownMod.Add("H", (int)Math.Round(averagines * averageH, 0));
                        unknownMod.Add("O", (int)Math.Round(averagines * averageO, 0));
                        unknownMod.Add("N", (int)Math.Round(averagines * averageN, 0));
                        unknownMod.Add("S", (int)Math.Round(averagines * averageS, 0));

                        // Because averagines can't perfectly represent every modification mass,
                        // an additional correction factor is needed.
                        peptideFormula.Add(unknownMod);
                        massCorrection -= unknownMod.MonoisotopicMass;
                    }
                }
            }

            List<IndexedMassSpectralPeak> imsPeaksOrderedByMz = peak.IsotopicEnvelopes
                .Select(e => e.IndexedPeak)
                .Where(p => Math.Abs(p.RetentionTime - peak.Apex.IndexedPeak.RetentionTime) < 0.0001)
                .OrderBy(p => p.Mz)
                .ToList();
            int z = peak.Apex.ChargeState;

            if (peptideFormula != null & imsPeaksOrderedByMz.IsNotNullOrEmpty())
            {
                IsotopicDistribution peptideDistribution = IsotopicDistribution
                    .GetDistribution(peptideFormula, fineResolution, minimumProbability);
                SpectralSimilarity isotopeSimilarity = new(
                    imsPeaksOrderedByMz.Select(p => p.Mz).ToArray(),
                    imsPeaksOrderedByMz.Select(p => p.Intensity).ToArray(),
                    peptideDistribution.Masses.Select(m => (m + massCorrection).ToMz(z)).ToArray(),
                    peptideDistribution.Intensities.ToArray(),
                    SpectralSimilarity.SpectrumNormalizationScheme.spectrumSum,
                    toleranceInPpm: 20.0,
                    allPeaks: true,
                    filterOutBelowThisMz: 0);
                isotopeAngle = isotopeSimilarity.SpectralContrastAngle();
            }

            return isotopeAngle;
        }
    }
}
