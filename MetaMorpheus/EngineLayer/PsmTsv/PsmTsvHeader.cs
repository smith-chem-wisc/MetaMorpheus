using System;
using System.Globalization;

namespace EngineLayer
{
    public static class PsmTsvHeader
    {
        // File and scan information
        public const string FileName = "File Name";
        public const string Ms2ScanNumber = "Scan Number";
        public const string Ms2ScanRetentionTime = "Scan Retention Time";
        public const string NumExperimentalPeaks = "Num Experimental Peaks";
        public const string TotalIonCurrent = "Total Ion Current";
        public const string PrecursorScanNum = "Precursor Scan Number";
        public const string PrecursorCharge = "Precursor Charge";
        public const string PrecursorMz = "Precursor MZ";
        public const string PrecursorMass = "Precursor Mass";
        public const string Score = "Score";
        public const string DeltaScore = "Delta Score";
        public const string Notch = "Notch";

        // Sequence information
        public const string BaseSequence = "Base Sequence";
        public const string FullSequence = "Full Sequence";
        public const string EssentialSequence = "Essential Sequence";
        public const string AmbiguityLevel = "Ambiguity Level";
        public const string PsmCount = "PSM Count (unambiguous, <0.01 q-value)";
        public const string Mods = "Mods";
        public const string ModsChemicalFormulas = "Mods Chemical Formulas";
        public const string ModsCombinedChemicalFormula = "Mods Combined Chemical Formula";
        public const string NumVariableMods = "Num Variable Mods";
        public const string MissedCleavages = "Missed Cleavages";
        public const string PeptideMonoMass = "Peptide Monoisotopic Mass";
        public const string MassDiffDa = "Mass Diff (Da)";
        public const string MassDiffPpm = "Mass Diff (ppm)";
        public const string ProteinAccession = "Protein Accession";
        public const string ProteinName = "Protein Name";
        public const string GeneName = "Gene Name";
        public const string OrganismName = "Organism Name";
        public const string IntersectingSequenceVariations = "Intersecting Sequence Variations";
        public const string IdentifiedSequenceVariations = "Identified Sequence Variations";
        public const string SpliceSites = "Splice Sites";
        public const string Contaminant = "Contaminant";
        public const string Decoy = "Decoy";
        public const string PeptideDesicription = "Peptide Description";
        public const string StartAndEndResiduesInProtein = "Start and End Residues In Protein";
        public const string PreviousAminoAcid = "Previous Amino Acid";
        public const string NextAminoAcid = "Next Amino Acid";
        public const string TheoreticalsSearched = "Theoreticals Searched";
        public const string DecoyContaminantTarget = "Decoy/Contaminant/Target";
        public const string MatchedIonSeries = "Matched Ion Series";
        public const string MatchedIonMzRatios = "Matched Ion Mass-To-Charge Ratios";
        public const string MatchedIonMassDiffDa = "Matched Ion Mass Diff (Da)";
        public const string MatchedIonMassDiffPpm = "Matched Ion Mass Diff (Ppm)";
        public const string MatchedIonIntensities = "Matched Ion Intensities";
        public const string MatchedIonCounts = "Matched Ion Counts";

        // Scoring
        public const string LocalizedScores = "Localized Scores";
        public const string ImprovementPossible = "Improvement Possible";
        public const string CumulativeTarget = "Cumulative Target";
        public const string CumulativeDecoy = "Cumulative Decoy";
        public const string CumulativeTargetNotch = "Cumulative Target Notch";
        public const string CumulativeDecoyNotch = "Cumulative Decoy Notch";
        public const string QValue = "QValue";
        public const string QValueNotch = "QValue Notch";
        public const string PEP = "PEP";
        public const string PEP_QValue = "PEP_QValue";
        public const string SpectralAngle = "Normalized Spectral Angle";

        // Crosslinks
        public const string CrossTypeLabel = "Cross Type";
        public const string LinkResiduesLabel = "Link Residues";
        public const string ProteinLinkSiteLabel = "Protein Link Site";
        public const string RankLabel = "Rank";
        public const string BetaPeptideProteinAccessionLabel = "Beta Peptide Protein Accession";
        public const string BetaPeptideProteinLinkSiteLabel = "Beta Peptide Protein LinkSite";
        public const string BetaPeptideBaseSequenceLabel = "Beta Peptide Base Sequence";
        public const string BetaPeptideFullSequenceLabel = "Beta Peptide Full Sequence";
        public const string BetaPeptideTheoreticalMassLabel = "Beta Peptide Theoretical Mass";
        public const string BetaPeptideScoreLabel = "Beta Peptide Score";
        public const string BetaPeptideRankLabel = "Beta Peptide Rank";
        public const string BetaPeptideMatchedIonsLabel = "Beta Peptide Matched Ion Mass-To-Charge Ratios";
        public const string BetaPeptideMatchedIonIntensitiesLabel = "Beta Peptide Matched Ion Intensities";
        public const string XLTotalScoreLabel = "XL Total Score";
        public const string ParentIonsLabel = "Parent Ions";
    }

    public static class PsmTsvHeader_Cross
    {

        // Crosslinks
        public const string CrossTypeLabel = "Cross Type";
        public const string LinkResiduesLabel = "Link Residues";
        public const string ProteinLinkSiteLabel = "Protein Link Site";
        public const string RankLabel = "Rank";
        public const string BetaPeptideProteinAccessionLabel = "Beta Peptide Protein Accession";
        public const string BetaPeptideProteinLinkSiteLabel = "Beta Peptide Protein LinkSite";
        public const string BetaPeptideBaseSequenceLabel = "Beta Peptide Base Sequence";
        public const string BetaPeptideFullSequenceLabel = "Beta Peptide Full Sequence";
        public const string BetaPeptideTheoreticalMassLabel = "Beta Peptide Theoretical Mass";
        public const string BetaPeptideScoreLabel = "Beta Peptide Score";
        public const string BetaPeptideRankLabel = "Beta Peptide Rank";
        public const string BetaPeptideMatchedIonsLabel = "Beta Peptide Matched Ion Mass-To-Charge Ratios";
        public const string XLTotalScoreLabel = "XL Total Score";
        public const string ParentIonsLabel = "Parent Ions";
    }

    public static class PsmTsvHeader_Glyco
    {
        public const string GlycanMass = "GlycanMass";
        public const string GlycanComposition = "Plausible GlycanComposition";
        public const string GlycanStructure = "Plausible GlycanStructure";
        public const string GlycanLocalizationLevel = "GlycanLocalizationLevel";
        public const string LocalizedGlycan = "Localized Glycans with Peptide Site Specific Probability";
    }

    public static class PsmTsvHeader_SpectralRecovery
    {
        public const string PeakApexRt = "Peak Apex RT (min)";
        public const string PeakShift = "RT Shift"; // Acceptor peak apex RT - predicted RT
        public const string PrecursorDeconvolutedBool = "Deconvolutable Precursor";
        public const string PrecursorIsotopicEnvelopeAngle = "Precursor Isotopic Envelope Score"; // Spectral Contrast angle exp vs theoretical
        public const string IsolationWindowCenter = "Isolation Window Center (Th)";
        public const string PrecursorOffset = "Precursor m/z - Isolation Center Distance (Th)";
        public const string IsolationWindowWidth = "Isolation Window Width (Th)";
        public const string OriginalPsmQ = "Original Psm QValue";
        public const string OriginalPsmPEP = "Original Psm PEP";
        public const string OriginalPsmPEP_QValue = "Original Psm PEP_QValue";

        public static string NullableToString<T>(this Nullable<T> value, CultureInfo cultureInfo)
            where T : struct, IConvertible
        {
            return value == null ? "" : ((IConvertible)value).ToString(cultureInfo);
        }
    }

    public class MaxQuantMsmsHeader
    {
        // File and scan information
        public const string FileName = "Raw file";
        public const string Ms2ScanNumber = "Scan number";
        public const string Ms2ScanRetentionTime = "Retention time";
        public const string PrecursorScanNum = "Precursor full scan number";
        public const string PrecursorCharge = "Charge";
        public const string PrecursorMz = "m/z";
        public const string PrecursorMass = "Mass";
        public const string Score = "Score";
        public const string DeltaScore = "Delta score";

        // Sequence information
        public const string BaseSequence = "Sequence";
        public const string FullSequence = "Modified sequence";
        public const string Mods = "Modifications";
        public const string MassDiffDa = "Mass error [Da]";
        public const string MassDiffPpm = "Mass error [ppm]";
        public const string ProteinAccession = "Proteins";
        public const string ProteinName = "Protein Names";
        public const string GeneName = "Gene Names";
        public const string Decoy = "Reverse";
        public const string MatchedIonSeries = "Matches";
        public const string MatchedIonMzRatios = "Masses";
        public const string MatchedIonMassDiffDa = "Mass deviations [Da]";
        public const string MatchedIonMassDiffPpm = "Mass deviations [ppm]";
        public const string MatchedIonIntensities = "Intensities";
        public const string MatchedIonCounts = "Number of matches";

        // Scoring
        public const string PEP = "PEP";
    }
}