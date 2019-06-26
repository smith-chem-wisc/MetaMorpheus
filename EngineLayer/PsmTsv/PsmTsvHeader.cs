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
}