using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace EngineLayer
{
    public class MaxQuantMsmsHeader
    {
        // File and scan information
        public const string FileName = "Raw file";
        public const string Ms2ScanNumber = "Scan number";
        public const string Ms2ScanRetentionTime = "Retention time";
        //public const string NumExperimentalPeaks = "Number of Matches";
        //public const string TotalIonCurrent = "Total Ion Current";
        public const string PrecursorScanNum = "Precursor full scan number";
        public const string PrecursorCharge = "Charge";
        public const string PrecursorMz = "m/z";
        public const string PrecursorMass = "Mass";
        public const string Score = "Score";
        public const string DeltaScore = "Delta score";
        //public const string Notch = "Notch";

        // Sequence information
        public const string BaseSequence = "Sequence";
        public const string FullSequence = "Modified sequence";
        //public const string EssentialSequence = "Essential Sequence";
        //public const string AmbiguityLevel = "Ambiguity Level";
        //public const string PsmCount = "PSM Count (unambiguous, <0.01 q-value)";
        public const string Mods = "Modifications";
        //public const string ModsChemicalFormulas = "Mods Chemical Formulas";
        //public const string ModsCombinedChemicalFormula = "Mods Combined Chemical Formula";
        //public const string NumVariableMods = "Num Variable Mods";
        //public const string MissedCleavages = "Missed Cleavages";
        //public const string PeptideMonoMass = "Peptide Monoisotopic Mass";
        public const string MassDiffDa = "Mass error [Da]";
        public const string MassDiffPpm = "Mass error [ppm]";
        public const string ProteinAccession = "Proteins";
        public const string ProteinName = "Protein Names";
        public const string GeneName = "Gene Names";
        //public const string OrganismName = "Organism Name";
        //public const string IntersectingSequenceVariations = "Intersecting Sequence Variations";
        //public const string IdentifiedSequenceVariations = "Identified Sequence Variations";
        //public const string SpliceSites = "Splice Sites";
        //public const string Contaminant = "Contaminant";
        public const string Decoy = "Reverse";
        //public const string PeptideDesicription = "Peptide Description";
        //public const string StartAndEndResiduesInProtein = "Start and End Residues In Protein";
        //public const string PreviousAminoAcid = "Previous Amino Acid";
        //public const string NextAminoAcid = "Next Amino Acid";
        //public const string TheoreticalsSearched = "Theoreticals Searched";
        //public const string DecoyContaminantTarget = "Decoy/Contaminant/Target";
        public const string MatchedIonSeries = "Matches";
        public const string MatchedIonMzRatios = "Masses";
        public const string MatchedIonMassDiffDa = "Mass deviations [Da]";
        public const string MatchedIonMassDiffPpm = "Mass deviations [ppm]";
        public const string MatchedIonIntensities = "Intensities";
        public const string MatchedIonCounts = "Number of matches";

        // Scoring
        public const string PEP = "PEP";
        //public const string LocalizedScores = "Localized Scores";
        //public const string ImprovementPossible = "Improvement Possible";
        //public const string CumulativeTarget = "Cumulative Target";
        //public const string CumulativeDecoy = "Cumulative Decoy";
        //public const string CumulativeTargetNotch = "Cumulative Target Notch";
        //public const string CumulativeDecoyNotch = "Cumulative Decoy Notch";
        //public const string QValue = "QValue";
        //public const string QValueNotch = "QValue Notch";
        //public const string PEP_QValue = "PEP_QValue";
        //public const string SpectralAngle = "Normalized Spectral Angle";


    }
}
