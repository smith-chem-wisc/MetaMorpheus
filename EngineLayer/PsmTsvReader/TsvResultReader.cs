using System;
using System.Collections.Generic;
using System.IO;

namespace EngineLayer
{
    public class PsmTsvReader
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
        public const string DifferentPeakMatches = "Different Peak Matches";

        // Sequence information

        public const string BaseSequence = "Base Sequence";
        public const string FullSequence = "Full Sequence";
        public const string EssentialSequence = "Essential Sequence";
        public const string PsmCount = "PSM Count";
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
        public const string AllScores = "All Scores";
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
        public const string EValue = "eValue";
        public const string EScore = "eScore";

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

        private static readonly char[] Split = { '\t' };

        public static List<PsmFromTsv> ReadTsv(string filePath, out List<string> warnings)
        {
            List<PsmFromTsv> psms = new List<PsmFromTsv>();
            warnings = new List<string>();

            StreamReader reader = null;
            try
            {
                reader = new StreamReader(filePath);
            }
            catch (Exception e)
            {
                throw new MetaMorpheusException("Could not read file: " + e.Message);
            }

            int lineCount = 0;

            string line;
            Dictionary<string, int> parsedHeader = null;

            while (reader.Peek() > 0)
            {
                lineCount++;

                line = reader.ReadLine();

                if (lineCount == 1)
                {
                    parsedHeader = ParseHeader(line);
                    continue;
                }

                try
                {
                    psms.Add(new PsmFromTsv(line, Split, parsedHeader));
                }
                catch (Exception)
                {
                    warnings.Add("Could not read line: " + lineCount);
                }
            }

            reader.Close();

            if ((lineCount - 1) != psms.Count)
            {
                warnings.Add("Warning: " + ((lineCount - 1) - psms.Count) + " PSMs were not read.");
            }

            return psms;
        }

        private static Dictionary<string, int> ParseHeader(string header)
        {
            var parsedHeader = new Dictionary<string, int>();
            var spl = header.Split(Split);

            parsedHeader.Add(FullSequence, Array.IndexOf(spl, FullSequence));
            parsedHeader.Add(Ms2ScanNumber, Array.IndexOf(spl, Ms2ScanNumber));
            parsedHeader.Add(FileName, Array.IndexOf(spl, FileName));
            parsedHeader.Add(TotalIonCurrent, Array.IndexOf(spl, TotalIonCurrent));
            parsedHeader.Add(PrecursorScanNum, Array.IndexOf(spl, PrecursorScanNum));
            parsedHeader.Add(PrecursorCharge, Array.IndexOf(spl, PrecursorCharge));
            parsedHeader.Add(PrecursorMz, Array.IndexOf(spl, PrecursorMz));
            parsedHeader.Add(PrecursorMass, Array.IndexOf(spl, PrecursorMass));
            parsedHeader.Add(Score, Array.IndexOf(spl, Score));
            parsedHeader.Add(DeltaScore, Array.IndexOf(spl, DeltaScore));
            parsedHeader.Add(Notch, Array.IndexOf(spl, Notch));
            parsedHeader.Add(BaseSequence, Array.IndexOf(spl, BaseSequence));
            parsedHeader.Add(EssentialSequence, Array.IndexOf(spl, EssentialSequence));
            parsedHeader.Add(MissedCleavages, Array.IndexOf(spl, MissedCleavages));
            parsedHeader.Add(PeptideMonoMass, Array.IndexOf(spl, PeptideMonoMass));
            parsedHeader.Add(MassDiffDa, Array.IndexOf(spl, MassDiffDa));
            parsedHeader.Add(MassDiffPpm, Array.IndexOf(spl, MassDiffPpm));
            parsedHeader.Add(ProteinAccession, Array.IndexOf(spl, ProteinAccession));
            parsedHeader.Add(ProteinName, Array.IndexOf(spl, ProteinName));
            parsedHeader.Add(GeneName, Array.IndexOf(spl, GeneName));
            parsedHeader.Add(OrganismName, Array.IndexOf(spl, OrganismName));
            parsedHeader.Add(PeptideDesicription, Array.IndexOf(spl, PeptideDesicription));
            parsedHeader.Add(StartAndEndResiduesInProtein, Array.IndexOf(spl, StartAndEndResiduesInProtein));
            parsedHeader.Add(PreviousAminoAcid, Array.IndexOf(spl, PreviousAminoAcid));
            parsedHeader.Add(NextAminoAcid, Array.IndexOf(spl, NextAminoAcid));
            parsedHeader.Add(DecoyContaminantTarget, Array.IndexOf(spl, DecoyContaminantTarget));
            parsedHeader.Add(MatchedIonMzRatios, Array.IndexOf(spl, MatchedIonMzRatios));
            parsedHeader.Add(QValue, Array.IndexOf(spl, QValue));
            parsedHeader.Add(QValueNotch, Array.IndexOf(spl, QValueNotch));

            parsedHeader.Add(CrossTypeLabel, Array.IndexOf(spl, CrossTypeLabel));
            parsedHeader.Add(LinkResiduesLabel, Array.IndexOf(spl, LinkResiduesLabel));
            parsedHeader.Add(ProteinLinkSiteLabel, Array.IndexOf(spl, ProteinLinkSiteLabel));
            parsedHeader.Add(RankLabel, Array.IndexOf(spl, RankLabel));
            parsedHeader.Add(BetaPeptideProteinAccessionLabel, Array.IndexOf(spl, BetaPeptideProteinAccessionLabel));
            parsedHeader.Add(BetaPeptideProteinLinkSiteLabel, Array.IndexOf(spl, BetaPeptideProteinLinkSiteLabel));
            parsedHeader.Add(BetaPeptideBaseSequenceLabel, Array.IndexOf(spl, BetaPeptideBaseSequenceLabel));
            parsedHeader.Add(BetaPeptideFullSequenceLabel, Array.IndexOf(spl, BetaPeptideFullSequenceLabel));
            parsedHeader.Add(BetaPeptideTheoreticalMassLabel, Array.IndexOf(spl, BetaPeptideTheoreticalMassLabel));
            parsedHeader.Add(BetaPeptideScoreLabel, Array.IndexOf(spl, BetaPeptideScoreLabel));
            parsedHeader.Add(BetaPeptideRankLabel, Array.IndexOf(spl, BetaPeptideRankLabel));
            parsedHeader.Add(BetaPeptideMatchedIonsLabel, Array.IndexOf(spl, BetaPeptideMatchedIonsLabel));
            parsedHeader.Add(XLTotalScoreLabel, Array.IndexOf(spl, XLTotalScoreLabel));
            parsedHeader.Add(ParentIonsLabel, Array.IndexOf(spl, ParentIonsLabel));

            return parsedHeader;
        }
    }
}