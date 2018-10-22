using System;
using System.Collections.Generic;
using System.Linq;
using System.IO;
using System.Text.RegularExpressions;
using System.Windows;
using MassSpectrometry;
using Chemistry;
using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using EngineLayer;

namespace MetaMorpheusGUI
{
    public class TsvResultReader
    {
        public const string FullSequenceLabel = "Full Sequence";
        public const string Ms2ScanNumberLabel = "Scan Number";
        public const string FilenameLabel = "File Name";
        public const string TotalIonCurrentLabel = "Total Ion Current";
        public const string PrecursorScanNumLabel = "Precursor Scan Number";
        public const string PrecursorChargeLabel = "Precursor Charge";
        public const string PrecursorMzLabel = "Precursor MZ";
        public const string PrecursorMassLabel = "Precursor Mass";
        public const string ScoreLabel = "Score";
        public const string DeltaScoreLabel = "Delta Score";
        public const string NotchLabel = "Notch";
        public const string BaseSeqLabel = "Base Sequence";
        public const string EssentialSeqLabel = "Essential Sequence";
        public const string MissedCleavageLabel = "Missed Cleavages";
        public const string PeptideMonoMassLabel = "Peptide Monoisotopic Mass";
        public const string MassDiffDaLabel = "Mass Diff (Da)";
        public const string MassDiffPpmLabel = "Mass Diff (ppm)";
        public const string ProteinAccessionLabel = "Protein Accession";
        public const string ProteinNameLabel = "Protein Name";
        public const string GeneNameLabel = "Gene Name";
        public const string SequenceVariationsLabel = "Sequence Variations";
        public const string OrganismNameLabel = "Organism Name";
        public const string PeptideDesicriptionLabel = "Peptide Description";
        public const string StartAndEndResiduesInProteinLabel = "Start and End Residues In Protein";
        public const string PreviousAminoAcidLabel = "Previous Amino Acid";
        public const string NextAminoAcidLabel = "Next Amino Acid";
        public const string DecoyContamTargetLabel = "Decoy/Contaminant/Target";
        public const string MatchedIonsLabel = "Matched Ion Mass-To-Charge Ratios";
        public const string QValueLabel = "QValue";
        public const string QValueNotchLabel = "QValue Notch";
        
        private static readonly char[] Split = { '\t' };
        
        public static List<MetaDrawPsm> ReadTsv(string filePath)
        {
            List<MetaDrawPsm> psms = new List<MetaDrawPsm>();

            StreamReader reader = null;
            try
            {
                reader = new StreamReader(filePath);
            }
            catch (Exception e)
            {
                MessageBox.Show("Could not read file: " + e.Message);
                return psms;
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

                if (parsedHeader.Values.Any(i => i < 0))
                {
                    MessageBox.Show("Could not read PSMs file. Is it from an older version of MetaMorpheus?");
                    return psms;
                }

                try
                {
                    psms.Add(new MetaDrawPsm(line, Split, parsedHeader));
                }
                catch (Exception)
                {
                    MessageBox.Show("Could not read line: " + lineCount);
                }
            }

            reader.Close();

            if ((lineCount - 1) != psms.Count)
            {
                MessageBox.Show("Warning: " + ((lineCount - 1) - psms.Count) + " PSMs were not read.");
            }

            return psms;
        }
        
        private static Dictionary<string, int> ParseHeader(string header)
        {
            var parsedHeader = new Dictionary<string, int>();
            var spl = header.Split(Split);

            parsedHeader.Add(FullSequenceLabel, Array.IndexOf(spl, FullSequenceLabel));
            parsedHeader.Add(Ms2ScanNumberLabel, Array.IndexOf(spl, Ms2ScanNumberLabel));
            parsedHeader.Add(FilenameLabel, Array.IndexOf(spl, FilenameLabel));
            parsedHeader.Add(TotalIonCurrentLabel, Array.IndexOf(spl, TotalIonCurrentLabel));
            parsedHeader.Add(PrecursorScanNumLabel, Array.IndexOf(spl, PrecursorScanNumLabel));
            parsedHeader.Add(PrecursorChargeLabel, Array.IndexOf(spl, PrecursorChargeLabel));
            parsedHeader.Add(PrecursorMzLabel, Array.IndexOf(spl, PrecursorMzLabel));
            parsedHeader.Add(PrecursorMassLabel, Array.IndexOf(spl, PrecursorMassLabel));
            parsedHeader.Add(ScoreLabel, Array.IndexOf(spl, ScoreLabel));
            parsedHeader.Add(DeltaScoreLabel, Array.IndexOf(spl, DeltaScoreLabel));
            parsedHeader.Add(NotchLabel, Array.IndexOf(spl, NotchLabel));
            parsedHeader.Add(BaseSeqLabel, Array.IndexOf(spl, BaseSeqLabel));
            parsedHeader.Add(EssentialSeqLabel, Array.IndexOf(spl, EssentialSeqLabel));
            parsedHeader.Add(MissedCleavageLabel, Array.IndexOf(spl, MissedCleavageLabel));
            parsedHeader.Add(PeptideMonoMassLabel, Array.IndexOf(spl, PeptideMonoMassLabel));
            parsedHeader.Add(MassDiffDaLabel, Array.IndexOf(spl, MassDiffDaLabel));
            parsedHeader.Add(MassDiffPpmLabel, Array.IndexOf(spl, MassDiffPpmLabel));
            parsedHeader.Add(ProteinAccessionLabel, Array.IndexOf(spl, ProteinAccessionLabel));
            parsedHeader.Add(ProteinNameLabel, Array.IndexOf(spl, ProteinNameLabel));
            parsedHeader.Add(GeneNameLabel, Array.IndexOf(spl, GeneNameLabel));
            parsedHeader.Add(SequenceVariationsLabel, Array.IndexOf(spl, SequenceVariationsLabel));
            parsedHeader.Add(OrganismNameLabel, Array.IndexOf(spl, OrganismNameLabel));
            parsedHeader.Add(PeptideDesicriptionLabel, Array.IndexOf(spl, PeptideDesicriptionLabel));
            parsedHeader.Add(StartAndEndResiduesInProteinLabel, Array.IndexOf(spl, StartAndEndResiduesInProteinLabel));
            parsedHeader.Add(PreviousAminoAcidLabel, Array.IndexOf(spl, PreviousAminoAcidLabel));
            parsedHeader.Add(NextAminoAcidLabel, Array.IndexOf(spl, NextAminoAcidLabel));
            parsedHeader.Add(DecoyContamTargetLabel, Array.IndexOf(spl, DecoyContamTargetLabel));
            parsedHeader.Add(MatchedIonsLabel, Array.IndexOf(spl, MatchedIonsLabel));
            parsedHeader.Add(QValueLabel, Array.IndexOf(spl, QValueLabel));
            parsedHeader.Add(QValueNotchLabel, Array.IndexOf(spl, QValueNotchLabel));

            return parsedHeader;
        }
    }
}