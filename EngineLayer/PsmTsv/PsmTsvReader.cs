using System;
using System.Collections.Generic;
using System.IO;

namespace EngineLayer
{
    public class PsmTsvReader
    {
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

            parsedHeader.Add(PsmTsvHeader.FullSequence, Array.IndexOf(spl, PsmTsvHeader.FullSequence));
            parsedHeader.Add(PsmTsvHeader.Ms2ScanNumber, Array.IndexOf(spl, PsmTsvHeader.Ms2ScanNumber));
            parsedHeader.Add(PsmTsvHeader.FileName, Array.IndexOf(spl, PsmTsvHeader.FileName));
            parsedHeader.Add(PsmTsvHeader.TotalIonCurrent, Array.IndexOf(spl, PsmTsvHeader.TotalIonCurrent));
            parsedHeader.Add(PsmTsvHeader.PrecursorScanNum, Array.IndexOf(spl, PsmTsvHeader.PrecursorScanNum));
            parsedHeader.Add(PsmTsvHeader.PrecursorCharge, Array.IndexOf(spl, PsmTsvHeader.PrecursorCharge));
            parsedHeader.Add(PsmTsvHeader.PrecursorMz, Array.IndexOf(spl, PsmTsvHeader.PrecursorMz));
            parsedHeader.Add(PsmTsvHeader.PrecursorMass, Array.IndexOf(spl, PsmTsvHeader.PrecursorMass));
            parsedHeader.Add(PsmTsvHeader.Score, Array.IndexOf(spl, PsmTsvHeader.Score));
            parsedHeader.Add(PsmTsvHeader.DeltaScore, Array.IndexOf(spl, PsmTsvHeader.DeltaScore));
            parsedHeader.Add(PsmTsvHeader.Notch, Array.IndexOf(spl, PsmTsvHeader.Notch));
            parsedHeader.Add(PsmTsvHeader.BaseSequence, Array.IndexOf(spl, PsmTsvHeader.BaseSequence));
            parsedHeader.Add(PsmTsvHeader.EssentialSequence, Array.IndexOf(spl, PsmTsvHeader.EssentialSequence));
            parsedHeader.Add(PsmTsvHeader.MissedCleavages, Array.IndexOf(spl, PsmTsvHeader.MissedCleavages));
            parsedHeader.Add(PsmTsvHeader.PeptideMonoMass, Array.IndexOf(spl, PsmTsvHeader.PeptideMonoMass));
            parsedHeader.Add(PsmTsvHeader.MassDiffDa, Array.IndexOf(spl, PsmTsvHeader.MassDiffDa));
            parsedHeader.Add(PsmTsvHeader.MassDiffPpm, Array.IndexOf(spl, PsmTsvHeader.MassDiffPpm));
            parsedHeader.Add(PsmTsvHeader.ProteinAccession, Array.IndexOf(spl, PsmTsvHeader.ProteinAccession));
            parsedHeader.Add(PsmTsvHeader.ProteinName, Array.IndexOf(spl, PsmTsvHeader.ProteinName));
            parsedHeader.Add(PsmTsvHeader.GeneName, Array.IndexOf(spl, PsmTsvHeader.GeneName));
            parsedHeader.Add(PsmTsvHeader.OrganismName, Array.IndexOf(spl, PsmTsvHeader.OrganismName));
            parsedHeader.Add(PsmTsvHeader.IntersectingSequenceVariations, Array.IndexOf(spl, PsmTsvHeader.IntersectingSequenceVariations));
            parsedHeader.Add(PsmTsvHeader.IdentifiedSequenceVariations, Array.IndexOf(spl, PsmTsvHeader.IdentifiedSequenceVariations));
            parsedHeader.Add(PsmTsvHeader.SpliceSites, Array.IndexOf(spl, PsmTsvHeader.SpliceSites));
            parsedHeader.Add(PsmTsvHeader.PeptideDesicription, Array.IndexOf(spl, PsmTsvHeader.PeptideDesicription));
            parsedHeader.Add(PsmTsvHeader.StartAndEndResiduesInProtein, Array.IndexOf(spl, PsmTsvHeader.StartAndEndResiduesInProtein));
            parsedHeader.Add(PsmTsvHeader.PreviousAminoAcid, Array.IndexOf(spl, PsmTsvHeader.PreviousAminoAcid));
            parsedHeader.Add(PsmTsvHeader.NextAminoAcid, Array.IndexOf(spl, PsmTsvHeader.NextAminoAcid));
            parsedHeader.Add(PsmTsvHeader.DecoyContaminantTarget, Array.IndexOf(spl, PsmTsvHeader.DecoyContaminantTarget));
            parsedHeader.Add(PsmTsvHeader.MatchedIonMzRatios, Array.IndexOf(spl, PsmTsvHeader.MatchedIonMzRatios));
            parsedHeader.Add(PsmTsvHeader.QValue, Array.IndexOf(spl, PsmTsvHeader.QValue));
            parsedHeader.Add(PsmTsvHeader.QValueNotch, Array.IndexOf(spl, PsmTsvHeader.QValueNotch));
            parsedHeader.Add(PsmTsvHeader.PEP, Array.IndexOf(spl, PsmTsvHeader.PEP));
            parsedHeader.Add(PsmTsvHeader.PEP_QValue, Array.IndexOf(spl, PsmTsvHeader.PEP_QValue));

            parsedHeader.Add(PsmTsvHeader.CrossTypeLabel, Array.IndexOf(spl, PsmTsvHeader.CrossTypeLabel));
            parsedHeader.Add(PsmTsvHeader.LinkResiduesLabel, Array.IndexOf(spl, PsmTsvHeader.LinkResiduesLabel));
            parsedHeader.Add(PsmTsvHeader.ProteinLinkSiteLabel, Array.IndexOf(spl, PsmTsvHeader.ProteinLinkSiteLabel));
            parsedHeader.Add(PsmTsvHeader.RankLabel, Array.IndexOf(spl, PsmTsvHeader.RankLabel));
            parsedHeader.Add(PsmTsvHeader.BetaPeptideProteinAccessionLabel, Array.IndexOf(spl, PsmTsvHeader.BetaPeptideProteinAccessionLabel));
            parsedHeader.Add(PsmTsvHeader.BetaPeptideProteinLinkSiteLabel, Array.IndexOf(spl, PsmTsvHeader.BetaPeptideProteinLinkSiteLabel));
            parsedHeader.Add(PsmTsvHeader.BetaPeptideBaseSequenceLabel, Array.IndexOf(spl, PsmTsvHeader.BetaPeptideBaseSequenceLabel));
            parsedHeader.Add(PsmTsvHeader.BetaPeptideFullSequenceLabel, Array.IndexOf(spl, PsmTsvHeader.BetaPeptideFullSequenceLabel));
            parsedHeader.Add(PsmTsvHeader.BetaPeptideTheoreticalMassLabel, Array.IndexOf(spl, PsmTsvHeader.BetaPeptideTheoreticalMassLabel));
            parsedHeader.Add(PsmTsvHeader.BetaPeptideScoreLabel, Array.IndexOf(spl, PsmTsvHeader.BetaPeptideScoreLabel));
            parsedHeader.Add(PsmTsvHeader.BetaPeptideRankLabel, Array.IndexOf(spl, PsmTsvHeader.BetaPeptideRankLabel));
            parsedHeader.Add(PsmTsvHeader.BetaPeptideMatchedIonsLabel, Array.IndexOf(spl, PsmTsvHeader.BetaPeptideMatchedIonsLabel));
            parsedHeader.Add(PsmTsvHeader.XLTotalScoreLabel, Array.IndexOf(spl, PsmTsvHeader.XLTotalScoreLabel));
            parsedHeader.Add(PsmTsvHeader.ParentIonsLabel, Array.IndexOf(spl, PsmTsvHeader.ParentIonsLabel));
            parsedHeader.Add(PsmTsvHeader.Ms2ScanRetentionTime, Array.IndexOf(spl, PsmTsvHeader.Ms2ScanRetentionTime));

            parsedHeader.Add(PsmTsvHeader_Glyco.GlycanIDs, Array.IndexOf(spl, PsmTsvHeader_Glyco.GlycanIDs));
            parsedHeader.Add(PsmTsvHeader_Glyco.GlycanMass, Array.IndexOf(spl, PsmTsvHeader_Glyco.GlycanMass));
            parsedHeader.Add(PsmTsvHeader_Glyco.GlycanStructure, Array.IndexOf(spl, PsmTsvHeader_Glyco.GlycanStructure));
            parsedHeader.Add(PsmTsvHeader_Glyco.GlycanComposition, Array.IndexOf(spl, PsmTsvHeader_Glyco.GlycanComposition));

            return parsedHeader;
        }
    }
}