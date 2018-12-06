using Chemistry;
using Proteomics.Fragmentation;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text.RegularExpressions;

namespace EngineLayer
{
    public class PsmFromTsv
    {
        private static readonly Regex IonParser = new Regex(@"([a-zA-Z]+)(\d+)");
        private static readonly char[] MzSplit = { '[', ',', ']', ';' };

        public string FullSequence { get; }
        public int Ms2ScanNumber { get; }
        public string Filename { get; }
        public int PrecursorScanNum { get; }
        public int PrecursorCharge { get; }
        public double PrecursorMz { get; }
        public double PrecursorMass { get; }
        public double Score { get; }
        public string ProteinAccession { get; }
        public List<MatchedFragmentIon> MatchedIons { get; }
        public double QValue { get; }

        public double? TotalIonCurrent { get; }
        public double? DeltaScore { get; }
        public string Notch { get; }
        public string BaseSeq { get; }
        public string EssentialSeq { get; }
        public string MissedCleavage { get; }
        public string PeptideMonoMass { get; }
        public string MassDiffDa { get; }
        public string MassDiffPpm { get; }
        public string ProteinName { get; }
        public string GeneName { get; }
        public string OrganismName { get; }
        public string PeptideDesicription { get; }
        public string StartAndEndResiduesInProtein { get; }
        public string PreviousAminoAcid { get; }
        public string NextAminoAcid { get; }
        public string DecoyContamTarget { get; }
        public double? QValueNotch { get; }

        //For crosslink

        public string CrossType { get; }
        public string LinkResidues { get; }
        public int? ProteinLinkSite { get; }
        public int? Rank { get; }
        public string BetaPeptideProteinAccession { get; }
        public int? BetaPeptideProteinLinkSite { get; }
        public string BetaPeptideBaseSequence { get; }
        public string BetaPeptideFullSequence { get; }
        public string BetaPeptideTheoreticalMass { get; }
        public double? BetaPeptideScore { get; }
        public int? BetaPeptideRank { get; }
        public List<MatchedFragmentIon> BetaPeptideMatchedIons { get; }
        public double? XLTotalScore { get; }
        public string ParentIons { get; }

        public PsmFromTsv(string line, char[] split, Dictionary<string, int> parsedHeader)
        {
            var spl = line.Split(split);

            //Required properties
            Filename = spl[parsedHeader[PsmTsvReader.FileName]].Trim();
            Ms2ScanNumber = int.Parse(spl[parsedHeader[PsmTsvReader.Ms2ScanNumber]]);
            PrecursorScanNum = int.Parse(spl[parsedHeader[PsmTsvReader.PrecursorScanNum]].Trim());
            PrecursorCharge = (int)double.Parse(spl[parsedHeader[PsmTsvReader.PrecursorCharge]].Trim());
            PrecursorMz = double.Parse(spl[parsedHeader[PsmTsvReader.PrecursorMz]].Trim());
            PrecursorMass = double.Parse(spl[parsedHeader[PsmTsvReader.PrecursorMass]].Trim());
            BaseSeq = spl[parsedHeader[PsmTsvReader.BaseSequence]].Trim();
            FullSequence = spl[parsedHeader[PsmTsvReader.FullSequence]];
            PeptideMonoMass = spl[parsedHeader[PsmTsvReader.PeptideMonoMass]].Trim();
            Score = double.Parse(spl[parsedHeader[PsmTsvReader.Score]].Trim());
            DecoyContamTarget = spl[parsedHeader[PsmTsvReader.DecoyContaminantTarget]].Trim();
            QValue = double.Parse(spl[parsedHeader[PsmTsvReader.QValue]].Trim());
            MatchedIons = ReadFragmentIonsFromString(spl[parsedHeader[PsmTsvReader.MatchedIonMzRatios]].Trim(), BaseSeq);

            //For general psms
            TotalIonCurrent = (parsedHeader[PsmTsvReader.TotalIonCurrent] < 0) ? null : (double?)double.Parse(spl[parsedHeader[PsmTsvReader.TotalIonCurrent]].Trim());
            DeltaScore = (parsedHeader[PsmTsvReader.DeltaScore] < 0) ? null : (double?)double.Parse(spl[parsedHeader[PsmTsvReader.DeltaScore]].Trim());
            Notch = (parsedHeader[PsmTsvReader.Notch] < 0) ? null : spl[parsedHeader[PsmTsvReader.Notch]].Trim();
            EssentialSeq = (parsedHeader[PsmTsvReader.EssentialSequence] < 0) ? null : spl[parsedHeader[PsmTsvReader.EssentialSequence]].Trim();
            MissedCleavage = (parsedHeader[PsmTsvReader.MissedCleavages] < 0) ? null : spl[parsedHeader[PsmTsvReader.MissedCleavages]].Trim();
            MassDiffDa = (parsedHeader[PsmTsvReader.MassDiffDa] < 0) ? null : spl[parsedHeader[PsmTsvReader.MassDiffDa]].Trim();
            MassDiffPpm = (parsedHeader[PsmTsvReader.MassDiffPpm] < 0) ? null : spl[parsedHeader[PsmTsvReader.MassDiffPpm]].Trim();
            ProteinAccession = (parsedHeader[PsmTsvReader.ProteinAccession] < 0) ? null : spl[parsedHeader[PsmTsvReader.ProteinAccession]].Trim();
            ProteinName = (parsedHeader[PsmTsvReader.ProteinName] < 0) ? null : spl[parsedHeader[PsmTsvReader.ProteinName]].Trim();
            GeneName = (parsedHeader[PsmTsvReader.GeneName] < 0) ? null : spl[parsedHeader[PsmTsvReader.GeneName]].Trim();
            OrganismName = (parsedHeader[PsmTsvReader.OrganismName] < 0) ? null : spl[parsedHeader[PsmTsvReader.OrganismName]].Trim();
            PeptideDesicription = (parsedHeader[PsmTsvReader.PeptideDesicription] < 0) ? null : spl[parsedHeader[PsmTsvReader.PeptideDesicription]].Trim();
            StartAndEndResiduesInProtein = (parsedHeader[PsmTsvReader.StartAndEndResiduesInProtein] < 0) ? null : spl[parsedHeader[PsmTsvReader.StartAndEndResiduesInProtein]].Trim();
            PreviousAminoAcid = (parsedHeader[PsmTsvReader.PreviousAminoAcid] < 0) ? null : spl[parsedHeader[PsmTsvReader.PreviousAminoAcid]].Trim();
            NextAminoAcid = (parsedHeader[PsmTsvReader.NextAminoAcid] < 0) ? null : spl[parsedHeader[PsmTsvReader.NextAminoAcid]].Trim();
            QValueNotch = (parsedHeader[PsmTsvReader.QValueNotch] < 0) ? null : (double?)double.Parse(spl[parsedHeader[PsmTsvReader.QValueNotch]].Trim());

            //For crosslinks
            CrossType = (parsedHeader[PsmTsvReader.CrossTypeLabel] < 0) ? null : spl[parsedHeader[PsmTsvReader.CrossTypeLabel]].Trim();
            LinkResidues = (parsedHeader[PsmTsvReader.LinkResiduesLabel] < 0) ? null : spl[parsedHeader[PsmTsvReader.LinkResiduesLabel]].Trim();
            ProteinLinkSite = (parsedHeader[PsmTsvReader.ProteinLinkSiteLabel] < 0) ? null : (int?)int.Parse(spl[parsedHeader[PsmTsvReader.ProteinLinkSiteLabel]].Trim());
            Rank = (parsedHeader[PsmTsvReader.RankLabel] < 0) ? null : (int?)int.Parse(spl[parsedHeader[PsmTsvReader.RankLabel]].Trim());
            BetaPeptideProteinAccession = (parsedHeader[PsmTsvReader.BetaPeptideProteinAccessionLabel] < 0) ? null : spl[parsedHeader[PsmTsvReader.BetaPeptideProteinAccessionLabel]].Trim();
            BetaPeptideProteinLinkSite = (parsedHeader[PsmTsvReader.BetaPeptideProteinLinkSiteLabel] < 0) ? null : (int?)int.Parse(spl[parsedHeader[PsmTsvReader.BetaPeptideProteinLinkSiteLabel]].Trim());
            BetaPeptideBaseSequence = (parsedHeader[PsmTsvReader.BetaPeptideBaseSequenceLabel] < 0) ? null : spl[parsedHeader[PsmTsvReader.BetaPeptideBaseSequenceLabel]].Trim();
            BetaPeptideFullSequence = (parsedHeader[PsmTsvReader.BetaPeptideFullSequenceLabel] < 0) ? null : spl[parsedHeader[PsmTsvReader.BetaPeptideFullSequenceLabel]].Trim();
            BetaPeptideTheoreticalMass = (parsedHeader[PsmTsvReader.BetaPeptideTheoreticalMassLabel] < 0) ? null : spl[parsedHeader[PsmTsvReader.BetaPeptideTheoreticalMassLabel]].Trim();
            BetaPeptideScore = (parsedHeader[PsmTsvReader.BetaPeptideScoreLabel] < 0) ? null : (double?)double.Parse(spl[parsedHeader[PsmTsvReader.BetaPeptideScoreLabel]].Trim());
            BetaPeptideRank = (parsedHeader[PsmTsvReader.BetaPeptideRankLabel] < 0) ? null : (int?)int.Parse(spl[parsedHeader[PsmTsvReader.BetaPeptideRankLabel]].Trim());
            BetaPeptideMatchedIons = (parsedHeader[PsmTsvReader.BetaPeptideMatchedIonsLabel] < 0) ? null : ReadFragmentIonsFromString(spl[parsedHeader[PsmTsvReader.BetaPeptideMatchedIonsLabel]].Trim(), BetaPeptideBaseSequence);
            XLTotalScore = (parsedHeader[PsmTsvReader.XLTotalScoreLabel] < 0) ? null : (double?)double.Parse(spl[parsedHeader[PsmTsvReader.XLTotalScoreLabel]].Trim());
            ParentIons = (parsedHeader[PsmTsvReader.ParentIonsLabel] < 0) ? null : spl[parsedHeader[PsmTsvReader.ParentIonsLabel]].Trim();
        }

        private static List<MatchedFragmentIon> ReadFragmentIonsFromString(string matchedMzString, string peptideBaseSequence)
        {
            var peaks = matchedMzString.Split(MzSplit, StringSplitOptions.RemoveEmptyEntries).Select(v => v.Trim())
                .ToList();
            peaks.RemoveAll(p => p.Contains("\""));

            List<MatchedFragmentIon> matchedIons = new List<MatchedFragmentIon>();

            foreach (var peak in peaks)
            {
                var split = peak.Split(new char[] { '+', ':' });

                string ionTypeAndNumber = split[0];
                Match result = IonParser.Match(ionTypeAndNumber);

                ProductType productType = (ProductType)Enum.Parse(typeof(ProductType), result.Groups[1].Value);

                int fragmentNumber = int.Parse(result.Groups[2].Value);
                int z = int.Parse(split[1]);
                double mz = double.Parse(split[2]);
                double neutralLoss = 0;

                // check for neutral loss
                if (ionTypeAndNumber.Contains("-"))
                {
                    string temp = ionTypeAndNumber.Replace("(", "");
                    temp = temp.Replace(")", "");
                    var split2 = temp.Split('-');
                    neutralLoss = double.Parse(split2[1]);
                }

                FragmentationTerminus terminus = FragmentationTerminus.None;
                if (TerminusSpecificProductTypes.ProductTypeToFragmentationTerminus.ContainsKey(productType))
                {
                    terminus = TerminusSpecificProductTypes.ProductTypeToFragmentationTerminus[productType];
                }

                int aminoAcidPosition = fragmentNumber;
                if (terminus == FragmentationTerminus.C)
                {
                    aminoAcidPosition = peptideBaseSequence.Length - fragmentNumber;
                }

                var t = new NeutralTerminusFragment(terminus, mz.ToMass(z) - DissociationTypeCollection.GetMassShiftFromProductType(productType), fragmentNumber, aminoAcidPosition);
                Product p = new Product(productType, t, neutralLoss);
                matchedIons.Add(new MatchedFragmentIon(p, mz, 1.0, z));
            }

            return matchedIons;
        }
    }
}