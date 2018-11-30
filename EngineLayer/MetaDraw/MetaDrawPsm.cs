using Proteomics.Fragmentation;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text.RegularExpressions;
using Chemistry;

namespace EngineLayer
{
    public class MetaDrawPsm
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
        public string BetaPeptideProteinAccession{get;}
        public int? BetaPeptideProteinLinkSite { get; }
        public string BetaPeptideBaseSequence { get; }
        public string BetaPeptideFullSequence { get; }
        public string BetaPeptideTheoreticalMass { get; }
        public double? BetaPeptideScore { get; }
        public int? BetaPeptideRank { get; }
        public List<MatchedFragmentIon> BetaPeptideMatchedIons { get; }
        public double? XLTotalScore { get; }
        public string ParentIons { get; }

        public MetaDrawPsm(string line, char[] split, Dictionary<string, int> parsedHeader)
        {
            var spl = line.Split(split);

            //Required properties
            Filename = spl[parsedHeader[TsvResultReader.FilenameLabel]].Trim();
            Ms2ScanNumber = int.Parse(spl[parsedHeader[TsvResultReader.Ms2ScanNumberLabel]]);
            PrecursorScanNum = int.Parse(spl[parsedHeader[TsvResultReader.PrecursorScanNumLabel]].Trim());
            PrecursorCharge = (int)double.Parse(spl[parsedHeader[TsvResultReader.PrecursorChargeLabel]].Trim());
            PrecursorMz = double.Parse(spl[parsedHeader[TsvResultReader.PrecursorMzLabel]].Trim());
            PrecursorMass = double.Parse(spl[parsedHeader[TsvResultReader.PrecursorMassLabel]].Trim());
            BaseSeq = spl[parsedHeader[TsvResultReader.BaseSeqLabel]].Trim();
            FullSequence = spl[parsedHeader[TsvResultReader.FullSequenceLabel]];           
            PeptideMonoMass = spl[parsedHeader[TsvResultReader.PeptideMonoMassLabel]].Trim();
            Score = double.Parse(spl[parsedHeader[TsvResultReader.ScoreLabel]].Trim());
            DecoyContamTarget = spl[parsedHeader[TsvResultReader.DecoyContamTargetLabel]].Trim();
            QValue = double.Parse(spl[parsedHeader[TsvResultReader.QValueLabel]].Trim());
            MatchedIons = ReadFragmentIonsFromString(spl[parsedHeader[TsvResultReader.MatchedIonsLabel]].Trim(), BaseSeq);

            //For general psms
            TotalIonCurrent = (parsedHeader[TsvResultReader.TotalIonCurrentLabel] < 0) ? null : (double?)double.Parse(spl[parsedHeader[TsvResultReader.TotalIonCurrentLabel]].Trim());            
            DeltaScore = (parsedHeader[TsvResultReader.DeltaScoreLabel] < 0) ? null : (double?)double.Parse(spl[parsedHeader[TsvResultReader.DeltaScoreLabel]].Trim());
            Notch = (parsedHeader[TsvResultReader.NotchLabel] < 0) ? null : spl[parsedHeader[TsvResultReader.NotchLabel]].Trim();          
            EssentialSeq = (parsedHeader[TsvResultReader.EssentialSeqLabel] < 0) ? null : spl[parsedHeader[TsvResultReader.EssentialSeqLabel]].Trim();
            MissedCleavage = (parsedHeader[TsvResultReader.MissedCleavageLabel] < 0) ? null : spl[parsedHeader[TsvResultReader.MissedCleavageLabel]].Trim();        
            MassDiffDa = (parsedHeader[TsvResultReader.MassDiffDaLabel] < 0) ? null : spl[parsedHeader[TsvResultReader.MassDiffDaLabel]].Trim();
            MassDiffPpm = (parsedHeader[TsvResultReader.MassDiffPpmLabel] < 0) ? null : spl[parsedHeader[TsvResultReader.MassDiffPpmLabel]].Trim();
            ProteinAccession = (parsedHeader[TsvResultReader.ProteinAccessionLabel] < 0) ? null : spl[parsedHeader[TsvResultReader.ProteinAccessionLabel]].Trim();
            ProteinName = (parsedHeader[TsvResultReader.ProteinNameLabel] < 0) ? null : spl[parsedHeader[TsvResultReader.ProteinNameLabel]].Trim();
            GeneName = (parsedHeader[TsvResultReader.GeneNameLabel] < 0) ? null : spl[parsedHeader[TsvResultReader.GeneNameLabel]].Trim();
            OrganismName = (parsedHeader[TsvResultReader.OrganismNameLabel] < 0) ? null : spl[parsedHeader[TsvResultReader.OrganismNameLabel]].Trim();
            PeptideDesicription = (parsedHeader[TsvResultReader.PeptideDesicriptionLabel] < 0) ? null : spl[parsedHeader[TsvResultReader.PeptideDesicriptionLabel]].Trim();
            StartAndEndResiduesInProtein = (parsedHeader[TsvResultReader.StartAndEndResiduesInProteinLabel] < 0) ? null : spl[parsedHeader[TsvResultReader.StartAndEndResiduesInProteinLabel]].Trim();
            PreviousAminoAcid = (parsedHeader[TsvResultReader.PreviousAminoAcidLabel] < 0) ? null : spl[parsedHeader[TsvResultReader.PreviousAminoAcidLabel]].Trim();
            NextAminoAcid = (parsedHeader[TsvResultReader.NextAminoAcidLabel] < 0) ? null : spl[parsedHeader[TsvResultReader.NextAminoAcidLabel]].Trim();      
            QValueNotch = (parsedHeader[TsvResultReader.QValueNotchLabel] < 0) ? null : (double?)double.Parse(spl[parsedHeader[TsvResultReader.QValueNotchLabel]].Trim());         

            //For crosslinks
            CrossType = (parsedHeader[TsvResultReader.CrossTypeLabel] < 0) ? null : spl[parsedHeader[TsvResultReader.CrossTypeLabel]].Trim();
            LinkResidues = (parsedHeader[TsvResultReader.LinkResiduesLabel] < 0) ? null : spl[parsedHeader[TsvResultReader.LinkResiduesLabel]].Trim();
            ProteinLinkSite = (parsedHeader[TsvResultReader.ProteinLinkSiteLabel] < 0) ? null : (int?)int.Parse(spl[parsedHeader[TsvResultReader.ProteinLinkSiteLabel]].Trim());
            Rank = (parsedHeader[TsvResultReader.RankLabel] < 0) ? null : (int?)int.Parse(spl[parsedHeader[TsvResultReader.RankLabel]].Trim());
            BetaPeptideProteinAccession = (parsedHeader[TsvResultReader.BetaPeptideProteinAccessionLabel] < 0) ? null : spl[parsedHeader[TsvResultReader.BetaPeptideProteinAccessionLabel]].Trim();
            BetaPeptideProteinLinkSite = (parsedHeader[TsvResultReader.BetaPeptideProteinLinkSiteLabel] < 0) ? null : (int?)int.Parse(spl[parsedHeader[TsvResultReader.BetaPeptideProteinLinkSiteLabel]].Trim());
            BetaPeptideBaseSequence = (parsedHeader[TsvResultReader.BetaPeptideBaseSequenceLabel] < 0) ? null : spl[parsedHeader[TsvResultReader.BetaPeptideBaseSequenceLabel]].Trim();
            BetaPeptideFullSequence = (parsedHeader[TsvResultReader.BetaPeptideFullSequenceLabel] < 0) ? null : spl[parsedHeader[TsvResultReader.BetaPeptideFullSequenceLabel]].Trim();
            BetaPeptideTheoreticalMass = (parsedHeader[TsvResultReader.BetaPeptideTheoreticalMassLabel] < 0) ? null : spl[parsedHeader[TsvResultReader.BetaPeptideTheoreticalMassLabel]].Trim();
            BetaPeptideScore = (parsedHeader[TsvResultReader.BetaPeptideScoreLabel] < 0) ? null : (double?)double.Parse(spl[parsedHeader[TsvResultReader.BetaPeptideScoreLabel]].Trim());
            BetaPeptideRank = (parsedHeader[TsvResultReader.BetaPeptideRankLabel] < 0) ? null : (int?)int.Parse(spl[parsedHeader[TsvResultReader.BetaPeptideRankLabel]].Trim());
            BetaPeptideMatchedIons = (parsedHeader[TsvResultReader.BetaPeptideMatchedIonsLabel] < 0) ? null : ReadFragmentIonsFromString(spl[parsedHeader[TsvResultReader.BetaPeptideMatchedIonsLabel]].Trim(), BetaPeptideBaseSequence);
            XLTotalScore = (parsedHeader[TsvResultReader.XLTotalScoreLabel] < 0) ? null : (double?)double.Parse(spl[parsedHeader[TsvResultReader.XLTotalScoreLabel]].Trim());
            ParentIons= (parsedHeader[TsvResultReader.ParentIonsLabel] < 0) ? null : spl[parsedHeader[TsvResultReader.ParentIonsLabel]].Trim();
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
