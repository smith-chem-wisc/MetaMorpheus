using Proteomics.Fragmentation;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text.RegularExpressions;
using Chemistry;

namespace MetaMorpheusGUI
{
    public class MetaDrawPsm
    {
        private static readonly Regex IonParser = new Regex(@"([a-zA-Z]+)(\d+)");
        private static readonly char[] MzSplit = { '[', ',', ']', ';' };

        public string FullSequence { get; }
        public int Ms2ScanNumber { get; }
        public string Filename { get; }
        public double TotalIonCurrent { get; }
        public int PrecursorScanNum { get; }
        public int PrecursorCharge { get; }
        public double PrecursorMz { get; }
        public double PrecursorMass { get; }
        public double Score { get; }
        public double DeltaScore { get; }
        public string Notch { get; }
        public string BaseSeq { get; }
        public string EssentialSeq { get; }
        public string MissedCleavage { get; }
        public string PeptideMonoMass { get; }
        public string MassDiffDa { get; }
        public string MassDiffPpm { get; }
        public string ProteinAccession { get; }
        public string ProteinName { get; }
        public string GeneName { get; }
        public string SequenceVariations { get; }
        public string OrganismName { get; }
        public string PeptideDesicription { get; }
        public string StartAndEndResiduesInProtein { get; }
        public string PreviousAminoAcid { get; }
        public string NextAminoAcid { get; }
        public string DecoyContamTarget { get; }
        public List<MatchedFragmentIon> MatchedIons { get; }
        public double QValue { get; }
        public double QValueNotch { get; }

        public MetaDrawPsm(string line, char[] split, Dictionary<string, int> parsedHeader)
        {
            var spl = line.Split(split);

            FullSequence = spl[parsedHeader[TsvResultReader.FullSequenceLabel]];
            Ms2ScanNumber = int.Parse(spl[parsedHeader[TsvResultReader.Ms2ScanNumberLabel]]);
            Filename = spl[parsedHeader[TsvResultReader.FilenameLabel]].Trim();
            TotalIonCurrent = double.Parse(spl[parsedHeader[TsvResultReader.TotalIonCurrentLabel]].Trim());
            PrecursorScanNum = int.Parse(spl[parsedHeader[TsvResultReader.PrecursorScanNumLabel]].Trim());
            PrecursorCharge = (int)double.Parse(spl[parsedHeader[TsvResultReader.PrecursorChargeLabel]].Trim());
            PrecursorMz = double.Parse(spl[parsedHeader[TsvResultReader.PrecursorMzLabel]].Trim());
            PrecursorMass = double.Parse(spl[parsedHeader[TsvResultReader.PrecursorMassLabel]].Trim());
            Score = double.Parse(spl[parsedHeader[TsvResultReader.ScoreLabel]].Trim());
            DeltaScore = double.Parse(spl[parsedHeader[TsvResultReader.DeltaScoreLabel]].Trim());
            Notch = spl[parsedHeader[TsvResultReader.NotchLabel]].Trim();
            BaseSeq = spl[parsedHeader[TsvResultReader.BaseSeqLabel]].Trim();
            EssentialSeq = spl[parsedHeader[TsvResultReader.EssentialSeqLabel]].Trim();
            MissedCleavage = spl[parsedHeader[TsvResultReader.MissedCleavageLabel]].Trim();
            PeptideMonoMass = spl[parsedHeader[TsvResultReader.PeptideMonoMassLabel]].Trim();
            MassDiffDa = spl[parsedHeader[TsvResultReader.MassDiffDaLabel]].Trim();
            MassDiffPpm = spl[parsedHeader[TsvResultReader.MassDiffPpmLabel]].Trim();
            ProteinAccession = spl[parsedHeader[TsvResultReader.ProteinAccessionLabel]].Trim();
            ProteinName = spl[parsedHeader[TsvResultReader.ProteinNameLabel]].Trim();
            GeneName = spl[parsedHeader[TsvResultReader.GeneNameLabel]].Trim();
            SequenceVariations = spl[parsedHeader[TsvResultReader.SequenceVariationsLabel]].Trim();
            OrganismName = spl[parsedHeader[TsvResultReader.OrganismNameLabel]].Trim();
            PeptideDesicription = spl[parsedHeader[TsvResultReader.PeptideDesicriptionLabel]].Trim();
            StartAndEndResiduesInProtein = spl[parsedHeader[TsvResultReader.StartAndEndResiduesInProteinLabel]].Trim();
            PreviousAminoAcid = spl[parsedHeader[TsvResultReader.PreviousAminoAcidLabel]].Trim();
            NextAminoAcid = spl[parsedHeader[TsvResultReader.NextAminoAcidLabel]].Trim();
            DecoyContamTarget = spl[parsedHeader[TsvResultReader.DecoyContamTargetLabel]].Trim();
            QValue = double.Parse(spl[parsedHeader[TsvResultReader.QValueLabel]].Trim());
            QValueNotch = double.Parse(spl[parsedHeader[TsvResultReader.QValueNotchLabel]].Trim());

            MatchedIons = ReadFragmentIonsFromString(spl[parsedHeader[TsvResultReader.MatchedIonsLabel]].Trim(), BaseSeq);
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
