using Chemistry;
using Proteomics.Fragmentation;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text.RegularExpressions;

namespace EngineLayer
{
    public class PsmFromTsv
    {
        private static readonly Regex PositionParser = new Regex(@"(\d+)\s+to\s+(\d+)");
        private static readonly Regex VariantParser = new Regex(@"[a-zA-Z]+(\d+)([a-zA-Z]+)");
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
        public Dictionary<int, List<MatchedFragmentIon>> ChildScanMatchedIons { get; } // this is only used in crosslink for now, but in the future will be used for other experiment types
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
        public string IntersectingSequenceVariations { get; }
        public string IdentifiedSequenceVariations { get; }
        public string SpliceSites { get; }
        public string PeptideDesicription { get; }
        public string StartAndEndResiduesInProtein { get; }
        public string PreviousAminoAcid { get; }
        public string NextAminoAcid { get; }
        public string DecoyContamTarget { get; }
        public double? QValueNotch { get; }

        public List<MatchedFragmentIon> VariantCrossingIons { get; }

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
        public Dictionary<int, List<MatchedFragmentIon>> BetaPeptideChildScanMatchedIons { get; }
        public double? XLTotalScore { get; }
        public string ParentIons { get; }
        public double? RetentionTime { get; }

        public PsmFromTsv(string line, char[] split, Dictionary<string, int> parsedHeader)
        {
            var spl = line.Split(split);

            //Required properties
            Filename = spl[parsedHeader[PsmTsvHeader.FileName]].Trim();
            Ms2ScanNumber = int.Parse(spl[parsedHeader[PsmTsvHeader.Ms2ScanNumber]]);
            PrecursorScanNum = int.Parse(spl[parsedHeader[PsmTsvHeader.PrecursorScanNum]].Trim());
            PrecursorCharge = (int)double.Parse(spl[parsedHeader[PsmTsvHeader.PrecursorCharge]].Trim(), CultureInfo.InvariantCulture);
            PrecursorMz = double.Parse(spl[parsedHeader[PsmTsvHeader.PrecursorMz]].Trim(), CultureInfo.InvariantCulture);
            PrecursorMass = double.Parse(spl[parsedHeader[PsmTsvHeader.PrecursorMass]].Trim(), CultureInfo.InvariantCulture);
            BaseSeq = RemoveParentheses(spl[parsedHeader[PsmTsvHeader.BaseSequence]].Trim());
            FullSequence = spl[parsedHeader[PsmTsvHeader.FullSequence]];
            PeptideMonoMass = spl[parsedHeader[PsmTsvHeader.PeptideMonoMass]].Trim();
            Score = double.Parse(spl[parsedHeader[PsmTsvHeader.Score]].Trim(), CultureInfo.InvariantCulture);
            DecoyContamTarget = spl[parsedHeader[PsmTsvHeader.DecoyContaminantTarget]].Trim();
            QValue = double.Parse(spl[parsedHeader[PsmTsvHeader.QValue]].Trim(), CultureInfo.InvariantCulture);
            MatchedIons = (spl[parsedHeader[PsmTsvHeader.MatchedIonMzRatios]].StartsWith("{")) ? ReadChildScanMatchedIons(spl[parsedHeader[PsmTsvHeader.MatchedIonMzRatios]].Trim(), BaseSeq).First().Value : ReadFragmentIonsFromString(spl[parsedHeader[PsmTsvHeader.MatchedIonMzRatios]].Trim(), BaseSeq);

            //For general psms
            TotalIonCurrent = (parsedHeader[PsmTsvHeader.TotalIonCurrent] < 0) ? null : (double?)double.Parse(spl[parsedHeader[PsmTsvHeader.TotalIonCurrent]].Trim(), CultureInfo.InvariantCulture);
            DeltaScore = (parsedHeader[PsmTsvHeader.DeltaScore] < 0) ? null : (double?)double.Parse(spl[parsedHeader[PsmTsvHeader.DeltaScore]].Trim(), CultureInfo.InvariantCulture);
            Notch = (parsedHeader[PsmTsvHeader.Notch] < 0) ? null : spl[parsedHeader[PsmTsvHeader.Notch]].Trim();
            EssentialSeq = (parsedHeader[PsmTsvHeader.EssentialSequence] < 0) ? null : spl[parsedHeader[PsmTsvHeader.EssentialSequence]].Trim();
            MissedCleavage = (parsedHeader[PsmTsvHeader.MissedCleavages] < 0) ? null : spl[parsedHeader[PsmTsvHeader.MissedCleavages]].Trim();
            MassDiffDa = (parsedHeader[PsmTsvHeader.MassDiffDa] < 0) ? null : spl[parsedHeader[PsmTsvHeader.MassDiffDa]].Trim();
            MassDiffPpm = (parsedHeader[PsmTsvHeader.MassDiffPpm] < 0) ? null : spl[parsedHeader[PsmTsvHeader.MassDiffPpm]].Trim();
            ProteinAccession = (parsedHeader[PsmTsvHeader.ProteinAccession] < 0) ? null : spl[parsedHeader[PsmTsvHeader.ProteinAccession]].Trim();
            ProteinName = (parsedHeader[PsmTsvHeader.ProteinName] < 0) ? null : spl[parsedHeader[PsmTsvHeader.ProteinName]].Trim();
            GeneName = (parsedHeader[PsmTsvHeader.GeneName] < 0) ? null : spl[parsedHeader[PsmTsvHeader.GeneName]].Trim();
            OrganismName = (parsedHeader[PsmTsvHeader.OrganismName] < 0) ? null : spl[parsedHeader[PsmTsvHeader.OrganismName]].Trim();
            IntersectingSequenceVariations = (parsedHeader[PsmTsvHeader.IntersectingSequenceVariations] < 0) ? null : spl[parsedHeader[PsmTsvHeader.IntersectingSequenceVariations]].Trim();
            IdentifiedSequenceVariations = (parsedHeader[PsmTsvHeader.IdentifiedSequenceVariations] < 0) ? null : spl[parsedHeader[PsmTsvHeader.IdentifiedSequenceVariations]].Trim();
            SpliceSites = (parsedHeader[PsmTsvHeader.SpliceSites] < 0) ? null : spl[parsedHeader[PsmTsvHeader.SpliceSites]].Trim();
            PeptideDesicription = (parsedHeader[PsmTsvHeader.PeptideDesicription] < 0) ? null : spl[parsedHeader[PsmTsvHeader.PeptideDesicription]].Trim();
            StartAndEndResiduesInProtein = (parsedHeader[PsmTsvHeader.StartAndEndResiduesInProtein] < 0) ? null : spl[parsedHeader[PsmTsvHeader.StartAndEndResiduesInProtein]].Trim();
            PreviousAminoAcid = (parsedHeader[PsmTsvHeader.PreviousAminoAcid] < 0) ? null : spl[parsedHeader[PsmTsvHeader.PreviousAminoAcid]].Trim();
            NextAminoAcid = (parsedHeader[PsmTsvHeader.NextAminoAcid] < 0) ? null : spl[parsedHeader[PsmTsvHeader.NextAminoAcid]].Trim();
            QValueNotch = (parsedHeader[PsmTsvHeader.QValueNotch] < 0) ? null : (double?)double.Parse(spl[parsedHeader[PsmTsvHeader.QValueNotch]].Trim(), CultureInfo.InvariantCulture);
            RetentionTime = (parsedHeader[PsmTsvHeader.Ms2ScanRetentionTime] < 0) ? null : (double?)double.Parse(spl[parsedHeader[PsmTsvHeader.Ms2ScanRetentionTime]].Trim(), CultureInfo.InvariantCulture);

            VariantCrossingIons = findVariantCrossingIons();

            //For crosslinks
            CrossType = (parsedHeader[PsmTsvHeader.CrossTypeLabel] < 0) ? null : spl[parsedHeader[PsmTsvHeader.CrossTypeLabel]].Trim();
            LinkResidues = (parsedHeader[PsmTsvHeader.LinkResiduesLabel] < 0) ? null : spl[parsedHeader[PsmTsvHeader.LinkResiduesLabel]].Trim();
            ProteinLinkSite = (parsedHeader[PsmTsvHeader.ProteinLinkSiteLabel] < 0) ? null : (spl[parsedHeader[PsmTsvHeader.ProteinLinkSiteLabel]] == "" ? null : (int?)int.Parse(spl[parsedHeader[PsmTsvHeader.ProteinLinkSiteLabel]].Trim()));
            Rank = (parsedHeader[PsmTsvHeader.RankLabel] < 0) ? null : (int?)int.Parse(spl[parsedHeader[PsmTsvHeader.RankLabel]].Trim());
            BetaPeptideProteinAccession = (parsedHeader[PsmTsvHeader.BetaPeptideProteinAccessionLabel] < 0) ? null : spl[parsedHeader[PsmTsvHeader.BetaPeptideProteinAccessionLabel]].Trim();
            BetaPeptideProteinLinkSite = (parsedHeader[PsmTsvHeader.BetaPeptideProteinLinkSiteLabel] < 0) ? null : (spl[parsedHeader[PsmTsvHeader.BetaPeptideProteinLinkSiteLabel]] == "" ? null : (int?)int.Parse(spl[parsedHeader[PsmTsvHeader.BetaPeptideProteinLinkSiteLabel]].Trim()));
            BetaPeptideBaseSequence = (parsedHeader[PsmTsvHeader.BetaPeptideBaseSequenceLabel] < 0) ? null : spl[parsedHeader[PsmTsvHeader.BetaPeptideBaseSequenceLabel]].Trim();
            BetaPeptideFullSequence = (parsedHeader[PsmTsvHeader.BetaPeptideFullSequenceLabel] < 0) ? null : spl[parsedHeader[PsmTsvHeader.BetaPeptideFullSequenceLabel]].Trim();
            BetaPeptideTheoreticalMass = (parsedHeader[PsmTsvHeader.BetaPeptideTheoreticalMassLabel] < 0) ? null : spl[parsedHeader[PsmTsvHeader.BetaPeptideTheoreticalMassLabel]].Trim();
            BetaPeptideScore = (parsedHeader[PsmTsvHeader.BetaPeptideScoreLabel] < 0) ? null : (double?)double.Parse(spl[parsedHeader[PsmTsvHeader.BetaPeptideScoreLabel]].Trim(), CultureInfo.InvariantCulture);
            BetaPeptideRank = (parsedHeader[PsmTsvHeader.BetaPeptideRankLabel] < 0) ? null : (int?)int.Parse(spl[parsedHeader[PsmTsvHeader.BetaPeptideRankLabel]].Trim());
            BetaPeptideMatchedIons = (parsedHeader[PsmTsvHeader.BetaPeptideMatchedIonsLabel] < 0) ? null :
                ((spl[parsedHeader[PsmTsvHeader.BetaPeptideMatchedIonsLabel]].StartsWith("{")) ? ReadChildScanMatchedIons(spl[parsedHeader[PsmTsvHeader.BetaPeptideMatchedIonsLabel]].Trim(), BetaPeptideBaseSequence).First().Value : ReadFragmentIonsFromString(spl[parsedHeader[PsmTsvHeader.BetaPeptideMatchedIonsLabel]].Trim(), BetaPeptideBaseSequence));
            XLTotalScore = (parsedHeader[PsmTsvHeader.XLTotalScoreLabel] < 0) ? null : (double?)double.Parse(spl[parsedHeader[PsmTsvHeader.XLTotalScoreLabel]].Trim(), CultureInfo.InvariantCulture);
            ParentIons = (parsedHeader[PsmTsvHeader.ParentIonsLabel] < 0) ? null : spl[parsedHeader[PsmTsvHeader.ParentIonsLabel]].Trim();

            // child scan matched ions (only for crosslinks for now, but in the future this will change) 
            ChildScanMatchedIons = (!spl[parsedHeader[PsmTsvHeader.MatchedIonMzRatios]].StartsWith("{")) ? null : ReadChildScanMatchedIons(spl[parsedHeader[PsmTsvHeader.MatchedIonMzRatios]].Trim(), BaseSeq);
            if (ChildScanMatchedIons != null && ChildScanMatchedIons.ContainsKey(Ms2ScanNumber))
            {
                ChildScanMatchedIons.Remove(Ms2ScanNumber);
            }

            // beta peptide child scan matched ions (for crosslinks)
            BetaPeptideChildScanMatchedIons = (parsedHeader[PsmTsvHeader.BetaPeptideMatchedIonsLabel] < 0) ? null :
                ((!spl[parsedHeader[PsmTsvHeader.BetaPeptideMatchedIonsLabel]].StartsWith("{")) ? null : ReadChildScanMatchedIons(spl[parsedHeader[PsmTsvHeader.BetaPeptideMatchedIonsLabel]].Trim(), BetaPeptideBaseSequence));
            if (BetaPeptideChildScanMatchedIons != null && BetaPeptideChildScanMatchedIons.ContainsKey(Ms2ScanNumber))
            {
                BetaPeptideChildScanMatchedIons.Remove(Ms2ScanNumber);
            }
        }

        //Used to remove Silac labels for proper annotation
        public static string RemoveParentheses(string baseSequence)
        {
            if (baseSequence.Contains("("))
            {
                string updatedBaseSequence = "";
                bool withinParentheses = false;
                foreach (char c in baseSequence)
                {
                    if (c == ')') //leaving the parentheses
                    {
                        withinParentheses = false;
                    }
                    else if (c == '(') //entering the parentheses
                    {
                        withinParentheses = true;
                    }
                    else if (!withinParentheses) //if outside the parentheses, preserve this amino acid
                    {
                        updatedBaseSequence += c;
                    }
                    //else do nothing
                }
                return updatedBaseSequence;
            }
            return baseSequence;
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
                double mz = double.Parse(split[2], CultureInfo.InvariantCulture);
                double neutralLoss = 0;

                // check for neutral loss
                if (ionTypeAndNumber.Contains("-"))
                {
                    string temp = ionTypeAndNumber.Replace("(", "");
                    temp = temp.Replace(")", "");
                    var split2 = temp.Split('-');
                    neutralLoss = double.Parse(split2[1], CultureInfo.InvariantCulture);
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
                
                Product p = new Product(productType, 
                    terminus, 
                    mz.ToMass(z) - DissociationTypeCollection.GetMassShiftFromProductType(productType), 
                    fragmentNumber, 
                    aminoAcidPosition, 
                    neutralLoss);

                matchedIons.Add(new MatchedFragmentIon(ref p, mz, 1.0, z));
            }

            return matchedIons;
        }

        private static Dictionary<int, List<MatchedFragmentIon>> ReadChildScanMatchedIons(string childScanMatchedMzString, string peptideBaseSequence)
        {
            var childScanMatchedIons = new Dictionary<int, List<MatchedFragmentIon>>();

            foreach (var childScan in childScanMatchedMzString.Split(new char[] { '}' }).Where(p => !string.IsNullOrWhiteSpace(p)))
            {
                var split1 = childScan.Split(new char[] { '@' });
                int scanNumber = int.Parse(split1[0].Trim(new char[] { '{' }));
                string matchedIonsString = split1[1];
                var childMatchedIons = ReadFragmentIonsFromString(matchedIonsString, peptideBaseSequence);
                childScanMatchedIons.Add(scanNumber, childMatchedIons);
            }

            return childScanMatchedIons;
        }

        // finds the ions that contain variant residues using the position in IdentifiedSequenceVariations. When the variation spans 
        // multiple residues, if any part is contained in an ion, the ion is marked as variant crossing.
        private List<MatchedFragmentIon> findVariantCrossingIons()
        {
            List<MatchedFragmentIon> variantCrossingIons = new List<MatchedFragmentIon>();


            if (StartAndEndResiduesInProtein != null && IdentifiedSequenceVariations != null)
            {
                Match positionMatch = PositionParser.Match(StartAndEndResiduesInProtein);
                Match variantMatch = VariantParser.Match(IdentifiedSequenceVariations);
                if (positionMatch.Success && variantMatch.Success)
                {
                    List<ProductType> abcProductTypes = new List<ProductType>() { ProductType.a, ProductType.aDegree, ProductType.aStar,
                                                                    ProductType.b, ProductType.bDegree, ProductType.bStar, ProductType.c };
                    List<ProductType> xyzProductTypes = new List<ProductType>() { ProductType.x, ProductType.y, ProductType.yDegree,
                                                                    ProductType.yStar, ProductType.zDot, ProductType.zPlusOne};
                    int peptideStart = int.Parse(positionMatch.Groups[1].Value);
                    int peptideEnd = int.Parse(positionMatch.Groups[2].Value);
                    int variantResidueStart = int.Parse(variantMatch.Groups[1].Value);
                    int variantResidueEnd = variantResidueStart + variantMatch.Groups[2].Value.Length - 1;

                    foreach (MatchedFragmentIon ion in MatchedIons)
                    {
                        Match ionMatch = IonParser.Match(ion.Annotation);
                        if (ionMatch.Success &&
                            (variantResidueEnd >= peptideStart && variantResidueStart <= peptideEnd) &&     // variant is within peptide
                            ((abcProductTypes.Contains(ion.NeutralTheoreticalProduct.ProductType) &&        // type a, b, or c
                              peptideStart + int.Parse(ionMatch.Groups[2].Value) > variantResidueStart) ||  // crosses variant
                             (xyzProductTypes.Contains(ion.NeutralTheoreticalProduct.ProductType) &&        // type x, y, or z
                              peptideEnd - int.Parse(ionMatch.Groups[2].Value) < variantResidueEnd)))       // crosses variant
                        {
                            variantCrossingIons.Add(ion);
                        }
                    }
                }
            }
            return variantCrossingIons;
        }
    }
}