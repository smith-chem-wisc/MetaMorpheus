using Chemistry;
using Proteomics.Fragmentation;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text.RegularExpressions;
using EngineLayer.GlycoSearch;
using System.IO;

namespace EngineLayer
{
    public class PsmFromTsv
    {
        private static readonly Regex PositionParser = new Regex(@"(\d+)\s+to\s+(\d+)");
        private static readonly Regex VariantParser = new Regex(@"[a-zA-Z]+(\d+)([a-zA-Z]+)");
        private static readonly Regex IonParser = new Regex(@"([a-zA-Z]+)(\d+)");

        public string FullSequence { get; }
        public int Ms2ScanNumber { get; }
        public string FileNameWithoutExtension { get; }
        public int PrecursorScanNum { get; }
        public int PrecursorCharge { get; }
        public double PrecursorMz { get; }
        public double PrecursorMass { get; }
        public double Score { get; }
        public string ProteinAccession { get; }
        public List<MatchedFragmentIon> MatchedIons { get; }
        public Dictionary<int, List<MatchedFragmentIon>> ChildScanMatchedIons { get; } // this is only used in crosslink for now, but in the future will be used for other experiment types
        public double QValue { get; }

        public double PEP { get; }

        public double PEP_QValue { get; }

        public double? TotalIonCurrent { get; }
        public double? DeltaScore { get; }
        public string Notch { get; }
        public string BaseSeq { get; }
        public string EssentialSeq { get; }
        public string AmbiguityLevel { get; }
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
        public string PeptideDescription { get; }
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

        //For Glyco
        public string GlycanStructure { get; set; }
        public double? GlycanMass { get; set; }
        public string GlycanComposition { get; set; }
        public LocalizationLevel? GlycanLocalizationLevel { get; set; }
        public string LocalizedGlycan { get; set; }

        public PsmFromTsv(string line, char[] split, Dictionary<string, int> parsedHeader)
        {
            var spl = line.Split(split);

            //Required properties
            FileNameWithoutExtension = spl[parsedHeader[PsmTsvHeader.FileName]].Trim();

            // remove file format, e.g., .raw, .mzML, .mgf
            // this is more robust but slower than Path.GetFileNameWithoutExtension
            if (FileNameWithoutExtension.Contains('.'))
            {
                foreach (var knownSpectraFileExtension in GlobalVariables.AcceptedSpectraFormats)
                {
                    FileNameWithoutExtension = Path.GetFileName(FileNameWithoutExtension.Replace(knownSpectraFileExtension, string.Empty, StringComparison.InvariantCultureIgnoreCase));
                }
            }

            Ms2ScanNumber = int.Parse(spl[parsedHeader[PsmTsvHeader.Ms2ScanNumber]]);

            // this will probably not be known in an .mgf data file
            if (int.TryParse(spl[parsedHeader[PsmTsvHeader.PrecursorScanNum]].Trim(), out int result))
            {
                PrecursorScanNum = result;
            }
            else
            {
                PrecursorScanNum = 0;
            }

            PrecursorCharge = (int)double.Parse(spl[parsedHeader[PsmTsvHeader.PrecursorCharge]].Trim(), CultureInfo.InvariantCulture);
            PrecursorMz = double.Parse(spl[parsedHeader[PsmTsvHeader.PrecursorMz]].Trim(), CultureInfo.InvariantCulture);
            PrecursorMass = double.Parse(spl[parsedHeader[PsmTsvHeader.PrecursorMass]].Trim(), CultureInfo.InvariantCulture);
            BaseSeq = RemoveParentheses(spl[parsedHeader[PsmTsvHeader.BaseSequence]].Trim());
            FullSequence = spl[parsedHeader[PsmTsvHeader.FullSequence]];
            PeptideMonoMass = spl[parsedHeader[PsmTsvHeader.PeptideMonoMass]].Trim();
            Score = double.Parse(spl[parsedHeader[PsmTsvHeader.Score]].Trim(), CultureInfo.InvariantCulture);
            DecoyContamTarget = spl[parsedHeader[PsmTsvHeader.DecoyContaminantTarget]].Trim();
            QValue = double.Parse(spl[parsedHeader[PsmTsvHeader.QValue]].Trim(), CultureInfo.InvariantCulture);
            MatchedIons = (spl[parsedHeader[PsmTsvHeader.MatchedIonMzRatios]].StartsWith("{")) ? ReadChildScanMatchedIons(spl[parsedHeader[PsmTsvHeader.MatchedIonMzRatios]].Trim(), spl[parsedHeader[PsmTsvHeader.MatchedIonIntensities]].Trim(), BaseSeq).First().Value : ReadFragmentIonsFromString(spl[parsedHeader[PsmTsvHeader.MatchedIonMzRatios]].Trim(), spl[parsedHeader[PsmTsvHeader.MatchedIonIntensities]].Trim(), BaseSeq);
            AmbiguityLevel = (parsedHeader[PsmTsvHeader.AmbiguityLevel] < 0) ? null : spl[parsedHeader[PsmTsvHeader.AmbiguityLevel]].Trim();

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
            PeptideDescription = (parsedHeader[PsmTsvHeader.PeptideDesicription] < 0) ? null : spl[parsedHeader[PsmTsvHeader.PeptideDesicription]].Trim();
            StartAndEndResiduesInProtein = (parsedHeader[PsmTsvHeader.StartAndEndResiduesInProtein] < 0) ? null : spl[parsedHeader[PsmTsvHeader.StartAndEndResiduesInProtein]].Trim();
            PreviousAminoAcid = (parsedHeader[PsmTsvHeader.PreviousAminoAcid] < 0) ? null : spl[parsedHeader[PsmTsvHeader.PreviousAminoAcid]].Trim();
            NextAminoAcid = (parsedHeader[PsmTsvHeader.NextAminoAcid] < 0) ? null : spl[parsedHeader[PsmTsvHeader.NextAminoAcid]].Trim();
            QValueNotch = (parsedHeader[PsmTsvHeader.QValueNotch] < 0) ? null : (double?)double.Parse(spl[parsedHeader[PsmTsvHeader.QValueNotch]].Trim(), CultureInfo.InvariantCulture);
            RetentionTime = (parsedHeader[PsmTsvHeader.Ms2ScanRetentionTime] < 0) ? null : (double?)double.Parse(spl[parsedHeader[PsmTsvHeader.Ms2ScanRetentionTime]].Trim(), CultureInfo.InvariantCulture);
            PEP = double.Parse(spl[parsedHeader[PsmTsvHeader.PEP]].Trim(), CultureInfo.InvariantCulture);
            PEP_QValue = double.Parse(spl[parsedHeader[PsmTsvHeader.PEP_QValue]].Trim(), CultureInfo.InvariantCulture);
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
                ((spl[parsedHeader[PsmTsvHeader.BetaPeptideMatchedIonsLabel]].StartsWith("{")) ? ReadChildScanMatchedIons(spl[parsedHeader[PsmTsvHeader.BetaPeptideMatchedIonsLabel]].Trim(), spl[parsedHeader[PsmTsvHeader.BetaPeptideMatchedIonIntensitiesLabel]].Trim(), BetaPeptideBaseSequence).First().Value : ReadFragmentIonsFromString(spl[parsedHeader[PsmTsvHeader.BetaPeptideMatchedIonsLabel]].Trim(), spl[parsedHeader[PsmTsvHeader.BetaPeptideMatchedIonIntensitiesLabel]].Trim(), BetaPeptideBaseSequence));
            XLTotalScore = (parsedHeader[PsmTsvHeader.XLTotalScoreLabel] < 0) ? null : (double?)double.Parse(spl[parsedHeader[PsmTsvHeader.XLTotalScoreLabel]].Trim(), CultureInfo.InvariantCulture);
            ParentIons = (parsedHeader[PsmTsvHeader.ParentIonsLabel] < 0) ? null : spl[parsedHeader[PsmTsvHeader.ParentIonsLabel]].Trim();

            // child scan matched ions (only for crosslinks for now, but in the future this will change) 
            ChildScanMatchedIons = (!spl[parsedHeader[PsmTsvHeader.MatchedIonMzRatios]].StartsWith("{")) ? null : ReadChildScanMatchedIons(spl[parsedHeader[PsmTsvHeader.MatchedIonMzRatios]].Trim(), spl[parsedHeader[PsmTsvHeader.MatchedIonIntensities]].Trim(), BaseSeq);
            if (ChildScanMatchedIons != null && ChildScanMatchedIons.ContainsKey(Ms2ScanNumber))
            {
                ChildScanMatchedIons.Remove(Ms2ScanNumber);
            }

            // beta peptide child scan matched ions (for crosslinks)
            BetaPeptideChildScanMatchedIons = (parsedHeader[PsmTsvHeader.BetaPeptideMatchedIonsLabel] < 0) ? null :
                ((!spl[parsedHeader[PsmTsvHeader.BetaPeptideMatchedIonsLabel]].StartsWith("{")) ? null : ReadChildScanMatchedIons(spl[parsedHeader[PsmTsvHeader.BetaPeptideMatchedIonsLabel]].Trim(), spl[parsedHeader[PsmTsvHeader.BetaPeptideMatchedIonIntensitiesLabel]].Trim(), BetaPeptideBaseSequence));
            if (BetaPeptideChildScanMatchedIons != null && BetaPeptideChildScanMatchedIons.ContainsKey(Ms2ScanNumber))
            {
                BetaPeptideChildScanMatchedIons.Remove(Ms2ScanNumber);
            }

            //For Glyco            
            GlycanMass = (parsedHeader[PsmTsvHeader_Glyco.GlycanMass] < 0) ? null : (double?)double.Parse(spl[parsedHeader[PsmTsvHeader_Glyco.GlycanMass]], CultureInfo.InvariantCulture);
            GlycanComposition = (parsedHeader[PsmTsvHeader_Glyco.GlycanComposition] < 0) ? null : spl[parsedHeader[PsmTsvHeader_Glyco.GlycanComposition]];
            GlycanStructure = (parsedHeader[PsmTsvHeader_Glyco.GlycanStructure] < 0) ? null : spl[parsedHeader[PsmTsvHeader_Glyco.GlycanStructure]];
            var localizationLevel = (parsedHeader[PsmTsvHeader_Glyco.GlycanLocalizationLevel] < 0) ? null : spl[parsedHeader[PsmTsvHeader_Glyco.GlycanLocalizationLevel]];
            if (localizationLevel != null)
            {
                GlycanLocalizationLevel = (LocalizationLevel)Enum.Parse(typeof(LocalizationLevel), localizationLevel);
            }
            LocalizedGlycan = (parsedHeader[PsmTsvHeader_Glyco.LocalizedGlycan] < 0) ? null : spl[parsedHeader[PsmTsvHeader_Glyco.LocalizedGlycan]];
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

        private static List<MatchedFragmentIon> ReadFragmentIonsFromString(string matchedMzString, string matchedIntensityString, string peptideBaseSequence)
        {
            List<MatchedFragmentIon> matchedIons = new List<MatchedFragmentIon>();

            if (matchedMzString.Length > 2) //check if there's an ion
            {
                List<string> peakMzs = CleanMatchedIonString(matchedMzString);
                List<string> peakIntensities = CleanMatchedIonString(matchedIntensityString);

                for (int index = 0; index < peakMzs.Count; index++)
                {
                    string peak = peakMzs[index];
                    string[] split = peak.Split(new char[] { '+', ':' }); //TODO: needs update for negative charges that doesn't break internal fragment ions or neutral losses
                    double intensity = peakMzs.Count == peakIntensities.Count ? //TODO: needs update for negative charges that doesn't break internal fragment ions or neutral losses
                        double.Parse(peakIntensities[index].Split(new char[] { '+', ':' })[2], CultureInfo.InvariantCulture) :
                        1.0;

                    int fragmentNumber = 0;
                    int secondaryFragmentNumber = 0;
                    ProductType productType;
                    ProductType? secondaryProductType = null;
                    FragmentationTerminus terminus = FragmentationTerminus.None; //default for internal fragments
                    int aminoAcidPosition;
                    double neutralLoss = 0;

                    //get theoretical fragment
                    string ionTypeAndNumber = split[0];

                    //if an internal fragment
                    if (ionTypeAndNumber.Contains("["))
                    {
                        string[] internalSplit = split[0].Split('[');
                        string[] productSplit = internalSplit[0].Split("I");
                        string[] positionSplit = internalSplit[1].Replace("]", "").Split('-');
                        productType = (ProductType)Enum.Parse(typeof(ProductType), productSplit[0]);
                        secondaryProductType = (ProductType)Enum.Parse(typeof(ProductType), productSplit[1]);
                        fragmentNumber = int.Parse(positionSplit[0]);
                        secondaryFragmentNumber = int.Parse(positionSplit[1]);
                        aminoAcidPosition = secondaryFragmentNumber - fragmentNumber;
                    }
                    else //terminal fragment
                    {
                        Match result = IonParser.Match(ionTypeAndNumber);
                        productType = (ProductType)Enum.Parse(typeof(ProductType), result.Groups[1].Value);
                        fragmentNumber = int.Parse(result.Groups[2].Value);
                        // check for neutral loss  
                        if (ionTypeAndNumber.Contains("("))
                        {
                            string temp = ionTypeAndNumber.Replace("(", "");
                            temp = temp.Replace(")", "");
                            var split2 = temp.Split('-');
                            neutralLoss = double.Parse(split2[1], CultureInfo.InvariantCulture);
                        }

                        //get terminus
                        if (TerminusSpecificProductTypes.ProductTypeToFragmentationTerminus.ContainsKey(productType))
                        {
                            terminus = TerminusSpecificProductTypes.ProductTypeToFragmentationTerminus[productType];
                        }

                        //get amino acid position
                        aminoAcidPosition = terminus == FragmentationTerminus.C ?
                            peptideBaseSequence.Length - fragmentNumber :
                            fragmentNumber;
                    }

                    //get charge and mz
                    int z = int.Parse(split[1]);
                    double mz = double.Parse(split[2], CultureInfo.InvariantCulture);

                    Product p = new Product(productType,
                      terminus,
                      mz.ToMass(z),
                      fragmentNumber,
                      aminoAcidPosition,
                      neutralLoss,
                      secondaryProductType,
                      secondaryFragmentNumber);

                    matchedIons.Add(new MatchedFragmentIon(ref p, mz, intensity, z));
                }
            }
            return matchedIons;
        }

        private static List<string> CleanMatchedIonString(string input)
        {
            List<string> ionProperty = input.Substring(1, input.Length - 2) //remove the brackets on the ends
                    .Replace("];[", ", ") //replace delimiter between ion series with the delimiter used between ions
                    .Split(", ") //split by delimiter between ions
                    .ToList();
            ionProperty.RemoveAll(p => p.Contains("\"") || p.Equals(""));
            return ionProperty;
        }

        private static Dictionary<int, List<MatchedFragmentIon>> ReadChildScanMatchedIons(string childScanMatchedMzString, string childScanMatchedIntensitiesString, string peptideBaseSequence)
        {
            var childScanMatchedIons = new Dictionary<int, List<MatchedFragmentIon>>();

            foreach (var childScan in childScanMatchedMzString.Split(new char[] { '}' }).Where(p => !string.IsNullOrWhiteSpace(p)))
            {
                var split1 = childScan.Split(new char[] { '@' });
                int scanNumber = int.Parse(split1[0].Trim(new char[] { '{' }));
                string matchedIonsString = split1[1];
                var childMatchedIons = ReadFragmentIonsFromString(matchedIonsString, childScanMatchedIntensitiesString, peptideBaseSequence);
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

        public static List<Tuple<int, string, double>> ReadLocalizedGlycan(string localizedGlycan)
        {
            List<Tuple<int, string, double>> tuples = new List<Tuple<int, string, double>>();
            if (localizedGlycan == null)
            {
                return tuples;
            }
            var lgs = localizedGlycan.Split(new string[] { "[", "]" }, StringSplitOptions.RemoveEmptyEntries);
            foreach (var lg in lgs)
            {
                var g = lg.Split(',', StringSplitOptions.RemoveEmptyEntries);

                Tuple<int, string, double> tuple = new Tuple<int, string, double>(int.Parse(g[0], CultureInfo.InvariantCulture), g[1], double.Parse(g[2], CultureInfo.InvariantCulture));
                tuples.Add(tuple);
            }

            return tuples;
        }
    }
}