using Chemistry;
using Proteomics.Fragmentation;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text.RegularExpressions;
using EngineLayer.GlycoSearch;
using System.IO;
using System.Text;
using Easy.Common.Extensions;
using EngineLayer.PsmTsv;
using Proteomics.ProteolyticDigestion;
using FlashLFQ;
using Proteomics;

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
        public double? SpectralAngle { get; }
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
        
        // First amino acid in protein is amino acid number 1, which differs from internal code numbering with N-terminus as 1
        // This numbering is for the peptide location within the protein
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
        public PeptideWithSetModifications PeptideWithSetModifications { get; }

        public PsmFromTsv(string line, char[] split, Dictionary<string, int> parsedHeader)
        {
            var spl = line.Split(split).Select(p => p.Trim('\"')).ToArray();

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
            MatchedIons = (spl[parsedHeader[PsmTsvHeader.MatchedIonMzRatios]].StartsWith("{")) ?
                ReadChildScanMatchedIons(spl[parsedHeader[PsmTsvHeader.MatchedIonMzRatios]].Trim(), spl[parsedHeader[PsmTsvHeader.MatchedIonIntensities]].Trim(), BaseSeq).First().Value : 
                ReadFragmentIonsFromString(spl[parsedHeader[PsmTsvHeader.MatchedIonMzRatios]].Trim(), spl[parsedHeader[PsmTsvHeader.MatchedIonIntensities]].Trim(), BaseSeq, spl[parsedHeader[PsmTsvHeader.MatchedIonMassDiffDa]].Trim());
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
            SpectralAngle = (parsedHeader[PsmTsvHeader.SpectralAngle] < 0) ? null : (double?)double.Parse(spl[parsedHeader[PsmTsvHeader.SpectralAngle]].Trim(), CultureInfo.InvariantCulture);

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
                if (localizationLevel.Equals("NA"))
                    GlycanLocalizationLevel = null;
                else
                    GlycanLocalizationLevel = (LocalizationLevel)Enum.Parse(typeof(LocalizationLevel), localizationLevel);
            }
            LocalizedGlycan = (parsedHeader[PsmTsvHeader_Glyco.LocalizedGlycan] < 0) ? null : spl[parsedHeader[PsmTsvHeader_Glyco.LocalizedGlycan]];
        }

        /// <summary>
        /// Used for creating PSMs from the outputs of different search engines. Currently only support MaxQuant
        /// </summary>
        /// <param name="line"></param>
        /// <param name="split"></param>
        /// <param name="parsedHeader"></param>
        /// <param name="fileType"></param>
        public PsmFromTsv(string[] splitLine, Dictionary<string, int> parsedHeader, PsmFileType fileType, 
            Dictionary<string, Modification> allKnownMods = null, int numFixedMods = 0, bool ignoreArtifactIons = false)
        {
            if (fileType != PsmFileType.MaxQuant)
                throw new NotImplementedException("Only MaxQuant msms.txt is currently supported");

            //Required properties
            FileNameWithoutExtension = splitLine[parsedHeader[MaxQuantMsmsHeader.FileName]].Trim();

            // remove file format, e.g., .raw, .mzML, .mgf
            // this is more robust but slower than Path.GetFileNameWithoutExtension
            if (FileNameWithoutExtension.Contains('.'))
            {
                foreach (var knownSpectraFileExtension in GlobalVariables.AcceptedSpectraFormats)
                {
                    FileNameWithoutExtension = Path.GetFileName(FileNameWithoutExtension.Replace(knownSpectraFileExtension, string.Empty, StringComparison.InvariantCultureIgnoreCase));
                }
            }

            Ms2ScanNumber = int.Parse(splitLine[parsedHeader[MaxQuantMsmsHeader.Ms2ScanNumber]]);

            // this will probably not be known in an .mgf data file
            if (int.TryParse(splitLine[parsedHeader[MaxQuantMsmsHeader.PrecursorScanNum]].Trim(), out int result))
            {
                PrecursorScanNum = result;
            }
            else
            {
                PrecursorScanNum = 0;
            }

            PrecursorCharge = (int)double.Parse(splitLine[parsedHeader[MaxQuantMsmsHeader.PrecursorCharge]].Trim(), CultureInfo.InvariantCulture);
            PrecursorMz = double.Parse(splitLine[parsedHeader[MaxQuantMsmsHeader.PrecursorMz]].Trim(), CultureInfo.InvariantCulture);
            PrecursorMass = double.Parse(splitLine[parsedHeader[MaxQuantMsmsHeader.PrecursorMass]].Trim(), CultureInfo.InvariantCulture);
            BaseSeq = RemoveParentheses(splitLine[parsedHeader[MaxQuantMsmsHeader.BaseSequence]].Trim());
            FullSequence = splitLine[parsedHeader[MaxQuantMsmsHeader.FullSequence]];
            Score = double.Parse(splitLine[parsedHeader[MaxQuantMsmsHeader.Score]].Trim(), CultureInfo.InvariantCulture);

            // MaxQuant msms sometimes contains empty fragment series, which causes problems down the line
            if (splitLine[parsedHeader[MaxQuantMsmsHeader.MatchedIonSeries]].IsNullOrEmptyOrWhiteSpace())
            {
                Score = 0;
                MatchedIons = null;
            }
            else
            {
                MatchedIons = ReadFragmentIonsFromList(
                    BuildFragmentStringsFromMaxQuant(
                        splitLine[parsedHeader[MaxQuantMsmsHeader.MatchedIonSeries]],
                        splitLine[parsedHeader[MaxQuantMsmsHeader.MatchedIonMzRatios]]
                    ),
                    BuildFragmentStringsFromMaxQuant(
                        splitLine[parsedHeader[MaxQuantMsmsHeader.MatchedIonSeries]],
                        splitLine[parsedHeader[MaxQuantMsmsHeader.MatchedIonIntensities]]
                    ),
                    BaseSeq,
                    ignoreArtifactIons: ignoreArtifactIons
                );
            }
            

            DeltaScore = (parsedHeader[MaxQuantMsmsHeader.DeltaScore] < 0) ? null : (double?)double.Parse(splitLine[parsedHeader[MaxQuantMsmsHeader.DeltaScore]].Trim(), CultureInfo.InvariantCulture);
            MassDiffDa = (parsedHeader[MaxQuantMsmsHeader.MassDiffDa] < 0) ? null : splitLine[parsedHeader[MaxQuantMsmsHeader.MassDiffDa]].Trim();
            MassDiffPpm = (parsedHeader[MaxQuantMsmsHeader.MassDiffPpm] < 0) ? null : splitLine[parsedHeader[MaxQuantMsmsHeader.MassDiffPpm]].Trim();
            ProteinAccession = (parsedHeader[MaxQuantMsmsHeader.ProteinAccession] < 0) ? null : splitLine[parsedHeader[MaxQuantMsmsHeader.ProteinAccession]].Trim();
            ProteinName = (parsedHeader[MaxQuantMsmsHeader.ProteinName] < 0) ? null : splitLine[parsedHeader[MaxQuantMsmsHeader.ProteinName]].Trim();
            GeneName = (parsedHeader[MaxQuantMsmsHeader.GeneName] < 0) ? null : splitLine[parsedHeader[MaxQuantMsmsHeader.GeneName]].Trim();
            RetentionTime = (parsedHeader[MaxQuantMsmsHeader.Ms2ScanRetentionTime] < 0) ? null : (double?)double.Parse(splitLine[parsedHeader[MaxQuantMsmsHeader.Ms2ScanRetentionTime]].Trim(), CultureInfo.InvariantCulture);
            PEP = double.Parse(splitLine[parsedHeader[MaxQuantMsmsHeader.PEP]].Trim(), CultureInfo.InvariantCulture);
            DecoyContamTarget = splitLine[parsedHeader[MaxQuantMsmsHeader.Decoy]].IsNullOrEmptyOrWhiteSpace() ? "T" : "D";

            if (!MassDiffDa.Equals("NaN") && Double.TryParse(MassDiffDa, out double massDiff)) PeptideMonoMass = (PrecursorMass + massDiff).ToString();
            else PeptideMonoMass = PrecursorMass.ToString();

            if (allKnownMods != null)
            {
                Protein protein = new Protein(sequence: BaseSeq, accession: ProteinAccession,
                    geneNames: new List<Tuple<string, string>> { new Tuple<string, string>("Unknown", GeneName) },
                    isDecoy: DecoyContamTarget == "D");
                PeptideWithSetModifications = new PeptideWithSetModifications(FullSequence, allKnownMods, numFixedMods, p: protein);
            }
        }

        /// <summary>
        /// Constructor used to disambiguate PsmFromTsv to a single psm object
        /// </summary>
        /// <param name="psm">psm to disambiguate</param>
        /// <param name="fullSequence">sequence of ambiguous psm to use</param>
        public PsmFromTsv(PsmFromTsv psm, string fullSequence, int index = 0, string baseSequence = "")
        {
            // psm is not ambiguous
            if (!psm.FullSequence.Contains("|"))
            {
                FullSequence = fullSequence;
                EssentialSeq = psm.EssentialSeq;
                BaseSeq = baseSequence == "" ? psm.BaseSeq : baseSequence;
                StartAndEndResiduesInProtein = psm.StartAndEndResiduesInProtein;
                ProteinAccession = psm.ProteinAccession;
                ProteinName = psm.ProteinName;
                GeneName = psm.GeneName;
                PeptideMonoMass = psm.PeptideMonoMass;
                MassDiffDa = psm.MassDiffDa;
                MassDiffPpm = psm.MassDiffPpm;
            }
            // potentially ambiguous fields
            else
            {
                FullSequence = fullSequence;
                EssentialSeq = psm.EssentialSeq.Split("|")[index];
                BaseSeq = baseSequence == "" ? psm.BaseSeq.Split("|")[index] : baseSequence;
                StartAndEndResiduesInProtein = psm.StartAndEndResiduesInProtein.Split("|")[index];
                ProteinAccession = psm.ProteinAccession.Split("|")[index];
                ProteinName = psm.ProteinName.Split("|")[index];
                GeneName = psm.GeneName.Split("|")[index];

                if (psm.PeptideMonoMass.Split("|").Count() == 1)
                {
                    PeptideMonoMass = psm.PeptideMonoMass.Split("|")[0];
                    MassDiffDa = psm.MassDiffDa.Split("|")[0];
                    MassDiffPpm = psm.MassDiffPpm.Split("|")[0];
                }
                else
                {
                    PeptideMonoMass = psm.PeptideMonoMass.Split("|")[index];
                    MassDiffDa = psm.MassDiffDa.Split("|")[index];
                    MassDiffPpm = psm.MassDiffPpm.Split("|")[index];
                }
            }
                       
            // non ambiguous fields
            Ms2ScanNumber = psm.Ms2ScanNumber;
            FileNameWithoutExtension = psm.FileNameWithoutExtension;
            PrecursorScanNum = psm.PrecursorScanNum;
            PrecursorCharge = psm.PrecursorCharge;
            Score = psm.Score;
            MatchedIons = psm.MatchedIons.ToList();
            ChildScanMatchedIons = psm.ChildScanMatchedIons;
            QValue = psm.QValue;
            PEP = psm.PEP;
            PEP_QValue = psm.PEP_QValue;
            TotalIonCurrent = psm.TotalIonCurrent;
            DeltaScore = psm.DeltaScore;
            AmbiguityLevel = psm.AmbiguityLevel;
            MissedCleavage = psm.MissedCleavage;
            OrganismName = psm.OrganismName;
            IntersectingSequenceVariations = psm.IntersectingSequenceVariations;
            SpliceSites = psm.SpliceSites;
            PeptideDescription = psm.PeptideDescription;
            PreviousAminoAcid = psm.PreviousAminoAcid;
            NextAminoAcid = psm.NextAminoAcid;
            DecoyContamTarget = psm.DecoyContamTarget;
            QValueNotch = psm.QValueNotch;
            VariantCrossingIons = psm.VariantCrossingIons?.ToList();
            CrossType = psm.CrossType;
            LinkResidues = psm.LinkResidues;
            ProteinLinkSite = psm.ProteinLinkSite;
            Rank = psm.Rank;
            BetaPeptideProteinAccession = psm.BetaPeptideProteinAccession;
            BetaPeptideProteinLinkSite = psm.BetaPeptideProteinLinkSite;
            BetaPeptideBaseSequence = psm.BetaPeptideBaseSequence;
            BetaPeptideFullSequence = psm.BetaPeptideFullSequence;
            BetaPeptideTheoreticalMass = psm.BetaPeptideTheoreticalMass;
            BetaPeptideScore = psm.BetaPeptideScore;
            BetaPeptideRank = psm.BetaPeptideRank;
            BetaPeptideMatchedIons = psm.BetaPeptideMatchedIons?.ToList();
            BetaPeptideChildScanMatchedIons = psm.BetaPeptideChildScanMatchedIons;
            XLTotalScore = psm.XLTotalScore;
            ParentIons = psm.ParentIons;
            RetentionTime = psm.RetentionTime;
            GlycanStructure = psm.GlycanStructure;
            GlycanMass = psm.GlycanMass;
            GlycanComposition = psm.GlycanComposition;
            GlycanLocalizationLevel = psm.GlycanLocalizationLevel;
            LocalizedGlycan = psm.LocalizedGlycan;
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

        /// <summary>
        /// Parses the full sequence to identify mods
        /// </summary>
        /// <param name="fullSequence"> Full sequence of the peptide in question</param>
        /// <returns> Dictionary with the key being the amino acid position of the mod and the value being the string representing the mod</returns>
        public static Dictionary<int, List<string>> ParseModifications(string fullSeq)
        {
            // use a regex to get all modifications
            string pattern = @"\[(.+?)\]";
            Regex regex = new(pattern);

            // remove each match after adding to the dict. Otherwise, getting positions
            // of the modifications will be rather difficult.
            //int patternMatches = regex.Matches(fullSeq).Count;
            Dictionary<int, List<string>> modDict = new();
            
            RemoveSpecialCharacters(ref fullSeq);
            MatchCollection matches = regex.Matches(fullSeq);
            int currentPosition = 0;
            foreach (Match match in matches)
            {
                GroupCollection group = match.Groups;
                string val = group[1].Value;
                int startIndex = group[0].Index;
                int captureLength = group[0].Length;
                int position = group["(.+?)"].Index;

                List<string> modList = new List<string>();
                modList.Add(val);
                // check to see if key already exist
                // if there is a missed cleavage, then there will be a label on K and a Label on X modification.
                // And, it'll be like [label]|[label] which complicates the positional stuff a little bit.
                // if the already key exists, update the current position with the capture length + 1.
                // otherwise, add the modification to the dict.

                // int to add is startIndex - current position
                int positionToAddToDict = startIndex - currentPosition;
                if (modDict.ContainsKey(positionToAddToDict))
                {
                    modDict[positionToAddToDict].Add(val);
                }
                else
                {
                    modDict.Add(positionToAddToDict, modList);
                }
                currentPosition += startIndex + captureLength;
            }
            return modDict;
        }

        /// <summary>
        /// Fixes an issue where the | appears and throws off the numbering if there are multiple mods on a single amino acid.
        /// </summary>
        /// <param name="fullSeq"></param>
        /// <param name="replacement"></param>
        /// <param name="specialCharacter"></param>
        /// <returns></returns>
        public static void RemoveSpecialCharacters(ref string fullSeq, string replacement = @"", string specialCharacter = @"\|")
        {
            // next regex is used in the event that multiple modifications are on a missed cleavage Lysine (K)
            Regex regexSpecialChar = new(specialCharacter);
            fullSeq = regexSpecialChar.Replace(fullSeq, replacement);
        }

        public static List<MatchedFragmentIon> ReadFragmentIonsFromList(List<string> peakMzs, List<string> peakIntensities,
            string peptideBaseSequence, List<string> peakMassErrorDa = null, bool ignoreArtifactIons = false)
        {
            List<MatchedFragmentIon> matchedIons = new List<MatchedFragmentIon>();

            for (int index = 0; index < peakMzs.Count; index++)
            {
                string peak = peakMzs[index];
                string[] split = peak.Split(new char[] { '+', ':' }); //TODO: needs update for negative charges that doesn't break internal fragment ions or neutral losses

                // if there is a mismatch between the number of peaks and number of intensities from the psmtsv, the intensity will be set to 1
                double intensity = peakMzs.Count == peakIntensities.Count ? //TODO: needs update for negative charges that doesn't break internal fragment ions or neutral losses
                    double.Parse(peakIntensities[index].Split(new char[] { '+', ':', ']' })[2], CultureInfo.InvariantCulture) :
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
                    // if there is no mismatch between intensity and peak counts from the psmtsv
                    if (!intensity.Equals(1.0))
                    {
                        intensity = double.Parse(peakIntensities[index].Split(new char[] { '+', ':', ']' })[3],
                            CultureInfo.InvariantCulture);
                    }
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
                        if (ignoreArtifactIons) continue;
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
                        peptideBaseSequence.Split('|')[0].Length - fragmentNumber :
                        fragmentNumber;
                }

                //get mass error in Daltons
                double errorDa = 0;
                if (peakMassErrorDa.IsNotNullOrEmpty() && peakMassErrorDa[index].IsNotNullOrEmpty())
                {
                    string peakError = peakMassErrorDa[index];
                    string[] errorSplit = peakError.Split(new char[] { '+', ':', ']' });
                    errorDa = double.Parse(errorSplit[2], CultureInfo.InvariantCulture);
                }

                //get charge and mz
                int z = int.Parse(split[1]);
                double mz = double.Parse(split[2], CultureInfo.InvariantCulture);
                double neutralExperimentalMass = mz.ToMass(z); //read in m/z converted to mass
                double neutralTheoreticalMass = neutralExperimentalMass - errorDa; //theoretical mass is measured mass - measured error

                //The product created here is the theoretical product, with the mass back-calculated from the measured mass and measured error
                Product theoreticalProduct = new Product(productType,
                  terminus,
                  neutralTheoreticalMass,
                  fragmentNumber,
                  aminoAcidPosition,
                  neutralLoss,
                  secondaryProductType,
                  secondaryFragmentNumber);

                matchedIons.Add(new MatchedFragmentIon(ref theoreticalProduct, mz, intensity, z));
            }
        
            return matchedIons;
        }

        public static List<string> BuildFragmentStringsFromMaxQuant(string fragmentString, string informationString)
        {
            List<string> fragments = new List<string>();
            string[] fragmentName = fragmentString.Split(';');
            string[] fragmentInfo = informationString.Split(';');
            Regex fragmentRegex = new Regex(@"[^\(\-]*");
            
            for (int i = 0; i < fragmentName.Length; i++)
            {
                string charge = "1";
                string fragment = fragmentName[i];
                if (fragmentName[i].Contains('('))
                {
                    Regex chargeRegex = new Regex(@"\((\d*)\+\)");
                    charge = chargeRegex.Match(fragmentName[i]).Groups[1].Value;
                    fragment = fragmentRegex.Match(fragmentName[i]).Groups[0].Value;
                } 
                if (fragmentName[i].Contains('-'))
                {
                    Regex neutralLossRegex = new Regex(@"\-([A-Z\d]+)");
                    ChemicalFormula neutralLoss = ChemicalFormula.ParseFormula(
                        neutralLossRegex.Match(fragmentName[i]).Groups[1].Value);
                    fragment = '(' + fragmentRegex.Match(fragmentName[i]).Groups[0].Value + 
                               '-' + neutralLoss.MonoisotopicMass.ToString() + ')';
                }
                fragments.Add(fragment + '+' + charge + ':' + fragmentInfo[i]);
            }

            return fragments;
        }

        private static List<MatchedFragmentIon> ReadFragmentIonsFromString(string matchedMzString, string matchedIntensityString, string peptideBaseSequence, string matchedMassErrorDaString = null)
        {
            List<MatchedFragmentIon> matchedIons = new List<MatchedFragmentIon>();

            if (matchedMzString.Length > 2) //check if there's an ion
            {
                List<string> peakMzs = CleanMatchedIonString(matchedMzString);
                List<string> peakIntensities = CleanMatchedIonString(matchedIntensityString);
                List<string> peakMassErrorDa = null;

                if (matchedMassErrorDaString.IsNotNullOrEmpty())
                {
                    peakMassErrorDa = CleanMatchedIonString(matchedMassErrorDaString);
                }

                for (int index = 0; index < peakMzs.Count; index++)
                {
                    string peak = peakMzs[index];
                    string[] split = peak.Split(new char[] { '+', ':' }); //TODO: needs update for negative charges that doesn't break internal fragment ions or neutral losses
                    
                    // if there is a mismatch between the number of peaks and number of intensities from the psmtsv, the intensity will be set to 1
                    double intensity = peakMzs.Count == peakIntensities.Count ? //TODO: needs update for negative charges that doesn't break internal fragment ions or neutral losses
                        double.Parse(peakIntensities[index].Split(new char[] { '+', ':', ']' })[2], CultureInfo.InvariantCulture) :
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
                        // if there is no mismatch between intensity and peak counts from the psmtsv
                        if (!intensity.Equals(1.0))
                        {
                            intensity = double.Parse(peakIntensities[index].Split(new char[] { '+', ':', ']' })[3],
                                CultureInfo.InvariantCulture);
                        }
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
                            peptideBaseSequence.Split('|')[0].Length - fragmentNumber :
                            fragmentNumber;
                    }

                    //get mass error in Daltons
                    double errorDa = 0; 
                    if (matchedMassErrorDaString.IsNotNullOrEmpty() && peakMassErrorDa[index].IsNotNullOrEmpty())
                    {
                        string peakError = peakMassErrorDa[index];
                        string[] errorSplit = peakError.Split(new char[] { '+', ':', ']' });
                        errorDa = double.Parse(errorSplit[2], CultureInfo.InvariantCulture);
                    }

                    //get charge and mz
                    int z = int.Parse(split[1]);
                    double mz = double.Parse(split[2], CultureInfo.InvariantCulture);
                    double neutralExperimentalMass = mz.ToMass(z); //read in m/z converted to mass
                    double neutralTheoreticalMass = neutralExperimentalMass - errorDa; //theoretical mass is measured mass - measured error

                    //The product created here is the theoretical product, with the mass back-calculated from the measured mass and measured error
                    Product theoreticalProduct = new Product(productType,
                      terminus,
                      neutralTheoreticalMass,
                      fragmentNumber,
                      aminoAcidPosition,
                      neutralLoss,
                      secondaryProductType,
                      secondaryFragmentNumber);

                    matchedIons.Add(new MatchedFragmentIon(ref theoreticalProduct, mz, intensity, z));
                }
            }
            return matchedIons;
        }
        /// <summary>
        /// Removes enclosing brackets and
        /// replaces delimimiters between ion series with comma
        /// then splits on comma
        /// </summary>
        /// <param name="input"> String containing ion series from .psmtsv </param>
        /// <returns> List of strings, with each entry containing one ion and associated property </returns>
        private static List<string> CleanMatchedIonString(string input)
        {
            List<string> ionProperty = input.Substring(1, input.Length - 2) 
                    .Replace("];[", ", ") 
                    .Split(", ") 
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
                                                                    ProductType.b, ProductType.bWaterLoss, ProductType.bAmmoniaLoss, ProductType.c };
                    List<ProductType> xyzProductTypes = new List<ProductType>() { ProductType.x, ProductType.y, ProductType.yAmmoniaLoss,
                                                                    ProductType.yWaterLoss, ProductType.zDot, ProductType.zPlusOne};
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

        public override string ToString()
        {
            return FullSequence;
        }
    }
}