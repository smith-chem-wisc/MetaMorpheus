using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace EngineLayer.CrosslinkSearch
{
    public class CrosslinkSpectralMatch : PeptideSpectralMatch
    {
        public CrosslinkSpectralMatch(PeptideWithSetModifications theBestPeptide, int notch, double score, int scanIndex, Ms2ScanWithSpecificMass scan, DigestionParams digestionParams, List<MatchedFragmentIon> matchedFragmentIons)
            : base(theBestPeptide, notch, score, scanIndex, scan, digestionParams, matchedFragmentIons)
        {
            this.XLTotalScore = score;
        }

        public CrosslinkSpectralMatch BetaPeptide { get; set; }
        public List<int> LinkPositions { get; set; }
        public double DeltaScore { get; set; }
        public double XLTotalScore { get; set; } //alpha + beta psmCross.
        public List<int> XlRank { get; set; } //only contain 2 intger, consider change to Tuple.
        public string ParentIonExist { get; set; }
        public int ParentIonExistNum { get; set; }
        public List<int> ParentIonMaxIntensityRanks { get; set; }
        public PsmCrossType CrossType { get; set; }
        public Dictionary<int, List<MatchedFragmentIon>> ChildMatchedFragmentIons { get; set; }
        public int? XlProteinPos { get; private set; }
        // loop crosslink protein position 2.
        public int? XlProteinPosLoop { get; private set; }

        public bool IsIntraCsm()
        {
            
            if (this.ProteinAccession != null && this.BetaPeptide.ProteinAccession != null)
            {
                if (this.ProteinAccession == this.BetaPeptide.ProteinAccession)
                {
                    return true;
                }
            }

            if (this.ProteinAccession == null)
            {
                var alphaProteins = BestMatchingPeptides.Select(p => p.Peptide.Protein.Accession).ToList();
                var betaProteins = BetaPeptide.BestMatchingPeptides.Select(p => p.Peptide.Protein.Accession).ToList();

                foreach (var alpha in alphaProteins)
                {
                    foreach (var beta in betaProteins)
                    {
                        if (alpha == beta)
                        {
                            return true;
                        }
                    }

                }
            }
            return false;
        }

        public void ResolveProteinPosAmbiguitiesForXl()
        {
            if (CrossType == PsmCrossType.Cross)
            {
                // alpha peptide crosslink residue in the protein
                XlProteinPos = OneBasedStartResidueInProtein == null ? (int?)null : OneBasedStartResidueInProtein.Value + LinkPositions[0] - 1;

                // beta crosslink residue in protein
                BetaPeptide.XlProteinPos = BetaPeptide.OneBasedStartResidueInProtein == null ? (int?)null : BetaPeptide.OneBasedStartResidueInProtein.Value + BetaPeptide.LinkPositions[0] - 1;
            }
            else if (CrossType == PsmCrossType.DeadEnd || CrossType == PsmCrossType.DeadEndH2O || CrossType == PsmCrossType.DeadEndNH2 || CrossType == PsmCrossType.DeadEndTris)
            {
                XlProteinPos = OneBasedStartResidueInProtein == null ? (int?)null : OneBasedStartResidueInProtein.Value + LinkPositions[0] - 1;
            }
            else if (CrossType == PsmCrossType.Loop)
            {
                XlProteinPos = OneBasedStartResidueInProtein == null ? (int?)null : OneBasedStartResidueInProtein.Value + LinkPositions[0] - 1;

                XlProteinPosLoop = OneBasedStartResidueInProtein == null ? (int?)null : OneBasedStartResidueInProtein.Value + LinkPositions[1] - 1;
            }
        }

        public static List<int> GetPossibleCrosslinkerModSites(char[] crosslinkerModSites, PeptideWithSetModifications peptide)
        {
            List<int> possibleXlPositions = new List<int>();
            bool wildcard = crosslinkerModSites.Any(p => p == 'X');

            for (int r = 0; r < peptide.BaseSequence.Length; r++)
            {
                if (crosslinkerModSites.Contains(peptide.BaseSequence[r]) || wildcard)
                {
                    //Try to eliminate those site with mod on it. Consider the possibility that the site is at Protein N terminal.       
                    if (!peptide.AllModsOneIsNterminus.Keys.Contains(r + 2))
                    {
                        possibleXlPositions.Add(r + 1);
                    }
                    else if (peptide.OneBasedStartResidueInProtein == 0 && r == 0 && !peptide.AllModsOneIsNterminus.Keys.Contains(1))
                    {
                        possibleXlPositions.Add(r + 1);
                    }
                }
            }
            return possibleXlPositions;
        }

        /// <summary>
        /// Rank experimental mass spectral peaks by intensity
        /// </summary>
        public static int[] GenerateIntensityRanks(double[] experimental_intensities)
        {
            var y = experimental_intensities.ToArray();
            var x = Enumerable.Range(1, y.Length).OrderBy(p => p).ToArray();
            Array.Sort(y, x);
            var experimental_intensities_rank = Enumerable.Range(1, y.Length).OrderByDescending(p => p).ToArray();
            Array.Sort(x, experimental_intensities_rank);
            return experimental_intensities_rank;
        }

        public static string GetTabSepHeaderCross()
        {
            var sb = new StringBuilder();
            sb.Append(PsmTsvHeader.FileName + '\t');
            sb.Append(PsmTsvHeader.Ms2ScanNumber + '\t');
            sb.Append(PsmTsvHeader.PrecursorScanNum + '\t');
            sb.Append(PsmTsvHeader.PrecursorMz + '\t');
            sb.Append(PsmTsvHeader.PrecursorCharge + '\t');
            sb.Append(PsmTsvHeader.PrecursorMass + '\t');
            sb.Append(PsmTsvHeader.CrossTypeLabel + '\t');
            sb.Append(PsmTsvHeader.LinkResiduesLabel + "\t");

            sb.Append("Peptide" + '\t');
            sb.Append(PsmTsvHeader.ProteinAccession + '\t');
            sb.Append(PsmTsvHeader.ProteinLinkSiteLabel + '\t');
            sb.Append(PsmTsvHeader.BaseSequence + '\t');
            sb.Append(PsmTsvHeader.FullSequence + '\t');
            sb.Append(PsmTsvHeader.PeptideMonoMass + '\t');
            sb.Append(PsmTsvHeader.Score + '\t');
            sb.Append(PsmTsvHeader.RankLabel + '\t');

            sb.Append(PsmTsvHeader.MatchedIonSeries + '\t');
            sb.Append(PsmTsvHeader.MatchedIonMzRatios + '\t');
            sb.Append(PsmTsvHeader.MatchedIonMassDiffDa + '\t');
            sb.Append(PsmTsvHeader.MatchedIonMassDiffPpm + '\t');
            sb.Append(PsmTsvHeader.MatchedIonIntensities + '\t');
            sb.Append(PsmTsvHeader.MatchedIonCounts + '\t');

            sb.Append("Beta Peptide" + '\t');
            sb.Append(PsmTsvHeader.BetaPeptideProteinAccessionLabel + '\t');
            sb.Append(PsmTsvHeader.BetaPeptideProteinLinkSiteLabel + '\t');
            sb.Append(PsmTsvHeader.BetaPeptideBaseSequenceLabel + '\t');
            sb.Append(PsmTsvHeader.BetaPeptideFullSequenceLabel + '\t');
            sb.Append(PsmTsvHeader.BetaPeptideTheoreticalMassLabel + '\t');
            sb.Append(PsmTsvHeader.BetaPeptideScoreLabel + '\t');
            sb.Append(PsmTsvHeader.BetaPeptideRankLabel + '\t');

            sb.Append("Beta Peptide Matched Ions" + '\t');
            sb.Append(PsmTsvHeader.BetaPeptideMatchedIonsLabel + '\t');
            sb.Append("Beta Peptide Matched Ion Mass Diff (Da)" + '\t');
            sb.Append("Beta Peptide Matched Ion Mass Diff (Ppm)" + '\t');
            sb.Append("Beta Peptide Matched Ion Intensities" + '\t');
            sb.Append("Beta Peptide Matched Ion Counts" + '\t');

            sb.Append("Summary" + '\t');
            sb.Append(PsmTsvHeader.XLTotalScoreLabel + '\t');
            sb.Append(PsmTsvHeader.MassDiffDa + '\t');
            sb.Append(PsmTsvHeader.ParentIonsLabel + '\t');
            sb.Append("ParentIonsNum" + '\t');
            sb.Append("ParentIonMaxIntensityRank" + '\t');
            sb.Append(PsmTsvHeader.DecoyContaminantTarget + '\t');
            sb.Append(PsmTsvHeader.QValue + '\t');
            
            return sb.ToString();
        }

        public static string GetTabSepHeaderSingle()
        {
            var sb = new StringBuilder();
            sb.Append("File Name" + '\t');
            sb.Append("Scan Number" + '\t');
            sb.Append("Precursor Scan Number" + '\t');
            sb.Append("Precursor MZ" + '\t');
            sb.Append("Precursor Charge" + '\t');
            sb.Append("Precursor Mass" + '\t');
            sb.Append("Cross Type" + '\t');
            sb.Append("Link Residues" + "\t");

            sb.Append("Peptide" + '\t');
            sb.Append("Protein Accession" + '\t');
            sb.Append("Protein Link Site" + '\t');
            sb.Append("Base Sequence" + '\t');
            sb.Append("Full Sequence" + '\t');
            sb.Append("Peptide Monoisotopic Mass" + '\t');
            sb.Append("Score" + '\t');
            sb.Append("Rank" + '\t');

            sb.Append("Matched Ion Series" + '\t');
            sb.Append("Matched Ion Mass-To-Charge Ratios" + '\t');
            sb.Append("Matched Ion Mass Diff (Da)" + '\t');
            sb.Append("Matched Ion Mass Diff (Ppm)" + '\t');
            sb.Append("Matched Ion Intensities" + '\t');
            sb.Append("Matched Ion Counts" + '\t');
            sb.Append("Child Scans Matched Ion Series" + '\t');
            sb.Append("Decoy/Contaminant/Target" + '\t');
            sb.Append("QValue" + '\t');

            return sb.ToString();
        }

        public override string ToString()
        {
            string position = "";
            switch (CrossType)
            {
                case PsmCrossType.Single:
                    break;

                case PsmCrossType.Loop:
                    position = "(" + LinkPositions[0].ToString() + "-" + LinkPositions[1].ToString() + ")";
                    break;

                default:
                    position = "(" + LinkPositions[0].ToString() + ")";
                    break;
            }

            var sb = new StringBuilder();
            sb.Append(FullFilePath + "\t");
            sb.Append(ScanNumber + "\t");
            sb.Append(PrecursorScanNumber + "\t");
            sb.Append(ScanPrecursorMonoisotopicPeakMz + "\t");
            sb.Append(ScanPrecursorCharge + "\t");
            sb.Append(ScanPrecursorMass + "\t");
            sb.Append(CrossType.ToString() + "\t");

            if (LinkPositions != null)
            {
                if (CrossType == PsmCrossType.Loop)
                {
                    sb.Append(BaseSequence[LinkPositions[0] - 1] + ";" + BaseSequence[LinkPositions[1] - 1] + "\t");
                }
                else if (CrossType == PsmCrossType.Inter || CrossType == PsmCrossType.Intra)
                {
                    sb.Append(BaseSequence[LinkPositions[0] - 1] + ";" + BetaPeptide.BaseSequence[BetaPeptide.LinkPositions[0] - 1] + "\t");
                }
                else
                {
                    // deadend
                    sb.Append(BaseSequence[LinkPositions[0] - 1] + "\t");
                }
            }
            else
            {
                sb.Append("\t");
            }

            sb.Append("\t");
            List<PeptideWithSetModifications> pepsWithMods = BestMatchingPeptides.Select(p => p.Peptide).ToList();
            var proteinAccessionString = ProteinAccession != null ? ProteinAccession : PsmTsvWriter.Resolve(pepsWithMods.Select(b => b.Protein.Accession), FullSequence).ResolvedString;
            sb.Append(proteinAccessionString + "\t");           
            sb.Append(XlProteinPos + (XlProteinPosLoop.HasValue? "~"+ XlProteinPosLoop.Value : null) + "\t");
            sb.Append(BaseSequence + "\t");
            sb.Append(FullSequence + position + "\t");
            sb.Append((PeptideMonisotopicMass.HasValue ? PeptideMonisotopicMass.Value.ToString() : "---")); sb.Append("\t");
            sb.Append(Score + "\t");
            sb.Append(XlRank[0] + "\t");

            if (ChildMatchedFragmentIons == null)
            {
                foreach (var mid in MatchedIonDataDictionary(this.MatchedFragmentIons))
                {
                    sb.Append(mid.Value);
                    sb.Append("\t");
                }
            }
            else
            {
                StringBuilder[] scanFragmentStringbuilder = new StringBuilder[6];
                int i =0;
                foreach (var mid in MatchedIonDataDictionary(this.MatchedFragmentIons))
                {
                    scanFragmentStringbuilder[i] = new StringBuilder();
                    scanFragmentStringbuilder[i].Append("{" + ScanNumber + "@" + mid.Value + "}");
                    i++;
                }
                foreach (var childScan in ChildMatchedFragmentIons)
                {
                    int j = 0;
                    int oneBasedScan = childScan.Key;
                    foreach (var mid in MatchedIonDataDictionary(childScan.Value))
                    {
                        scanFragmentStringbuilder[j].Append("{" + oneBasedScan + "@" + mid.Value + "}");
                        j++;
                    }

                }
                foreach (var s in scanFragmentStringbuilder)
                {
                    sb.Append(s.ToString() + "\t");
                }       
            }


            if (BetaPeptide != null)
            {
                sb.Append("\t");
                List<PeptideWithSetModifications> betaPepsWithMods = BetaPeptide.BestMatchingPeptides.Select(p => p.Peptide).ToList();
                var betaProteinAccessionString = BetaPeptide.ProteinAccession != null ? BetaPeptide.ProteinAccession : PsmTsvWriter.Resolve(betaPepsWithMods.Select(b => b.Protein.Accession), FullSequence).ResolvedString;
                sb.Append(betaProteinAccessionString + "\t");
                sb.Append(BetaPeptide.XlProteinPos + "\t");
                sb.Append(BetaPeptide.BaseSequence + "\t");
                sb.Append(BetaPeptide.FullSequence + "(" + BetaPeptide.LinkPositions[0].ToString() + ")" + "\t");
                sb.Append(BetaPeptide.PeptideMonisotopicMass.ToString() + "\t");
                sb.Append(BetaPeptide.Score + "\t");
                sb.Append(XlRank[1] + "\t");

                if (BetaPeptide.ChildMatchedFragmentIons == null)
                {
                    foreach (var betamid in MatchedIonDataDictionary(this.BetaPeptide.MatchedFragmentIons))
                    {
                        sb.Append(betamid.Value);
                        sb.Append("\t");
                    }
                }
                else
                {
                    StringBuilder[] scanFragmentStringbuilder = new StringBuilder[6];
                    int i = 0;
                    foreach (var betamid in MatchedIonDataDictionary(this.BetaPeptide.MatchedFragmentIons))
                    {
                        scanFragmentStringbuilder[i] = new StringBuilder();
                        scanFragmentStringbuilder[i].Append("{" + ScanNumber + "@" + betamid.Value + "}");
                        i++;
                    }
                    foreach (var betaChildScan in BetaPeptide.ChildMatchedFragmentIons)
                    {
                        int j = 0;
                        int betaOneBasedScan = betaChildScan.Key;
                        foreach (var betamid in MatchedIonDataDictionary(betaChildScan.Value))
                        {
                            scanFragmentStringbuilder[j].Append("{" + betaOneBasedScan + "@" + betamid.Value + "}");
                            j++;
                        }
                    }
                    foreach (var s in scanFragmentStringbuilder)
                    {
                        sb.Append(s.ToString() + "\t");
                    }
                }

                sb.Append("\t");
                sb.Append(XLTotalScore + "\t");

                // mass of crosslinker
                sb.Append(((PeptideMonisotopicMass.HasValue) ? (ScanPrecursorMass - BetaPeptide.PeptideMonisotopicMass - PeptideMonisotopicMass.Value).ToString() : "---")); sb.Append("\t");

                int alphaNumParentIons = MatchedFragmentIons.Count(p => p.NeutralTheoreticalProduct.ProductType == ProductType.M);
                int betaNumParentIons = BetaPeptide.MatchedFragmentIons.Count(p => p.NeutralTheoreticalProduct.ProductType == ProductType.M);

                sb.Append(alphaNumParentIons + ";" + betaNumParentIons + "\t");
                sb.Append(alphaNumParentIons + betaNumParentIons + "\t");
                sb.Append(((ParentIonMaxIntensityRanks != null) && (ParentIonMaxIntensityRanks.Any()) ? ParentIonMaxIntensityRanks.Min().ToString() : "-")); sb.Append("\t");
            }

            if (BetaPeptide == null)
            {
                sb.Append((IsDecoy) ? "D" : (IsContaminant) ? "C" : "T");
                sb.Append("\t");
            }
            else
            {
                sb.Append((IsDecoy || BetaPeptide.IsDecoy) ? "D" : (IsContaminant || BetaPeptide.IsContaminant) ? "C" : "T");
                sb.Append("\t");
            }

            sb.Append(FdrInfo.QValue.ToString());
            sb.Append("\t");
            
            return sb.ToString();
        }

        public static Dictionary<string, string> MatchedIonDataDictionary(List<MatchedFragmentIon> matchedFragmentIons)
        {
            Dictionary<string, string> s = new Dictionary<string, string>();
            PsmTsvWriter.AddMatchedIonsData(s, matchedFragmentIons);
            return s;
        }
    }
}