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
        public double XLTotalScore { get; set; } //alpha + beta psmCross
        public int XlProteinPos { get; set; }
        public List<int> XlRank { get; set; } //only contain 2 intger, consider change to Tuple
        public string ParentIonExist { get; set; }
        public int ParentIonExistNum { get; set; }
        public List<int> ParentIonMaxIntensityRanks { get; set; }
        public PsmCrossType CrossType { get; set; }

        public static List<int> GetPossibleCrosslinkerModSites(char[] crosslinkerModSites, PeptideWithSetModifications peptide)
        {
            List<int> possibleXlPositions = new List<int>();
            bool wildcard = crosslinkerModSites.Any(p => p == 'X');

            for (int r = 0; r < peptide.BaseSequence.Length; r++)
            {
                if (crosslinkerModSites.Contains(peptide.BaseSequence[r]) || wildcard)
                {
                    possibleXlPositions.Add(r + 1);
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

            sb.Append("Beta Peptide" + '\t');
            sb.Append("Beta Peptide Protein Accession" + '\t');
            sb.Append("Beta Peptide Protein LinkSite" + '\t');
            sb.Append("Beta Peptide Base Sequence" + '\t');
            sb.Append("Beta Peptide Full Sequence" + '\t');
            sb.Append("Beta Peptide Theoretical Mass" + '\t');
            sb.Append("Beta Peptide Score" + '\t');
            sb.Append("Beta Peptide Rank" + '\t');

            sb.Append("Beta Peptide Matched Ion Series" + '\t');
            sb.Append("Beta Peptide Matched Ion Mass-To-Charge Ratios" + '\t');
            sb.Append("Beta Peptide Matched Ion Mass Diff (Da)" + '\t');
            sb.Append("Beta Peptide Matched Ion Mass Diff (Ppm)" + '\t');
            sb.Append("Beta Peptide Matched Ion Intensities" + '\t');
            sb.Append("Beta Peptide Matched Ion Counts" + '\t');

            sb.Append("Summary" + '\t');
            sb.Append("XL Total Score" + '\t');
            sb.Append("Mass Diff (Da)" + '\t');
            sb.Append("Parent Ions" + '\t');
            sb.Append("ParentIonsNum" + '\t');
            sb.Append("ParentIonMaxIntensityRank" + '\t');
            sb.Append("Decoy/Contaminant/Target" + '\t');
            sb.Append("QValue" + '\t');
       

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
            sb.Append("Decoy/Contaminant/Target" + '\t');
            sb.Append("QValue" + '\t');

            return sb.ToString();
        }

        public static string GetTabSepHeaderGlyco()
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

            sb.Append("Decoy/Contaminant/Target" + '\t');
            sb.Append("QValue" + '\t');

            sb.Append("GlyID" + '\t');
            sb.Append("GlyMass" + '\t');
            sb.Append("GlyStruct(H,N,A,G,F)" + '\t');
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
            sb.Append(ProteinAccession + "\t");
            sb.Append(XlProteinPos + "\t");
            sb.Append(BaseSequence + "\t");
            sb.Append(FullSequence + position + "\t");
            sb.Append((PeptideMonisotopicMass.HasValue ? PeptideMonisotopicMass.Value.ToString() : "---")); sb.Append("\t");
            sb.Append(Score + "\t");
            sb.Append(XlRank[0] + "\t");

            foreach (var mid in MatchedIonDataDictionary(this))
            {
                sb.Append(mid.Value);
                sb.Append("\t");
            }

            if (BetaPeptide != null)
            {
                sb.Append("\t");
                sb.Append(BetaPeptide.ProteinAccession + "\t");
                sb.Append(BetaPeptide.XlProteinPos + "\t");
                sb.Append(BetaPeptide.BaseSequence + "\t");
                sb.Append(BetaPeptide.FullSequence + "(" + BetaPeptide.LinkPositions[0].ToString() + ")" + "\t");
                sb.Append(BetaPeptide.PeptideMonisotopicMass.ToString() + "\t");
                sb.Append(BetaPeptide.Score + "\t");
                sb.Append(XlRank[1] + "\t");

                foreach (var betamid in MatchedIonDataDictionary(this.BetaPeptide))
                {
                    sb.Append(betamid.Value);
                    sb.Append("\t");
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

        public static Dictionary<string, string> MatchedIonDataDictionary(PeptideSpectralMatch psm)
        {
            Dictionary<string, string> s = new Dictionary<string, string>();
            PsmTsvWriter.AddMatchedIonsData(s, psm);
            return s;
        }
    }
}