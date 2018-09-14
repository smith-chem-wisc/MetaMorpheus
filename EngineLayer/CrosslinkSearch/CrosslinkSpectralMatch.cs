using Chemistry;
using MassSpectrometry;
using MzLibUtil;
using Proteomics.AminoAcidPolymer;
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
        private static readonly double waterMonoisotopicMass = PeriodicTable.GetElement("H").PrincipalIsotope.AtomicMass * 2 + PeriodicTable.GetElement("O").PrincipalIsotope.AtomicMass;
        private static readonly double nitrogenAtomMonoisotopicMass = PeriodicTable.GetElement("N").PrincipalIsotope.AtomicMass;
        private static readonly double oxygenAtomMonoisotopicMass = PeriodicTable.GetElement("O").PrincipalIsotope.AtomicMass;
        private static readonly double hydrogenAtomMonoisotopicMass = PeriodicTable.GetElement("H").PrincipalIsotope.AtomicMass;

        public CrosslinkSpectralMatch(PeptideWithSetModifications theBestPeptide, int notch, double score, int scanIndex, Ms2ScanWithSpecificMass scan, DigestionParams digestionParams, List<MatchedFragmentIon> matchedFragmentIons)
            : base(theBestPeptide, notch, score, scanIndex, scan, digestionParams, matchedFragmentIons)
        {
            this.BestScore = score;
            this.XLTotalScore = score;
        }

        public CrosslinkSpectralMatch BetaPeptide { get; set; }
        public double BestScore { get; set; } //For the current psmCross
        public List<int> ModPositions { get; set; }
        public double DeltaScore { get; set; }
        public double XLTotalScore { get; set; } //alpha + beta psmCross
        public double XLQvalueTotalScore { get; set; } //Calc based on XLtotalScore for Qvalue
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

        public void GetBestMatch(Ms2ScanWithSpecificMass theScan, Dictionary<List<int>, List<Product>> productsForEachModPosition, CommonParameters commonParameters)
        {
            foreach (var listOfProducts in productsForEachModPosition)
            {
                List<MatchedFragmentIon> matchedIons = MetaMorpheusEngine.MatchFragmentIons(theScan.TheScan.MassSpectrum, listOfProducts.Value, commonParameters);
                double score = MetaMorpheusEngine.CalculatePeptideScore(theScan.TheScan, matchedIons, 0);

                if (score > BestScore)
                {
                    BestScore = score;
                    MatchedFragmentIons = matchedIons;
                    ModPositions = listOfProducts.Key;
                }
            }

            if (MatchedFragmentIons != null)
            {
                double[] experimental_intensities = theScan.TheScan.MassSpectrum.YArray;
                int[] experimental_intensities_rank = GenerateIntensityRanks(experimental_intensities);
                foreach (var matchedIon in MatchedFragmentIons)
                {
                    // get the closest peak in the spectrum to the theoretical peak
                    int matchedPeakIndex = theScan.TheScan.MassSpectrum.GetClosestPeakIndex(matchedIon.Mz).Value;

                    // TODO: add intensity rank in mzLib
                    //matchedIon.IntensityRank = experimental_intensities_rank[matchedPeakIndex];
                }
            }
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
            sb.Append("Scan Numer" + '\t');
            sb.Append("Precusor MZ" + '\t');
            sb.Append("Precusor charge" + '\t');
            sb.Append("Precusor mass" + '\t');
            sb.Append("CrossType" + '\t');

            sb.Append("Pep1" + '\t');
            sb.Append("Pep1 Protein Access" + '\t');
            sb.Append("Protein link site" + '\t');
            sb.Append("Pep1 Base sequence(crosslink site)" + '\t');
            sb.Append("Pep1 Full sequence" + '\t');
            sb.Append("Pep1 mass" + '\t');
            sb.Append("Pep1 BestScore" + '\t');
            sb.Append("Pep1 Rank" + '\t');

            sb.Append("Pep2" + '\t');
            sb.Append("Pep2 Protein Access" + '\t');
            sb.Append("Protein link site" + '\t');
            sb.Append("Pep2 Base sequence(crosslink site)" + '\t');
            sb.Append("Pep2 Full sequence" + '\t');
            sb.Append("Pep2 mass" + '\t');
            sb.Append("Pep2 BestScore" + '\t');
            sb.Append("Pep2 Rank" + '\t');

            sb.Append("Summary" + '\t');
            sb.Append("QvalueTotalScore" + '\t');
            sb.Append("Mass diff" + '\t');
            sb.Append("ParentIons" + '\t');
            sb.Append("ParentIonsNum" + '\t');
            sb.Append("ParentIonMaxIntensityRank" + '\t');
            sb.Append("Charge2/3Number" + '\t');
            sb.Append("Target/Decoy" + '\t');
            sb.Append("QValue" + '\t');
            return sb.ToString();
        }

        public static string GetTabSepHeaderSingle()
        {
            var sb = new StringBuilder();
            sb.Append("File Name" + '\t');
            sb.Append("Scan Numer" + '\t');
            sb.Append("Precusor MZ" + '\t');
            sb.Append("Precusor charge" + '\t');
            sb.Append("Precusor mass" + '\t');
            sb.Append("CrossType" + '\t');

            sb.Append("Pep" + '\t');
            sb.Append("Pep Protein Access" + '\t');
            sb.Append("Protein site" + '\t');
            sb.Append("Pep Base sequence" + '\t');
            sb.Append("Pep Full sequence" + '\t');
            sb.Append("Pep mass" + '\t');
            sb.Append("Pep BestScore" + '\t');
            sb.Append("Pep Rank" + '\t');
            sb.Append("Charge2/3Number" + '\t');
            sb.Append("Target/Decoy" + '\t');
            sb.Append("QValue" + '\t');
            return sb.ToString();
        }

        public static string GetTabSepHeaderGlyco()
        {
            var sb = new StringBuilder();
            sb.Append("File Name" + '\t');
            sb.Append("Scan Numer" + '\t');
            sb.Append("Precusor MZ" + '\t');
            sb.Append("Precusor charge" + '\t');
            sb.Append("Precusor mass" + '\t');
            sb.Append("CrossType" + '\t');

            sb.Append("Pep" + '\t');
            sb.Append("Pep Protein Access" + '\t');
            sb.Append("Protein mod site" + '\t');
            sb.Append("Pep Base sequence" + '\t');
            sb.Append("Pep Full sequence" + '\t');
            sb.Append("Pep mass" + '\t');
            sb.Append("Pep BestScore" + '\t');
            sb.Append("Pep Rank" + '\t');
            sb.Append("Charge2/3Number" + '\t');
            sb.Append("Target/Decoy" + '\t');
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
                    position = "(" + ModPositions[0].ToString() + "-" + ModPositions[1].ToString() + ")";
                    break;

                default:
                    position = "(" + ModPositions[0].ToString() + ")";
                    break;
            }

            var sb = new StringBuilder();
            sb.Append(FullFilePath); sb.Append("\t");
            sb.Append(ScanNumber); sb.Append("\t");
            sb.Append(ScanPrecursorMonoisotopicPeakMz); sb.Append("\t");
            sb.Append(ScanPrecursorCharge); sb.Append("\t");
            sb.Append(ScanPrecursorMass); sb.Append("\t");
            sb.Append(CrossType.ToString()); sb.Append("\t");

            sb.Append(""); sb.Append("\t");
            sb.Append(BestMatchingPeptideWithSetMods.First().Pwsm.Protein.Accession); sb.Append("\t");
            sb.Append(XlProteinPos); sb.Append("\t");
            sb.Append(BaseSequence); sb.Append("\t");
            sb.Append(FullSequence + position); sb.Append("\t");
            sb.Append((PeptideMonisotopicMass.HasValue ? PeptideMonisotopicMass.Value.ToString() : "---")); sb.Append("\t");
            sb.Append(BestScore); sb.Append("\t");
            sb.Append(XlRank[0]); sb.Append("\t");

            if (BetaPeptide != null)
            {
                sb.Append(""); sb.Append("\t");
                sb.Append(BetaPeptide.ProteinAccession); sb.Append("\t");
                sb.Append(BetaPeptide.XlProteinPos); sb.Append("\t");
                sb.Append(BetaPeptide.BaseSequence); sb.Append("\t");
                sb.Append(BetaPeptide.FullSequence + "(" + ModPositions[0].ToString() + ")"); sb.Append("\t");
                sb.Append(BetaPeptide.PeptideMonisotopicMass.ToString()); sb.Append("\t");
                sb.Append(BetaPeptide.BestScore); sb.Append("\t");
                sb.Append(XlRank[1]); sb.Append("\t");

                sb.Append(""); sb.Append("\t");
                sb.Append(XLQvalueTotalScore); sb.Append("\t");

                // mass of crosslinker
                sb.Append(((PeptideMonisotopicMass.HasValue) ? (ScanPrecursorMass - BetaPeptide.PeptideMonisotopicMass - PeptideMonisotopicMass.Value).ToString() : "---")); sb.Append("\t");

                sb.Append(ParentIonExist + "." + BetaPeptide.ParentIonExist); sb.Append("\t");
                sb.Append(ParentIonExistNum); sb.Append("\t");
                sb.Append(((ParentIonMaxIntensityRanks != null) && (ParentIonMaxIntensityRanks.Any()) ? ParentIonMaxIntensityRanks.Min().ToString() : "-")); sb.Append("\t");
            }

            sb.Append(MatchedFragmentIons.Count(p => p.Charge > 1));
            sb.Append("\t");

            if (BetaPeptide == null)
            {
                sb.Append((IsDecoy) ? "-1" : "1");
                sb.Append("\t");
            }
            else
            {
                sb.Append((IsDecoy || BetaPeptide.IsDecoy) ? "-1" : "1");
                sb.Append("\t");
            }

            sb.Append((FdrInfo != null) ? FdrInfo.QValue.ToString() : "-");
            sb.Append("\t");

            return sb.ToString();
        }


    }
}