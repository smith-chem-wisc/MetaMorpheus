using Easy.Common.Extensions;
using Omics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Readers;

namespace EngineLayer.CrosslinkSearch
{
    public class CrosslinkSpectralMatch : PeptideSpectralMatch
    {
        public CrosslinkSpectralMatch(
            PeptideWithSetModifications theBestPeptide,
            int notch,
            double score,
            int scanIndex,
            Ms2ScanWithSpecificMass scan,
            CommonParameters commonParameters,
            List<MatchedFragmentIon> matchedFragmentIons)
            : base(theBestPeptide, notch, score, scanIndex, scan, commonParameters, matchedFragmentIons)
        {
            //The XLTotalScore is set here because some CSMs are not crosslinks and we need this score to be non-zero.
            XLTotalScore = score;
            _BestMatchingBioPolymersWithSetMods.Clear();

            _BestMatchingBioPolymersWithSetMods.Add(new SpectrumMatch.SpectralMatchHypothesis(0, theBestPeptide, matchedFragmentIons, score));

        }

        public CrosslinkSpectralMatch BetaPeptide { get; set; }
        /// <summary>
        /// UniqueSequence uniquely identifies crosslinked peptides to enable library lookup.
        /// For crosslinks, the unique sequence is identical to the spectral library name ( PEPTIDEK(8)EDITPEPK(8) ),
        /// where the numbers in parentheses represent the one-based residue where the crosslink occurs on each peptide.
        /// For non crosslinked peptides, UniqueSequence is identical to the full sequence.
        /// </summary>
        public string UniqueSequence
        {
            get
            {
                if (_uniqueSequence.IsNullOrEmpty())
                {
                    string position = "";
                    if (!LinkPositions.IsNotNullOrEmpty())
                    {
                        _uniqueSequence = FullSequence;
                        return _uniqueSequence;
                    }
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
                    if (BetaPeptide != null && BetaPeptide.LinkPositions.IsNotNullOrEmpty())
                    {
                        _uniqueSequence = FullSequence + position + BetaPeptide.FullSequence + "(" + BetaPeptide.LinkPositions[0].ToString() + ")";
                    }
                    else
                    {
                        _uniqueSequence = FullSequence + position;
                    }
                }
                return _uniqueSequence;
            }
        }
        private string _uniqueSequence;
        public List<int> LinkPositions { get; set; }
        public double XLTotalScore { get; set; }    //alpha + beta psmCross.

        public double SecondBestXlScore { get; set; } // score of the second-best CSM; this is used to calculate delta score
        public int XlRank { get; set; }   //Rank after indexing score. Could be used for PEP
        public List<int> ParentIonMaxIntensityRanks { get; set; }
        public PsmCrossType CrossType { get; set; }
        public double MS3ChildScore { get; set; }
        public Dictionary<int, List<MatchedFragmentIon>> ChildMatchedFragmentIons { get; set; }
        public int? XlProteinPos { get; private set; }

        // loop crosslink protein position 2.
        public int? XlProteinPosLoop { get; private set; }
        public new double DeltaScore { get { return (XLTotalScore - SecondBestXlScore); } }


        public void ResolveProteinPosAmbiguitiesForXl()
        {
            if (CrossType == PsmCrossType.Cross)
            {
                // alpha peptide crosslink residue in the protein
                XlProteinPos = OneBasedStartResidue == null ? (int?)null : OneBasedStartResidue.Value + LinkPositions[0] - 1;

                // beta crosslink residue in protein
                BetaPeptide.XlProteinPos = BetaPeptide.OneBasedStartResidue == null ? (int?)null : BetaPeptide.OneBasedStartResidue.Value + BetaPeptide.LinkPositions[0] - 1;
            }
            else if (CrossType == PsmCrossType.DeadEnd || CrossType == PsmCrossType.DeadEndH2O || CrossType == PsmCrossType.DeadEndNH2 || CrossType == PsmCrossType.DeadEndTris)
            {
                XlProteinPos = OneBasedStartResidue == null ? (int?)null : OneBasedStartResidue.Value + LinkPositions[0] - 1;
            }
            else if (CrossType == PsmCrossType.Loop)
            {
                XlProteinPos = OneBasedStartResidue == null ? (int?)null : OneBasedStartResidue.Value + LinkPositions[0] - 1;

                XlProteinPosLoop = OneBasedStartResidue == null ? (int?)null : OneBasedStartResidue.Value + LinkPositions[1] - 1;
            }
        }

        public static bool IsIntraCsm(CrosslinkSpectralMatch csm)
        {
            //The pair "ProteinA and Decoy_ProteinA" is count for intra-crosslink. 
            if (csm.Accession != null && csm.BetaPeptide.Accession != null)
            {
                if (csm.Accession == csm.BetaPeptide.Accession ||
                    csm.Accession == "DECOY_"+ csm.BetaPeptide.Accession ||
                    csm.BetaPeptide.Accession == "DECOY_" + csm.Accession)
                {
                    return true;
                }
            }

            if (csm.Accession == null)
            {
                var alphaProteins = csm.BestMatchingBioPolymersWithSetMods.Select(p => p.SpecificBioPolymer.Parent.Accession).ToList();
                var betaProteins = csm.BetaPeptide.BestMatchingBioPolymersWithSetMods.Select(p => p.SpecificBioPolymer.Parent.Accession).ToList();

                foreach (var alpha in alphaProteins)
                {
                    foreach (var beta in betaProteins)
                    {
                        if (alpha == beta || alpha == "DECOY_" + beta || beta == "DECOY_" + alpha)
                        {
                            return true;
                        }
                    }
                }
            }

            return false;
        }

        public static void ResolveProteinPosAmbiguitiesForXl(CrosslinkSpectralMatch csm)
        {
            csm.ResolveAllAmbiguities();
            if (csm.BetaPeptide != null)
            {
                csm.BetaPeptide.ResolveAllAmbiguities();
            }

            //Assign PsmCrossType.Cross to Intra or Inter.
            if (csm.CrossType == PsmCrossType.Cross)
            {
                if (IsIntraCsm(csm))
                {
                    csm.CrossType = PsmCrossType.Intra;
                }
                else
                {
                    csm.CrossType = PsmCrossType.Inter;
                }
            }

            if (csm.CrossType == PsmCrossType.Cross || csm.CrossType == PsmCrossType.Intra || csm.CrossType == PsmCrossType.Inter)
            {
                // alpha peptide crosslink residue in the protein
                csm.XlProteinPos = csm.OneBasedStartResidue == null ? (int?)null : csm.OneBasedStartResidue.Value + csm.LinkPositions[0] - 1;

                // beta crosslink residue in protein
                csm.BetaPeptide.XlProteinPos = csm.BetaPeptide.OneBasedStartResidue == null ? (int?)null : csm.BetaPeptide.OneBasedStartResidue.Value + csm.BetaPeptide.LinkPositions[0] - 1;
            }
            else if (csm.CrossType == PsmCrossType.DeadEnd || csm.CrossType == PsmCrossType.DeadEndH2O || csm.CrossType == PsmCrossType.DeadEndNH2 || csm.CrossType == PsmCrossType.DeadEndTris)
            {
                csm.XlProteinPos = csm.OneBasedStartResidue == null ? (int?)null : csm.OneBasedStartResidue.Value + csm.LinkPositions[0] - 1;
            }
            else if (csm.CrossType == PsmCrossType.Loop)
            {
                csm.XlProteinPos = csm.OneBasedStartResidue == null ? (int?)null : csm.OneBasedStartResidue.Value + csm.LinkPositions[0] - 1;

                csm.XlProteinPosLoop = csm.OneBasedStartResidue == null ? (int?)null : csm.OneBasedStartResidue.Value + csm.LinkPositions[1] - 1;
            }
        }

        public static List<int> GetPossibleCrosslinkerModSites(char[] crosslinkerModSites, PeptideWithSetModifications peptide, InitiatorMethionineBehavior initiatorMethionineBehavior, bool CrosslinkAtCleavageSite)
        {
            List<int> possibleXlPositions = null;

            bool wildcard = crosslinkerModSites.Any(p => p == 'X');

            var range = Enumerable.Range(0, peptide.BaseSequence.Length);
            if (!CrosslinkAtCleavageSite && peptide.OneBasedEndResidue != peptide.Protein.Length 
                && !peptide.Protein.TruncationProducts.Any(x => x.OneBasedEndPosition == peptide.OneBasedEndResidue))
            {
                //The C termial cannot be crosslinked and cleaved.
                range = Enumerable.Range(0, peptide.BaseSequence.Length - 1);
            }
            foreach (var r in range)
            {
                if (crosslinkerModSites.Contains(peptide.BaseSequence[r]) || wildcard)
                {
                    //Try to eliminate those site with mod on it. Consider the possibility that the site is at Protein N terminal.       
                    if (!peptide.AllModsOneIsNterminus.Keys.Contains(r + 2))
                    {
                        if (possibleXlPositions == null)
                        {
                            possibleXlPositions = new List<int>();
                        }
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
            sb.Append(SpectrumMatchFromTsvHeader.FileName + '\t');
            sb.Append(SpectrumMatchFromTsvHeader.UniqueSequence + '\t');
            sb.Append(SpectrumMatchFromTsvHeader.Ms2ScanNumber + '\t');
            sb.Append(SpectrumMatchFromTsvHeader.PrecursorScanNum + '\t');
            sb.Append(SpectrumMatchFromTsvHeader.PrecursorMz + '\t');
            sb.Append(SpectrumMatchFromTsvHeader.PrecursorCharge + '\t');
            sb.Append(SpectrumMatchFromTsvHeader.PrecursorMass + '\t');
            sb.Append(SpectrumMatchFromTsvHeader.CrossTypeLabel + '\t');
            sb.Append(SpectrumMatchFromTsvHeader.LinkResiduesLabel + "\t");

            sb.Append("Peptide Info -->" + '\t');
            sb.Append(SpectrumMatchFromTsvHeader.ProteinAccession + '\t');
            sb.Append(SpectrumMatchFromTsvHeader.ProteinLinkSiteLabel + '\t');
            sb.Append(SpectrumMatchFromTsvHeader.BaseSequence + '\t');
            sb.Append(SpectrumMatchFromTsvHeader.FullSequence + '\t');
            sb.Append(SpectrumMatchFromTsvHeader.PeptideMonoMass + '\t');
            sb.Append(SpectrumMatchFromTsvHeader.Score + '\t');

            sb.Append(SpectrumMatchFromTsvHeader.MatchedIonSeries + '\t');
            sb.Append(SpectrumMatchFromTsvHeader.MatchedIonMzRatios + '\t');
            sb.Append(SpectrumMatchFromTsvHeader.MatchedIonMassDiffDa + '\t');
            sb.Append(SpectrumMatchFromTsvHeader.MatchedIonMassDiffPpm + '\t');
            sb.Append(SpectrumMatchFromTsvHeader.MatchedIonIntensities + '\t');
            sb.Append(SpectrumMatchFromTsvHeader.MatchedIonCounts + '\t');

            sb.Append("Beta Peptide Info -->" + '\t');
            sb.Append(SpectrumMatchFromTsvHeader.BetaPeptideProteinAccessionLabel + '\t');
            sb.Append(SpectrumMatchFromTsvHeader.BetaPeptideProteinLinkSiteLabel + '\t');
            sb.Append(SpectrumMatchFromTsvHeader.BetaPeptideBaseSequenceLabel + '\t');
            sb.Append(SpectrumMatchFromTsvHeader.BetaPeptideFullSequenceLabel + '\t');
            sb.Append(SpectrumMatchFromTsvHeader.BetaPeptideTheoreticalMassLabel + '\t');
            sb.Append(SpectrumMatchFromTsvHeader.BetaPeptideScoreLabel + '\t');

            sb.Append("Beta Peptide Matched Ions" + '\t');
            sb.Append(SpectrumMatchFromTsvHeader.BetaPeptideMatchedIonsLabel + '\t');
            sb.Append("Beta Peptide Matched Ion Mass Diff (Da)" + '\t');
            sb.Append("Beta Peptide Matched Ion Mass Diff (Ppm)" + '\t');
            sb.Append("Beta Peptide Matched Ion Intensities" + '\t');
            sb.Append("Beta Peptide Matched Ion Counts" + '\t');

            sb.Append("Summary Info -->" + '\t');
            sb.Append(SpectrumMatchFromTsvHeader.XLTotalScoreLabel + '\t');
            sb.Append(SpectrumMatchFromTsvHeader.MassDiffDa + '\t');
            sb.Append("AlphaIndexingRank" + '\t');
            sb.Append(SpectrumMatchFromTsvHeader.ParentIonsLabel + '\t');
            sb.Append("ParentIonsNum" + '\t');
            sb.Append("AlphaParentIonMaxIntensityRank" + '\t');
            sb.Append("BetaParentIonMaxIntensityRank" + '\t');
            sb.Append(SpectrumMatchFromTsvHeader.DecoyContaminantTarget + '\t');
            sb.Append(SpectrumMatchFromTsvHeader.QValue + '\t');
            sb.Append(SpectrumMatchFromTsvHeader.PEP + '\t');
            sb.Append(SpectrumMatchFromTsvHeader.PEP_QValue + '\t');

            return sb.ToString();
        }

        public static string GetTabSepHeaderSingle()
        {
            var sb = new StringBuilder();
            sb.Append("File Name" + '\t');
            sb.Append(SpectrumMatchFromTsvHeader.UniqueSequence + '\t');
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

            sb.Append("Matched Ion Series" + '\t');
            sb.Append("Matched Ion Mass-To-Charge Ratios" + '\t');
            sb.Append("Matched Ion Mass Diff (Da)" + '\t');
            sb.Append("Matched Ion Mass Diff (Ppm)" + '\t');
            sb.Append("Matched Ion Intensities" + '\t');
            sb.Append("Matched Ion Counts" + '\t');
            sb.Append("Decoy/Contaminant/Target" + '\t');
            sb.Append("QValue" + '\t');
            sb.Append(SpectrumMatchFromTsvHeader.PEP + '\t');
            sb.Append(SpectrumMatchFromTsvHeader.PEP_QValue + '\t');

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
            sb.Append(UniqueSequence + "\t");
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

            sb.Append("\t"); //Intentionally left empty for readability in the tsv file.
            List<PeptideWithSetModifications> pepsWithMods = BestMatchingBioPolymersWithSetMods.Select(p => p.SpecificBioPolymer as PeptideWithSetModifications).ToList();
            var proteinAccessionString = Accession ?? PsmTsvWriter.Resolve(pepsWithMods.Select(b => b.Protein.Accession), FullSequence).ResolvedString;
            sb.Append(proteinAccessionString + "\t");
            sb.Append(XlProteinPos + (XlProteinPosLoop.HasValue ? "~" + XlProteinPosLoop.Value : null) + "\t");
            sb.Append(BaseSequence + "\t");
            sb.Append(FullSequence + position + "\t");
            sb.Append((BioPolymerWithSetModsMonoisotopicMass.HasValue ? BioPolymerWithSetModsMonoisotopicMass.Value.ToString() : "---"));
            sb.Append("\t");
            sb.Append(Score + "\t");

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
                int i = 0;
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
                sb.Append("\t"); //Intentionally left empty for readability in the tsv file.
                List<PeptideWithSetModifications> betaPepsWithMods = BetaPeptide.BestMatchingBioPolymersWithSetMods.Select(p => p.SpecificBioPolymer as PeptideWithSetModifications).ToList();
                var betaProteinAccessionString = BetaPeptide.Accession ?? PsmTsvWriter.Resolve(betaPepsWithMods.Select(b => b.Protein.Accession), FullSequence).ResolvedString;
                sb.Append(betaProteinAccessionString + "\t");
                sb.Append(BetaPeptide.XlProteinPos + "\t");
                sb.Append(BetaPeptide.BaseSequence + "\t");
                sb.Append(BetaPeptide.FullSequence + "(" + BetaPeptide.LinkPositions[0].ToString() + ")" + "\t");
                sb.Append(BetaPeptide.BioPolymerWithSetModsMonoisotopicMass.ToString() + "\t");
                sb.Append(BetaPeptide.Score + "\t");

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
                sb.Append(((BioPolymerWithSetModsMonoisotopicMass.HasValue) ? (ScanPrecursorMass - BetaPeptide.BioPolymerWithSetModsMonoisotopicMass - BioPolymerWithSetModsMonoisotopicMass.Value).ToString() : "---")); sb.Append("\t");

                sb.Append(XlRank.ToString() + "\t");

                int alphaNumParentIons = MatchedFragmentIons.Count(p => p.NeutralTheoreticalProduct.ProductType == ProductType.M);
                int betaNumParentIons = BetaPeptide.MatchedFragmentIons.Count(p => p.NeutralTheoreticalProduct.ProductType == ProductType.M);

                sb.Append(alphaNumParentIons + ";" + betaNumParentIons + "\t");
                sb.Append(alphaNumParentIons + betaNumParentIons + "\t");
                sb.Append(((ParentIonMaxIntensityRanks != null) && (ParentIonMaxIntensityRanks.Any()) ? ParentIonMaxIntensityRanks.Min().ToString() : "-")); sb.Append("\t");
                sb.Append(((BetaPeptide.ParentIonMaxIntensityRanks != null) && (BetaPeptide.ParentIonMaxIntensityRanks.Any()) ? BetaPeptide.ParentIonMaxIntensityRanks.Min().ToString() : "-")); sb.Append("\t");
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
            if (FdrInfo == null)
            {
                sb.Append("0\t0\t0\t");
            }
            else
            {
                sb.Append(FdrInfo.QValue.ToString());
                sb.Append("\t");

                sb.Append(FdrInfo.PEP.ToString()); sb.Append("\t");

                sb.Append(FdrInfo.PEP_QValue.ToString()); sb.Append("\t");
            }
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