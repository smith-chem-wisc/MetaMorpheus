using Chemistry;
using MassSpectrometry;
using MzLibUtil;
using Proteomics.AminoAcidPolymer;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace EngineLayer.CrosslinkSearch
{
    public class PsmCross : PeptideSpectralMatch
    {
        public CompactPeptide compactPeptide;

        private static readonly double waterMonoisotopicMass = PeriodicTable.GetElement("H").PrincipalIsotope.AtomicMass * 2 + PeriodicTable.GetElement("O").PrincipalIsotope.AtomicMass;
        private static readonly double nitrogenAtomMonoisotopicMass = PeriodicTable.GetElement("N").PrincipalIsotope.AtomicMass;
        private static readonly double oxygenAtomMonoisotopicMass = PeriodicTable.GetElement("O").PrincipalIsotope.AtomicMass;
        private static readonly double hydrogenAtomMonoisotopicMass = PeriodicTable.GetElement("H").PrincipalIsotope.AtomicMass;

        public PsmCross(CompactPeptide theBestPeptide, int notch, double score, int scanIndex, Ms2ScanWithSpecificMass scan, DigestionParams digestionParams) : base(theBestPeptide, notch, score, scanIndex, scan, digestionParams)
        {
            compactPeptide = theBestPeptide;
        }

        public double BestScore { get; set; } //For the current psmCross
        public List<int> ModPositions { get; set; }
        public PsmCross BetaPsmCross { get; set; }
        public double DScore { get; set; }
        public List<MatchedFragmentIon> MatchedIons { get; set; }

        public double XLTotalScore { get; set; } //alpha + beta psmCross
        public double XLQvalueTotalScore { get; set; } //Calc based on XLtotalScore for Qvalue
        public int XlProteinPos { get; set; }
        public List<int> XlRank { get; set; } //only contain 2 intger, consider change to Tuple
        public string ParentIonExist { get; set; }
        public int ParentIonExistNum { get; set; }
        public List<int> ParentIonMaxIntensityRanks { get; set; }
        public PsmCrossType CrossType { get; set; }

        public List<TheoreticalFragmentIon> GetTheoreticalFragmentIons(List<ProductType> productTypes)
        {
            List<TheoreticalFragmentIon> theoreticalFragmentIons = new List<TheoreticalFragmentIon>();
            foreach (var pt in productTypes)
            {
                int ionNumberAdd = 1;
                if (pt == ProductType.BnoB1ions)
                {
                    // first generated b ion is b2, not b1, if we're skipping b1 ions
                    ionNumberAdd++;
                }
                List<ProductType> temp = new List<ProductType> { pt };
                var productMasses = compactPeptide.ProductMassesMightHaveDuplicatesAndNaNs(temp);

                for (int i = 0; i < productMasses.Length; i++)
                {
                    theoreticalFragmentIons.Add(new TheoreticalFragmentIon(productMasses[i], double.NaN, 1, pt, i + ionNumberAdd));
                }
            }


            return theoreticalFragmentIons;
        }

        public Dictionary<List<int>, List<TheoreticalFragmentIon>> XlGetTheoreticalFramentIons(List<ProductType> productTypes, bool Charge_2_3, CrosslinkerTypeClass crosslinker,List<int> modPos, double modMass)
        {
            Dictionary<List<int>, List<TheoreticalFragmentIon>> AllTheoreticalFragmentIonsLists = new Dictionary<List<int>, List<TheoreticalFragmentIon>>();

            List<TheoreticalFragmentIon> baseTheoreticalFragmentIons = GetTheoreticalFragmentIons(productTypes);

            foreach (var iPos in modPos)
            {
                List<TheoreticalFragmentIon> currentIons = new List<TheoreticalFragmentIon>();

                if (crosslinker.Cleavable)
                {
                    currentIons.Add(new TheoreticalFragmentIon(compactPeptide.MonoisotopicMassIncludingFixedMods + crosslinker.CleaveMassShort, double.NaN, 1, ProductType.None, 0));
                    currentIons.Add(new TheoreticalFragmentIon(compactPeptide.MonoisotopicMassIncludingFixedMods + crosslinker.CleaveMassLong, double.NaN, 1, ProductType.None, 0));
                }

                foreach (var iIon in baseTheoreticalFragmentIons)
                {
                    var iType = iIon.ProductType;
                    switch (iType)
                    {
                        case ProductType.BnoB1ions:
                            if (iIon.IonNumber < iPos)
                            { currentIons.Add(iIon); }
                            else
                            {
                                if (crosslinker.Cleavable)
                                {
                                    currentIons.Add(new TheoreticalFragmentIon(iIon.Mass + crosslinker.CleaveMassShort, double.NaN, 1, iIon.ProductType, iIon.IonNumber));
                                    currentIons.Add(new TheoreticalFragmentIon(iIon.Mass + crosslinker.CleaveMassLong, double.NaN, 1, iIon.ProductType, iIon.IonNumber));
                                }
                                else
                                {
                                    currentIons.Add(new TheoreticalFragmentIon(iIon.Mass + modMass + crosslinker.TotalMass, double.NaN, 1, iIon.ProductType, iIon.IonNumber));
                                }
                            }
                            break;
                        case ProductType.C:
                            if (iIon.IonNumber < iPos)
                            {
                                currentIons.Add(iIon);
                            }
                            else { currentIons.Add(new TheoreticalFragmentIon(iIon.Mass + modMass + crosslinker.TotalMass, double.NaN, 1, iIon.ProductType, iIon.IonNumber)); }
                            break;
                        case ProductType.Y:
                            if (iIon.IonNumber < compactPeptide.CTerminalMasses.Length - iPos + 2)
                            { currentIons.Add(iIon); }
                            else
                            {
                                if (crosslinker.Cleavable)
                                {
                                    currentIons.Add(new TheoreticalFragmentIon(iIon.Mass + crosslinker.CleaveMassShort, double.NaN, 1, iIon.ProductType, iIon.IonNumber));
                                    currentIons.Add(new TheoreticalFragmentIon(iIon.Mass + crosslinker.CleaveMassLong, double.NaN, 1, iIon.ProductType, iIon.IonNumber));
                                }
                                else
                                { currentIons.Add(new TheoreticalFragmentIon(iIon.Mass + modMass + crosslinker.TotalMass, double.NaN, 1, iIon.ProductType, iIon.IonNumber)); }
                            }
                            break;
                        case ProductType.Zdot:
                            if (iIon.IonNumber < compactPeptide.CTerminalMasses.Length - iPos + 2)
                            { currentIons.Add(iIon); }
                            else { currentIons.Add(new TheoreticalFragmentIon(iIon.Mass + modMass + crosslinker.TotalMass, double.NaN, 1, iIon.ProductType, iIon.IonNumber)); }
                            break;
                        default:
                            break;
                    }
                }

                if (Charge_2_3)
                {
                    var length = currentIons.Count;
                    for (int i = 0; i < length; i++)
                    {
                        currentIons.Add(new TheoreticalFragmentIon(currentIons[i].Mass, double.NaN, 2, currentIons[i].ProductType, currentIons[i].IonNumber));
                        currentIons.Add(new TheoreticalFragmentIon(currentIons[i].Mass, double.NaN, 3, currentIons[i].ProductType, currentIons[i].IonNumber));
                    }
                }

                AllTheoreticalFragmentIonsLists.Add(new List<int> { iPos }, currentIons);
            }

            return AllTheoreticalFragmentIonsLists;
        }

        //TO DO: the second ModPostion jPos is not recorded. 
        public Dictionary<List<int>, List<TheoreticalFragmentIon>> XlLoopGetTheoreticalFramentIons(List<ProductType> productTypes, bool Charge_2_3, CrosslinkerTypeClass crosslinker, List<int> modPos, double modMass)
        {
            Dictionary<List<int>, List<TheoreticalFragmentIon>> AllTheoreticalFragmentIonsLists = new Dictionary<List<int>, List<TheoreticalFragmentIon>>();

            List<TheoreticalFragmentIon> baseTheoreticalFragmentIons = GetTheoreticalFragmentIons(productTypes);

            if (modPos.Count() >= 2)
            {
                for (int iPos = 0; iPos < modPos.Count() - 1; iPos++)
                {
                    for (int jPos = iPos + 1; jPos < modPos.Count(); jPos++)
                    {
                        List<TheoreticalFragmentIon> currentIons = new List<TheoreticalFragmentIon>();

                        foreach (var iIon in baseTheoreticalFragmentIons)
                        {
                            var iType = iIon.ProductType;
                            switch (iType)
                            {
                                case ProductType.BnoB1ions:
                                    if (iIon.IonNumber < modPos[iPos])
                                    {
                                        currentIons.Add(iIon);
                                    }
                                    else if(iIon.IonNumber >= modPos[jPos])
                                    {
                                        currentIons.Add(new TheoreticalFragmentIon(iIon.Mass + modMass + crosslinker.LoopMass, double.NaN, 1, iIon.ProductType, iIon.IonNumber));
                                    }
                                    break;
                                case ProductType.C:
                                    if (iIon.IonNumber < modPos[iPos])
                                    {
                                        currentIons.Add(iIon);
                                    }
                                    else if(iIon.IonNumber >= modPos[jPos]) { currentIons.Add(new TheoreticalFragmentIon(iIon.Mass + modMass + crosslinker.LoopMass, double.NaN, 1, iIon.ProductType, iIon.IonNumber)); }
                                    break;
                                case ProductType.Y:
                                    if (iIon.IonNumber < compactPeptide.CTerminalMasses.Length -modPos[jPos] + 2)
                                    {
                                        currentIons.Add(iIon);
                                    }
                                    else if(iIon.IonNumber >= compactPeptide.CTerminalMasses.Length - modPos[iPos] + 2)
                                    {
                                        currentIons.Add(new TheoreticalFragmentIon(iIon.Mass + modMass + crosslinker.LoopMass, double.NaN, 1, iIon.ProductType, iIon.IonNumber));          
                                    }
                                    break;
                                case ProductType.Zdot:
                                    if (iIon.IonNumber < compactPeptide.CTerminalMasses.Length - modPos[jPos] + 2)
                                    {
                                        currentIons.Add(iIon);
                                    }
                                    else if (iIon.IonNumber >= compactPeptide.CTerminalMasses.Length - modPos[iPos] + 2) { currentIons.Add(new TheoreticalFragmentIon(iIon.Mass + modMass + crosslinker.LoopMass, double.NaN, 1, iIon.ProductType, iIon.IonNumber)); }
                                    break;
                                default:
                                    break;
                            }
                        }

                        if (Charge_2_3)
                        {
                            var length = currentIons.Count;
                            for (int i = 0; i < length; i++)
                            {
                                currentIons.Add(new TheoreticalFragmentIon(currentIons[i].Mass, double.NaN, 2, currentIons[i].ProductType, currentIons[i].IonNumber));
                                currentIons.Add(new TheoreticalFragmentIon(currentIons[i].Mass, double.NaN, 3, currentIons[i].ProductType, currentIons[i].IonNumber));
                            }
                        }
                        AllTheoreticalFragmentIonsLists.Add(new List<int> { modPos[iPos], modPos[jPos] }, currentIons);
                    }
                }
            }

            return AllTheoreticalFragmentIonsLists;
        }

        public static int[] GenerateIntensityRanks(double[] experimental_intensities)
        {
            var y = experimental_intensities.ToArray();
            var x = Enumerable.Range(1, y.Length).OrderBy(p => p).ToArray();
            Array.Sort(y, x);
            var experimental_intensities_rank = Enumerable.Range(1, y.Length).OrderByDescending(p => p).ToArray();
            Array.Sort(x, experimental_intensities_rank);
            return experimental_intensities_rank;
        }

        public List<int> XlPosCal(string crosslinkerModSites)
        {
            Tolerance tolerance = new PpmTolerance(1);
            List<int> xlpos = new List<int>();
            foreach (char item in crosslinkerModSites)
            {
                if (tolerance.Within(compactPeptide.NTerminalMasses[0], Residue.GetResidue(item).MonoisotopicMass))
                {
                    xlpos.Add(1);
                }
                for (int i = 1; i < compactPeptide.NTerminalMasses.Length; i++)
                {
                    if (tolerance.Within(compactPeptide.NTerminalMasses[i] - compactPeptide.NTerminalMasses[i - 1], Residue.GetResidue(item).MonoisotopicMass))
                    {
                        xlpos.Add(i+1);
                    }
                }
                if (tolerance.Within(compactPeptide.CTerminalMasses[0], Residue.GetResidue(item).MonoisotopicMass))
                {
                    xlpos.Add(compactPeptide.NTerminalMasses.Length);
                }
            }
            xlpos.Sort();
            return xlpos;
        }

        public void GetBestMatch(Ms2ScanWithSpecificMass theScan, Dictionary<List<int>, List<TheoreticalFragmentIon>> pmmhList, CommonParameters commonParameters)
        {
            BestScore = 0;
            foreach (var pmmh in pmmhList)
            {
                var matchedIons = MetaMorpheusEngine.MatchFragmentIons( theScan.TheScan.MassSpectrum, pmmh.Value, commonParameters);
                var score = MetaMorpheusEngine.CalculatePeptideScore(theScan.TheScan, matchedIons, 0);
                if (score > BestScore)
                {
                    BestScore = score;
                    MatchedIons = matchedIons;
                    ModPositions = pmmh.Key;          
                }
            }

            if (MatchedIons != null)
            {
                double[] experimental_intensities = theScan.TheScan.MassSpectrum.YArray;
                int[] experimental_intensities_rank = GenerateIntensityRanks(experimental_intensities);
                foreach (var tIon in MatchedIons)
                {
                    // get the closest peak in the spectrum to the theoretical peak
                    int matchedPeakIndex = theScan.TheScan.MassSpectrum.GetClosestPeakIndex(tIon.Mz).Value;
                    tIon.IntensityRank = experimental_intensities_rank[matchedPeakIndex];
                }
            }
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
                case PsmCrossType.Singe:
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
            sb.Append(CompactPeptides.First().Value.Item2.Select(p => p.Protein.Accession).First().ToString()); sb.Append("\t");
            sb.Append(XlProteinPos); sb.Append("\t");
            sb.Append(BaseSequence); sb.Append("\t");
            sb.Append(FullSequence + position); sb.Append("\t");
            sb.Append((PeptideMonisotopicMass.HasValue ? PeptideMonisotopicMass.Value.ToString() : "---")); sb.Append("\t");
            sb.Append(BestScore); sb.Append("\t");
            sb.Append(XlRank[0]); sb.Append("\t");

            if (BetaPsmCross!= null)
            {             
                sb.Append(""); sb.Append("\t");
                sb.Append(BetaPsmCross.CompactPeptides.First().Value.Item2.Select(p => p.Protein.Accession).First().ToString()); sb.Append("\t");
                sb.Append(BetaPsmCross.XlProteinPos); sb.Append("\t");
                sb.Append(BetaPsmCross.BaseSequence); sb.Append("\t");
                sb.Append(BetaPsmCross.FullSequence + "(" + ModPositions[0].ToString() + ")"); sb.Append("\t");
                sb.Append((BetaPsmCross.PeptideMonisotopicMass.HasValue ? BetaPsmCross.PeptideMonisotopicMass.Value.ToString() : "---")); sb.Append("\t");
                sb.Append(BetaPsmCross.BestScore); sb.Append("\t");
                sb.Append(XlRank[1]); sb.Append("\t");

                sb.Append(""); sb.Append("\t");
                sb.Append(XLQvalueTotalScore); sb.Append("\t");
                sb.Append(((PeptideMonisotopicMass.HasValue && BetaPsmCross.PeptideMonisotopicMass.HasValue) ? (BetaPsmCross.ScanPrecursorMass - BetaPsmCross.PeptideMonisotopicMass.Value - PeptideMonisotopicMass.Value).ToString() : "---")); sb.Append("\t");
                
                sb.Append(ParentIonExist + "." + BetaPsmCross.ParentIonExist); sb.Append("\t");
                sb.Append(ParentIonExistNum); sb.Append("\t");
                sb.Append(((ParentIonMaxIntensityRanks != null) && (ParentIonMaxIntensityRanks.Any()) ? ParentIonMaxIntensityRanks.Min().ToString() : "-")); sb.Append("\t");               
            }

            sb.Append(MatchedIons.Count(p=>p.TheoreticalFragmentIon.Charge > 1)); sb.Append("\t");
            if (BetaPsmCross == null) { sb.Append((IsDecoy) ? "-1" : "1"); sb.Append("\t"); }
            else { sb.Append((IsDecoy || BetaPsmCross.IsDecoy) ? "-1" : "1"); sb.Append("\t"); }
            
            sb.Append((FdrInfo != null) ? FdrInfo.QValue.ToString() : "-"); sb.Append("\t");

            return sb.ToString();
        }
    }
}