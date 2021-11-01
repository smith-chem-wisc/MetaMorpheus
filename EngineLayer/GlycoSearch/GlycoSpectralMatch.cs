﻿using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Proteomics;

namespace EngineLayer.GlycoSearch
{
    public class GlycoSpectralMatch : PeptideSpectralMatch
    {
        public GlycoSpectralMatch(PeptideWithSetModifications theBestPeptide, int notch, double score, int scanIndex, Ms2ScanWithSpecificMass scan, CommonParameters commonParameters, List<MatchedFragmentIon> matchedFragmentIons)
            : base(theBestPeptide, notch, score, scanIndex, scan, commonParameters, matchedFragmentIons)
        {

        }

        #region Proterties

        public int Rank { get; set; }
        public Dictionary<int, List<MatchedFragmentIon>> ChildMatchedFragmentIons { get; set; }

        public double PredictedHydrophobicity { get; set; }

        public double PredictedRT { get; set; }

        //Glyco properties
        public List<int> ModPos { get; set; }
        public LocalizationLevel LocalizationLevel { get; set; }

        public double MScore { get; set; }  //The original Morpheus score.
        public double PeptideScore { get; set; } //Scores from only mathced peptide fragments.
        public double GlycanScore { get; set; } //Scores from only matched Y ions. 
        public double DiagnosticIonScore { get; set; } //Since every glycopeptide generate DiagnosticIon, it is important to seperate the score. 
        public static GlycanBox GetFirstGraphGlycanBox(GlycoSpectralMatch gsm)
        {

            if (gsm.GlycanType == GlycoType.OGlycoPep)
            {
                return GlycanBox.OGlycanBoxes[gsm.LocalizationGraphs.First().ModBoxId];
            }
            else if (gsm.GlycanType == GlycoType.NGlycoPep)
            {
                return GlycanBox.NGlycanBoxes[gsm.LocalizationGraphs.First().ModBoxId];
            }
            else
            {
                return GlycanBox.MixedModBoxes[gsm.LocalizationGraphs.First().ModBoxId];
            }

        }

        public static Glycan[] GetFirstGraphGlycans(GlycoSpectralMatch gsm, GlycanBox glycanBox)
        {
            var glycans = new Glycan[glycanBox.ModCount];
            if (gsm.GlycanType == GlycoType.OGlycoPep)
            {
                for (int i = 0; i < glycanBox.ModCount; i++)
                {
                    glycans[i] = GlycanBox.GlobalOGlycans[glycanBox.ModIds[i]];
                }
                return glycans;
            }
            else if (gsm.GlycanType == GlycoType.NGlycoPep)
            {
                for (int i = 0; i < glycanBox.ModCount; i++)
                {
                    glycans[i] = GlycanBox.GlobalNGlycans[glycanBox.ModIds[i]];
                }
                return glycans;
            }
            else
            {
                for (int i = 0; i < glycanBox.ModCount; i++)
                {
                    glycans[i] = GlycanBox.GlobalMixedGlycans[glycanBox.ModIds[i]];
                }
                return glycans;
            }
        }

        //Glycan type indicator
        public GlycoType GlycanType { get; set; } 

        public bool NGlycanMotifExist { get; set; } //NGlycan Motif exist. 

        public double R138to144 { get; set; } //The intensity ratio of this 138 and 144 could be a signature for O-glycan or N-glycan.

        public double[] OxoniumIonIntensity { get; set; } //Please Check Glycan.AllOxoniumIons

        public bool PepNN { get; set; } //Yion Pep + 2HexNAc, supposed to be exist more in N-Glycopeptide

        public bool PepNH { get; set; } //Yion Pep + Hex + HexNAc, supposed to be only exist in O-Glycopeptide

        public int LongestconcatenatedYion { get; set; } //Currently Simplified with GlycanScore.

        //Glyco Info
        public List<LocalizationGraph> LocalizationGraphs { get; set; }  //Graph-based Localization information.
        public List<Route> Routes { get; set; } //Localized modification sites and modfication ID.

        public double ScanInfo_p { get; set; }  //Scan P value, Used for Localization probability calculation. Ref PhosphoRS paper.

        public int Thero_n { get; set; } //Scan n value. Used for Localization probability calculation. Ref PhosphoRS paper.

        public Dictionary<int, List<Tuple<int, double>>> SiteSpeciLocalProb { get; set; } // Data <modPos, List<glycanId, site probability>>

        public List<GlycoSite> LocalizedGlycan { get; set; } 

        #endregion

        public static string GetTabSepHeaderGlyco(bool IsGlycopepitde)
        {
            var sb = new StringBuilder();
            sb.Append("File Name" + '\t');
            sb.Append("Scan Number" + '\t');
            sb.Append("Scan Retention Time" + '\t');
            sb.Append("Precursor Scan Number" + '\t');
            sb.Append("Precursor MZ" + '\t');
            sb.Append("Precursor Charge" + '\t');
            sb.Append("Precursor Mass" + '\t');

            sb.Append("Protein Accession" + '\t');
            sb.Append("Organism" + '\t');
            sb.Append("Protein Name" + '\t');
            sb.Append("Start and End Residues In Protein" + '\t');
            sb.Append("Base Sequence" + '\t');
            sb.Append("FlankingResidues" + '\t');
            sb.Append("Full Sequence" + '\t');
            sb.Append("Number of Mods" + '\t');
            sb.Append("Peptide Monoisotopic Mass" + '\t');
            //sb.Append("Predicted Hydrophobicity" + '\t');
            sb.Append("Score" + '\t');
            sb.Append("MScore" + '\t');
            sb.Append("Rank" + '\t');

            sb.Append("Matched Ion Series" + '\t');
            sb.Append("Matched Ion Mass-To-Charge Ratios" + '\t');
            sb.Append("Matched Ion Mass Diff (Da)" + '\t');
            sb.Append("Matched Ion Mass Diff (Ppm)" + '\t');
            sb.Append("Matched Ion Intensities" + '\t');
            sb.Append("Matched Ion Counts" + '\t');

            sb.Append("Decoy/Contaminant/Target" + '\t');
            sb.Append("QValue" + '\t');
            sb.Append("PEP" + '\t');
            sb.Append("PEP_QValue" + '\t');

            if (IsGlycopepitde)
            {
                sb.Append("Localization Score" + '\t');
                sb.Append("Yion Score" + '\t');
                sb.Append("DiagonosticIon Score" + '\t');
                sb.Append("Plausible Number Of Glycans" + '\t');
                sb.Append("Total Glycosylation sites" + '\t');
                sb.Append("GlycanMass" + '\t');
                sb.Append("Plausible GlycanComposition" + '\t');               
                sb.Append("Plausible GlycanStructure" + '\t');
                sb.Append("GlycanLocalizationLevel" + '\t');
                sb.Append("Localized Glycans with Peptide Site Specific Probability" + '\t');
                sb.Append("Localized Glycans with Protein Site Specific Probability" + '\t');


                sb.Append("All potential glycan localizations" + '\t');
                sb.Append("AllSiteSpecificLocalizationProbability" + '\t');


                //Glycan Type indicator
                sb.Append("GlycanType" + '\t');
                sb.Append("N-Glycan motif Check" + '\t');
                sb.Append("R138/144" + '\t');
                sb.Append("I_168" + '\t');
                sb.Append("I_186" + '\t');
                sb.Append("I_366" + '\t');
                sb.Append("I_274" + '\t');
                sb.Append("I_292" + '\t');
                sb.Append("PepNN" + '\t');
                sb.Append("PepNH" + '\t');
                sb.Append("NumOfContateYion" + '\t');
            }
           
            return sb.ToString();
        }

        public override string ToString()
        {
            var sb = new StringBuilder();
            sb.Append(FullFilePath + "\t");
            sb.Append(ScanNumber + "\t");
            sb.Append(ScanRetentionTime + "\t");
            sb.Append(PrecursorScanNumber + "\t");
            sb.Append(ScanPrecursorMonoisotopicPeakMz + "\t");
            sb.Append(ScanPrecursorCharge + "\t");
            sb.Append(ScanPrecursorMass + "\t");

            var proteinAccessionString = ProteinAccession ?? PsmTsvWriter.Resolve(BestMatchingPeptides.Select(p => p.Peptide.Protein.Accession), FullSequence).ResolvedString;
            sb.Append(proteinAccessionString + "\t");
            sb.Append(Organism + "\t");
            sb.Append(PsmTsvWriter.Resolve(BestMatchingPeptides.Select(b => b.Peptide.Protein.FullName), FullSequence).ResolvedString + "\t");
            int _FirstOneBasedStartResidueInProtein = OneBasedStartResidueInProtein.HasValue ? OneBasedStartResidueInProtein.Value : BestMatchingPeptides.First().Peptide.OneBasedStartResidueInProtein;
            int _FirstOneBasedEndResidueInProtein = OneBasedEndResidueInProtein.HasValue ? OneBasedEndResidueInProtein.Value : BestMatchingPeptides.First().Peptide.OneBasedEndResidueInProtein; ;

            if (OneBasedStartResidueInProtein.HasValue)
            {
                sb.Append("[" + OneBasedStartResidueInProtein.Value.ToString() + " to " +  OneBasedEndResidueInProtein.Value.ToString() + "]" + '\t');
            }
            else
            {
                sb.Append("\t");
            }
            
            sb.Append(BaseSequence + "\t");
            sb.Append(BestMatchingPeptides.First().Peptide.PreviousAminoAcid + "," + BestMatchingPeptides.First().Peptide.NextAminoAcid + "\t");
            sb.Append(FullSequence + "\t");
            sb.Append(BestMatchingPeptides.First().Peptide.AllModsOneIsNterminus.Count + "\t");

            sb.Append((PeptideMonisotopicMass.HasValue ? PeptideMonisotopicMass.Value.ToString() : "---")); sb.Append("\t");
            //sb.Append(PredictedHydrophobicity); sb.Append("\t");
            //sb.Append(PredictedRT); sb.Append("\t");
            sb.Append(Score + "\t");
            sb.Append(MScore + "\t");
            sb.Append(Rank + "\t");

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

            sb.Append((IsDecoy) ? "D" : (IsContaminant) ? "C" : "T"); sb.Append("\t");


            sb.Append(FdrInfo!=null? FdrInfo.QValue.ToString() : "-1" );  sb.Append("\t");

            sb.Append("0" + "\t"); //PEP

            sb.Append("0" + "\t");  //PEP_QValue

            if (LocalizationGraphs != null)
            {
                sb.Append(LocalizationGraphs.First().TotalScore + "\t");

                sb.Append(GlycanScore + "\t");

                sb.Append(DiagnosticIonScore + "\t");

                var glycanBox = GetFirstGraphGlycanBox(this);

                sb.Append(glycanBox.ModCount + "\t");

                sb.Append(LocalizationGraphs.First().ModPos.Length + "\t");

                sb.Append(glycanBox.Mass + "\t");

                sb.Append(Glycan.GetKindString(glycanBox.Kind)); sb.Append("\t");

                //Get glycans
                var glycans = GetFirstGraphGlycans(this, glycanBox);

                //if (glycans.First().Struc!=null)
                //{
                //    sb.Append(string.Join(",", glycans.Select(p => p.Struc.ToString()).ToArray())); 
                //}
                sb.Append(string.Join(",", glycans.Select(p => p.Composition).ToArray()));
                sb.Append("\t");

                if (Routes!=null)
                {
                    sb.Append(LocalizationLevel); sb.Append("\t");

                    string local_peptide = "";
                    string local_protein = "";
                    if (GlycanType == GlycoType.OGlycoPep)
                    {
                        LocalizedSiteSpeciLocalInfo(SiteSpeciLocalProb, LocalizedGlycan, OneBasedStartResidueInProtein, GlycanBox.GlobalOGlycans, ref local_peptide, ref local_protein);
                    }
                    else if (GlycanType == GlycoType.NGlycoPep)
                    {
                        LocalizedSiteSpeciLocalInfo(SiteSpeciLocalProb, LocalizedGlycan, OneBasedStartResidueInProtein, GlycanBox.GlobalNGlycans, ref local_peptide, ref local_protein);
                    }
                    else
                    {
                        LocalizedSiteSpeciLocalInfo(SiteSpeciLocalProb, LocalizedGlycan, OneBasedStartResidueInProtein, GlycanBox.GlobalMixedGlycans, ref local_peptide, ref local_protein);
                    }

                    sb.Append(local_peptide); sb.Append("\t");
                    sb.Append(local_protein); sb.Append("\t");

                    sb.Append(AllLocalizationInfo(Routes)); sb.Append("\t");

                    sb.Append(SiteSpeciLocalInfo(SiteSpeciLocalProb)); sb.Append("\t");
                }
                else
                {
                    sb.Append(LocalizationLevel); sb.Append("\t");
                    sb.Append("\t");
                    sb.Append("\t");
                    sb.Append("\t");
                    sb.Append("\t");
                }
            }

            //Output for Glycan type indicator
            if (LocalizationGraphs != null)
            {
                sb.Append(GlycanType); sb.Append("\t");
                var NSiteExist = GlycoPeptides.MotifExist(BaseSequence, new string[] { "Nxt", "Nxs" });
                sb.Append(NSiteExist); sb.Append("\t");

                var R138vs144 = 1.0;
                if (OxoniumIonIntensity[5] <= 0.00000001)
                {
                    R138vs144 = 100000000;
                }
                else
                {
                    R138vs144 = OxoniumIonIntensity[4] / OxoniumIonIntensity[5];
                }

                sb.Append(R138vs144.ToString()); sb.Append("\t");
                //Intensity of oxonium ion 168. Please Check Glycan.AllOxoniumIons
                sb.Append(OxoniumIonIntensity[7].ToString()); sb.Append("\t");
                //Intensity of oxonium ion 186
                sb.Append(OxoniumIonIntensity[8].ToString()); sb.Append("\t");
                //Intensity of oxonium ion 366
                sb.Append(OxoniumIonIntensity[14].ToString()); sb.Append("\t");
                //Intensity of oxonium ion 274
                sb.Append(OxoniumIonIntensity[10].ToString()); sb.Append("\t");
                //Intensity of oxonium ion 292
                sb.Append(OxoniumIonIntensity[12].ToString()); sb.Append("\t");
                sb.Append(PepNN); sb.Append("\t");
                sb.Append(PepNH); sb.Append("\t");
                //This here should be LongestconcatenatedYion
                sb.Append(GlycanScore.ToString()); sb.Append("\t");
            }

            return sb.ToString();
        }

        public static Dictionary<string, string> MatchedIonDataDictionary(List<MatchedFragmentIon> matchedFragmentIons)
        {
            Dictionary<string, string> s = new Dictionary<string, string>();
            PsmTsvWriter.AddMatchedIonsData(s, matchedFragmentIons);
            return s;
        }

        #region Glycopeptide Localization Output
        public static string AllLocalizationInfo(List<Route> routes)
        {
            string local = "";

            if (routes == null || routes.Count == 0)
            {
                return local;
            }
            //Some GSP have a lot paths, in which case only output first 10 paths and the total number of the paths.
            int maxOutputPath = 10;
            if (routes.Count <= maxOutputPath)
            {
                maxOutputPath = routes.Count;
            }

            int i = 0;
            while (i < maxOutputPath)
            {
                var ogl = routes[i];
                local += "{@" + ogl.ModBoxId.ToString() + "[";
                var g = string.Join(",", ogl.Mods.Select(p => (p.ModSite - 1).ToString() + "-" + p.GlycanID.ToString()));
                local += g + "]}";
                i++;
            }

            if (routes.Count > maxOutputPath)
            {
                local += "... In Total:" + routes.Count.ToString() + " Paths";
            }

            return local;
        }

        public static void LocalizedSiteSpeciLocalInfo(Dictionary<int, List<Tuple<int, double>>> siteSpeciLocalProb, List<GlycoSite> localizedGlycan, int? OneBasedStartResidueInProtein, Glycan[] globalGlycans, ref string local, ref string local_protein)
        {
            if (siteSpeciLocalProb == null)
            {
                return;
            }

            foreach (var loc in localizedGlycan.Where(p => p.IsLocalized))
            {
                var x = siteSpeciLocalProb[loc.ModSite].Where(p => p.Item1 == loc.GlycanID).First().Item2;
                var peptide_site = loc.ModSite - 1;
                local += "[" + peptide_site + "," + globalGlycans[loc.GlycanID].Composition + "," + x.ToString("0.000") + "]";

                var protein_site = OneBasedStartResidueInProtein.HasValue ? OneBasedStartResidueInProtein.Value + loc.ModSite - 2 : -1;
                local_protein += "[" + protein_site + "," + globalGlycans[loc.GlycanID].Composition + "," + x.ToString("0.000") + "]";
            }

        }
        
        public static string SiteSpeciLocalInfo(Dictionary<int, List<Tuple<int, double>>> siteSpeciLocalProb)
        {
            string local = "";

            if (siteSpeciLocalProb == null)
            {
                return local;
            }

            foreach (var sitep in siteSpeciLocalProb)
            {
                var site_1 = sitep.Key - 1;
                local += "{@" + site_1;
                foreach (var s in sitep.Value)
                {
                    local += "[" + s.Item1 + "," + s.Item2.ToString("0.000") + "]";
                }
                local += "}";
            }

            return local;
        }

        #endregion
    }
}