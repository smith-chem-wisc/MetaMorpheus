using Proteomics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Proteomics;

namespace EngineLayer.GlycoSearch
{
    //Localization of multiple glycans on one peptides can be divide into the following groups based on the quanlity of the localization. Similar to Proteomform Level.
    public enum LocalizationLevel
    {
        Level1,
        Level1b,
        Level2,
        Level3
    }

    public class GlycoSpectralMatch : PeptideSpectralMatch
    {
        public GlycoSpectralMatch(PeptideWithSetModifications theBestPeptide, int notch, double score, int scanIndex, Ms2ScanWithSpecificMass scan, CommonParameters commonParameters, List<MatchedFragmentIon> matchedFragmentIons)
            : base(theBestPeptide, notch, score, scanIndex, scan, commonParameters, matchedFragmentIons)
        {

        }

        public int Rank { get; set; }
        public Dictionary<int, List<MatchedFragmentIon>> ChildMatchedFragmentIons { get; set; }
        //Glyco properties
        public List<Glycan> NGlycan { get; set; }  //Identified NGlycan
        public List<int> NGlycanLocalizations { get; set; }


        public List<LocalizationGraph> LocalizationGraphs { get; set; }  //Graph-based Localization information.
        public List<Route> Routes { get; set; } //Localized modification sites and modfication ID.

        public double ScanInfo_p { get; set; }  //Scan P value, Used for Localization probability calculation. Ref PhosphoRS paper.

        public int Thero_n { get; set; } //Scan n value. Used for Localization probability calculation. Ref PhosphoRS paper.

        public Dictionary<int, List<Tuple<int, double>>> SiteSpeciLocalProb { get; set; } // Data <modPos, List<glycanId, site probability>>
        public double PeptideScore { get; set; } //Scores from only mathced peptide fragments.
        public double GlycanScore { get; set; } //Scores from only matched Y ions. 
        public double DiagnosticIonScore { get; set; } //Since every glycopeptide generate DiagnosticIon, it is important to seperate the score. 

        public double R138vs144 { get; set; } // The intensity ratio of this 138 and 144 could be a signature for O-glycan or N-glycan.
        public List<Tuple<int, int, bool>> LocalizedGlycan { get; set; } //<mod site, glycanID, isLocalized> All seen glycans identified.
        public LocalizationLevel LocalizationLevel { get; set; }  

        //Motif should be writen with required form
        public static List<int> GetPossibleModSites(PeptideWithSetModifications peptide, string[] motifs)
        {
            List<int> possibleModSites = new List<int>();

            List<Modification> modifications = new List<Modification>();

            foreach (var mtf in motifs)
            {
                if (ModificationMotif.TryGetMotif(mtf, out ModificationMotif aMotif))
                {
                    Modification modWithMotif = new Modification(_target: aMotif, _locationRestriction: "Anywhere.");
                    modifications.Add(modWithMotif);
                }
            }

            foreach (var modWithMotif in modifications)
            {
                for (int r = 0; r < peptide.Length; r++)
                {
                    if (peptide.AllModsOneIsNterminus.Keys.Contains(r+2))
                    {
                        continue;
                    }
                    
                    //FullSequence is used here to avoid duplicated modification on same sites?
                    if (ModificationLocalization.ModFits(modWithMotif, peptide.BaseSequence, r + 1, peptide.Length, r + 1))
                    {
                        possibleModSites.Add(r + 2);
                    }
                }
            }

            return possibleModSites;
        }

        public static bool MotifExist(string baseSeq, string[] motifs)
        {
            List<Modification> modifications = new List<Modification>();

            foreach (var mtf in motifs)
            {
                if (ModificationMotif.TryGetMotif(mtf, out ModificationMotif aMotif))
                {
                    Modification modWithMotif = new Modification(_target: aMotif, _locationRestriction: "Anywhere.");
                    modifications.Add(modWithMotif);
                }
            }

            foreach (var modWithMotif in modifications)
            {
                for (int r = 0; r < baseSeq.Length; r++)
                {
                    //Modification is not considered.                  
                    if (ModificationLocalization.ModFits(modWithMotif, baseSeq, r + 1, baseSeq.Length, r + 1))
                    {
                        return true;
                    }
                }
            }

            return false;
        }

        public static string GetTabSepHeaderSingle()
        {
            var sb = new StringBuilder();
            sb.Append("File Name" + '\t');
            sb.Append("Scan Number" + '\t');
            sb.Append("Retention Time" + '\t');
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
            sb.Append("PEP" + '\t');
            sb.Append("PEP_QValue" + '\t');

            return sb.ToString();
        }

        public static string GetTabSepHeaderOGlyco()
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
            sb.Append("PEP" + '\t');
            sb.Append("PEP_QValue" + '\t');

            sb.Append("Localization Score" + '\t');
            sb.Append("Yion Score" + '\t');
            sb.Append("DiagonosticIon Score" + '\t');
            sb.Append("Plausible Number Of Glycans" + '\t');
            sb.Append("Total Glycosylation sites" + '\t');
            sb.Append("GlycanMass" + '\t');
            sb.Append("Plausible GlycanComposition" + '\t');
            sb.Append("N-Glycan motif Check" + '\t');
            sb.Append("R138/144" + '\t');
            sb.Append("Plausible GlycanStructure" + '\t');
            sb.Append("GlycanLocalizationLevel" + '\t');
            sb.Append("Localized Glycans with Peptide Site Specific Probability" + '\t');
            sb.Append("Localized Glycans with Protein Site Specific Probability" + '\t');
            sb.Append("All potential glycan localizations" + '\t');
            sb.Append("AllSiteSpecificLocalizationProbability" + '\t');

            return sb.ToString();
        }

        public static string GetTabSepHeaderNGlyco()
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
            sb.Append("PEP" + '\t');
            sb.Append("PEP_QValue" + '\t');

            sb.Append("Localization Score" + '\t');
            sb.Append("Yion Score" + '\t');
            sb.Append("DiagonosticIon Score" + '\t');
            sb.Append("GlycanMass" + '\t');
            sb.Append("Plausible GlycanComposition" + '\t');
            sb.Append("R138/144" + '\t');
            sb.Append("Plausible GlycanStructure" + '\t');
            sb.Append("GlycanLocalizationLevel" + '\t');
            sb.Append("Localized Glycans with Peptide Site Specific Probability" + '\t');
            sb.Append("Localized Glycans with Protein Site Specific Probability" + '\t');

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
            sb.Append(Score + "\t");
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

            sb.Append("0" + "\t");

            sb.Append("0" + "\t");

            if (NGlycan != null)
            {            
                sb.Append(PeptideScore + "\t");
                sb.Append(GlycanScore + "\t");
                sb.Append(DiagnosticIonScore + "\t");
                sb.Append((double)NGlycan.First().Mass / 1E5); sb.Append("\t");
                sb.Append(Glycan.GetKindString(NGlycan.First().Kind)); sb.Append("\t");
                sb.Append(R138vs144.ToString()); sb.Append("\t");
                if (NGlycan.First().Struc!=null)
                {
                    sb.Append(NGlycan.First().Struc); sb.Append("\t");
                }
            }

            if (Routes != null)
            {
                sb.Append(LocalizationGraphs.First().TotalScore + "\t");

                sb.Append(GlycanScore + "\t");

                sb.Append(DiagnosticIonScore + "\t");              

                var glycanBox = GlycanBox.OGlycanBoxes[Routes.First().ModBoxId];

                sb.Append(glycanBox.NumberOfMods + "\t");

                sb.Append(LocalizationGraphs.First().ModPos.Length + "\t");

                sb.Append(glycanBox.Mass + "\t");

                sb.Append(Glycan.GetKindString(glycanBox.Kind)); sb.Append("\t");

                var NSiteExist = MotifExist(BaseSequence, new string[] { "Nxt", "Nxs" });

                sb.Append(NSiteExist); sb.Append("\t");

                sb.Append(R138vs144.ToString()); sb.Append("\t");

                //Get glycans
                var glycans = new Glycan[glycanBox.NumberOfMods];
                for (int i = 0; i < glycanBox.NumberOfMods; i++)
                {
                    glycans[i] = GlycanBox.GlobalOGlycans[glycanBox.ModIds[i]];
                }

                if (glycans.First().Struc!=null)
                {
                    sb.Append(string.Join(",", glycans.Select(p => p.Struc.ToString()).ToArray())); 
                }
                sb.Append("\t");

                sb.Append(CorrectLocalizationLevel(SiteSpeciLocalProb, LocalizationGraphs.First(), Routes.First(), LocalizedGlycan,  LocalizationLevel)) ; sb.Append("\t");

                //string localizedGlycan = LocalizedGlycan.Where(p=>p.Item3).Count() > 0 ? "[" + string.Join(",", LocalizedGlycan.Where(p => p.Item3).Select(p => p.Item1.ToString() + "-" + p.Item2.ToString())) + "]" : "";
                //sb.Append(localizedGlycan); sb.Append("\t");
                string local_peptide = "";
                string local_protein = "";
                LocalizedSiteSpeciLocalInfo(SiteSpeciLocalProb, LocalizedGlycan, OneBasedStartResidueInProtein, ref local_peptide, ref local_protein);
                sb.Append(local_peptide); sb.Append("\t");
                sb.Append(local_protein); sb.Append("\t");

                sb.Append(AllLocalizationInfo(Routes)); sb.Append("\t");

                sb.Append(SiteSpeciLocalInfo(SiteSpeciLocalProb));
            }

            return sb.ToString();
        }

        public static Dictionary<string, string> MatchedIonDataDictionary(List<MatchedFragmentIon> matchedFragmentIons)
        {
            Dictionary<string, string> s = new Dictionary<string, string>();
            PsmTsvWriter.AddMatchedIonsData(s, matchedFragmentIons);
            return s;
        }

        //Output: <int, int, string> <ModBoxId, ModPosition, is localized>; Input: List<Route>
        public static List<Tuple<int, int, bool>> GetLocalizedGlycan(List<Route> OGlycanBoxLocalization, out LocalizationLevel localizationLevel)
        {
            List<Tuple<int, int, bool>> localizedGlycan = new List<Tuple<int, int, bool>>();

            //Dictionary<string, int>: modsite-id, count
            Dictionary<string, int> seenModSite = new Dictionary<string, int>();

            foreach (var ogl in OGlycanBoxLocalization)
            {
                foreach (var og in ogl.Mods)
                {
                    var k = og.Item1.ToString() + "-" + og.Item2.ToString();
                    if (seenModSite.ContainsKey(k))
                    {
                        seenModSite[k] += 1;
                    }
                    else
                    {
                        seenModSite.Add(k, 1);
                    }
                }
            }

            localizationLevel = LocalizationLevel.Level3;
            if (OGlycanBoxLocalization.Count == 1)
            {
                localizationLevel = LocalizationLevel.Level1;
            }
            else if (OGlycanBoxLocalization.Count > 1)
            {
                if (seenModSite.Values.Where(p => p == OGlycanBoxLocalization.Count).Count() > 0)
                {
                    localizationLevel = LocalizationLevel.Level2;
                }
                else
                {
                    localizationLevel = LocalizationLevel.Level3;
                }
            }

            foreach (var seenMod in seenModSite)
            {
                if (seenMod.Value == OGlycanBoxLocalization.Count)
                {
                    localizedGlycan.Add(new Tuple<int, int, bool>(int.Parse(seenMod.Key.Split('-')[0]), int.Parse(seenMod.Key.Split('-')[1]), true));
                }
                else
                {
                    localizedGlycan.Add(new Tuple<int, int, bool>(int.Parse(seenMod.Key.Split('-')[0]), int.Parse(seenMod.Key.Split('-')[1]), false));
                }
            }

            return localizedGlycan;
        }

        public static string AllLocalizationInfo(List<Route> OGlycanBoxLocalization)
        {
            string local = "";

            if (OGlycanBoxLocalization == null || OGlycanBoxLocalization.Count == 0)
            {
                return local;
            }
            //Some GSP have a lot paths, in which case only output first 10 paths and the total number of the paths.
            int maxOutputPath = 10;
            if (OGlycanBoxLocalization.Count <= maxOutputPath)
            {
                maxOutputPath = OGlycanBoxLocalization.Count;
            }

            int i = 0;
            while (i < maxOutputPath)
            {
                var ogl = OGlycanBoxLocalization[i];
                local += "{@" + ogl.ModBoxId.ToString() + "[";
                var g = string.Join(",", ogl.Mods.Select(p => (p.Item1 - 1).ToString() + "-" + p.Item2.ToString()));
                local += g + "]}";
                i++;
            }

            if (OGlycanBoxLocalization.Count > maxOutputPath)
            {
                local += "... In Total:" + OGlycanBoxLocalization.Count.ToString() + " Paths";
            }

            return local;
        }

        //Correct Localization Level based on site specific probability. If LocalizationLevel = 1, and there are site probability lower than 0.75, Correct the level to 1b.
        public static LocalizationLevel CorrectLocalizationLevel(Dictionary<int, List<Tuple<int, double>>> siteSpeciLocalProb, LocalizationGraph localizationGraph, Route route, List<Tuple<int, int, bool>> localizedGlycan, LocalizationLevel localizationLevel)
        {
            if (siteSpeciLocalProb == null || localizationLevel!=LocalizationLevel.Level1)
            {
                return localizationLevel;
            }

            if (localizationGraph.ModPos.Length == 1 && localizationGraph.TotalScore == 0)
            {
                return LocalizationLevel.Level1b;
            }


            for (int i = 0; i < localizedGlycan.Count; i++)
            {
                var g = localizedGlycan[i];
                if (siteSpeciLocalProb[g.Item1].Where(p => p.Item1 == g.Item2).First().Item2 < 0.75)
                {
                    return LocalizationLevel.Level1b;
                }

                if (!route.Mods[i].Item3)
                {
                    return LocalizationLevel.Level1b;
                }
            }


            return localizationLevel;

        }
        public static void LocalizedSiteSpeciLocalInfo(Dictionary<int, List<Tuple<int, double>>> siteSpeciLocalProb, List<Tuple<int, int, bool>> localizedGlycan, int? OneBasedStartResidueInProtein, ref string local, ref string local_protein)
        {
            if (siteSpeciLocalProb == null)
            {
                return;
            }

            foreach (var loc in localizedGlycan.Where(p => p.Item3))
            {
                var x = siteSpeciLocalProb[loc.Item1].Where(p => p.Item1 == loc.Item2).First().Item2;
                var peptide_site = loc.Item1 - 1;
                local += "[" + peptide_site + "," + GlycanBox.GlobalOGlycans[loc.Item2].Composition + "," + x.ToString("0.000") + "]";

                var protein_site = OneBasedStartResidueInProtein.HasValue ? OneBasedStartResidueInProtein.Value + loc.Item1 - 2 : -1;
                local_protein += "[" + protein_site + "," + GlycanBox.GlobalOGlycans[loc.Item2].Composition + "," + x.ToString("0.000") + "]";
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
    }
}