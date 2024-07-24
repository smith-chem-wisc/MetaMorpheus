using Omics.Fragmentation;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using Omics.Modifications;
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

    public class GlycoSpectralMatch : SpectralMatch
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
        /// <summary>
        /// Try to get the ModSite in the right format.
        /// </summary>
        /// <param name="peptide"> full peptide sequence ex. "PTLFKNVSLYK" </param>
        /// <param name="motifs"> modificatino AA ex. "S","T"</param>
        /// <returns> int[], the Modpositon index list ex.[9,3] </returns>
        public static List<int> GetPossibleModSites(PeptideWithSetModifications peptide, string[] motifs)
        {
            List<int> possibleModSites = new List<int>();

            List<Modification> modifications = new List<Modification>();

            foreach (var mtf in motifs)
            {
                if (ModificationMotif.TryGetMotif(mtf, out ModificationMotif aMotif)) //Check if the motif is valid, and creat the motif object from the string.
                {
                    Modification modWithMotif = new Modification(_target: aMotif, _locationRestriction: "Anywhere."); 
                    modifications.Add(modWithMotif);
                }
            }

            foreach (var modWithMotif in modifications) //interate through all the modifications with motif.
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

        /// <summary>
        /// Generate the peptide header, ex File name, Precursor m/z, Score…
        /// </summary>
        /// <returns></returns>
        public static string GetTabSepHeaderSingle() //Most complicate part in this class
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
            sb.Append("PEP_QValue");

            return sb.ToString();
        }

        /// <summary>
        /// Generate the glyco header ex Localization Score, Yion Score…
        /// </summary>
        /// <returns></returns>
        public static string GetTabSeperatedHeaderGlyco()
        {
            var sb = new StringBuilder();
            sb.Append("\t");//provides the separation needed from GetTabSepHeaderSingle()
            sb.Append("Localization Score" + '\t');
            sb.Append("Yion Score" + '\t');
            sb.Append("DiagonosticIon Score" + '\t');
            sb.Append("Plausible Number Of Glycans" + '\t');//Not used for N-Glyco
            sb.Append("Total Glycosylation sites" + '\t');//Not used for N-Glyco
            sb.Append("GlycanMass" + '\t');
            sb.Append("Plausible GlycanComposition" + '\t');
            sb.Append("N-Glycan motif Check" + '\t');//Not used for N-Glyco
            sb.Append("R138/144" + '\t');
            sb.Append("Plausible GlycanStructure" + '\t');
            sb.Append("GlycanLocalizationLevel" + '\t');
            sb.Append("Localized Glycans with Peptide Site Specific Probability" + '\t');
            sb.Append("Localized Glycans with Protein Site Specific Probability" + '\t');
            sb.Append("All potential glycan localizations" + '\t');//Not used for N-Glyco
            sb.Append("AllSiteSpecificLocalizationProbability");//Not used for N-Glyco

            return sb.ToString();
        }

        /// <summary>
        /// Put the psm data into the corresponding columns.
        /// </summary>
        /// <returns></returns>
        public string SingleToString()
        {
            var sb = new StringBuilder();
            sb.Append(FullFilePath + "\t");
            sb.Append(ScanNumber + "\t");
            sb.Append(ScanRetentionTime + "\t");
            sb.Append(PrecursorScanNumber + "\t");
            sb.Append(ScanPrecursorMonoisotopicPeakMz + "\t");
            sb.Append(ScanPrecursorCharge + "\t");
            sb.Append(ScanPrecursorMass + "\t");

            var proteinAccessionString = Accession ?? PsmTsvWriter.Resolve(BestMatchingBioPolymersWithSetMods.Select(p => p.Peptide.Parent.Accession), FullSequence).ResolvedString;
            sb.Append(proteinAccessionString + "\t");
            sb.Append(Organism + "\t");
            sb.Append(PsmTsvWriter.Resolve(BestMatchingBioPolymersWithSetMods.Select(b => b.Peptide.Parent.FullName), FullSequence).ResolvedString + "\t"); //protein name
            int _FirstOneBasedStartResidueInProtein = OneBasedStartResidue.HasValue ? OneBasedStartResidue.Value : BestMatchingBioPolymersWithSetMods.First().Peptide.OneBasedStartResidue;
            int _FirstOneBasedEndResidueInProtein = OneBasedEndResidue.HasValue ? OneBasedEndResidue.Value : BestMatchingBioPolymersWithSetMods.First().Peptide.OneBasedEndResidue; ;

            if (OneBasedStartResidue.HasValue)
            {
                sb.Append("[" + OneBasedStartResidue.Value.ToString() + " to " +  OneBasedEndResidue.Value.ToString() + "]" + '\t');
            }
            else
            {
                sb.Append("\t");
            }
            
            sb.Append(BaseSequence + "\t");
            sb.Append(BestMatchingBioPolymersWithSetMods.First().Peptide.PreviousResidue + "," + BestMatchingBioPolymersWithSetMods.First().Peptide.NextResidue + "\t");
            sb.Append(FullSequence + "\t");
            sb.Append(BestMatchingBioPolymersWithSetMods.First().Peptide.AllModsOneIsNterminus.Count + "\t");

            sb.Append((BioPolymerWithSetModsMonoisotopicMass.HasValue ? BioPolymerWithSetModsMonoisotopicMass.Value.ToString() : "---")); sb.Append("\t");
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

            sb.Append("0" + "\t"); //This is space for PEP

            sb.Append("0"); //This is space for PEP Q-value

            return sb.ToString();
        }

        /// <summary>
        /// Put the glycan data into the corresponding columns.
        /// </summary>
        /// <returns></returns>
        public string GlycoToString()
        {
            var sb = new StringBuilder();

            if (Routes != null)//this gets the o-glyco
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
                } //Convert the glycanBox index into the real glycan object. ex. [H1N1, H2N2A1, H2N2A1F1]

                if (glycans.First().Struc != null)
                {
                    sb.Append(string.Join(",", glycans.Select(p => p.Struc.ToString()).ToArray())); //ex. (N(H)),(N(H(A))(N(H))),(N(H)(N(H(A))(F))
                }
                sb.Append("\t");

                sb.Append(CorrectLocalizationLevel(SiteSpeciLocalProb, LocalizationGraphs.First(), Routes.First(), LocalizedGlycan, LocalizationLevel)); sb.Append("\t");

                string local_peptide = "";
                string local_protein = "";
                LocalizedSiteSpeciLocalInfo(SiteSpeciLocalProb, LocalizedGlycan, OneBasedStartResidue, ref local_peptide, ref local_protein);
                sb.Append(local_peptide); sb.Append("\t");
                sb.Append(local_protein); sb.Append("\t");

                sb.Append(AllLocalizationInfo(Routes)); sb.Append("\t");

                sb.Append(SiteSpeciLocalInfo(SiteSpeciLocalProb));
            }
            else if (GlycanScore > 0)//this gets the N-glcyo that remain
            {
                sb.Append("\t"); //Localization score

                sb.Append(GlycanScore + "\t");

                sb.Append(DiagnosticIonScore + "\t");

                sb.Append("\t"); //number of mods

                sb.Append( "\t"); //mod pos length

                sb.Append((double)NGlycan.First().Mass / 1E5 + "\t");

                sb.Append(Glycan.GetKindString(NGlycan.First().Kind) + "\t");

                var NSiteExist = MotifExist(BaseSequence, new string[] { "Nxt", "Nxs" });

                sb.Append(NSiteExist); sb.Append("\t");

                sb.Append(R138vs144.ToString()); sb.Append("\t");

                if (NGlycan.First().Struc != null)
                {
                    sb.Append(NGlycan.First().Struc);
                }
                sb.Append("\t");

                sb.Append("\t");

                sb.Append("\t");
                sb.Append("\t");

                 sb.Append("\t");

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


        /// <summary>
        /// Two function included: 
        /// (1) Analysis all pair, and evaluate any site is occured in all cases, if yes set a true on that. If not, set a false.
        /// (2) Classify the localization level base on the localization.
        /// </summary>
        /// <param name="OGlycanBoxLocalization"> all case of the pair </param>
        /// <param name="localizationLevel"> level 1 to level 3 </param>
        /// <returns> A tuple, represent the pair and its confidience ex. [3,5,ture] means glycan 5 located on glycosite 3, and very confidience </returns>
        public static List<Tuple<int, int, bool>> GetLocalizedGlycan(List<Route> OGlycanBoxLocalization, out LocalizationLevel localizationLevel)
        {
            List<Tuple<int, int, bool>> localizedGlycan = new List<Tuple<int, int, bool>>();

            Dictionary<string, int> modSiteSeenCount = new Dictionary<string, int>(); // all possible glycan-sites pair, Dictionary<string, int>: site-glycan pair, count

            foreach (var ogl in OGlycanBoxLocalization) // ogl means one case, there are three glycan located on the same peptide: (5,1,False),(9,8,Flase),(10,9,Ture)
            {
                foreach (var og in ogl.Mods)            // og means one glycan locaization, like (5,1,False) -> glycan 1 attached on postion5.
                {
                    var k = og.Item1.ToString() + "-" + og.Item2.ToString(); // k = 5-1(glycosite-glycan) means the glycan-site pair
                    if (modSiteSeenCount.ContainsKey(k)) // accout the number of the same glycan-site pair
                    {
                        modSiteSeenCount[k] += 1;   // this pair cpunt +1
                    }
                    else
                    {
                        modSiteSeenCount.Add(k, 1); // If the pair is first time to seen, add it to the dictionary.
                    }
                }
            }

            localizationLevel = LocalizationLevel.Level3;
            if (OGlycanBoxLocalization.Count == 1) // we just have one situation(route), no other possibility
            {
                localizationLevel = LocalizationLevel.Level1;
            }
            else if (OGlycanBoxLocalization.Count > 1)
            {
                if (modSiteSeenCount.Values.Where(p => p == OGlycanBoxLocalization.Count).Count() > 0) //If anyone of the glycan-site pair is localized in all the cases, then the localization level is 2.
                {
                    localizationLevel = LocalizationLevel.Level2;
                }
                else
                {
                    localizationLevel = LocalizationLevel.Level3;
                }
            }

            foreach (var seenMod in modSiteSeenCount)
            {
                if (seenMod.Value == OGlycanBoxLocalization.Count) // Try to fine the glycan-site pair that always localized in all the cases.
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

        /// <summary>
        /// convert the Route information into the string format.
        /// </summary>
        /// <param name="OGlycanBoxLocalization"> Route collection ex. [(9,4),(8,4),(7,4)...], ModBoxId = 7 </param>
        /// <returns> string {@7[8-4]}{@7[7-4]}{@7[6-4]} means three case, glycan 4 located on glycosite 6,  glycan 4 located on glycosite 7, glycan 4 located on glycosite 8 </returns>
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
                var g = string.Join(",", ogl.Mods.Select(p => (p.Item1 - 1).ToString() + "-" + p.Item2.ToString())); //why we have to -1 here?
                local += g + "]}";
                i++;
            }

            if (OGlycanBoxLocalization.Count > maxOutputPath)
            {
                local += "... In Total:" + OGlycanBoxLocalization.Count.ToString() + " Paths";
            }

            return local;
        }

        /// <summary>
        /// Just for the case at Level1 and Level1b.
        /// </summary>
        /// <param name="siteSpeciLocalProb"></param>
        /// <param name="localizationGraph"></param>
        /// <param name="route"></param>
        /// <param name="localizedGlycan"></param>
        /// <param name="localizationLevel"></param>
        /// <returns> level 1 or level 1b</returns>
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

                if (!route.Mods[i].Item3) // if the peak is not exist.
                {
                    return LocalizationLevel.Level1b;
                }
            }


            return localizationLevel;

        }
        /// <summary>
        /// Output the special localization information. String store in Local_peptide and Local_protein. ex. [9,H2N2A1F1,0.589] means glycan H2N2A1F1 located on glycosite 9 with 0.589 probability.
        /// </summary>
        /// <param name="siteSpeciLocalProb"> site : (glycan, probility)[] ex. site2 : [(glycan1, 5%), (glycan2, 5%), (glycan3, 90%)] </param>
        /// <param name="localizedGlycan"> [(6,4,false),(7,4,false),(7,2,true)], glycosite,glycan,confidience respectively </param>
        /// <param name="OneBasedStartResidueInProtein"></param>
        /// <param name="local"></param>
        /// <param name="local_protein"></param>
        public static void LocalizedSiteSpeciLocalInfo(Dictionary<int, List<Tuple<int, double>>> siteSpeciLocalProb, List<Tuple<int, int, bool>> localizedGlycan, int? OneBasedStartResidueInProtein, ref string local_peptide, ref string local_protein)
        {
            if (siteSpeciLocalProb == null)
            {
                return;
            }

            foreach (var glycositePair in localizedGlycan.Where(p => p.Item3)) // get the most confidient glycosite-glycan pair, loc is a pair of glycosite and glycan. Item 1 is glycosite, Item 2 is glycanId.
            {
                var site_glycanProb = siteSpeciLocalProb[glycositePair.Item1].Where(p => p.Item1 == glycositePair.Item2).First().Item2; // get the probability of the specfic glycan on the specific site.
                var peptide_site = glycositePair.Item1 - 1;
                local_peptide += "[" + peptide_site + "," + GlycanBox.GlobalOGlycans[glycositePair.Item2].Composition + "," + site_glycanProb.ToString("0.000") + "]";

                var protein_site = OneBasedStartResidueInProtein.HasValue ? OneBasedStartResidueInProtein.Value + glycositePair.Item1 - 2 : -1;
                local_protein += "[" + protein_site + "," + GlycanBox.GlobalOGlycans[glycositePair.Item2].Composition + "," + site_glycanProb.ToString("0.000") + "]";
            }

        }

        /// <summary>
        /// Generate the site specific localization information.
        /// </summary>
        /// <param name="siteSpeciLocalProb"></param>
        /// <returns> Site specific localization information. ex. {1[1,0.2][2,0.8]} means glycan 1 and 2 are located on glycosite 1 and 2 with 20% and 80% probability. </returns>
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