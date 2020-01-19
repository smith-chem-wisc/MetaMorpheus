using Proteomics.Fragmentation;
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
        public GlycoSpectralMatch(PeptideWithSetModifications theBestPeptide, int notch, double score, int scanIndex, Ms2ScanWithSpecificMass scan, DigestionParams digestionParams, List<MatchedFragmentIon> matchedFragmentIons)
            : base(theBestPeptide, notch, score, scanIndex, scan, digestionParams, matchedFragmentIons)
        {
            this.TotalScore = score;
        }

        public double TotalScore { get; set; } //peptide + glycan psmCross
        public int Rank { get; set; } 
        public Dictionary<int, List<MatchedFragmentIon>> ChildMatchedFragmentIons { get; set; }
        //Glyco properties
        public List<Glycan> NGlycan { get; set; }
        public List<int> NGlycanLocalizations { get; set; }

        public List<Tuple<int, Tuple<int, int>[]>> OGlycanBoxLocalization;

        public double PeptideScore { get; set; }
        public double GlycanScore { get; set; }
        public double DiagnosticIonScore { get; set; }
        public double R138vs144 { get; set; }

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
            sb.Append("Protein Name" + '\t');
            sb.Append("Start and End Residues In Protein" + '\t');
            sb.Append("Base Sequence" + '\t');
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
            sb.Append("Child Scans Matched Ion Series" + '\t');
            sb.Append("Decoy/Contaminant/Target" + '\t');
            sb.Append("QValue" + '\t');

            return sb.ToString();
        }

        public static string GetTabSepHeaderGlyco()
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
            sb.Append("Protein Name" + '\t');
            sb.Append("Start and End Residues In Protein" + '\t');
            sb.Append("Base Sequence" + '\t');
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

            sb.Append("Total Score" + '\t');
            //sb.Append("Peptide Score" + '\t');
            //sb.Append("Glycan Score" + '\t');
            //sb.Append("DiagonosticIon Score" + '\t');
            sb.Append("GlycanMass" + '\t');
            //sb.Append("GlycanDecoy" + '\t');   
            sb.Append("Plausible GlycanComposition" + '\t');
            sb.Append("R138/144" + '\t');
            sb.Append("Plausible GlycanStructure" + '\t');                       
            sb.Append("GlycanLocalization" + '\t');
            sb.Append("GlycanLocalizationLevel" + '\t');
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

            sb.Append(ProteinAccession + "\t");
            sb.Append(BestMatchingPeptides.First().Peptide.Protein.FullName + "\t");
            sb.Append("[" + OneBasedStartResidueInProtein.Value.ToString() + " to " + OneBasedEndResidueInProtein.Value.ToString() + "]" + '\t');

            sb.Append(BaseSequence + "\t");
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

            sb.Append((IsDecoy) ? "D" : (IsContaminant) ? "C" : "T");
            sb.Append("\t");


            sb.Append(FdrInfo.QValue.ToString() + "\t");
            sb.Append("0" + "\t");
            sb.Append("0" + "\t");
            if (NGlycan != null)
            {
                sb.Append(TotalScore + "\t");             
                sb.Append(PeptideScore + "\t");
                sb.Append(GlycanScore + "\t");
                sb.Append(DiagnosticIonScore + "\t");
                sb.Append(string.Join("|", NGlycan.Select(p => p.GlyId.ToString()).ToArray())); sb.Append("\t");
                sb.Append(NGlycan.First().Decoy? "D": "T"); sb.Append("\t");
                sb.Append(NGlycan.First().Struc); sb.Append("\t");
                sb.Append((double)NGlycan.First().Mass/1E5); sb.Append("\t");
                sb.Append(string.Join(" ", NGlycan.First().Kind.Select(p => p.ToString()).ToArray())); sb.Append("\t");
            }

            if (OGlycanBoxLocalization != null)
            {
                sb.Append(TotalScore + "\t");

                var glycanBox = GlycanBox.OGlycanBoxes[OGlycanBoxLocalization.First().Item1];

                sb.Append(glycanBox.Mass); sb.Append("\t");

                //sb.Append( "T" + '\t');  

                sb.Append(Glycan.GetKindString(glycanBox.Kind)); sb.Append("\t");

                sb.Append(R138vs144.ToString()); sb.Append("\t");

                //Get glycans
                var glycans = new Glycan[glycanBox.NumberOfMods];
                for (int i = 0; i < glycanBox.NumberOfMods; i++)
                {
                    glycans[i] = GlycanBox.GlobalOGlycans[glycanBox.ModIds[i]];
                }
                sb.Append(string.Join(",", glycans.Select(p => p.Struc.ToString()).ToArray())); sb.Append("\t");
                                 
                string localizationLevel;
                sb.Append(localizationInfo(OGlycanBoxLocalization, out localizationLevel)); sb.Append("\t");

                sb.Append(localizationLevel); sb.Append("\t");
            }

            return sb.ToString();
        }

        public static Dictionary<string, string> MatchedIonDataDictionary(List<MatchedFragmentIon> matchedFragmentIons)
        {
            Dictionary<string, string> s = new Dictionary<string, string>();
            PsmTsvWriter.AddMatchedIonsData(s, matchedFragmentIons);
            return s;
        }

        public static string localizationInfo(List<Tuple<int, Tuple<int, int>[]>> OGlycanBoxLocalization, out string localizationLevel)
        {
            HashSet<int> allGlycanIds = new HashSet<int>(OGlycanBoxLocalization.Select(p=>p.Item2).SelectMany(p => p.Select(q => q.Item2)));

            Dictionary<int, int> seenGlycanIds = new Dictionary<int, int>();

            HashSet<int> seenGlycanBoxIds = new HashSet<int>(OGlycanBoxLocalization.Select(p=>p.Item1));

            //Dictionary<string, int>: mod-id, count
            Dictionary<string, int> seenModSite = new Dictionary<string, int>();
            foreach (var ogl in OGlycanBoxLocalization)
            {
                foreach (var og in ogl.Item2)
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

                    if (seenGlycanIds.ContainsKey(og.Item2))
                    {
                        seenGlycanIds[og.Item2] += 1;                   
                    }
                    else
                    {
                        seenGlycanIds.Add(og.Item2, 1);
                    }
                }
            }

            localizationLevel = "Level5";
            if (OGlycanBoxLocalization.Count == 1)
            {
                localizationLevel = "Level1";
            }
            else if (OGlycanBoxLocalization.Count > 1 && seenGlycanBoxIds.Count == 1)
            {
                if (seenModSite.Values.Where(p=>p == OGlycanBoxLocalization.Count).Count() > 0)
                {
                    localizationLevel = "Level2";
                }
                else
                {
                    localizationLevel = "Level3a";
                }              
            }
            else if (OGlycanBoxLocalization.Count > 1 && seenGlycanBoxIds.Count > 1)
            {
                if (seenModSite.Values.Where(p => p == OGlycanBoxLocalization.Count).Count() > 0)
                {
                    localizationLevel = "Level3b";
                }

                if (seenGlycanIds.Values.Where(p=>p== OGlycanBoxLocalization.Count).Count() >0)
                {
                    localizationLevel = "Level4";
                }

            }


            string local = "";
            //Some GSP have a lot paths, in which case only output first 10 paths and the total number of the paths.
            int maxOutputPath = 10;
            if (OGlycanBoxLocalization.Count <= maxOutputPath)
            {
                maxOutputPath = OGlycanBoxLocalization.Count;
            }

            int i = 0;
            while (i<maxOutputPath)
            {
                var ogl = OGlycanBoxLocalization[i];
                local += "{@" + ogl.Item1.ToString() + "[";
                var g = string.Join(",", ogl.Item2.Select(p => p.Item1.ToString() + "-" + p.Item2.ToString()));
                local += g + "]}";
                i++;
            }

            if (OGlycanBoxLocalization.Count > maxOutputPath)
            {
                local += "... In Total:" + OGlycanBoxLocalization.Count.ToString() + " Paths";
            }

            return local;
        }
    }
}