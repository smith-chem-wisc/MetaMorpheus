using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using EngineLayer;
using GuiFunctions.ViewModels.Legends;
using Proteomics;
using Proteomics.ProteolyticDigestion;

namespace GuiFunctions
{
    public class ChimeraLegendViewModel : LegendViewModel
    {
        #region Private Properties

        private Dictionary<string, List<ChimeraLegendItemViewModel>> chimeraLegendItems;

        #endregion

        #region Public Properties

        public Dictionary<string, List<ChimeraLegendItemViewModel>> ChimeraLegendItems
        {
            get => chimeraLegendItems;
            set { chimeraLegendItems = value; OnPropertyChanged(nameof(ChimeraLegendItems)); }
        }

        #endregion

        public ChimeraLegendViewModel(List<PsmFromTsv> chimericIDs, double offset = 0) : base()
        {
            TopOffset = offset;
            if (MetaDrawSettings.ShowLegend)
            {
                Visibility = Visibility.Visible;
            }
            else
            {
                Visibility = Visibility.Hidden;
            }
            ChimeraLegendItems = new();
            ParseLegendItemsFromPsms(chimericIDs);
        }

        public void ParseLegendItemsFromPsms(List<PsmFromTsv> chimericIDs)
        {
            var groupedByProtein = chimericIDs.GroupBy(p => p.BaseSeq).OrderByDescending(p => p.Count());
            int proteinIndex = 0;
            foreach (var protein in groupedByProtein)
            {
                ChimeraLegendItems.Add(protein.Key, new List<ChimeraLegendItemViewModel>());
                if (protein.Count() > 1)
                {
                    ChimeraLegendItems[protein.Key].Add(new("Shared Ions", ChimeraSpectrumMatchPlot.ColorByProteinDictionary[proteinIndex][0]));
                }
                for (int i = 0; i < protein.Count(); i++)
                {
                    PeptideWithSetModifications peptideWithSetMods =
                        new(protein.ToList()[i].FullSequence.Split("|")[0], GlobalVariables.AllModsKnownDictionary);
                    var modsString = String.Join(", ",
                        peptideWithSetMods.AllModsOneIsNterminus.Select(p => p.Key + " - " + p.Value.IdWithMotif));
                    ChimeraLegendItems[protein.Key].Add(new(modsString, ChimeraSpectrumMatchPlot.ColorByProteinDictionary[proteinIndex][i + 1]));
                }

            }
        }
    }

    /// <summary>
    /// For design time display of data
    /// </summary>
    public class ChimeraLegendModel : ChimeraLegendViewModel
    {
        private static List<PsmFromTsv> chimericPsms { get; set; }
        public static ChimeraLegendModel Instance => new ChimeraLegendModel(chimericPsms);
        public ChimeraLegendModel(List<PsmFromTsv> chimericPsms) : base(chimericPsms)
        {

        }

        static ChimeraLegendModel()
        {
            string header =
                "File Name\tScan Number\tScan Retention Time\tNum Experimental Peaks\tTotal Ion Current\tPrecursor Scan Number\tPrecursor Charge\tPrecursor MZ\tPrecursor Mass\tScore\tDelta Score\tNotch\tBase Sequence\tFull Sequence\tEssential Sequence\tAmbiguity Level\tPSM Count (unambiguous, <0.01 q-value)\tMods\tMods Chemical Formulas\tMods Combined Chemical Formula\tNum Variable Mods\tMissed Cleavages\tPeptide Monoisotopic Mass\tMass Diff (Da)\tMass Diff (ppm)\tProtein Accession\tProtein Name\tGene Name\tOrganism Name\tIdentified Sequence Variations\tSplice Sites\tContaminant\tDecoy\tPeptide Description\tStart and End Residues In Protein\tPrevious Amino Acid\tNext Amino Acid\tTheoreticals Searched\tDecoy/Contaminant/Target\tMatched Ion Series\tMatched Ion Mass-To-Charge Ratios\tMatched Ion Mass Diff (Da)\tMatched Ion Mass Diff (Ppm)\tMatched Ion Intensities\tMatched Ion Counts\tNormalized Spectral Angle\tLocalized Scores\tImprovement Possible\tCumulative Target\tCumulative Decoy\tQValue\tCumulative Target Notch\tCumulative Decoy Notch\tQValue Notch\tPEP\tPEP_QValue";
            string psm1 = "FXN3_tr1_032017-calib\t1558\t63.25946\t2758\t3770149.922\t1557\t22\t700.57742\t15390.54319\t15.012\t1\t0\tARTKQTARKSTGGKAPRKQLATKAARKSAPATGGVKKPHRYRPGTVALREIRRYQKSTELLIRKLPFQRLVREIAQDFKTDLRFQSSAVMALQEASEAYLVGLFEDTNLCAIHAKRVTIMPKDIQLARRIRGERA\tARTKQTARK[UniProt:N6-(2-hydroxyisobutyryl)lysine on K]STGGKAPRKQLATKAARKSAPATGGVKKPHRYRPGTVALREIRRYQKSTELLIRKLPFQRLVREIAQDFKTDLRFQSSAVMALQEASEAYLVGLFEDTNLC[Common Fixed:Carbamidomethyl on C]AIHAKRVTIMPKDIQLARRIRGERA|ARTKQTARKSTGGK[UniProt:N6-(2-hydroxyisobutyryl)lysine on K]APRKQLATKAARKSAPATGGVKKPHRYRPGTVALREIRRYQKSTELLIRKLPFQRLVREIAQDFKTDLRFQSSAVMALQEASEAYLVGLFEDTNLC[Common Fixed:Carbamidomethyl on C]AIHAKRVTIMPKDIQLARRIRGERA|ARTKQTARKSTGGKAPRK[UniProt:N6-(2-hydroxyisobutyryl)lysine on K]QLATKAARKSAPATGGVKKPHRYRPGTVALREIRRYQKSTELLIRKLPFQRLVREIAQDFKTDLRFQSSAVMALQEASEAYLVGLFEDTNLC[Common Fixed:Carbamidomethyl on C]AIHAKRVTIMPKDIQLARRIRGERA|ARTKQTARKSTGGKAPRKQLATK[UniProt:N6-(2-hydroxyisobutyryl)lysine on K]AARKSAPATGGVKKPHRYRPGTVALREIRRYQKSTELLIRKLPFQRLVREIAQDFKTDLRFQSSAVMALQEASEAYLVGLFEDTNLC[Common Fixed:Carbamidomethyl on C]AIHAKRVTIMPKDIQLARRIRGERA|ARTKQTARKSTGGKAPRKQLATKAARK[UniProt:N6-(2-hydroxyisobutyryl)lysine on K]SAPATGGVKKPHRYRPGTVALREIRRYQKSTELLIRKLPFQRLVREIAQDFKTDLRFQSSAVMALQEASEAYLVGLFEDTNLC[Common Fixed:Carbamidomethyl on C]AIHAKRVTIMPKDIQLARRIRGERA|ARTKQTARKSTGGKAPRKQLATKAARKSAPATGGVK[UniProt:N6-(2-hydroxyisobutyryl)lysine on K]KPHRYRPGTVALREIRRYQKSTELLIRKLPFQRLVREIAQDFKTDLRFQSSAVMALQEASEAYLVGLFEDTNLC[Common Fixed:Carbamidomethyl on C]AIHAKRVTIMPKDIQLARRIRGERA\tARTKQTARK[UniProt:N6-(2-hydroxyisobutyryl)lysine on K]STGGKAPRKQLATKAARKSAPATGGVKKPHRYRPGTVALREIRRYQKSTELLIRKLPFQRLVREIAQDFKTDLRFQSSAVMALQEASEAYLVGLFEDTNLCAIHAKRVTIMPKDIQLARRIRGERA|ARTKQTARKSTGGK[UniProt:N6-(2-hydroxyisobutyryl)lysine on K]APRKQLATKAARKSAPATGGVKKPHRYRPGTVALREIRRYQKSTELLIRKLPFQRLVREIAQDFKTDLRFQSSAVMALQEASEAYLVGLFEDTNLCAIHAKRVTIMPKDIQLARRIRGERA|ARTKQTARKSTGGKAPRK[UniProt:N6-(2-hydroxyisobutyryl)lysine on K]QLATKAARKSAPATGGVKKPHRYRPGTVALREIRRYQKSTELLIRKLPFQRLVREIAQDFKTDLRFQSSAVMALQEASEAYLVGLFEDTNLCAIHAKRVTIMPKDIQLARRIRGERA|ARTKQTARKSTGGKAPRKQLATK[UniProt:N6-(2-hydroxyisobutyryl)lysine on K]AARKSAPATGGVKKPHRYRPGTVALREIRRYQKSTELLIRKLPFQRLVREIAQDFKTDLRFQSSAVMALQEASEAYLVGLFEDTNLCAIHAKRVTIMPKDIQLARRIRGERA|ARTKQTARKSTGGKAPRKQLATKAARK[UniProt:N6-(2-hydroxyisobutyryl)lysine on K]SAPATGGVKKPHRYRPGTVALREIRRYQKSTELLIRKLPFQRLVREIAQDFKTDLRFQSSAVMALQEASEAYLVGLFEDTNLCAIHAKRVTIMPKDIQLARRIRGERA|ARTKQTARKSTGGKAPRKQLATKAARKSAPATGGVK[UniProt:N6-(2-hydroxyisobutyryl)lysine on K]KPHRYRPGTVALREIRRYQKSTELLIRKLPFQRLVREIAQDFKTDLRFQSSAVMALQEASEAYLVGLFEDTNLCAIHAKRVTIMPKDIQLARRIRGERA\t2A\t0\tCarbamidomethyl on C N6-(2-hydroxyisobutyryl)lysine on K\tC6H9NO3\tC6H9NO3\t1\t0\t15390.52395\t0.01924\t1.25\tQ71DI3\tHistone H3.2\tprimary:H3C15, synonym:HIST2H3A, primary:H3C14, synonym:H3F2, synonym:H3FM, synonym:HIST2H3C, primary:H3C13, synonym:HIST2H3D\tHomo sapiens\t\t\tN\tN\tfull:M cleaved\t[2 to 136]\tM\t-\t \tT\t[b8+2, b47+7];[y12+3, y15+3, y22+4, y23+4, y24+4, y25+3, y26+5, y27+3, y29+3, y30+5, y31+4, y32+4, y34+4]\t[b8+2:457.26943, b47+7:716.41503];[y12+3:480.29562, y15+3:593.68657, y22+4:645.38557, y23+4:679.65062, y24+4:707.92216, y25+3:967.23788, y26+5:612.75224, y27+3:1058.27618, y29+3:1129.97223, y30+5:701.39234, y31+4:908.74910, y32+4:945.51573, y34+4:988.04208]\t[b8+2:-0.00104, b47+7:-0.01167];[y12+3:0.00097, y15+3:-0.00085, y22+4:0.00062, y23+4:0.00189, y24+4:0.00401, y25+3:-0.00085, y26+5:0.00151, y27+3:-0.00064, y29+3:-0.00310, y30+5:0.00041, y31+4:-0.00022, y32+4:-0.00209, y34+4:-0.00224]\t[b8+2:-1.14, b47+7:-2.33];[y12+3:0.67, y15+3:-0.48, y22+4:0.24, y23+4:0.69, y24+4:1.42, y25+3:-0.29, y26+5:0.49, y27+3:-0.20, y29+3:-0.91, y30+5:0.12, y31+4:-0.06, y32+4:-0.55, y34+4:-0.57]\t[b8+2:1121, b47+7:1541];[y12+3:11850, y15+3:586, y22+4:2542, y23+4:4115, y24+4:3082, y25+3:1639, y26+5:2692, y27+3:1347, y29+3:903, y30+5:4251, y31+4:5211, y32+4:3481, y34+4:1486]\t15\t0\t \t \t13304\t4\t0.00028\t6454\t2\t0.000281\t1.79E-07\t0";
            string psm2 = "FXN3_tr1_032017-calib\t1558\t63.25946\t2758\t3770149.922\t1557\t22\t702.57761\t15434.54742\t15.012\t1\t2\tARTKQTARKSTGGKAPRKQLATKAARKSAPATGGVKKPHRYRPGTVALREIRRYQKSTELLIRKLPFQRLVREIAQDFKTDLRFQSSAVMALQEASEAYLVGLFEDTNLCAIHAKRVTIMPKDIQLARRIRGERA\tARTKQTARKSTGGKAPRKQLATKAARK[UniProt:N6-(2-hydroxyisobutyryl)lysine on K]SAPATGGVKKPHRYRPGTVALREIRRYQK[UniProt:N6-acetyllysine on K]STELLIRKLPFQRLVREIAQDFKTDLRFQSSAVMALQEASEAYLVGLFEDTNLC[Common Fixed:Carbamidomethyl on C]AIHAKRVTIMPKDIQLARRIRGERA|ARTKQTARKSTGGKAPRKQLATKAARK[UniProt:N6-(2-hydroxyisobutyryl)lysine on K]SAPATGGVKKPHRYRPGTVALREIRRYQK[UniProt:N6,N6,N6-trimethyllysine on K]STELLIRKLPFQRLVREIAQDFKTDLRFQSSAVMALQEASEAYLVGLFEDTNLC[Common Fixed:Carbamidomethyl on C]AIHAKRVTIMPKDIQLARRIRGERA|ARTKQTARKSTGGKAPRKQLATKAARK[UniProt:N6-(2-hydroxyisobutyryl)lysine on K]SAPATGGVKKPHRYRPGTVALREIRRYQKSTELLIRKLPFQRLVREIAQDFK[UniProt:N6-acetyllysine on K]TDLRFQSSAVMALQEASEAYLVGLFEDTNLC[Common Fixed:Carbamidomethyl on C]AIHAKRVTIMPKDIQLARRIRGERA|ARTKQTARKSTGGKAPRKQLATKAARK[UniProt:N6-(2-hydroxyisobutyryl)lysine on K]SAPATGGVKKPHRYRPGTVALREIRRYQKSTELLIRKLPFQRLVREIAQDFK[UniProt:N6,N6,N6-trimethyllysine on K]TDLRFQSSAVMALQEASEAYLVGLFEDTNLC[Common Fixed:Carbamidomethyl on C]AIHAKRVTIMPKDIQLARRIRGERA|ARTKQTARKSTGGKAPRKQLATKAARKSAPATGGVK[UniProt:N6-(2-hydroxyisobutyryl)lysine on K]KPHRYRPGTVALREIRRYQK[UniProt:N6-acetyllysine on K]STELLIRKLPFQRLVREIAQDFKTDLRFQSSAVMALQEASEAYLVGLFEDTNLC[Common Fixed:Carbamidomethyl on C]AIHAKRVTIMPKDIQLARRIRGERA|ARTKQTARKSTGGKAPRKQLATKAARKSAPATGGVK[UniProt:N6-(2-hydroxyisobutyryl)lysine on K]KPHRYRPGTVALREIRRYQK[UniProt:N6,N6,N6-trimethyllysine on K]STELLIRKLPFQRLVREIAQDFKTDLRFQSSAVMALQEASEAYLVGLFEDTNLC[Common Fixed:Carbamidomethyl on C]AIHAKRVTIMPKDIQLARRIRGERA|ARTKQTARKSTGGKAPRKQLATKAARKSAPATGGVK[UniProt:N6-(2-hydroxyisobutyryl)lysine on K]KPHRYRPGTVALREIRRYQKSTELLIRKLPFQRLVREIAQDFK[UniProt:N6-acetyllysine on K]TDLRFQSSAVMALQEASEAYLVGLFEDTNLC[Common Fixed:Carbamidomethyl on C]AIHAKRVTIMPKDIQLARRIRGERA|ARTKQTARKSTGGKAPRKQLATKAARKSAPATGGVK[UniProt:N6-(2-hydroxyisobutyryl)lysine on K]KPHRYRPGTVALREIRRYQKSTELLIRKLPFQRLVREIAQDFK[UniProt:N6,N6,N6-trimethyllysine on K]TDLRFQSSAVMALQEASEAYLVGLFEDTNLC[Common Fixed:Carbamidomethyl on C]AIHAKRVTIMPKDIQLARRIRGERA\tARTKQTARKSTGGKAPRKQLATKAARK[UniProt:N6-(2-hydroxyisobutyryl)lysine on K]SAPATGGVKKPHRYRPGTVALREIRRYQK[UniProt:N6-acetyllysine on K]STELLIRKLPFQRLVREIAQDFKTDLRFQSSAVMALQEASEAYLVGLFEDTNLCAIHAKRVTIMPKDIQLARRIRGERA|ARTKQTARKSTGGKAPRKQLATKAARK[UniProt:N6-(2-hydroxyisobutyryl)lysine on K]SAPATGGVKKPHRYRPGTVALREIRRYQK[UniProt:N6,N6,N6-trimethyllysine on K]STELLIRKLPFQRLVREIAQDFKTDLRFQSSAVMALQEASEAYLVGLFEDTNLCAIHAKRVTIMPKDIQLARRIRGERA|ARTKQTARKSTGGKAPRKQLATKAARK[UniProt:N6-(2-hydroxyisobutyryl)lysine on K]SAPATGGVKKPHRYRPGTVALREIRRYQKSTELLIRKLPFQRLVREIAQDFK[UniProt:N6-acetyllysine on K]TDLRFQSSAVMALQEASEAYLVGLFEDTNLCAIHAKRVTIMPKDIQLARRIRGERA|ARTKQTARKSTGGKAPRKQLATKAARK[UniProt:N6-(2-hydroxyisobutyryl)lysine on K]SAPATGGVKKPHRYRPGTVALREIRRYQKSTELLIRKLPFQRLVREIAQDFK[UniProt:N6,N6,N6-trimethyllysine on K]TDLRFQSSAVMALQEASEAYLVGLFEDTNLCAIHAKRVTIMPKDIQLARRIRGERA|ARTKQTARKSTGGKAPRKQLATKAARKSAPATGGVK[UniProt:N6-(2-hydroxyisobutyryl)lysine on K]KPHRYRPGTVALREIRRYQK[UniProt:N6-acetyllysine on K]STELLIRKLPFQRLVREIAQDFKTDLRFQSSAVMALQEASEAYLVGLFEDTNLCAIHAKRVTIMPKDIQLARRIRGERA|ARTKQTARKSTGGKAPRKQLATKAARKSAPATGGVK[UniProt:N6-(2-hydroxyisobutyryl)lysine on K]KPHRYRPGTVALREIRRYQK[UniProt:N6,N6,N6-trimethyllysine on K]STELLIRKLPFQRLVREIAQDFKTDLRFQSSAVMALQEASEAYLVGLFEDTNLCAIHAKRVTIMPKDIQLARRIRGERA|ARTKQTARKSTGGKAPRKQLATKAARKSAPATGGVK[UniProt:N6-(2-hydroxyisobutyryl)lysine on K]KPHRYRPGTVALREIRRYQKSTELLIRKLPFQRLVREIAQDFK[UniProt:N6-acetyllysine on K]TDLRFQSSAVMALQEASEAYLVGLFEDTNLCAIHAKRVTIMPKDIQLARRIRGERA|ARTKQTARKSTGGKAPRKQLATKAARKSAPATGGVK[UniProt:N6-(2-hydroxyisobutyryl)lysine on K]KPHRYRPGTVALREIRRYQKSTELLIRKLPFQRLVREIAQDFK[UniProt:N6,N6,N6-trimethyllysine on K]TDLRFQSSAVMALQEASEAYLVGLFEDTNLCAIHAKRVTIMPKDIQLARRIRGERA\t3\t0\tCarbamidomethyl on C N6-(2-hydroxyisobutyryl)lysine on K N6-acetyllysine on K|Carbamidomethyl on C N6-(2-hydroxyisobutyryl)lysine on K N6,N6,N6-trimethyllysine on K|Carbamidomethyl on C N6-(2-hydroxyisobutyryl)lysine on K N6-acetyllysine on K|Carbamidomethyl on C N6-(2-hydroxyisobutyryl)lysine on K N6,N6,N6-trimethyllysine on K|Carbamidomethyl on C N6-(2-hydroxyisobutyryl)lysine on K N6-acetyllysine on K|Carbamidomethyl on C N6-(2-hydroxyisobutyryl)lysine on K N6,N6,N6-trimethyllysine on K|Carbamidomethyl on C N6-(2-hydroxyisobutyryl)lysine on K N6-acetyllysine on K|Carbamidomethyl on C N6-(2-hydroxyisobutyryl)lysine on K N6,N6,N6-trimethyllysine on K\tC8H11NO4|C9H15NO3|C8H11NO4|C9H15NO3|C8H11NO4|C9H15NO3|C8H11NO4|C9H15NO3\tC8H11NO4|C9H15NO3|C8H11NO4|C9H15NO3|C8H11NO4|C9H15NO3|C8H11NO4|C9H15NO3\t2\t0\t15432.53451|15432.57090|15432.53451|15432.57090|15432.53451|15432.57090|15432.53451|15432.57090\t2.01291|1.97652|2.01291|1.97652|2.01291|1.97652|2.01291|1.97652\t130.43|128.07|130.43|128.07|130.43|128.07|130.43|128.07\tQ71DI3\tHistone H3.2\tprimary:H3C15, synonym:HIST2H3A, primary:H3C14, synonym:H3F2, synonym:H3FM, synonym:HIST2H3C, primary:H3C13, synonym:HIST2H3D\tHomo sapiens\t\t\tN\tN\tfull:M cleaved\t[2 to 136]\tM\t-\t \tT\t[b8+2, b47+7];[y12+3, y15+3, y22+4, y23+4, y24+4, y25+3, y26+5, y27+3, y29+3, y30+5, y31+4, y32+4, y34+4]\t[b8+2:457.26943, b47+7:716.41503];[y12+3:480.29562, y15+3:593.68657, y22+4:645.38557, y23+4:679.65062, y24+4:707.92216, y25+3:967.23788, y26+5:612.75224, y27+3:1058.27618, y29+3:1129.97223, y30+5:701.39234, y31+4:908.74910, y32+4:945.51573, y34+4:988.04208]\t[b8+2:-0.00104, b47+7:-0.01167];[y12+3:0.00097, y15+3:-0.00085, y22+4:0.00062, y23+4:0.00189, y24+4:0.00401, y25+3:-0.00085, y26+5:0.00151, y27+3:-0.00064, y29+3:-0.00310, y30+5:0.00041, y31+4:-0.00022, y32+4:-0.00209, y34+4:-0.00224]\t[b8+2:-1.14, b47+7:-2.33];[y12+3:0.67, y15+3:-0.48, y22+4:0.24, y23+4:0.69, y24+4:1.42, y25+3:-0.29, y26+5:0.49, y27+3:-0.20, y29+3:-0.91, y30+5:0.12, y31+4:-0.06, y32+4:-0.55, y34+4:-0.57]\t[b8+2:1121, b47+7:1541];[y12+3:11850, y15+3:586, y22+4:2542, y23+4:4115, y24+4:3082, y25+3:1639, y26+5:2692, y27+3:1347, y29+3:903, y30+5:4251, y31+4:5211, y32+4:3481, y34+4:1486]\t15\t0\t \t \t13305\t4\t0.00028\t1681\t0\t0\t2.98E-07\t0";
            string psm3 = "FXN3_tr1_032017-calib\t1558\t63.25946\t2758\t3770149.922\t1557\t22\t703.1695\t15447.56897\t15.012\t1\t3\tARTKQTARKSTGGKAPRKQLATKAARKSAPSTGGVKKPHRYRPGTVALREIRRYQKSTELLIRKLPFQRLVREIAQDFKTDLRFQSAAIGALQEASEAYLVGLFEDTNLCAIHAKRVTIMPKDIQLARRIRGERA\tARTKQTARKSTGGKAPRKQLATKAARK[UniProt:N6-(2-hydroxyisobutyryl)lysine on K]SAPSTGGVKKPHRYRPGTVALREIRRYQK[UniProt:N6-glutaryllysine on K]STELLIRKLPFQRLVREIAQDFKTDLRFQSAAIGALQEASEAYLVGLFEDTNLC[Common Fixed:Carbamidomethyl on C]AIHAKRVTIMPKDIQLARRIRGERA|ARTKQTARKSTGGKAPRKQLATKAARK[UniProt:N6-(2-hydroxyisobutyryl)lysine on K]SAPSTGGVKKPHRYRPGTVALREIRRYQK[UniProt:N6-glutaryllysine on K]STELLIRKLPFQRLVREIAQDFKTDLRFQSAAIGALQEASEAYLVGLFEDTNLC[Common Fixed:Carbamidomethyl on C]AIHAKRVTIMPKDIQLARRIRGERA|ARTKQTARKSTGGKAPRKQLATKAARK[UniProt:N6-(2-hydroxyisobutyryl)lysine on K]SAPSTGGVKKPHRYRPGTVALREIRRYQKSTELLIRKLPFQRLVREIAQDFK[UniProt:N6-glutaryllysine on K]TDLRFQSAAIGALQEASEAYLVGLFEDTNLC[Common Fixed:Carbamidomethyl on C]AIHAKRVTIMPKDIQLARRIRGERA|ARTKQTARKSTGGKAPRKQLATKAARK[UniProt:N6-(2-hydroxyisobutyryl)lysine on K]SAPSTGGVKKPHRYRPGTVALREIRRYQKSTELLIRKLPFQRLVREIAQDFK[UniProt:N6-glutaryllysine on K]TDLRFQSAAIGALQEASEAYLVGLFEDTNLC[Common Fixed:Carbamidomethyl on C]AIHAKRVTIMPKDIQLARRIRGERA|ARTKQTARKSTGGKAPRKQLATKAARKSAPSTGGVK[UniProt:N6-(2-hydroxyisobutyryl)lysine on K]KPHRYRPGTVALREIRRYQK[UniProt:N6-glutaryllysine on K]STELLIRKLPFQRLVREIAQDFKTDLRFQSAAIGALQEASEAYLVGLFEDTNLC[Common Fixed:Carbamidomethyl on C]AIHAKRVTIMPKDIQLARRIRGERA|ARTKQTARKSTGGKAPRKQLATKAARKSAPSTGGVK[UniProt:N6-(2-hydroxyisobutyryl)lysine on K]KPHRYRPGTVALREIRRYQK[UniProt:N6-glutaryllysine on K]STELLIRKLPFQRLVREIAQDFKTDLRFQSAAIGALQEASEAYLVGLFEDTNLC[Common Fixed:Carbamidomethyl on C]AIHAKRVTIMPKDIQLARRIRGERA|ARTKQTARKSTGGKAPRKQLATKAARKSAPSTGGVK[UniProt:N6-(2-hydroxyisobutyryl)lysine on K]KPHRYRPGTVALREIRRYQKSTELLIRKLPFQRLVREIAQDFK[UniProt:N6-glutaryllysine on K]TDLRFQSAAIGALQEASEAYLVGLFEDTNLC[Common Fixed:Carbamidomethyl on C]AIHAKRVTIMPKDIQLARRIRGERA|ARTKQTARKSTGGKAPRKQLATKAARKSAPSTGGVK[UniProt:N6-(2-hydroxyisobutyryl)lysine on K]KPHRYRPGTVALREIRRYQKSTELLIRKLPFQRLVREIAQDFK[UniProt:N6-glutaryllysine on K]TDLRFQSAAIGALQEASEAYLVGLFEDTNLC[Common Fixed:Carbamidomethyl on C]AIHAKRVTIMPKDIQLARRIRGERA\tARTKQTARKSTGGKAPRKQLATKAARK[UniProt:N6-(2-hydroxyisobutyryl)lysine on K]SAPSTGGVKKPHRYRPGTVALREIRRYQK[UniProt:N6-glutaryllysine on K]STELLIRKLPFQRLVREIAQDFKTDLRFQSAAIGALQEASEAYLVGLFEDTNLCAIHAKRVTIMPKDIQLARRIRGERA|ARTKQTARKSTGGKAPRKQLATKAARK[UniProt:N6-(2-hydroxyisobutyryl)lysine on K]SAPSTGGVKKPHRYRPGTVALREIRRYQK[UniProt:N6-glutaryllysine on K]STELLIRKLPFQRLVREIAQDFKTDLRFQSAAIGALQEASEAYLVGLFEDTNLCAIHAKRVTIMPKDIQLARRIRGERA|ARTKQTARKSTGGKAPRKQLATKAARK[UniProt:N6-(2-hydroxyisobutyryl)lysine on K]SAPSTGGVKKPHRYRPGTVALREIRRYQKSTELLIRKLPFQRLVREIAQDFK[UniProt:N6-glutaryllysine on K]TDLRFQSAAIGALQEASEAYLVGLFEDTNLCAIHAKRVTIMPKDIQLARRIRGERA|ARTKQTARKSTGGKAPRKQLATKAARK[UniProt:N6-(2-hydroxyisobutyryl)lysine on K]SAPSTGGVKKPHRYRPGTVALREIRRYQKSTELLIRKLPFQRLVREIAQDFK[UniProt:N6-glutaryllysine on K]TDLRFQSAAIGALQEASEAYLVGLFEDTNLCAIHAKRVTIMPKDIQLARRIRGERA|ARTKQTARKSTGGKAPRKQLATKAARKSAPSTGGVK[UniProt:N6-(2-hydroxyisobutyryl)lysine on K]KPHRYRPGTVALREIRRYQK[UniProt:N6-glutaryllysine on K]STELLIRKLPFQRLVREIAQDFKTDLRFQSAAIGALQEASEAYLVGLFEDTNLCAIHAKRVTIMPKDIQLARRIRGERA|ARTKQTARKSTGGKAPRKQLATKAARKSAPSTGGVK[UniProt:N6-(2-hydroxyisobutyryl)lysine on K]KPHRYRPGTVALREIRRYQK[UniProt:N6-glutaryllysine on K]STELLIRKLPFQRLVREIAQDFKTDLRFQSAAIGALQEASEAYLVGLFEDTNLCAIHAKRVTIMPKDIQLARRIRGERA|ARTKQTARKSTGGKAPRKQLATKAARKSAPSTGGVK[UniProt:N6-(2-hydroxyisobutyryl)lysine on K]KPHRYRPGTVALREIRRYQKSTELLIRKLPFQRLVREIAQDFK[UniProt:N6-glutaryllysine on K]TDLRFQSAAIGALQEASEAYLVGLFEDTNLCAIHAKRVTIMPKDIQLARRIRGERA|ARTKQTARKSTGGKAPRKQLATKAARKSAPSTGGVK[UniProt:N6-(2-hydroxyisobutyryl)lysine on K]KPHRYRPGTVALREIRRYQKSTELLIRKLPFQRLVREIAQDFK[UniProt:N6-glutaryllysine on K]TDLRFQSAAIGALQEASEAYLVGLFEDTNLCAIHAKRVTIMPKDIQLARRIRGERA\t2A\t0\tCarbamidomethyl on C N6-(2-hydroxyisobutyryl)lysine on K N6-glutaryllysine on K\tC11H15NO6\tC11H15NO6\t2\t0\t15444.55227\t3.01669\t195.32\tP84243\tHistone H3.3\tprimary:H3-3A, synonym:H3.3A, synonym:H3F3, synonym:H3F3A, ORF:PP781, primary:H3-3B, synonym:H3.3B, synonym:H3F3B\tHomo sapiens\t\t\tN\tN\tfull:M cleaved|chain|full:M cleaved|chain|full:M cleaved|chain|full:M cleaved|chain\t[2 to 136]\tM\t-\t \tT\t[b8+2, b47+5];[y12+3, y15+3, y22+4, y23+4, y24+4, y25+3, y26+5, y27+3, y29+3, y30+5, y31+4, y32+4, y34+4]\t[b8+2:457.26943, b47+5:1005.77868];[y12+3:480.29562, y15+3:593.68657, y22+4:645.38557, y23+4:679.65062, y24+4:707.92216, y25+3:967.23788, y26+5:612.75224, y27+3:1058.27618, y29+3:1129.97223, y30+5:701.39234, y31+4:908.74910, y32+4:945.51573, y34+4:988.04208]\t[b8+2:-0.00104, b47+5:-0.00384];[y12+3:0.00097, y15+3:-0.00085, y22+4:0.00062, y23+4:0.00189, y24+4:0.00401, y25+3:-0.00085, y26+5:0.00151, y27+3:-0.00064, y29+3:-0.00310, y30+5:0.00041, y31+4:-0.00022, y32+4:-0.00209, y34+4:-0.00224]\t[b8+2:-1.14, b47+5:-0.77];[y12+3:0.67, y15+3:-0.48, y22+4:0.24, y23+4:0.69, y24+4:1.42, y25+3:-0.29, y26+5:0.49, y27+3:-0.20, y29+3:-0.91, y30+5:0.12, y31+4:-0.06, y32+4:-0.55, y34+4:-0.57]\t[b8+2:1121, b47+5:1161];[y12+3:11850, y15+3:586, y22+4:2542, y23+4:4115, y24+4:3082, y25+3:1639, y26+5:2692, y27+3:1347, y29+3:903, y30+5:4251, y31+4:5211, y32+4:3481, y34+4:1486]\t15\t0\t \t \t13308\t4\t0.00028\t1765\t0\t0\t5.96E-08\t0";
            Dictionary<string, int> parsedHeader = ParseHeader(header);

            List<PsmFromTsv> psms = new();
            psms.Add(new PsmFromTsv(psm1, "\t".ToCharArray(), parsedHeader));
            psms.Add(new PsmFromTsv(psm2, "\t".ToCharArray(), parsedHeader));
            psms.Add(new PsmFromTsv(psm3, "\t".ToCharArray(), parsedHeader));
            chimericPsms = psms;
        }

        private static Dictionary<string, int> ParseHeader(string header)
        {
            var parsedHeader = new Dictionary<string, int>();
            var spl = header.Split("\t");

            parsedHeader.Add(PsmTsvHeader.FullSequence, Array.IndexOf(spl, PsmTsvHeader.FullSequence));
            parsedHeader.Add(PsmTsvHeader.Ms2ScanNumber, Array.IndexOf(spl, PsmTsvHeader.Ms2ScanNumber));
            parsedHeader.Add(PsmTsvHeader.FileName, Array.IndexOf(spl, PsmTsvHeader.FileName));
            parsedHeader.Add(PsmTsvHeader.TotalIonCurrent, Array.IndexOf(spl, PsmTsvHeader.TotalIonCurrent));
            parsedHeader.Add(PsmTsvHeader.PrecursorScanNum, Array.IndexOf(spl, PsmTsvHeader.PrecursorScanNum));
            parsedHeader.Add(PsmTsvHeader.PrecursorCharge, Array.IndexOf(spl, PsmTsvHeader.PrecursorCharge));
            parsedHeader.Add(PsmTsvHeader.PrecursorMz, Array.IndexOf(spl, PsmTsvHeader.PrecursorMz));
            parsedHeader.Add(PsmTsvHeader.PrecursorMass, Array.IndexOf(spl, PsmTsvHeader.PrecursorMass));
            parsedHeader.Add(PsmTsvHeader.Score, Array.IndexOf(spl, PsmTsvHeader.Score));
            parsedHeader.Add(PsmTsvHeader.DeltaScore, Array.IndexOf(spl, PsmTsvHeader.DeltaScore));
            parsedHeader.Add(PsmTsvHeader.Notch, Array.IndexOf(spl, PsmTsvHeader.Notch));
            parsedHeader.Add(PsmTsvHeader.BaseSequence, Array.IndexOf(spl, PsmTsvHeader.BaseSequence));
            parsedHeader.Add(PsmTsvHeader.EssentialSequence, Array.IndexOf(spl, PsmTsvHeader.EssentialSequence));
            parsedHeader.Add(PsmTsvHeader.AmbiguityLevel, Array.IndexOf(spl, PsmTsvHeader.AmbiguityLevel));
            parsedHeader.Add(PsmTsvHeader.MissedCleavages, Array.IndexOf(spl, PsmTsvHeader.MissedCleavages));
            parsedHeader.Add(PsmTsvHeader.PeptideMonoMass, Array.IndexOf(spl, PsmTsvHeader.PeptideMonoMass));
            parsedHeader.Add(PsmTsvHeader.MassDiffDa, Array.IndexOf(spl, PsmTsvHeader.MassDiffDa));
            parsedHeader.Add(PsmTsvHeader.MassDiffPpm, Array.IndexOf(spl, PsmTsvHeader.MassDiffPpm));
            parsedHeader.Add(PsmTsvHeader.ProteinAccession, Array.IndexOf(spl, PsmTsvHeader.ProteinAccession));
            parsedHeader.Add(PsmTsvHeader.ProteinName, Array.IndexOf(spl, PsmTsvHeader.ProteinName));
            parsedHeader.Add(PsmTsvHeader.GeneName, Array.IndexOf(spl, PsmTsvHeader.GeneName));
            parsedHeader.Add(PsmTsvHeader.OrganismName, Array.IndexOf(spl, PsmTsvHeader.OrganismName));
            parsedHeader.Add(PsmTsvHeader.IntersectingSequenceVariations, Array.IndexOf(spl, PsmTsvHeader.IntersectingSequenceVariations));
            parsedHeader.Add(PsmTsvHeader.IdentifiedSequenceVariations, Array.IndexOf(spl, PsmTsvHeader.IdentifiedSequenceVariations));
            parsedHeader.Add(PsmTsvHeader.SpliceSites, Array.IndexOf(spl, PsmTsvHeader.SpliceSites));
            parsedHeader.Add(PsmTsvHeader.PeptideDesicription, Array.IndexOf(spl, PsmTsvHeader.PeptideDesicription));
            parsedHeader.Add(PsmTsvHeader.StartAndEndResiduesInProtein, Array.IndexOf(spl, PsmTsvHeader.StartAndEndResiduesInProtein));
            parsedHeader.Add(PsmTsvHeader.PreviousAminoAcid, Array.IndexOf(spl, PsmTsvHeader.PreviousAminoAcid));
            parsedHeader.Add(PsmTsvHeader.NextAminoAcid, Array.IndexOf(spl, PsmTsvHeader.NextAminoAcid));
            parsedHeader.Add(PsmTsvHeader.DecoyContaminantTarget, Array.IndexOf(spl, PsmTsvHeader.DecoyContaminantTarget));
            parsedHeader.Add(PsmTsvHeader.MatchedIonMzRatios, Array.IndexOf(spl, PsmTsvHeader.MatchedIonMzRatios));
            parsedHeader.Add(PsmTsvHeader.MatchedIonIntensities, Array.IndexOf(spl, PsmTsvHeader.MatchedIonIntensities));
            parsedHeader.Add(PsmTsvHeader.SpectralAngle, Array.IndexOf(spl, PsmTsvHeader.SpectralAngle));
            parsedHeader.Add(PsmTsvHeader.QValue, Array.IndexOf(spl, PsmTsvHeader.QValue));
            parsedHeader.Add(PsmTsvHeader.QValueNotch, Array.IndexOf(spl, PsmTsvHeader.QValueNotch));
            parsedHeader.Add(PsmTsvHeader.PEP, Array.IndexOf(spl, PsmTsvHeader.PEP));
            parsedHeader.Add(PsmTsvHeader.PEP_QValue, Array.IndexOf(spl, PsmTsvHeader.PEP_QValue));

            parsedHeader.Add(PsmTsvHeader.CrossTypeLabel, Array.IndexOf(spl, PsmTsvHeader.CrossTypeLabel));
            parsedHeader.Add(PsmTsvHeader.LinkResiduesLabel, Array.IndexOf(spl, PsmTsvHeader.LinkResiduesLabel));
            parsedHeader.Add(PsmTsvHeader.ProteinLinkSiteLabel, Array.IndexOf(spl, PsmTsvHeader.ProteinLinkSiteLabel));
            parsedHeader.Add(PsmTsvHeader.RankLabel, Array.IndexOf(spl, PsmTsvHeader.RankLabel));
            parsedHeader.Add(PsmTsvHeader.BetaPeptideProteinAccessionLabel, Array.IndexOf(spl, PsmTsvHeader.BetaPeptideProteinAccessionLabel));
            parsedHeader.Add(PsmTsvHeader.BetaPeptideProteinLinkSiteLabel, Array.IndexOf(spl, PsmTsvHeader.BetaPeptideProteinLinkSiteLabel));
            parsedHeader.Add(PsmTsvHeader.BetaPeptideBaseSequenceLabel, Array.IndexOf(spl, PsmTsvHeader.BetaPeptideBaseSequenceLabel));
            parsedHeader.Add(PsmTsvHeader.BetaPeptideFullSequenceLabel, Array.IndexOf(spl, PsmTsvHeader.BetaPeptideFullSequenceLabel));
            parsedHeader.Add(PsmTsvHeader.BetaPeptideTheoreticalMassLabel, Array.IndexOf(spl, PsmTsvHeader.BetaPeptideTheoreticalMassLabel));
            parsedHeader.Add(PsmTsvHeader.BetaPeptideScoreLabel, Array.IndexOf(spl, PsmTsvHeader.BetaPeptideScoreLabel));
            parsedHeader.Add(PsmTsvHeader.BetaPeptideRankLabel, Array.IndexOf(spl, PsmTsvHeader.BetaPeptideRankLabel));
            parsedHeader.Add(PsmTsvHeader.BetaPeptideMatchedIonsLabel, Array.IndexOf(spl, PsmTsvHeader.BetaPeptideMatchedIonsLabel));
            parsedHeader.Add(PsmTsvHeader.BetaPeptideMatchedIonIntensitiesLabel, Array.IndexOf(spl, PsmTsvHeader.BetaPeptideMatchedIonIntensitiesLabel));
            parsedHeader.Add(PsmTsvHeader.XLTotalScoreLabel, Array.IndexOf(spl, PsmTsvHeader.XLTotalScoreLabel));
            parsedHeader.Add(PsmTsvHeader.ParentIonsLabel, Array.IndexOf(spl, PsmTsvHeader.ParentIonsLabel));
            parsedHeader.Add(PsmTsvHeader.Ms2ScanRetentionTime, Array.IndexOf(spl, PsmTsvHeader.Ms2ScanRetentionTime));


            parsedHeader.Add(PsmTsvHeader_Glyco.GlycanMass, Array.IndexOf(spl, PsmTsvHeader_Glyco.GlycanMass));
            parsedHeader.Add(PsmTsvHeader_Glyco.GlycanStructure, Array.IndexOf(spl, PsmTsvHeader_Glyco.GlycanStructure));
            parsedHeader.Add(PsmTsvHeader_Glyco.GlycanComposition, Array.IndexOf(spl, PsmTsvHeader_Glyco.GlycanComposition));
            parsedHeader.Add(PsmTsvHeader_Glyco.GlycanLocalizationLevel, Array.IndexOf(spl, PsmTsvHeader_Glyco.GlycanLocalizationLevel));
            parsedHeader.Add(PsmTsvHeader_Glyco.LocalizedGlycan, Array.IndexOf(spl, PsmTsvHeader_Glyco.LocalizedGlycan));
            return parsedHeader;
        }
    }

}
