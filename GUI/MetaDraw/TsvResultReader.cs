using System;
using System.Collections.Generic;
using System.Linq;
using System.IO;
using System.Text.RegularExpressions;
using System.Windows;
using MassSpectrometry;
using Chemistry;

namespace MetaMorpheusGUI
{
    public class TsvResultReader
    {
        //1.
        private const string FULL_SEQUENCE_HEADER = "Full Sequence";
        private const string SCAN_NUMBER_HEADER = "Scan Number";
        private const string FILENAME_HEADER = "File Name";
        private const string MATCHED_MZ_HEADER = "Matched Ion Mass-To-Charge Ratios";
        private const string NUM_PEAK = "Num Experimental Peaks";
        private const string NUM_ION = "Total Ion Current";
        private const string NUM_PreScan = "Precursor Scan Number";
        private const string NUM_PreChar = "Precursor Charge";
        private const string NUM_PreMZ = "Precursor MZ";
        private const string NUM_PreMass = "Precursor Mass"; 
        private const string SCORE = "Score";
        private const string DLSCORE = "Delta Score";
        private const string NOTCH = "Notch";
        private const string DIFFPEAKMATCH = "Different Peak Matches";
        private const string PEPSHARESAMEPEAK = "Peptides Sharing Same Peaks";
        private const string BSEQ = "Base Sequence";
        private const string ESSEQ = "Essential Sequence";
        private const string MOD = "Mods";
        private const string MODSCF = "Mods Chemical Formulas";
        private const string MODSCCF= "Mods Combined Chemical Formula";
        private const string NVM = "Num Variable Mods";
        private const string MISSCLEAV = "Missed Cleavages";
        private const string PEPMONOMASS = "Peptide Monoisotopic Mass";
        private const string MDDA = "Mass Diff (Da)";
        private const string MDPPM = "Mass Diff (ppm)";
        private const string PROTEINACC = "Protein Accession";
        private const string PROTEINNAME = "Protein Name";
        private const string GAMENAME = "Gene Name";
        private const string SEQVAR = "Sequence Variations";
        private const string ORGNAME = "Organism Name";
        private const string CON = "Contaminant";
        private const string DEC = "Decoy";
        private const string PEPDES = "Peptide Description";
        private const string STRTENDRESINPROTEIN = "Start and End Residues In Protein";
        private const string PREVAMINO = "Previous Amino Acid";
        private const string NXTAMINO = "Next Amino Acid";
        private const string ALLSCORES = "All Scores";
        private const string THEOSEARCH = "Theoreticals Searched";
        private const string DECCONTAMINATETARGET = "Decoy/Contaminant/Target";
        private const string MATCHIONCOUNT = "Matched Ion Counts";
        private const string MATCHIONSERIES = "Matched Ion Series";
        private const string MATCHIONMASSTOCHARGE = "Matched Ion Mass-To-Charge Ratios";
        private const string MATCHEDIONMDDA = "Matched Ion Mass Diff (Da)";
        private const string MATCHEDIONMDPPM = "Matched Ion Mass Diff (Ppm)";
        private const string MATCHEDIONINT = "Matched Ion Intensities";
        private const string LOCALSCORES = "Localized Scores";           
        private const string IMPROVPOSS = "Improvement Possible";
        private const string CUMUTARGET = "Cumulative Target";
        private const string CUMUDECO = "Cumulative Decoy";
        private const string QVALUE = "QValue";
        private const string CUMUTARNOTCH = "Cumulative Target Notch";
        private const string CUMUDECNOTCH = "Cumulative Decoy Notch";
        

        private static Regex ionParser = new Regex(@"([a-zA-Z]+)(\d+)");
        private static char[] split = new char[] { '\t' };
        private static char[] mzSplit = new char[] { '[', ',', ']', ';' };
        
        public static List<MetaDrawPsm> ReadTsv(string filePath)
        {
            List<MetaDrawPsm> psms = new List<MetaDrawPsm>();

            StreamReader reader = null;
            try
            {
                reader = new StreamReader(filePath);
            }
            catch (Exception e)
            {
                MessageBox.Show("Could not read the file: " + e.Message);
                return psms;
            }

            //2.
            int lineCount = 0, fullSeqIndex = -1, scanNumberIndex = -1, fileNameIndex = -1, matchedMzIndex = -1, peakExpIndex=-1, ionNumIndex=-1, preScanIndex = -1, preCharIndex=-1, preMZIndex=1, preMassIndex=-1, scoreIndex=-1, dlscoreIndex=-1, notchIndex=-1, diffPeakIndex=-1,pepShareIndex=-1, bseqIndex=-1, esseqIndex=-1, modsIndex=-1, modscfIndex = -1, modsccfIndex = -1, nvmIndex = -1, misscleavIndex = -1, pepmonomassIndex = -1, mddaIndex = -1, mdppmIndex = -1, proteinaccIndex = -1, proteinnameIndex = -1, gamenameIndex = -1, seqvarIndex=-1, organameIndex = -1, conIndex = -1, decIndex=-1, pepdesIndex = -1, strtendIndex = -1, prevamIndex = -1, nxtamIndex = -1, allscoreIndex = -1, theosearchIndex = -1, deccontarIndex = -1, matchcountIndex = -1, matchseriesIndex = -1, matchmasstochargeIndex = -1, matchmddaIndex = -1, matchmdppmIndex = -1, matchintIndex = -1, localscoreIndex = -1, improveIndex = -1, cumutarIndex = -1, cumudecIndex = -1, qvalueIndex = -1, cumutarnotchIndex = -1, cumudecnotchIndex = -1; //Index=-1

            string line;
            while (reader.Peek() > 0)
            {
                lineCount++;

                line = reader.ReadLine();
                var spl = line.Split(split);

                if (lineCount == 1)
                {
                    //3.
                    fullSeqIndex = Array.IndexOf(spl, FULL_SEQUENCE_HEADER);
                    scanNumberIndex = Array.IndexOf(spl, SCAN_NUMBER_HEADER);
                    fileNameIndex = Array.IndexOf(spl, FILENAME_HEADER);
                    matchedMzIndex = Array.IndexOf(spl, MATCHED_MZ_HEADER);
                    peakExpIndex = Array.IndexOf(spl,NUM_PEAK);
                    ionNumIndex = Array.IndexOf(spl, NUM_ION);
                    preScanIndex = Array.IndexOf(spl, NUM_PreScan);
                    preCharIndex = Array.IndexOf(spl, NUM_PreChar);
                    preMZIndex = Array.IndexOf(spl, NUM_PreMZ);
                    preMassIndex = Array.IndexOf(spl, NUM_PreMass);
                    scoreIndex = Array.IndexOf(spl, SCORE);
                    dlscoreIndex= Array.IndexOf(spl, DLSCORE);
                    notchIndex = Array.IndexOf(spl, NOTCH);
                    diffPeakIndex = Array.IndexOf(spl, DIFFPEAKMATCH);
                    pepShareIndex = Array.IndexOf(spl, PEPSHARESAMEPEAK);
                    bseqIndex = Array.IndexOf(spl, BSEQ);
                    esseqIndex = Array.IndexOf(spl, ESSEQ);
                    modsIndex = Array.IndexOf(spl, MOD);
                    modscfIndex = Array.IndexOf(spl,MODSCF);
                    modsccfIndex = Array.IndexOf(spl,MODSCCF);
                    nvmIndex = Array.IndexOf(spl,NVM);
                    misscleavIndex = Array.IndexOf(spl,MISSCLEAV);
                    pepmonomassIndex = Array.IndexOf(spl,PEPMONOMASS);
                    mddaIndex = Array.IndexOf(spl,MDDA);
                    mdppmIndex = Array.IndexOf(spl,MDPPM);
                    proteinaccIndex = Array.IndexOf(spl,PROTEINACC);
                    proteinnameIndex = Array.IndexOf(spl,PROTEINNAME);
                    gamenameIndex = Array.IndexOf(spl,GAMENAME);
                    seqvarIndex = Array.IndexOf(spl,SEQVAR);
                    organameIndex = Array.IndexOf(spl,ORGNAME);
                    conIndex = Array.IndexOf(spl,CON);
                    decIndex = Array.IndexOf(spl,DEC);
                    pepdesIndex = Array.IndexOf(spl,PEPDES);
                    strtendIndex = Array.IndexOf(spl,STRTENDRESINPROTEIN);
                    prevamIndex = Array.IndexOf(spl,PREVAMINO);
                    nxtamIndex = Array.IndexOf(spl,NXTAMINO);
                    allscoreIndex = Array.IndexOf(spl,ALLSCORES);
                    theosearchIndex = Array.IndexOf(spl,THEOSEARCH);
                    deccontarIndex = Array.IndexOf(spl,DECCONTAMINATETARGET);
                    matchcountIndex = Array.IndexOf(spl,MATCHIONCOUNT);
                    matchseriesIndex = Array.IndexOf(spl,MATCHIONSERIES);
                    matchmasstochargeIndex = Array.IndexOf(spl,MATCHIONMASSTOCHARGE);
                    matchmddaIndex = Array.IndexOf(spl,MATCHEDIONMDDA);
                    matchmdppmIndex = Array.IndexOf(spl,MATCHEDIONMDPPM);
                    matchintIndex = Array.IndexOf(spl,MATCHEDIONINT);
                    localscoreIndex = Array.IndexOf(spl,LOCALSCORES);
                    improveIndex = Array.IndexOf(spl,IMPROVPOSS);
                    cumutarIndex = Array.IndexOf(spl,CUMUTARGET);
                    cumudecIndex = Array.IndexOf(spl,CUMUDECO);
                    qvalueIndex = Array.IndexOf(spl,QVALUE);
                    cumutarnotchIndex = Array.IndexOf(spl,CUMUTARNOTCH);
                    cumudecnotchIndex = Array.IndexOf(spl,CUMUDECNOTCH);
                    continue;
                }

                //4.
                if (fullSeqIndex < 0 || scanNumberIndex < 0 || fileNameIndex < 0 || matchedMzIndex < 0 || peakExpIndex < 0 || ionNumIndex < 0 || preScanIndex<0 || preCharIndex<0 || preMZIndex<0 || preMassIndex<0 || scoreIndex<0 || dlscoreIndex<0 || diffPeakIndex<0 || pepShareIndex<0 || bseqIndex<0 || esseqIndex<0 || modsIndex<0 || modscfIndex < 0 || modsccfIndex < 0 || nvmIndex < 0 || misscleavIndex < 0 || pepmonomassIndex < 0 || mddaIndex < 0 || mdppmIndex < 0 || proteinaccIndex < 0 || proteinnameIndex < 0 || gamenameIndex < 0 || seqvarIndex < 0 || organameIndex<0 || pepdesIndex < 0 || strtendIndex < 0 || prevamIndex < 0 || nxtamIndex < 0 || allscoreIndex < 0 || theosearchIndex < 0 || deccontarIndex < 0 || matchcountIndex < 0 || matchseriesIndex < 0 || matchmasstochargeIndex < 0 || matchmddaIndex < 0 || matchmdppmIndex < 0 || matchintIndex < 0 || localscoreIndex < 0 || improveIndex < 0 || cumutarIndex < 0 || cumudecIndex < 0 || qvalueIndex < 0 || cumutarnotchIndex < 0 || cumudecnotchIndex < 0 || conIndex<0 || decIndex<0)//
                {
                    MessageBox.Show("Could not read PSMs file. Is it from an older version of MetaMorpheus?");
                    return psms;
                }

                try
                {
                    //5.
                    int oneBasedScanNumber = int.Parse(spl[scanNumberIndex]);
                    int preScan = int.Parse(spl[preScanIndex]);
                    int peakNum = (int)double.Parse(spl[peakExpIndex]);
                    int preChar = (int)double.Parse(spl[preCharIndex]);
                    int diffPeak = (int)double.Parse(spl[diffPeakIndex]);
                    string peptideSequence = spl[fullSeqIndex];
                    string fileName = spl[fileNameIndex];
                    string matchedPeaks = spl[matchedMzIndex];
                    string notch = spl[notchIndex].Trim();
                    string pepshare = spl[pepShareIndex].Trim();
                    string bseq = spl[bseqIndex].Trim();
                    string esseq = spl[esseqIndex].Trim();
                    string mods = spl[modsIndex].Trim();
                    string modscf = spl[modscfIndex].Trim();
                    string modsccf = spl[modsccfIndex].Trim();
                    string nvm = spl[nvmIndex].Trim();
                    string misscleav = spl[misscleavIndex].Trim();
                    string pepmonomass = spl[pepmonomassIndex].Trim();
                    string mdda = spl[mddaIndex].Trim();
                    string mdppm = spl[mdppmIndex].Trim();
                    string proteinacc = spl[proteinaccIndex].Trim();
                    string proteinname = spl[proteinnameIndex].Trim();
                    string gamename = spl[gamenameIndex].Trim();
                    string seqvar = spl[seqvarIndex].Trim();
                    string orgname = spl[organameIndex].Trim();
                    string con = spl[conIndex].Trim();
                    string dec = spl[decIndex].Trim();
                    string pepdes = spl[pepdesIndex].Trim();
                    string strtendresinprotein = spl[strtendIndex].Trim();
                    string prevamino = spl[prevamIndex].Trim();
                    string nxtamino = spl[nxtamIndex].Trim();
                    string allscores = spl[allscoreIndex].Trim();
                    string theosearch = spl[theosearchIndex].Trim();
                    string deccontaminatetarget = spl[deccontarIndex].Trim();
                    string matchioncount = spl[matchcountIndex].Trim();
                    string matchionseries = spl[matchseriesIndex].Trim();
                    string matchionmasstocharge = spl[matchmasstochargeIndex].Trim();
                    string matchedionmdda = spl[matchmddaIndex].Trim();
                    string matchedionmdppm = spl[matchmdppmIndex].Trim();
                    string matchedionint = spl[matchintIndex].Trim();
                    string localscores = spl[localscoreIndex].Trim();
                    string improvposs = spl[improveIndex].Trim();
                    string cumutarget = spl[cumutarIndex].Trim();
                    string cumudeco = spl[cumudecIndex].Trim();
                    string qvalue = spl[qvalueIndex].Trim();
                    string cumutarnotch = spl[cumutarnotchIndex].Trim();
                    string cumudecnotch = spl[cumudecnotchIndex].Trim();
                    List<TheoreticalFragmentIon> peaks = ReadFragmentIonsFromString(matchedPeaks);
                    double ionNum = double.Parse(spl[ionNumIndex].TrimEnd(new char[] { '0'}));
                    double preMZ=double.Parse(spl[preMZIndex].TrimEnd(new char[] { '0' }));
                    double preMass =double.Parse(spl[preMassIndex].TrimEnd(new char[] { '0' }));
                    double score = double.Parse(spl[scoreIndex].TrimEnd(new char[] { '0' }));
                    double dlScore=double.Parse(spl[dlscoreIndex].TrimEnd(new char[] { '0' }));
                    psms.Add(new MetaDrawPsm(oneBasedScanNumber, fileName, peptideSequence, peaks,peakNum, ionNum, preScan, preChar, preMZ, preMass, score, dlScore, notch,diffPeak, pepshare, bseq, esseq, mods, modscf, modsccf, nvm, misscleav, pepmonomass, mdda, mdppm, proteinacc, proteinname, gamename, seqvar, orgname, con, dec, pepdes, strtendresinprotein, prevamino, nxtamino, allscores, theosearch, deccontaminatetarget, matchioncount, matchionseries, matchionmasstocharge, matchedionmdda, matchedionmdppm, matchedionint, localscores, improvposs, cumutarget, cumudeco, qvalue, cumutarnotch, cumudecnotch));
                }
                catch (Exception)
                {
                    MessageBox.Show(lineCount + "- Exception in TSV Reader");
                }
            }

            reader.Close();

            if ((lineCount - 1) != psms.Count)
            {
                MessageBox.Show("Warning: " + ((lineCount - 1) - psms.Count) + " PSMs were not read.");
            }

            return psms;
        }

        private static List<TheoreticalFragmentIon> ReadFragmentIonsFromString(string matchedMzString)
        {
            var peaks = matchedMzString.Split(mzSplit, StringSplitOptions.RemoveEmptyEntries).Select(v => v.Trim());
            List<TheoreticalFragmentIon> matchedIons = new List<TheoreticalFragmentIon>();

            foreach (var peak in peaks)
            {
                var split = peak.Split(new char[] { '+', ':' });

                string ionTypeAndNumber = split[0];
                Match result = ionParser.Match(ionTypeAndNumber);
                string ionType = result.Groups[1].Value;
                int ionNumber = int.Parse(result.Groups[2].Value);

                ProductType p = ProductType.None;
                if (ionType.Contains("b"))
                    p = ProductType.B;
                else if (ionType.Contains("y"))
                    p = ProductType.Y;
                else if (ionType.Contains("c"))
                    p = ProductType.C;
                else if (ionType.Contains("z"))
                    p = ProductType.Zdot;

                int z = int.Parse(split[1]);
                double mz = double.Parse(split[2]);

                matchedIons.Add(new TheoreticalFragmentIon(mz.ToMass(z), double.NaN, z, p, ionNumber));
            }

            return matchedIons;
        }
    }
}
