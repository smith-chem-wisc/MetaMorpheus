using MassSpectrometry;
using System.Collections.Generic;

namespace MetaMorpheusGUI
{
    public class MetaDrawPsm
    {
        public int ScanNum { get; private set; }
        public string FullSequence { get; private set; }
        public string FileName { get; private set; }
        public List<TheoreticalFragmentIon> FragmentIons { get; private set; }
        public int NumExperimentalPeaks { get; private set; }
        public double TotalIonCurrent { get; private set; }
        public int PrecursorScanNumber { get; private set; }
        public int PrecursorCharge { get; private set; }
        public double PrecursorMZ { get; private set; }
        public double PrecursorMass { get; private set; }
        public double Score { get; private set; }
        public double DeltaScore { get; private set; }
        public string Notch { get; private set; }
        public int DifferentPeakMatches { get; private set; }
        public string PeptidesSharingSamePeaks { get; private set; }
        public string BaseSequence { get; private set; }
        public string EssentialSequence { get; private set; }
        public string Mods { get; private set; }
        public string ModsChemicalFormulas { get; private set; }
        public string ModsCombinedChemicalFormula { get; private set; }
        public string NumVariableMods { get; private set; }
        public string MissedCleavages { get; private set; }
        public string PeptideMonoisotopicMass { get; private set; }
        public string MassDiffDa { get; private set; }
        public string MassDiffppm { get; private set; }
        public string ProteinAccession { get; private set; }
        public string ProteinName { get; private set; }
        public string GeneName { get; private set; }
        public string SequenceVariations { get; private set; }
        public string OrganismName { get; private set; }
        public string Contaminant { get; private set; }
        public string Decoy { get; private set; }
        public string PeptideDescription { get; private set; }
        public string StartandEndResiduesInProtein { get; private set; }
        public string PreviousAminoAcid { get; private set; }
        public string NextAminoAcid { get; private set; }
        public string AllScores { get; private set; }
        public string TheoreticalsSearched { get; private set; }
        public string DecoyContaminantTarget      {get;private set;}
        public string MatchedIonCounts { get; private set; }
        public string MatchedIonSeries { get; private set; }
        public string MatchedIonMassToChargeRatios{get;private set;}
        public string MatchedIonMassDiffDa { get; private set; }
        public string MatchedIonMassDiffPpm { get; private set; }
        public string MatchedIonIntensities { get; private set; }
        public string LocalizedScores { get; private set; }
        public string ImprovementPossible { get; private set; }
        public string CumulativeTarget { get; private set; }
        public string CumulativeDecoy { get; private set; }
        public string QValue { get; private set; }
        public string CumulativeTargetNotch { get; private set; }
        public string CumulativeDecoyNotch { get; private set; }
        
        public MetaDrawPsm(int oneBasedScanNumber, string fileName, string fullSequence, List<TheoreticalFragmentIon> fragmentIons, int expPeaks, double ionNum, int preScan, int preChar, double preMZ, double preMass, double score, double dlScore, string notch, int diffpeak,string pepshare, string bseq, string esseq, string mods, string modscf, string modsccf, string nvm, string misscleav, string pepmonomass, string mdda, string mdppm, string proteinacc, string proteinname, string gamename, string seqvar, string orgname, string con, string dec, string pepdes, string strtendresinprotein, string prevamino, string nxtamino, string allscores, string theosearch, string deccontaminatetarget, string matchioncount, string matchionseries, string matchionmasstocharge, string matchedionmdda, string matchedionmdppm, string matchedionint, string localscores, string improvposs, string cumutarget, string cumudeco, string qvalue, string cumutarnotch, string cumudecnotch)
        {
            this.ScanNum = oneBasedScanNumber;
            this.FileName = fileName;
            this.FullSequence = fullSequence;
            this.FragmentIons = fragmentIons;
            this.NumExperimentalPeaks = expPeaks;
            this.TotalIonCurrent = ionNum;
            this.PrecursorScanNumber = preScan;
            this.PrecursorCharge = preChar;
            this.PrecursorMZ = preMZ;
            this.PrecursorMass = preMass;
            this.Score = score;
            this.DeltaScore = dlScore;
            this.Notch = notch;
            this.DifferentPeakMatches = diffpeak;
            this.PeptidesSharingSamePeaks = pepshare;
            this.BaseSequence = bseq;
            this.EssentialSequence = esseq;
            this.Mods = mods;
            this.ModsChemicalFormulas = modscf;
            this.ModsCombinedChemicalFormula = modsccf;
            this.NumVariableMods = nvm; 
            this.MissedCleavages = misscleav; 
            this.PeptideMonoisotopicMass = pepmonomass; 
            this.MassDiffDa = mdda; 
            this.MassDiffppm = mdppm; 
            this.ProteinAccession = proteinacc; 
            this.ProteinName = proteinname; 
            this.GeneName = gamename; 
            this.SequenceVariations = seqvar; 
            this.OrganismName = orgname; 
            this.Contaminant = con; 
            this.Decoy = dec;  
            this.PeptideDescription = pepdes; 
            this.StartandEndResiduesInProtein = strtendresinprotein; 
            this.PreviousAminoAcid = prevamino;                
            this.NextAminoAcid = nxtamino;
            this.AllScores = allscores;
            this.TheoreticalsSearched = theosearch;
            this.DecoyContaminantTarget = deccontaminatetarget;
            this.MatchedIonCounts = matchioncount;
            this.MatchedIonSeries = matchionseries;
            this.MatchedIonMassToChargeRatios = matchionmasstocharge;
            this.MatchedIonMassDiffDa = matchedionmdda;
            this.MatchedIonMassDiffPpm = matchedionmdppm;
            this.MatchedIonIntensities = matchedionint;
            this.LocalizedScores = localscores;
            this.ImprovementPossible = improvposs;
            this.CumulativeTarget = cumutarget;
            this.CumulativeDecoy = cumudeco;
            this.QValue = qvalue;
            this.CumulativeTargetNotch = cumutarnotch;
            this.CumulativeDecoyNotch = cumudecnotch;
        }
    }
}
