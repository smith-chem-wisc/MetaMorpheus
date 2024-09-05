using EngineLayer.CrosslinkSearch;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text.RegularExpressions;

namespace EngineLayer.FdrAnalysis
{
    public class FdrAnalysisEngine : MetaMorpheusEngine
    {
        private List<SpectralMatch> AllPsms;
        private readonly int MassDiffAcceptorNumNotches;
        private readonly string AnalysisType;
        private readonly string OutputFolder; // used for storing PEP training models  
        private readonly bool DoPEP;

        public FdrAnalysisEngine(List<SpectralMatch> psms, int massDiffAcceptorNumNotches, CommonParameters commonParameters,
            List<(string fileName, CommonParameters fileSpecificParameters)> fileSpecificParameters, List<string> nestedIds, string analysisType = "PSM", 
            bool doPEP = true, string outputFolder = null) : base(commonParameters, fileSpecificParameters, nestedIds)
        {
            AllPsms = psms.OrderByDescending(p => p).ToList();
            MassDiffAcceptorNumNotches = massDiffAcceptorNumNotches;
            AnalysisType = analysisType;
            OutputFolder = outputFolder;
            DoPEP = doPEP;
            if (AllPsms.Any())
                AddPsmAndPeptideFdrInfoIfNotPresent();
            if (fileSpecificParameters == null) throw new ArgumentNullException("file specific parameters cannot be null");
        }

        private void AddPsmAndPeptideFdrInfoIfNotPresent()
        {
            foreach (var psm in AllPsms.Where(p=> p.PsmFdrInfo == null))
            {
                psm.PsmFdrInfo = new FdrInfo();                
            }

            foreach (var psm in AllPsms.Where(p => p.PeptideFdrInfo == null))
            {
                psm.PeptideFdrInfo = new FdrInfo();
            }
        }

        protected override MetaMorpheusEngineResults RunSpecific()
        {
            FdrAnalysisResults myAnalysisResults = new FdrAnalysisResults(this, AnalysisType);

            Status("Running FDR analysis...");
            DoFalseDiscoveryRateAnalysis(myAnalysisResults);
            Status("Done.");
            myAnalysisResults.PsmsWithin1PercentFdr = AllPsms.Count(b => b.FdrInfo.QValue <= 0.01 && !b.IsDecoy);

            return myAnalysisResults;
        }

        private void DoFalseDiscoveryRateAnalysis(FdrAnalysisResults myAnalysisResults)
        {
            // Stop if canceled
            if (GlobalVariables.StopLoops) { return; }

            // calculate FDR on a per-protease basis (targets and decoys for a specific protease)
            var psmsGroupedByProtease = AllPsms.GroupBy(p => p.DigestionParams.Protease);

            foreach (var proteasePsms in psmsGroupedByProtease)
            {
                var psms = proteasePsms.OrderByDescending(p => p).ToList();
                var peptides = psms
                        .OrderByDescending(p => p)
                        .GroupBy(b => b.FullSequence)
                        .Select(b => b.FirstOrDefault())
                        .ToList();

                if (psms.Count > 100 & DoPEP)
                {
                    CalculateQValue(psms, peptideLevelCalculation: false, pepCalculation: false);
                    if (peptides.Count > 100 )
                    {
                        CalculateQValue(peptides, peptideLevelCalculation: true, pepCalculation: false);

                        //PEP will model will be developed using peptides and then applied to all PSMs. 
                        Compute_PEPValue(myAnalysisResults, psms);

                        // peptides are ordered by MM score from good to bad in order to select the best PSM for each peptide
                        peptides = psms
                            .OrderByDescending(p => p)
                            .GroupBy(p => p.FullSequence)
                            .Select(p => p.FirstOrDefault())
                            .OrderBy(p => p.FdrInfo.PEP) // Then order by PEP (PSM PEP and Peptide PEP are the same)
                            .ThenByDescending(p => p)
                            .ToList();
                        CalculateQValue(peptides, peptideLevelCalculation: true, pepCalculation: true);

                        psms = psms.OrderBy(p => p.PsmFdrInfo.PEP).ThenByDescending(p => p).ToList();
                        CalculateQValue(psms, peptideLevelCalculation: false, pepCalculation: true);
                        
                    }
                    else //we have more than 100 psms but less than 100 peptides so
                    {
                        //this will be done using PSMs because we dont' have enough peptides
                        Compute_PEPValue(myAnalysisResults, psms);
                        psms = psms.OrderBy(p => p.PsmFdrInfo.PEP).ThenByDescending(p => p).ToList();
                        CalculateQValue(psms, peptideLevelCalculation: false, pepCalculation: true);
                    }
                }
                else if(psms.Any(psm => psm.FdrInfo.PEP > 0)) 
                {
                    // If PEP's have been calculated, but doPEP = false, then we don't want to train another model,
                    // but we do want to calculate pep q-values
                    // really, in this case, we only need to run one or the other (i.e., only peptides or psms are passed in)
                    // but there's no mechanism to pass that info to the FDR analysis engine, so we'll do this for now
                    peptides = psms
                            .OrderByDescending(p => p)
                            .GroupBy(p => p.FullSequence)
                            .Select(p => p.FirstOrDefault()) // Get the best psm for each peptide based on MBR score
                            .OrderBy(p => p.FdrInfo.PEP) // Then order by PEP (PSM PEP and Peptide PEP are the same)
                            .ThenByDescending(p => p)
                            .ToList();
                    CalculateQValue(peptides, peptideLevelCalculation: true, pepCalculation: true);

                    psms = psms
                        .OrderBy(p => p.PsmFdrInfo.PEP)
                        .ThenByDescending(p => p)
                        .ToList();
                    CalculateQValue(psms, peptideLevelCalculation: false, pepCalculation: true);

                }

                //we do this section last so that target and decoy counts written in the psmtsv files are appropriate for the sort order which is by MM score
                peptides = psms
                    .OrderByDescending(p => p)
                    .GroupBy(b => b.FullSequence)
                    .Select(b => b.FirstOrDefault())
                    .ToList();
                CalculateQValue(peptides, peptideLevelCalculation: true, pepCalculation: false);

                psms = psms.OrderByDescending(p => p).ToList();
                CalculateQValue(psms, peptideLevelCalculation: false, pepCalculation: false);
                
                CountPsm(psms);
            }
        }

        /// <summary>
        /// This methods assumes that PSMs are already sorted appropriately for downstream usage
        /// Then, it counts the number of targets and (fractional) decoys, writes those values to the 
        /// appropriate FdrInfo (PSM or Peptide level), and calculates q-values
        /// </summary>
        public void CalculateQValue(List<SpectralMatch> psms, bool peptideLevelCalculation, bool pepCalculation = false)
        {
            double cumulativeTarget = 0;
            double cumulativeDecoy = 0;

            //set up arrays for local FDRs
            double[] cumulativeTargetPerNotch = new double[MassDiffAcceptorNumNotches + 1];
            double[] cumulativeDecoyPerNotch = new double[MassDiffAcceptorNumNotches + 1];

            //Assign FDR values to PSMs
            foreach (var psm in psms)
            {
                // Stop if canceled
                if (GlobalVariables.StopLoops) { break; }

                int notch = psm.Notch ?? MassDiffAcceptorNumNotches;
                if (psm.IsDecoy)
                {
                    // in that case, count it as the fraction of decoy hits
                    // e.g. if the PSM matched to 1 target and 2 decoys, it counts as 2/3 decoy
                    double decoyHits = 0;
                    double totalHits = 0;
                    var hits = psm.BestMatchingBioPolymersWithSetMods.GroupBy(p => p.Peptide.FullSequence);
                    foreach (var hit in hits)
                    {
                        if (hit.First().Peptide.Parent.IsDecoy)
                        {
                            decoyHits++;
                        }
                        totalHits++;
                    }
                    cumulativeDecoy += decoyHits / totalHits;
                    cumulativeDecoyPerNotch[notch] += decoyHits / totalHits;
                }
                else
                {
                    cumulativeTarget++;
                    cumulativeTargetPerNotch[notch]++;
                }

                psm.GetFdrInfo(peptideLevelCalculation).CumulativeDecoy = cumulativeDecoy;
                psm.GetFdrInfo(peptideLevelCalculation).CumulativeTarget = cumulativeTarget;
                psm.GetFdrInfo(peptideLevelCalculation).CumulativeDecoyNotch = cumulativeDecoyPerNotch[notch];
                psm.GetFdrInfo(peptideLevelCalculation).CumulativeTargetNotch = cumulativeTargetPerNotch[notch];

            }

            if (pepCalculation)
            {
                PepQValueInverted(psms, peptideLevelAnalysis: peptideLevelCalculation);
            }
            else
            {
                if(psms.Count < 100)
                {

                   QValueTraditional(psms, peptideLevelAnalysis: peptideLevelCalculation);
                }
                else
                {
                    QValueInverted(psms, peptideLevelAnalysis: peptideLevelCalculation);
                }
            }
        }

        /// <summary>
        /// This method is used only to calculate q-values for total PSM counts below 100
        /// </summary>
        private void QValueTraditional(List<SpectralMatch> psms, bool peptideLevelAnalysis)
        {
            double qValue = 0;
            double qValueNotch = 0;
            for (int i = 0; i < psms.Count; i++)
            {
                // Stop if canceled
                if (GlobalVariables.StopLoops) { break; }

                qValue = Math.Max(qValue, psms[i].GetFdrInfo(peptideLevelAnalysis).CumulativeDecoy / Math.Max(psms[i].GetFdrInfo(peptideLevelAnalysis).CumulativeTarget, 1));
                qValueNotch = Math.Max(qValueNotch, psms[i].GetFdrInfo(peptideLevelAnalysis).CumulativeDecoyNotch / Math.Max(psms[i].GetFdrInfo(peptideLevelAnalysis).CumulativeTargetNotch, 1));

                psms[i].GetFdrInfo(peptideLevelAnalysis).QValue = Math.Min(qValue, 1);
                psms[i].GetFdrInfo(peptideLevelAnalysis).QValueNotch = Math.Min(qValueNotch, 1);
            }
        }

        private static void QValueInverted(List<SpectralMatch> psms, bool peptideLevelAnalysis)
        {
            psms.Reverse();
            //this calculation is performed from bottom up. So, we begin the loop by computing qValue
            //and qValueNotch for the last/lowest scoring psm in the bunch
            double qValue = (psms[0].GetFdrInfo(peptideLevelAnalysis).CumulativeDecoy + 1) / psms[0].GetFdrInfo(peptideLevelAnalysis).CumulativeTarget;
            double qValueNotch = (psms[0].GetFdrInfo(peptideLevelAnalysis).CumulativeDecoyNotch + 1) / psms[0].GetFdrInfo(peptideLevelAnalysis).CumulativeTargetNotch;

            //Assign FDR values to PSMs
            for (int i = 0; i < psms.Count; i++)
            {
                // Stop if canceled
                if (GlobalVariables.StopLoops) { break; }

                qValue = Math.Min(qValue, (psms[i].GetFdrInfo(peptideLevelAnalysis).CumulativeDecoy + 1) / Math.Max(psms[i].GetFdrInfo(peptideLevelAnalysis).CumulativeTarget, 1));
                qValueNotch = Math.Min(qValueNotch, (psms[i].GetFdrInfo(peptideLevelAnalysis).CumulativeDecoyNotch + 1) / Math.Max(psms[i].GetFdrInfo(peptideLevelAnalysis).CumulativeTargetNotch, 1));

                psms[i].GetFdrInfo(peptideLevelAnalysis).QValue = Math.Min(qValue, 1);
                psms[i].GetFdrInfo(peptideLevelAnalysis).QValueNotch = Math.Min(qValueNotch, 1);
            }
            psms.Reverse(); //we inverted the psms for this calculation. now we need to put them back into the original order
        }

        public static void PepQValueInverted(List<SpectralMatch> psms, bool peptideLevelAnalysis)
        {
            psms.Reverse();
            //this calculation is performed from bottom up. So, we begin the loop by computing qValue
            //and qValueNotch for the last/lowest scoring psm in the bunch
            double qValue = (psms[0].GetFdrInfo(peptideLevelAnalysis).CumulativeDecoy + 1) / psms[0].GetFdrInfo(peptideLevelAnalysis).CumulativeTarget;

            //Assign FDR values to PSMs
            for (int i = 0; i < psms.Count; i++)
            {
                // Stop if canceled
                if (GlobalVariables.StopLoops) { break; }

                qValue = Math.Min(qValue, (psms[i].GetFdrInfo(peptideLevelAnalysis).CumulativeDecoy + 1) / psms[i].GetFdrInfo(peptideLevelAnalysis).CumulativeTarget);

                psms[i].GetFdrInfo(peptideLevelAnalysis).PEP_QValue = qValue;
            }
            psms.Reverse(); //we inverted the psms for this calculation. now we need to put them back into the original order
        }

        public void Compute_PEPValue(FdrAnalysisResults myAnalysisResults, List<SpectralMatch> psms)
        {
            string searchType;
            // Currently, searches of mixed data (bottom-up + top-down) are not supported
            // PEP will be calculated based on the search type of the first file/PSM in the list, which isn't ideal
            // This will be addressed in a future release
            switch(psms[0].DigestionParams.Protease.Name)
            {
               case "top-down":
                    searchType = "top-down";
                    break;
                default:
                    searchType = "standard";
                    break;
            }
            if (psms[0] is CrosslinkSpectralMatch)
            {
                searchType = "crosslink";
            }
            myAnalysisResults.BinarySearchTreeMetrics = new PepAnalysisEngine(psms, searchType, FileSpecificParameters, OutputFolder).ComputePEPValuesForAllPSMs();

        }

        /// <summary>
        /// This method gets the count of PSMs with the same full sequence (with q-value < 0.01) to include in the psmtsv output
        /// </summary>
        public void CountPsm(List<SpectralMatch> proteasePsms)
        {
            // exclude ambiguous psms and has a fdr cutoff = 0.01
            var allUnambiguousPsms = proteasePsms.Where(psm => psm.FullSequence != null).ToList();

            var unambiguousPsmsLessThanOnePercentFdr = allUnambiguousPsms.Where(psm =>
                    psm.FdrInfo.QValue <= 0.01
                    && psm.FdrInfo.QValueNotch <= 0.01)
                .GroupBy(p => p.FullSequence);

            Dictionary<string, int> sequenceToPsmCount = new Dictionary<string, int>();

            foreach (var sequenceGroup in unambiguousPsmsLessThanOnePercentFdr)
            {
                sequenceToPsmCount.TryAdd(sequenceGroup.First().FullSequence, sequenceGroup.Count());
            }

            foreach (SpectralMatch psm in allUnambiguousPsms)
            {
                if (sequenceToPsmCount.TryGetValue(psm.FullSequence, out int count))
                {
                    psm.PsmCount = count;
                }
            }
        }
    }
}