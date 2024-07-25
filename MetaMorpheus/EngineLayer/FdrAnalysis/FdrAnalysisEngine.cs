using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

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
            AllPsms = psms.ToList();
            MassDiffAcceptorNumNotches = massDiffAcceptorNumNotches;
            AnalysisType = analysisType;
            this.OutputFolder = outputFolder;
            this.DoPEP = doPEP;
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
                var psms = proteasePsms.OrderByDescending(p=>p).ToList();
                if (psms.Count > 100)
                {
                    ComputeCumulativeTargetAndDecoyCountsOnSortedPSMs(psms, true);
                    QValueInvertedPsms(psms);

                    if (DoPEP)
                    {
                        //PEP will model will be developed using peptides and then applied to all PSMs. 
                        Compute_PEPValue(myAnalysisResults, psms);
                        //some PSMs will be eliminated during the PEP calculation. So, we need to recompute the cumulative target and decoy counts
                        foreach (var psm in psms) psm.ResolveAllAmbiguities(); // Resolve ambiguities after PEP calculation before TDC

                        //peptiides are first ordered by PEP from good to bad and then by MM score from good to bad
                        // Peptide level and PSM level PEPs are identical
                        var peptides = psms
                            .OrderBy(p => p.PeptideFdrInfo.PEP)
                            .ThenByDescending(p => p)
                            .GroupBy(p => p.FullSequence)
                            .Select(p => p.FirstOrDefault())
                            .ToList();
                        ComputeCumulativeTargetAndDecoyCountsOnSortedPSMs(peptides, psmLevelCalculation: false, pepCalculation: true);

                        psms = psms.OrderBy(p => p.PsmFdrInfo.PEP).ThenByDescending(p => p).ToList();
                        ComputeCumulativeTargetAndDecoyCountsOnSortedPSMs(psms, psmLevelCalculation: true, pepCalculation: true);

                        //we do this section last so that target and decoy counts written in the psmtsv files are appropriate for the sort order which is by MM score
                        peptides = psms
                            .OrderByDescending(p => p)
                            .GroupBy(p => p.FullSequence)
                            .Select(p => p.FirstOrDefault())
                            .ToList();
                        ComputeCumulativeTargetAndDecoyCountsOnSortedPSMs(peptides, psmLevelCalculation: false, pepCalculation: false);

                        psms = psms.OrderByDescending(p => p).ToList();
                        ComputeCumulativeTargetAndDecoyCountsOnSortedPSMs(psms, psmLevelCalculation: true, pepCalculation: false);
                    }
                }
                else //psms .Count <= 100
                {
                    var peptides = psms
                        .OrderByDescending(p => p)
                        .GroupBy(b => b.FullSequence)
                        .Select(b => b.FirstOrDefault()).ToList();
                    ComputeCumulativeTargetAndDecoyCountsOnSortedPSMs(peptides, psmLevelCalculation: false, pepCalculation: false);
                    ComputeCumulativeTargetAndDecoyCountsOnSortedPSMs(psms.OrderByDescending(p => p).ToList(), psmLevelCalculation: true, pepCalculation: false);
                }
                CountPsm(psms);
            }
        }
        /// <summary>
        /// This methods assumes that PSMs are already sorted appropriately for downstream usage
        /// For traditional q-value calculation, the PSMs should be sorted from highest to lowest score
        /// For PEP q-value calculation, the PSMs should be sorted from lowest to highest PEP
        /// </summary>
        public void ComputeCumulativeTargetAndDecoyCountsOnSortedPSMs(List<SpectralMatch> psms, bool psmLevelCalculation, bool pepCalculation = false)
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
                    cumulativeDecoy ++;
                    cumulativeDecoyPerNotch[notch] ++;
                }
                else
                {
                    cumulativeTarget++;
                    cumulativeTargetPerNotch[notch]++;
                }

                if (psmLevelCalculation)
                {
                    psm.PsmFdrInfo.CumulativeDecoy = cumulativeDecoy;
                    psm.PsmFdrInfo.CumulativeTarget = cumulativeTarget;
                    psm.PsmFdrInfo.CumulativeDecoyNotch = cumulativeDecoyPerNotch[notch];
                    psm.PsmFdrInfo.CumulativeTargetNotch = cumulativeTargetPerNotch[notch];
                }
                else
                {
                    psm.PeptideFdrInfo.CumulativeDecoy = cumulativeDecoy;
                    psm.PeptideFdrInfo.CumulativeTarget = cumulativeTarget;
                    psm.PeptideFdrInfo.CumulativeDecoyNotch = cumulativeDecoyPerNotch[notch];
                    psm.PeptideFdrInfo.CumulativeTargetNotch = cumulativeTargetPerNotch[notch];
                }

            }

            if(psmLevelCalculation)
            {
                if (pepCalculation) PepQValueInvertedPsms(psms);
                else QValueInvertedPsms(psms);
            }
            else
            {
                if (pepCalculation) PepQValueInvertedPeptides(psms);
                else QValueInvertedPeptides(psms);
            }
        }
        
        private static void QValueInvertedPsms(List<SpectralMatch> psms)
        {
            psms.Reverse();
            //this calculation is performed from bottom up. So, we begin the loop by computing qValue
            //and qValueNotch for the last/lowest scoring psm in the bunch
            double qValue = (psms[0].PsmFdrInfo.CumulativeDecoy + 1) / psms[0].PsmFdrInfo.CumulativeTarget;
            double qValueNotch = (psms[0].PsmFdrInfo.CumulativeDecoyNotch + 1) / psms[0].PsmFdrInfo.CumulativeTargetNotch;

            //Assign FDR values to PSMs
            for (int i = 0; i < psms.Count; i++)
            {
                // Stop if canceled
                if (GlobalVariables.StopLoops) { break; }

                qValue = Math.Min(qValue, (psms[i].PsmFdrInfo.CumulativeDecoy + 1) / Math.Max(psms[i].PsmFdrInfo.CumulativeTarget,1));
                qValueNotch = Math.Min(qValueNotch, (psms[i].PsmFdrInfo.CumulativeDecoyNotch + 1) / Math.Max(psms[i].PsmFdrInfo.CumulativeTargetNotch,1));

                psms[i].PsmFdrInfo.QValue = Math.Min(qValue,1);
                psms[i].PsmFdrInfo.QValueNotch = Math.Min(qValueNotch,1);
            }
            psms.Reverse(); //we inverted the psms for this calculation. now we need to put them back into the original order
        }
        public static void QValueInvertedPeptides(List<SpectralMatch> psms)
        {
            psms.Reverse();
            //this calculation is performed from bottom up. So, we begin the loop by computing qValue
            //and qValueNotch for the last/lowest scoring psm in the bunch
            double qValue = (psms[0].PeptideFdrInfo.CumulativeDecoy + 1) / psms[0].PeptideFdrInfo.CumulativeTarget;
            double qValueNotch = (psms[0].PeptideFdrInfo.CumulativeDecoyNotch + 1) / psms[0].PeptideFdrInfo.CumulativeTargetNotch;

            //Assign FDR values to PSMs
            for (int i = 0; i < psms.Count; i++)
            {
                // Stop if canceled
                if (GlobalVariables.StopLoops) { break; }

                qValue = Math.Min(qValue, (psms[i].PeptideFdrInfo.CumulativeDecoy + 1) / psms[i].PeptideFdrInfo.CumulativeTarget);
                qValueNotch = Math.Min(qValueNotch, (psms[i].PeptideFdrInfo.CumulativeDecoyNotch + 1) / psms[i].PeptideFdrInfo.CumulativeTargetNotch);
                psms[i].PeptideFdrInfo.QValue = qValue;
                psms[i].PeptideFdrInfo.QValueNotch = qValueNotch;
            }
            psms.Reverse(); //we inverted the psms for this calculation. now we need to put them back into the original order
        }
        public static void PepQValueInvertedPsms(List<SpectralMatch> psms)
        {
            psms.Reverse();
            //this calculation is performed from bottom up. So, we begin the loop by computing qValue
            //and qValueNotch for the last/lowest scoring psm in the bunch
            double qValue = (psms[0].PsmFdrInfo.CumulativeDecoy + 1) / psms[0].PsmFdrInfo.CumulativeTarget;

            //Assign FDR values to PSMs
            for (int i = 0; i < psms.Count; i++)
            {
                // Stop if canceled
                if (GlobalVariables.StopLoops) { break; }

                psms[i].PsmFdrInfo.PEP_QValue = Math.Min(qValue, (psms[i].PsmFdrInfo.CumulativeDecoy + 1) / psms[i].PsmFdrInfo.CumulativeTarget);
            }
            psms.Reverse(); //we inverted the psms for this calculation. now we need to put them back into the original order
        }

        public static void PepQValueInvertedPeptides(List<SpectralMatch> psms)
        {
            psms.Reverse();
            //this calculation is performed from bottom up. So, we begin the loop by computing qValue
            //and qValueNotch for the last/lowest scoring psm in the bunch
            double qValue = (psms[0].PeptideFdrInfo.CumulativeDecoy + 1) / psms[0].PeptideFdrInfo.CumulativeTarget;

            //Assign FDR values to PSMs
            for (int i = 0; i < psms.Count; i++)
            {
                // Stop if canceled
                if (GlobalVariables.StopLoops) { break; }
                qValue = Math.Min(qValue, (psms[i].PeptideFdrInfo.CumulativeDecoy + 1) / psms[i].PeptideFdrInfo.CumulativeTarget);
                psms[i].PeptideFdrInfo.PEP_QValue = qValue;
            }
            psms.Reverse(); //we inverted the psms for this calculation. now we need to put them back into the original order
        }

        public void Compute_PEPValue(FdrAnalysisResults myAnalysisResults, List<SpectralMatch> psms)
        {
            if (psms[0].DigestionParams.Protease.Name == "top-down")
            {
                myAnalysisResults.BinarySearchTreeMetrics = PEP_Analysis_Cross_Validation.ComputePEPValuesForAllPSMsGeneric(psms, "top-down", this.FileSpecificParameters, this.OutputFolder, this);
            }
            else if (psms[0].DigestionParams.Protease.Name == "crosslink")
            {
                myAnalysisResults.BinarySearchTreeMetrics = PEP_Analysis_Cross_Validation.ComputePEPValuesForAllPSMsGeneric(psms, "crosslink", this.FileSpecificParameters, this.OutputFolder, this);
            }
            else
            {
                myAnalysisResults.BinarySearchTreeMetrics = PEP_Analysis_Cross_Validation.ComputePEPValuesForAllPSMsGeneric(psms, "standard", this.FileSpecificParameters, this.OutputFolder, this);
            }
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