using System;
using System.Collections.Generic;
using System.Linq;
using System.Text.RegularExpressions;
using EngineLayer;
using EngineLayer.FdrAnalysis;
using Newtonsoft.Json.Linq;

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
            List<(string fileName, CommonParameters fileSpecificParameters)> fileSpecificParameters, List<string> nestedIds, string analysisType = "PSM", bool doPEP = true, string outputFolder = null) : base(commonParameters, fileSpecificParameters, nestedIds)
        {
            AllPsms = psms.OrderBy(p=>p).ToList();
            MassDiffAcceptorNumNotches = massDiffAcceptorNumNotches;
            AnalysisType = analysisType;
            this.OutputFolder = outputFolder;
            this.DoPEP = doPEP;
            if (psms[0].FdrInfo == null)
                AddPsmAndPeptideFdrInfoIfNotPresent();
            if (fileSpecificParameters == null) throw new ArgumentNullException("file specific parameters cannot be null");
        }

        private void AddPsmAndPeptideFdrInfoIfNotPresent()
        {
            foreach (var psm in AllPsms)
            {
                psm.PsmFdrInfo = new FdrInfo();                
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
                var psms = proteasePsms.ToList();
                if (psms.Count > 100)
                {
                    var peptides = psms
                        .GroupBy(b => b.FullSequence)
                        .Select(b => b.FirstOrDefault()).ToList();
                    if (peptides.Count > 100)
                    {
                        ComputeCumulativeTargetAndDecoyCountsOnSortedPSMs(peptides, false);
                        QValueInvertedPeptides(peptides);
                        if (DoPEP)
                        {
                            //PEP will model will be developed using peptides and then applied to all PSMs. 
                            Compute_PEPValue(myAnalysisResults, psms);
                            //some PSMs will be eliminated during the PEP calculation. So, we need to recompute the cumulative target and decoy counts
                            //peptiides are first ordered by PEP from good to bad and then by MM score from good to bad
                            peptides = peptides.OrderBy(p => p.PeptideFdrInfo.PEP).ThenBy(p => p).ToList();
                            ComputeCumulativeTargetAndDecoyCountsOnSortedPSMs(peptides, false);
                            PepQValueInvertedPeptides(peptides);
                            psms = psms.OrderBy(p => p.PsmFdrInfo.PEP).ThenBy(p => p).ToList();
                            ComputeCumulativeTargetAndDecoyCountsOnSortedPSMs(psms, false);
                            PepQValueInvertedPsms(psms);

                            //we do this section last so that target and decoy counts written in the psmtsv files are appropriate for the sort order which is by MM score
                            peptides = peptides.OrderBy(p => p).ToList();
                            ComputeCumulativeTargetAndDecoyCountsOnSortedPSMs(peptides.OrderBy(p=>p).ToList(), false);
                            QValueInvertedPeptides(peptides);
                            psms = psms.OrderBy(p => p).ToList();
                            ComputeCumulativeTargetAndDecoyCountsOnSortedPSMs(psms, true);
                            QValueInvertedPsms(psms);
                        }
                    }
                    else //we have more than 100 psms but less than 100 peptides so
                    {
                        if (DoPEP)
                        {
                            //this will be done using PSMs because we dont' have enough peptides
                            Compute_PEPValue(myAnalysisResults, psms);
                            psms = psms.OrderBy(p => p.PsmFdrInfo.PEP).ThenBy(p => p).ToList();
                            ComputeCumulativeTargetAndDecoyCountsOnSortedPSMs(psms, false);
                            PepQValueInvertedPsms(psms);
                        }
                        //we do this section last so that target and decoy counts written in the psmtsv files are appropriate for the sort order which is by MM score
                        peptides = peptides.OrderBy(p => p).ToList();
                        ComputeCumulativeTargetAndDecoyCountsOnSortedPSMs(peptides.OrderBy(p => p).ToList(), false);
                        QValueTraditionalPeptides(peptides);
                        psms = psms.OrderBy(p => p).ToList();
                        ComputeCumulativeTargetAndDecoyCountsOnSortedPSMs(psms, true);
                        QValueInvertedPsms(psms);
                    }
                }
                else //psms .Count <= 100
                {
                    var peptides = psms
                        .GroupBy(b => b.FullSequence)
                        .Select(b => b.FirstOrDefault()).ToList();
                    ComputeCumulativeTargetAndDecoyCountsOnSortedPSMs(peptides.OrderBy(p => p).ToList(), false);
                    QValueTraditionalPeptides(peptides);
                    ComputeCumulativeTargetAndDecoyCountsOnSortedPSMs(psms.OrderBy(p => p).ToList(), true);
                    QValueTraditionalPsms(psms);
                }
                CountPsm(psms);
            }
        }
        /// <summary>
        /// This methods assumes that PSMs are already sorted appropriately for downstream usage
        /// For traditional q-value calculation, the PSMs should be sorted from highest to lowest score
        /// For PEP q-value calculation, the PSMs should be sorted from lowest to highest PEP
        /// </summary>
        private void ComputeCumulativeTargetAndDecoyCountsOnSortedPSMs(List<SpectralMatch> psms, bool isPsmNotPeptide)
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
                    // the PSM can be ambiguous between a target and a decoy sequence
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

                if (!isPsmNotPeptide)
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
        }
        /// <summary>
        /// This method is used only to calculate q-values for total PSM counts below 100
        /// </summary>
        private void QValueTraditionalPsms(List<SpectralMatch> psms)
        {
            double qValue = 0;
            double qValueNotch = 0;
            for (int i = 0; i < psms.Count; i++)
            {
                // Stop if canceled
                if (GlobalVariables.StopLoops) { break; }

                qValue = Math.Max(qValue, psms[i].PsmFdrInfo.CumulativeDecoy / psms[i].PsmFdrInfo.CumulativeTarget);
                qValueNotch = Math.Max(qValueNotch, psms[i].PsmFdrInfo.CumulativeDecoyNotch / psms[i].PsmFdrInfo.CumulativeTargetNotch);

                psms[i].PsmFdrInfo.QValue = qValue;
                psms[i].PsmFdrInfo.QValueNotch = qValueNotch;
            }
        }
        /// <summary>
        /// This method is used only to calculate q-values for total Peptide counts below 100
        /// </summary>
        /// <param name="psms"></param>
        private void QValueTraditionalPeptides(List<SpectralMatch> psms)
        {
            double qValue = 0;
            double qValueNotch = 0;
            for (int i = 0; i < psms.Count; i++)
            {
                // Stop if canceled
                if (GlobalVariables.StopLoops) { break; }

                qValue = Math.Max(qValue, psms[i].PeptideFdrInfo.CumulativeDecoy / psms[i].PeptideFdrInfo.CumulativeTarget);
                qValueNotch = Math.Max(qValueNotch, psms[i].PeptideFdrInfo.CumulativeDecoyNotch / psms[i].PeptideFdrInfo.CumulativeTargetNotch);

                psms[i].PeptideFdrInfo.QValue = qValue;
                psms[i].PeptideFdrInfo.QValueNotch = qValueNotch;
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

                psms[i].PsmFdrInfo.QValue = Math.Min(qValue, (psms[i].PsmFdrInfo.CumulativeDecoy + 1) / psms[i].PsmFdrInfo.CumulativeTarget);
                psms[i].PsmFdrInfo.QValueNotch = Math.Min(qValueNotch, (psms[i].PsmFdrInfo.CumulativeDecoyNotch + 1) / psms[i].PsmFdrInfo.CumulativeTargetNotch);
            }
            psms.Reverse(); //we inverted the psms for this calculation. now we need to put them back into the original order
        }

        private static void QValueInvertedPeptides(List<SpectralMatch> psms)
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

                psms[i].PeptideFdrInfo.QValue = Math.Min(qValue, (psms[i].PeptideFdrInfo.CumulativeDecoy + 1) / psms[i].PeptideFdrInfo.CumulativeTarget);
                psms[i].PeptideFdrInfo.QValueNotch = Math.Min(qValueNotch, (psms[i].PeptideFdrInfo.CumulativeDecoyNotch + 1) / psms[i].PeptideFdrInfo.CumulativeTargetNotch);
            }
            psms.Reverse(); //we inverted the psms for this calculation. now we need to put them back into the original order
        }
        private static void PepQValueInvertedPsms(List<SpectralMatch> psms)
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

        private static void PepQValueInvertedPeptides(List<SpectralMatch> psms)
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

                psms[i].PeptideFdrInfo.PEP_QValue = Math.Min(qValue, (psms[i].PeptideFdrInfo.CumulativeDecoy + 1) / psms[i].PeptideFdrInfo.CumulativeTarget);
            }
            psms.Reverse(); //we inverted the psms for this calculation. now we need to put them back into the original order
        }

        public void Compute_PEPValue(FdrAnalysisResults myAnalysisResults, List<SpectralMatch> psms)
        {
            if (psms[0].DigestionParams.Protease.Name != "top-down")
            {
                myAnalysisResults.BinarySearchTreeMetrics = PEP_Analysis_Cross_Validation.ComputePEPValuesForAllPSMsGeneric(psms, AnalysisType, this.FileSpecificParameters, this.OutputFolder);
            }
            else
            {
                myAnalysisResults.BinarySearchTreeMetrics = PEP_Analysis_Cross_Validation.ComputePEPValuesForAllPSMsGeneric(psms, "top-down", this.FileSpecificParameters, this.OutputFolder);
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