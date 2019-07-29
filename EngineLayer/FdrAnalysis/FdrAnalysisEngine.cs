using System;
using System.Collections.Generic;
using System.Linq;

namespace EngineLayer.FdrAnalysis
{
    public class FdrAnalysisEngine : MetaMorpheusEngine
    {
        private List<PeptideSpectralMatch> AllPsms;
        private readonly int MassDiffAcceptorNumNotches;
        private readonly bool UseDeltaScore;
        private readonly double ScoreCutoff;
        private readonly string AnalysisType;

        public FdrAnalysisEngine(List<PeptideSpectralMatch> psms, int massDiffAcceptorNumNotches, CommonParameters commonParameters, List<string> nestedIds, string analysisType = "PSM") : base(commonParameters, nestedIds)
        {
            AllPsms = psms;
            MassDiffAcceptorNumNotches = massDiffAcceptorNumNotches;
            UseDeltaScore = commonParameters.UseDeltaScore;
            ScoreCutoff = commonParameters.ScoreCutoff;
            AnalysisType = analysisType;
        }

        protected override MetaMorpheusEngineResults RunSpecific()
        {
            FdrAnalysisResults myAnalysisResults = new FdrAnalysisResults(this, AnalysisType);

            Status("Running FDR analysis...");
            DoFalseDiscoveryRateAnalysis(myAnalysisResults);

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

                //determine if Score or DeltaScore performs better
                if (UseDeltaScore)
                {
                    const double qValueCutoff = 0.01; //optimize to get the most PSMs at a 1% FDR

                    List<PeptideSpectralMatch> scoreSorted = psms.OrderByDescending(b => b.Score).ThenBy(b => b.PeptideMonisotopicMass.HasValue ? Math.Abs(b.ScanPrecursorMass - b.PeptideMonisotopicMass.Value) : double.MaxValue).GroupBy(b => new Tuple<string, int, double?>(b.FullFilePath, b.ScanNumber, b.PeptideMonisotopicMass)).Select(b => b.First()).ToList();
                    int ScorePSMs = GetNumPSMsAtqValueCutoff(scoreSorted, qValueCutoff);
                    scoreSorted = psms.OrderByDescending(b => b.DeltaScore).ThenBy(b => b.PeptideMonisotopicMass.HasValue ? Math.Abs(b.ScanPrecursorMass - b.PeptideMonisotopicMass.Value) : double.MaxValue).GroupBy(b => new Tuple<string, int, double?>(b.FullFilePath, b.ScanNumber, b.PeptideMonisotopicMass)).Select(b => b.First()).ToList();
                    int DeltaScorePSMs = GetNumPSMsAtqValueCutoff(scoreSorted, qValueCutoff);

                    //sort by best method
                    myAnalysisResults.DeltaScoreImprovement = DeltaScorePSMs > ScorePSMs;
                    psms = myAnalysisResults.DeltaScoreImprovement ?
                        psms.OrderByDescending(b => b.DeltaScore).ThenBy(b => b.PeptideMonisotopicMass.HasValue ? Math.Abs(b.ScanPrecursorMass - b.PeptideMonisotopicMass.Value) : double.MaxValue).ToList() :
                        psms.OrderByDescending(b => b.Score).ThenBy(b => b.PeptideMonisotopicMass.HasValue ? Math.Abs(b.ScanPrecursorMass - b.PeptideMonisotopicMass.Value) : double.MaxValue).ToList();
                }
                else //sort by score
                {
                    psms = psms.OrderByDescending(b => b.Score).ThenBy(b => b.PeptideMonisotopicMass.HasValue ? Math.Abs(b.ScanPrecursorMass - b.PeptideMonisotopicMass.Value) : double.MaxValue).ToList();
                }

                double cumulativeTarget = 0;
                double cumulativeDecoy = 0;

                //set up arrays for local FDRs
                double[] cumulativeTargetPerNotch = new double[MassDiffAcceptorNumNotches + 1];
                double[] cumulativeDecoyPerNotch = new double[MassDiffAcceptorNumNotches + 1];

                //Assign FDR values to PSMs
                for (int i = 0; i < psms.Count; i++)
                {
                    // Stop if canceled
                    if (GlobalVariables.StopLoops) { break; }

                    PeptideSpectralMatch psm = psms[i];
                    int notch = psm.Notch ?? MassDiffAcceptorNumNotches;
                    if (psm.IsDecoy)
                    {
                        // the PSM can be ambiguous between a target and a decoy sequence
                        // in that case, count it as the fraction of decoy hits
                        // e.g. if the PSM matched to 1 target and 2 decoys, it counts as 2/3 decoy
                        double decoyHits = 0;
                        double totalHits = 0;
                        var hits = psm.BestMatchingPeptides.GroupBy(p => p.Peptide.FullSequence);
                        foreach (var hit in hits)
                        {
                            if (hit.First().Peptide.Protein.IsDecoy)
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

                    double qValue = Math.Min(1, cumulativeDecoy / cumulativeTarget);
                    double qValueNotch = Math.Min(1, cumulativeDecoyPerNotch[notch] / cumulativeTargetPerNotch[notch]);

                    double pep = psm.FdrInfo == null ? double.NaN : psm.FdrInfo.PEP;
                    double pepQValue = psm.FdrInfo == null ? double.NaN : psm.FdrInfo.PEP_QValue;

                    psm.SetFdrValues(cumulativeTarget, cumulativeDecoy, qValue, cumulativeTargetPerNotch[notch], cumulativeDecoyPerNotch[notch], qValueNotch, pep, pepQValue);
                }

                // set q-value thresholds such that a lower scoring PSM can't have
                // a higher confidence than a higher scoring PSM
                //Populate min qValues
                double qValueThreshold = 1.0;
                double[] qValueNotchThreshold = new double[MassDiffAcceptorNumNotches + 1];
                for (int i = 0; i < qValueNotchThreshold.Length; i++)
                {
                    qValueNotchThreshold[i] = 1.0;
                }

                for (int i = psms.Count - 1; i >= 0; i--)
                {
                    PeptideSpectralMatch psm = psms[i];

                    // threshold q-values
                    if (psm.FdrInfo.QValue > qValueThreshold)
                    {
                        psm.FdrInfo.QValue = qValueThreshold;
                    }
                    else if (psm.FdrInfo.QValue < qValueThreshold)
                    {
                        qValueThreshold = psm.FdrInfo.QValue;
                    }

                    // threshold notch q-values
                    int notch = psm.Notch ?? MassDiffAcceptorNumNotches;
                    if (psm.FdrInfo.QValueNotch > qValueNotchThreshold[notch])
                    {
                        psm.FdrInfo.QValueNotch = qValueNotchThreshold[notch];
                    }
                    else if (psm.FdrInfo.QValueNotch < qValueNotchThreshold[notch])
                    {
                        qValueNotchThreshold[notch] = psm.FdrInfo.QValueNotch;
                    }
                }
            }

            if (AnalysisType == "PSM")
            {
                CountPsm();
                //Need some reasonable number of PSMs to train on to get a reasonable estimation of the PEP
                if (AllPsms.Count > 100)
                {
                    myAnalysisResults.BinarySearchTreeMetrics = PEP_Analysis.ComputePEPValuesForAllPSMsGeneric(AllPsms);
                    Compute_PEPValue_Based_QValue(AllPsms);
                }
            }

            if (AnalysisType == "Peptide")
            {
                Compute_PEPValue_Based_QValue(AllPsms);
            }
        }

        public static void Compute_PEPValue_Based_QValue(List<PeptideSpectralMatch> psms)
        {
            double[] allPEPValues = psms.Select(p => p.FdrInfo.PEP).ToArray();
            int[] psmsArrayIndicies = Enumerable.Range(0, psms.Count).ToArray();
            Array.Sort(allPEPValues, psmsArrayIndicies);//sort the second thing by the first

            double runningSum = 0;
            for (int i = 0; i < allPEPValues.Length; i++)
            {
                runningSum += allPEPValues[i];
                double qValue = runningSum / (i + 1);
                psms[psmsArrayIndicies[i]].FdrInfo.PEP_QValue = Math.Round(qValue, 6);
            }
        }

        private static int GetNumPSMsAtqValueCutoff(List<PeptideSpectralMatch> psms, double qValueCutoff)
        {
            int cumulative_target = 0;
            int cumulative_decoy = 0;
            foreach (PeptideSpectralMatch psm in psms)
            {
                if (psm.IsDecoy)
                {
                    cumulative_decoy++;
                    if ((double)cumulative_decoy / cumulative_target >= qValueCutoff)
                    {
                        return cumulative_target;
                    }
                }
                else
                    cumulative_target++;
            }
            return cumulative_target;
        }

        public void CountPsm()
        {
            // exclude ambiguous psms and has a fdr cutoff = 0.01
            var allUnambiguousPsms = AllPsms.Where(psm => psm.FullSequence != null);

            var unambiguousPsmsLessThanOnePercentFdr = allUnambiguousPsms.Where(psm =>
                psm.FdrInfo.QValue <= 0.01
                && psm.FdrInfo.QValueNotch <= 0.01)
                .GroupBy(p => p.FullSequence);

            Dictionary<string, int> sequenceToPsmCount = new Dictionary<string, int>();

            foreach (var sequenceGroup in unambiguousPsmsLessThanOnePercentFdr)
            {
                if (!sequenceToPsmCount.ContainsKey(sequenceGroup.First().FullSequence))
                {
                    sequenceToPsmCount.Add(sequenceGroup.First().FullSequence, sequenceGroup.Count());
                }
            }

            foreach (PeptideSpectralMatch psm in allUnambiguousPsms)
            {
                if (sequenceToPsmCount.ContainsKey(psm.FullSequence))
                {
                    psm.PsmCount = sequenceToPsmCount[psm.FullSequence];
                }
            }
        }
    }
}