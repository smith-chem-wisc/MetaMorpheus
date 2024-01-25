using EngineLayer.CrosslinkSearch;
using Newtonsoft.Json.Linq;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text.RegularExpressions;
using System.Threading;
using Easy.Common.Extensions;

namespace EngineLayer.FdrAnalysis
{
    public class FdrAnalysisEngine : MetaMorpheusEngine
    {
        private List<PeptideSpectralMatch> AllPsms;
        private readonly int MassDiffAcceptorNumNotches;
        private readonly bool UseDeltaScore;
        private readonly double ScoreCutoff;
        private readonly string AnalysisType;
        private readonly string OutputFolder; // used for storing PEP training models
        private readonly bool DoPEP;

        public FdrAnalysisEngine(List<PeptideSpectralMatch> psms, int massDiffAcceptorNumNotches, CommonParameters commonParameters,
            List<(string fileName, CommonParameters fileSpecificParameters)> fileSpecificParameters, List<string> nestedIds, string analysisType = "PSM", bool doPEP = true, string outputFolder = null) : base(commonParameters, fileSpecificParameters, nestedIds)
        {
            AllPsms = psms;
            MassDiffAcceptorNumNotches = massDiffAcceptorNumNotches;
            UseDeltaScore = commonParameters.UseDeltaScore;
            ScoreCutoff = commonParameters.ScoreCutoff;
            AnalysisType = analysisType;
            this.OutputFolder = outputFolder;
            this.DoPEP = doPEP;
            if (fileSpecificParameters == null) throw new ArgumentNullException("file specific parameters cannot be null");
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

                ComputeTargetAndDecoyCounts(psms);
                QValueTraditional(psms);
                if (psms.Count > 100)
                {
                    QValueInverted(psms);
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
            if (DoPEP)
            {
                Compute_PEPValue(myAnalysisResults);
            }
        }

        private void ComputeTargetAndDecoyCounts(List<PeptideSpectralMatch> psms)
        {
            double cumulativeTarget = 0;
            double cumulativeDecoy = 0;

            double[] cumulativeTargetPerNotch = new double[MassDiffAcceptorNumNotches + 1];
            double[] cumulativeDecoyPerNotch = new double[MassDiffAcceptorNumNotches + 1];

            //Assign FDR values to PSMs
            for (int i = 0; i < psms.Count; i++)
            {
                // Stop if canceled
                if (GlobalVariables.StopLoops) { break; }

                // the PSM can be ambiguous between a target and a decoy sequence
                // in that case, count it as the fraction of decoy hits
                // e.g. if the PSM matched to 1 target and 2 decoys, it counts as 2/3 decoy

                double[] targetPerNotch = new double[MassDiffAcceptorNumNotches + 1];
                double[] decoyPerNotch = new double[MassDiffAcceptorNumNotches + 1];

                foreach (var bestMatchingPeptide in psms[i].BestMatchingPeptides)
                {
                    if (bestMatchingPeptide.Peptide == null)
                    {
                        continue;
                    }

                    if (bestMatchingPeptide.Peptide.Protein.IsDecoy)
                    {
                        decoyPerNotch[bestMatchingPeptide.Notch]++;
                    }
                    else
                    {
                        targetPerNotch[bestMatchingPeptide.Notch]++;
                    }
                }

                cumulativeTarget += targetPerNotch.Sum() / (targetPerNotch.Sum() + decoyPerNotch.Sum());
                cumulativeDecoy += decoyPerNotch.Sum() / (targetPerNotch.Sum() + decoyPerNotch.Sum());

                List<(double, int)> myPairs = new List<(double, int)>();

                for (int j = 0; j < MassDiffAcceptorNumNotches; j++)
                {
                    double denominator = targetPerNotch[j] + decoyPerNotch[j];
                    if (denominator != 0)
                    {
                        cumulativeTargetPerNotch[j] += (targetPerNotch[j]) / denominator;
                        cumulativeDecoyPerNotch[j] += (decoyPerNotch[j]) / denominator;
                        myPairs.Add((cumulativeTargetPerNotch[j],j));
                    }
                }

                double worstQvalue = myPairs.Select(k => k.Item1).Max();
                int worstNotch = myPairs.Where(k => k.Item1 == worstQvalue).Select(k => k.Item2).First();

                psms[i].SetFdrTargetAndDecoyCounts(cumulativeTarget,cumulativeDecoy, cumulativeTargetPerNotch[worstNotch], cumulativeDecoyPerNotch[worstNotch]);
            }
        }

        private void QValueInverted(List<PeptideSpectralMatch> psms)
        {
            psms.Reverse();
            double qValue = 0;
            double qValueNotch = 0;

            //Assign FDR values to PSMs
            for (int i = 0; i < psms.Count; i++)
            {
                // Stop if canceled
                if (GlobalVariables.StopLoops) { break; }

                qValue = Math.Min(qValue, (psms[i].FdrInfo.CumulativeDecoy + 1) / psms[i].FdrInfo.CumulativeTarget);
                qValueNotch = Math.Min(qValueNotch, (psms[i].FdrInfo.CumulativeDecoyNotch + 1) / psms[i].FdrInfo.CumulativeTargetNotch);

                double pep = psms[i].FdrInfo == null ? double.NaN : psms[i].FdrInfo.PEP;
                double pepQValue = psms[i].FdrInfo == null ? double.NaN : psms[i].FdrInfo.PEP_QValue;

                psms[i].SetQandPEPvalues(qValue, qValueNotch, pep, pepQValue);
            }
            psms.Reverse(); //we inverted the psms for this calculation. now we need to put them back into the original order
        }

        private void QValueTraditional(List<PeptideSpectralMatch> psms)
        {
            double qValue = 0;
            double qValueNotch = 0;

            //Assign FDR values to PSMs
            for (int i = 0; i < psms.Count; i++)
            {
                // Stop if canceled
                if (GlobalVariables.StopLoops) { break; }

                qValue = Math.Max(qValue, psms[i].FdrInfo.CumulativeDecoy / psms[i].FdrInfo.CumulativeTarget);
                qValueNotch = Math.Max(qValueNotch, psms[i].FdrInfo.CumulativeDecoyNotch / psms[i].FdrInfo.CumulativeTargetNotch);

                double pep = psms[i].FdrInfo == null ? double.NaN : psms[i].FdrInfo.PEP;
                double pepQValue = psms[i].FdrInfo == null ? double.NaN : psms[i].FdrInfo.PEP_QValue;

                psms[i].SetQandPEPvalues(qValue, qValueNotch, pep, pepQValue);
            }
        }

        public void Compute_PEPValue(FdrAnalysisResults myAnalysisResults)
        {
            if (AnalysisType == "PSM")
            {
                CountPsm();
                //Need some reasonable number of PSMs to train on to get a reasonable estimation of the PEP
                if (AllPsms.Count > 100)
                {
                    string searchType = "standard";
                    if (AllPsms[0].DigestionParams.Protease.Name == "top-down")
                    {
                        searchType = "top-down";
                    }

                    myAnalysisResults.BinarySearchTreeMetrics = PEP_Analysis_Cross_Validation.ComputePEPValuesForAllPSMsGeneric(AllPsms, searchType, this.FileSpecificParameters, this.OutputFolder);

                    Compute_PEPValue_Based_QValue(AllPsms);
                }
                CountPsm(); // recounting Psm's after PEP based disambiguation
            }

            if (AnalysisType == "Peptide")
            {
                Compute_PEPValue_Based_QValue(AllPsms);
            }

            if (AnalysisType == "crosslink" && AllPsms.Count > 100)
            {
                myAnalysisResults.BinarySearchTreeMetrics = PEP_Analysis_Cross_Validation.ComputePEPValuesForAllPSMsGeneric(AllPsms, "crosslink", this.FileSpecificParameters, this.OutputFolder);
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