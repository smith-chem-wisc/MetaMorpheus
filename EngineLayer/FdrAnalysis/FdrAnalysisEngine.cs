using MathNet.Numerics;
using System;
using System.Collections.Generic;
using System.Linq;

namespace EngineLayer.FdrAnalysis
{
    public class FdrAnalysisEngine : MetaMorpheusEngine
    {
        private List<PeptideSpectralMatch> Psms;
        private readonly int MassDiffAcceptorNumNotches;
        private readonly bool UseDeltaScore;
        private readonly bool CalculateEValue;
        private readonly double ScoreCutoff;

        public FdrAnalysisEngine(List<PeptideSpectralMatch> psms, int massDiffAcceptorNumNotches, CommonParameters commonParameters, List<string> nestedIds) : base(commonParameters, nestedIds)
        {
            Psms = psms;
            MassDiffAcceptorNumNotches = massDiffAcceptorNumNotches;
            UseDeltaScore = commonParameters.UseDeltaScore;
            ScoreCutoff = commonParameters.ScoreCutoff;
            CalculateEValue = commonParameters.CalculateEValue;
        }

        protected override MetaMorpheusEngineResults RunSpecific()
        {
            FdrAnalysisResults myAnalysisResults = new FdrAnalysisResults(this);

            Status("Running FDR analysis...");
            DoFalseDiscoveryRateAnalysis(myAnalysisResults);

            myAnalysisResults.PsmsWithin1PercentFdr = Psms.Count(b => b.FdrInfo.QValue < 0.01);

            return myAnalysisResults;
        }

        private void DoFalseDiscoveryRateAnalysis(FdrAnalysisResults myAnalysisResults)
        {
            // Stop if canceled
            if (GlobalVariables.StopLoops) { return; }

            // generate the null distribution for e-value calculations
            double globalMeanScore = 0;
            int globalMeanCount = 0;

            if (CalculateEValue && Psms.Any())
            {
                List<double> combinedScores = new List<double>();

                foreach (PeptideSpectralMatch psm in Psms)
                {
                    psm.AllScores.Sort();
                    combinedScores.AddRange(psm.AllScores);

                    //remove top scoring peptide
                    if (combinedScores.Any())
                    {
                        combinedScores.RemoveAt(combinedScores.Count - 1);
                    }
                }

                if (combinedScores.Any())
                {
                    globalMeanScore = combinedScores.Average();
                    globalMeanCount = (int)((double)combinedScores.Count / Psms.Count);
                }
                else
                {
                    // should be a very rare case... if there are PSMs but each PSM only has one hit
                    globalMeanScore = 0;
                    globalMeanCount = 0;
                }
            }

            int cumulativeTarget = 0;
            int cumulativeDecoy = 0;

            //Calculate delta scores for the psms (regardless of if we are using them)
            foreach (PeptideSpectralMatch psm in Psms)
            {
                if (psm != null)
                {
                    psm.CalculateDeltaScore(ScoreCutoff);
                }
            }

            //determine if Score or DeltaScore performs better
            if (UseDeltaScore)
            {
                const double qValueCutoff = 0.01; //optimize to get the most PSMs at a 1% FDR

                List<PeptideSpectralMatch> scoreSorted = Psms.OrderByDescending(b => b.Score).ThenBy(b => b.PeptideMonisotopicMass.HasValue ? Math.Abs(b.ScanPrecursorMass - b.PeptideMonisotopicMass.Value) : double.MaxValue).GroupBy(b => new Tuple<string, int, double?>(b.FullFilePath, b.ScanNumber, b.PeptideMonisotopicMass)).Select(b => b.First()).ToList();
                int ScorePSMs = GetNumPSMsAtqValueCutoff(scoreSorted, qValueCutoff);
                scoreSorted = Psms.OrderByDescending(b => b.DeltaScore).ThenBy(b => b.PeptideMonisotopicMass.HasValue ? Math.Abs(b.ScanPrecursorMass - b.PeptideMonisotopicMass.Value) : double.MaxValue).GroupBy(b => new Tuple<string, int, double?>(b.FullFilePath, b.ScanNumber, b.PeptideMonisotopicMass)).Select(b => b.First()).ToList();
                int DeltaScorePSMs = GetNumPSMsAtqValueCutoff(scoreSorted, qValueCutoff);

                //sort by best method BUT DON'T GROUP BECAUSE WE NEED THAT FOR LOCALIZATION
                Psms = (DeltaScorePSMs > ScorePSMs) ?
                    Psms.OrderByDescending(b => b.DeltaScore).ThenBy(b => b.PeptideMonisotopicMass.HasValue ? Math.Abs(b.ScanPrecursorMass - b.PeptideMonisotopicMass.Value) : double.MaxValue).ToList() :
                    Psms.OrderByDescending(b => b.Score).ThenBy(b => b.PeptideMonisotopicMass.HasValue ? Math.Abs(b.ScanPrecursorMass - b.PeptideMonisotopicMass.Value) : double.MaxValue).ToList();
                myAnalysisResults.DeltaScoreImprovement = DeltaScorePSMs > ScorePSMs;
            }
            else //sort by score
            {
                Psms = Psms.OrderByDescending(b => b.Score).ThenBy(b => b.PeptideMonisotopicMass.HasValue ? Math.Abs(b.ScanPrecursorMass - b.PeptideMonisotopicMass.Value) : double.MaxValue).ToList();
            }

            //set up arrays for local FDRs
            int[] cumulativeTargetPerNotch = new int[MassDiffAcceptorNumNotches + 1];
            int[] cumulativeDecoyPerNotch = new int[MassDiffAcceptorNumNotches + 1];

            //Assign FDR values to PSMs
            for (int i = 0; i < Psms.Count; i++)
            {
                // Stop if canceled
                if (GlobalVariables.StopLoops) { break; }

                var psm = Psms[i];
                int notch = psm.Notch ?? MassDiffAcceptorNumNotches;
                if (psm.IsDecoy)
                {
                    cumulativeDecoy++;
                    cumulativeDecoyPerNotch[notch]++;
                }
                else
                {
                    cumulativeTarget++;
                    cumulativeTargetPerNotch[notch]++;
                }

                double qValue = (double)cumulativeDecoy / cumulativeTarget;
                double qValueNotch = (double)cumulativeDecoyPerNotch[notch] / cumulativeTargetPerNotch[notch];

                double maximumLikelihood = 0;
                double eValue = 0;
                double eScore = 0;
                if (CalculateEValue)
                {
                    eValue = GetEValue(psm, globalMeanCount, globalMeanScore, out maximumLikelihood);
                    eScore = -Math.Log(eValue, 10);
                }

                if (qValue > 1)
                {
                    qValue = 1;
                }
                if (qValueNotch > 1)
                {
                    qValueNotch = 1;
                }

                psm.SetFdrValues(cumulativeTarget, cumulativeDecoy, qValue, cumulativeTargetPerNotch[notch], cumulativeDecoyPerNotch[notch], qValueNotch, maximumLikelihood, eValue, eScore, CalculateEValue);
            }

            //Populate min qValues
            double min_q_value = double.PositiveInfinity;
            double[] min_q_value_notch = new double[MassDiffAcceptorNumNotches + 1];
            for (int i = 0; i < MassDiffAcceptorNumNotches + 1; i++)
            {
                min_q_value_notch[i] = double.PositiveInfinity;
            }

            //The idea here is to set previous qValues as thresholds,
            //such that a lower scoring PSM can't have a higher confidence than a higher scoring PSM
            for (int i = Psms.Count - 1; i >= 0; i--)
            {
                PeptideSpectralMatch psm = Psms[i];
                if (psm.FdrInfo.QValue > min_q_value)
                {
                    psm.FdrInfo.QValue = min_q_value;
                }
                else if (psm.FdrInfo.QValue < min_q_value)
                {
                    min_q_value = psm.FdrInfo.QValue;
                }
                int notch = psm.Notch ?? MassDiffAcceptorNumNotches;
                if (psm.FdrInfo.QValueNotch > min_q_value_notch[notch])
                {
                    psm.FdrInfo.QValueNotch = min_q_value_notch[notch];
                }
                else if (psm.FdrInfo.QValueNotch < min_q_value_notch[notch])
                {
                    min_q_value_notch[notch] = psm.FdrInfo.QValueNotch;
                }
            }
        }

        private static double GetEValue(PeptideSpectralMatch psm, int globalMeanCount, double globalMeanScore, out double maximumLikelihood)
        {
            // get all of the PSM's scores for all hits, sort them, then remove the last value (the best score)
            List<double> scoresWithoutBestHit = new List<double>();
            scoresWithoutBestHit.AddRange(psm.AllScores);
            scoresWithoutBestHit.Sort();

            if (scoresWithoutBestHit.Any())
            {
                scoresWithoutBestHit.RemoveAt(scoresWithoutBestHit.Count - 1);
            }

            // this is the "default" case for when there are no scores except the best hit
            // it uses a global mean score (all scores across all PSMs) to generate the null Poisson distribution
            // this will be overriden by the next few lines if there are enough scores in this PSM to estimate a null distribution
            double preValue = SpecialFunctions.GammaLowerRegularized(globalMeanScore, psm.Score);
            maximumLikelihood = globalMeanScore;

            // calculate single-spectrum evalue if there are enough hits besides the best scoring peptide
            if (psm.Score == 0)
            {
                preValue = 1;
                maximumLikelihood = 0;
            }
            else if (scoresWithoutBestHit.Any())
            {
                maximumLikelihood = scoresWithoutBestHit.Average();

                // this is the cumulative distribution for the poisson at each score up to but not including the score of the winner.
                // This is the probability that the winner has of getting that score at random by matching against a SINGLE spectrum
                if (maximumLikelihood > 0)
                {
                    preValue = SpecialFunctions.GammaLowerRegularized(maximumLikelihood, psm.Score);
                }
            }

            // Now the probability of getting the winner's score goes up for each spectrum match.
            // We multiply the preValue by the number of theoretical spectrum within the tolerance to get this new probability.
            int count = scoresWithoutBestHit.Count;
            if (count == 0)
            {
                count = globalMeanCount;
            }

            double probabilityOfScore = 1 - Math.Pow(preValue, count);
            return count * probabilityOfScore;
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
    }
}