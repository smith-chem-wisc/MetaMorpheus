using MathNet.Numerics;
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
        private readonly bool CalculateEValue;
        private readonly double ScoreCutoff;

        public FdrAnalysisEngine(List<PeptideSpectralMatch> psms, int massDiffAcceptorNumNotches, CommonParameters commonParameters, List<string> nestedIds) : base(commonParameters, nestedIds)
        {
            AllPsms = psms;
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

            myAnalysisResults.PsmsWithin1PercentFdr = AllPsms.Count(b => b.FdrInfo.QValue < 0.01);

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
                
                // generate the null distribution for e-value calculations
                double globalMeanScore = 0;
                int globalMeanCount = 0;

                if (CalculateEValue && psms.Any())
                {
                    List<double> combinedScores = new List<double>();

                    foreach (PeptideSpectralMatch psm in psms)
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
                        globalMeanCount = (int)((double)combinedScores.Count / psms.Count);
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
                foreach (PeptideSpectralMatch psm in psms)
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

                //set up arrays for local FDRs
                int[] cumulativeTargetPerNotch = new int[MassDiffAcceptorNumNotches + 1];
                int[] cumulativeDecoyPerNotch = new int[MassDiffAcceptorNumNotches + 1];

                //Assign FDR values to PSMs
                for (int i = 0; i < psms.Count; i++)
                {
                    // Stop if canceled
                    if (GlobalVariables.StopLoops) { break; }

                    var psm = psms[i];
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

                // set q-value thresholds such that a lower scoring PSM can't have 
                // a higher confidence than a higher scoring PSM
                double qThrehold = 0;
                double[] notchQThresholds = new double[MassDiffAcceptorNumNotches + 1];
                for(int i = 0; i < psms.Count; i++)
                {
                    PeptideSpectralMatch psm = psms[i];

                    if (psm.FdrInfo.QValue > qThrehold)
                    {
                        qThrehold = psm.FdrInfo.QValue;
                    }
                    else
                    {
                        psm.FdrInfo.QValue = qThrehold;
                    }

                    int notch = psm.Notch ?? MassDiffAcceptorNumNotches;
                    if (psm.FdrInfo.QValueNotch > notchQThresholds[notch])
                    {
                        notchQThresholds[notch] = psm.FdrInfo.QValueNotch;
                    }
                    else
                    {
                        psm.FdrInfo.QValueNotch = notchQThresholds[notch];
                    }
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
