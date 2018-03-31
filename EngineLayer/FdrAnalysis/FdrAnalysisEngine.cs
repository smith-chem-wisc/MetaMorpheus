using MathNet.Numerics;
using System;
using System.Collections.Generic;
using System.Linq;

namespace EngineLayer.FdrAnalysis
{
    public class FdrAnalysisEngine : MetaMorpheusEngine
    {
        #region Private Fields

        private List<PeptideSpectralMatch> psms;
        private readonly int massDiffAcceptorNumNotches;
        private readonly bool calculateDeltaScore;
        private readonly bool calculateEValue;
        private readonly double scoreCutoff;

        #endregion Private Fields

        #region Public Constructors

        public FdrAnalysisEngine(List<PeptideSpectralMatch> psms, int massDiffAcceptorNumNotches, ICommonParameters commonParameters, List<string> nestedIds) : base(nestedIds)
        {
            this.psms = psms;
            this.massDiffAcceptorNumNotches = massDiffAcceptorNumNotches;
            this.calculateDeltaScore = commonParameters.CalculateDeltaScore;
            this.scoreCutoff = commonParameters.ScoreCutoff;
            this.calculateEValue = commonParameters.CalculateEValue;
        }

        #endregion Public Constructors

        #region Protected Methods

        protected override MetaMorpheusEngineResults RunSpecific()
        {
            FdrAnalysisResults myAnalysisResults = new FdrAnalysisResults(this);

            Status("Running FDR analysis...");
            DoFalseDiscoveryRateAnalysis(myAnalysisResults);

            myAnalysisResults.PsmsWithin1PercentFdr = psms.Count(b => b.FdrInfo.QValue < 0.01);

            return myAnalysisResults;
        }

        #endregion Protected Methods

        #region Private Methods

        private void DoFalseDiscoveryRateAnalysis(FdrAnalysisResults myAnalysisResults)
        {
            double globalMeanScore = 0;
            double globalMeanCount = 0;

            if (calculateEValue)
            {
                List<Int64> combinedAllScores = new List<Int64>();
                Int64 totalCount = 0;
                Int64 totalSum = 0;

                foreach (PeptideSpectralMatch psm in psms)
                {
                    psm.AllScores[psm.AllScores.Count - 1]--; //remove top scoring peptide
                    while (combinedAllScores.Count < psm.AllScores.Count) //expand array if neccessary
                    {
                        combinedAllScores.Add(0);
                    }
                    for (int score = 0; score < psm.AllScores.Count; score++) //add scores
                    {
                        combinedAllScores[score] += psm.AllScores[score];
                    }
                }
                for (int i = 0; i < combinedAllScores.Count; i++)
                {
                    totalCount += combinedAllScores[i];
                    totalSum += combinedAllScores[i] * i;
                }
                globalMeanScore = totalCount != 0 ? totalSum * 1.0d / totalCount : 0;
                globalMeanCount = totalCount * 1.0d / psms.Count;
            }

            int cumulative_target = 0;
            int cumulative_decoy = 0;

            //determine if Score or DeltaScore performs better
            if (calculateDeltaScore)
            {
                //Calculate delta scores for the psms
                foreach (PeptideSpectralMatch psm in psms)
                {
                    if (psm != null)
                    {
                        psm.CalculateDeltaScore(scoreCutoff);
                    }
                }
                const double qValueCutoff = 0.01; //optimize to get the most PSMs at a 1% FDR

                List<PeptideSpectralMatch> scoreSorted = psms.Where(b => b != null).OrderByDescending(b => b.Score).ThenBy(b => b.PeptideMonisotopicMass.HasValue ? Math.Abs(b.ScanPrecursorMass - b.PeptideMonisotopicMass.Value) : double.MaxValue).GroupBy(b => new Tuple<string, int, double?>(b.FullFilePath, b.ScanNumber, b.PeptideMonisotopicMass)).Select(b => b.First()).ToList();
                int ScorePSMs = GetNumPSMsAtqValueCutoff(scoreSorted, qValueCutoff);
                scoreSorted = psms.Where(b => b != null).OrderByDescending(b => b.DeltaScore).ThenBy(b => b.PeptideMonisotopicMass.HasValue ? Math.Abs(b.ScanPrecursorMass - b.PeptideMonisotopicMass.Value) : double.MaxValue).GroupBy(b => new Tuple<string, int, double?>(b.FullFilePath, b.ScanNumber, b.PeptideMonisotopicMass)).Select(b => b.First()).ToList();
                int DeltaScorePSMs = GetNumPSMsAtqValueCutoff(scoreSorted, qValueCutoff);

                //sort by best method BUT DON'T GROUP BECAUSE WE NEED THAT FOR LOCALIZATION
                psms = (DeltaScorePSMs > ScorePSMs) ?
                    psms.Where(b => b != null).OrderByDescending(b => b.DeltaScore).ThenBy(b => b.PeptideMonisotopicMass.HasValue ? Math.Abs(b.ScanPrecursorMass - b.PeptideMonisotopicMass.Value) : double.MaxValue).ToList() :
                    psms.Where(b => b != null).OrderByDescending(b => b.Score).ThenBy(b => b.PeptideMonisotopicMass.HasValue ? Math.Abs(b.ScanPrecursorMass - b.PeptideMonisotopicMass.Value) : double.MaxValue).ToList();
                myAnalysisResults.DeltaScoreImprovement = DeltaScorePSMs > ScorePSMs;
            }
            else //sort by score
            {
                psms = psms.Where(b => b != null).OrderByDescending(b => b.Score).ThenBy(b => b.PeptideMonisotopicMass.HasValue ? Math.Abs(b.ScanPrecursorMass - b.PeptideMonisotopicMass.Value) : double.MaxValue).ToList();
            }

            //set up arrays for local FDRs
            int[] cumulative_target_per_notch = new int[massDiffAcceptorNumNotches + 1];
            int[] cumulative_decoy_per_notch = new int[massDiffAcceptorNumNotches + 1];

            //Assign FDR values to PSMs
            for (int i = 0; i < psms.Count; i++)
            {
                var psm = psms[i];
                var isDecoy = psm.IsDecoy;
                int notch = psm.Notch ?? massDiffAcceptorNumNotches;
                if (isDecoy)
                {
                    cumulative_decoy++;
                    cumulative_decoy_per_notch[notch]++;
                }
                else
                {
                    cumulative_target++;
                    cumulative_target_per_notch[notch]++;
                }

                double temp_q_value = (double)cumulative_decoy / cumulative_target;
                double temp_q_value_for_notch = (double)cumulative_decoy_per_notch[notch] / cumulative_target_per_notch[notch];

                double maximumLikelihood = 0;
                decimal eValue = 0;
                double eScore = 0;
                if (calculateEValue)
                {
                    (int sum, int count) sumAndCount = GetSumAndCount(psm.AllScores);
                    maximumLikelihood = GetMaximumLikelihood(sumAndCount.sum, sumAndCount.count);
                    eValue = GetEValue(psm.Score, sumAndCount.count, maximumLikelihood, globalMeanCount, globalMeanScore);
                    eScore = GetEScore(eValue);
                }
                psm.SetFdrValues(cumulative_target, cumulative_decoy, temp_q_value, cumulative_target_per_notch[notch], cumulative_decoy_per_notch[notch], temp_q_value_for_notch, maximumLikelihood, eValue, eScore, calculateEValue);
            }

            //Populate min qValues
            double min_q_value = double.PositiveInfinity;
            double[] min_q_value_notch = new double[massDiffAcceptorNumNotches + 1];
            for (int i = 0; i < massDiffAcceptorNumNotches + 1; i++)
            {
                min_q_value_notch[i] = double.PositiveInfinity;
            }
            //The idea here is to set previous qValues as thresholds, 
            //such that a lower scoring PSM can't have a higher confidence than a higher scoring PSM
            for (int i = psms.Count - 1; i >= 0; i--)
            {
                PeptideSpectralMatch psm = psms[i];
                if (psm.FdrInfo.QValue > min_q_value)
                {
                    psm.FdrInfo.QValue = min_q_value;
                }
                else if (psm.FdrInfo.QValue < min_q_value)
                {
                    min_q_value = psm.FdrInfo.QValue;
                }
                int notch = psm.Notch ?? massDiffAcceptorNumNotches;
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

        private static (int sum, int count) GetSumAndCount(List<int> allScores)
        {
            if (!allScores.Any(x => x != 0))
            {
                return (0, 0);
            }
            else
            {
                int count = 0;
                int sum = 0;
                for (int i = 0; i < allScores.Count; i++)
                {
                    count += allScores[i];
                    sum += i * allScores[i];
                }
                return (sum, count);
            }
        }

        private static double GetMaximumLikelihood(int sum, int count)
        {
            return count == 0 ? 0 : (1.0d / count) * sum;
        }

        private static decimal GetEValue(double dScore, double count, double maximumLikelihood, double globalMeanCount, double globalMeanScore)
        {
            int score = (int)Math.Floor(dScore);
            double preValue; // this is the cumulative distribution for the poisson at each score up to but not including the score of the winner. This is the probability that the winner has of getting that score at random by matching against a SINGLE spectrum

            if (score == 0)
            {
                preValue = 1;
            }
            else if (maximumLikelihood == 0)
            {
                preValue = SpecialFunctions.GammaLowerRegularized(globalMeanScore, (score - 1));
            }
            else
            {
                preValue = SpecialFunctions.GammaLowerRegularized(maximumLikelihood, (score - 1));
            }
            // Now the probability of getting the winner's score goes up for each spectrum match. We multiply the preValue by the number of theoretical spectrum within the tolerance to get this new probability.
            if (count > 0)
            {
                double four = (1 - Math.Pow(preValue, (count)));
                return (Convert.ToDecimal((count) * four));
            }
            else
            {
                double four = (1 - Math.Pow(preValue, (globalMeanCount)));
                return (Convert.ToDecimal((count) * four));
            }
        }

        private static double GetEScore(decimal eValue)
        {
            return (-10d * Math.Log10((double)eValue));
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

        #endregion Private Methods
    }
}