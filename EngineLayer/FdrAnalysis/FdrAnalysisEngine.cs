using MathNet.Numerics;
using System;
using System.Collections.Generic;
using System.Linq;

namespace EngineLayer.FdrAnalysis
{
    public class FdrAnalysisEngine : MetaMorpheusEngine
    {
        #region Private Fields

        private readonly IEnumerable<PeptideSpectralMatch> newPsms;
        private readonly int massDiffAcceptorNumNotches;
        private readonly bool calculateEValue;

        #endregion Private Fields

        #region Public Constructors

        public FdrAnalysisEngine(List<PeptideSpectralMatch> newPsms, int massDiffAcceptorNumNotches, bool calculateEValue, List<string> nestedIds) : base(nestedIds)
        {
            this.newPsms = newPsms;
            this.massDiffAcceptorNumNotches = massDiffAcceptorNumNotches;
            this.calculateEValue = calculateEValue;
        }

        #endregion Public Constructors

        #region Protected Methods

        protected override MetaMorpheusEngineResults RunSpecific()
        {
            FdrAnalysisResults myAnalysisResults = new FdrAnalysisResults(this);

            Status("Running FDR analysis...");
            DoFalseDiscoveryRateAnalysis(newPsms, massDiffAcceptorNumNotches, calculateEValue);

            myAnalysisResults.PsmsWithin1PercentFdr = newPsms.Count(b => b.FdrInfo.QValue < 0.01);

            return myAnalysisResults;
        }

        #endregion Protected Methods

        #region Private Methods

        private static List<PeptideSpectralMatch> DoFalseDiscoveryRateAnalysis(IEnumerable<PeptideSpectralMatch> items, int massDiffAcceptorNumNotches, bool calculateEValue)
        {
            double globalMeanScore = 0;
            double globalMeanCount = 0;

            if (calculateEValue)
            {
                List<Int64> combinedAllScores = new List<Int64>();
                Int64 totalCount = 0;
                Int64 totalSum = 0;

                foreach (PeptideSpectralMatch item in items)
                {
                    item.AllScores[item.AllScores.Count - 1]--; //remove top scoring peptide
                    while (combinedAllScores.Count < item.AllScores.Count) //expand array if neccessary
                        combinedAllScores.Add(0);
                    for (int score = 0; score < item.AllScores.Count; score++) //add scores
                        combinedAllScores[score] += item.AllScores[score];
                }
                for (int i = 0; i < combinedAllScores.Count; i++)
                {
                    totalCount += combinedAllScores[i];
                    totalSum += combinedAllScores[i] * i;
                }
                globalMeanScore = totalCount != 0 ? totalSum * 1.0d / totalCount : 0;
                globalMeanCount = totalCount * 1.0d / items.Count();
            }
            var ids = new List<PeptideSpectralMatch>();
            foreach (PeptideSpectralMatch item in items)
                ids.Add(item);

            int cumulative_target = 0;
            int cumulative_decoy = 0;

            int[] cumulative_target_per_notch = new int[massDiffAcceptorNumNotches + 1];
            int[] cumulative_decoy_per_notch = new int[massDiffAcceptorNumNotches + 1];

            for (int i = 0; i < ids.Count; i++)
            {
                var item = ids[i];
                var isDecoy = item.IsDecoy;
                int notch = item.Notch ?? massDiffAcceptorNumNotches;
                if (isDecoy)
                    cumulative_decoy++;
                else
                    cumulative_target++;

                if (isDecoy)
                    cumulative_decoy_per_notch[notch]++;
                else
                    cumulative_target_per_notch[notch]++;

                double temp_q_value = (double)cumulative_decoy / cumulative_target;
                double temp_q_value_for_notch = (double)cumulative_decoy_per_notch[notch] / cumulative_target_per_notch[notch];

                double maximumLikelihood = 0;
                decimal eValue = 0;
                double eScore = 0;
                if (calculateEValue)
                {
                    (int sum, int count) sumAndCount = GetSumAndCount(item.AllScores);
                    maximumLikelihood = GetMaximumLikelihood(sumAndCount.sum, sumAndCount.count);
                    eValue = GetEValue(item.Score, sumAndCount.count, maximumLikelihood, globalMeanCount, globalMeanScore);
                    eScore = GetEScore(eValue);
                }
                item.SetFdrValues(cumulative_target, cumulative_decoy, temp_q_value, cumulative_target_per_notch[notch], cumulative_decoy_per_notch[notch], temp_q_value_for_notch, maximumLikelihood, eValue, eScore, calculateEValue);
            }

            double min_q_value = double.PositiveInfinity;
            double[] min_q_value_notch = new double[massDiffAcceptorNumNotches + 1];
            for (int i = 0; i < massDiffAcceptorNumNotches + 1; i++)
                min_q_value_notch[i] = double.PositiveInfinity;

            for (int i = ids.Count - 1; i >= 0; i--)
            {
                PeptideSpectralMatch id = ids[i];
                if (id.FdrInfo.QValue > min_q_value)
                    id.FdrInfo.QValue = min_q_value;
                else if (id.FdrInfo.QValue < min_q_value)
                    min_q_value = id.FdrInfo.QValue;

                int notch = id.Notch ?? massDiffAcceptorNumNotches;
                if (id.FdrInfo.QValueNotch > min_q_value_notch[notch])
                    id.FdrInfo.QValueNotch = min_q_value_notch[notch];
                else if (id.FdrInfo.QValueNotch < min_q_value_notch[notch])
                    min_q_value_notch[notch] = id.FdrInfo.QValueNotch;
            }

            return ids;
        }

        private static (int sum, int count) GetSumAndCount(List<int> allScores)
        {
            if (!allScores.Any(x => x != 0))
                return (0, 0);
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
                preValue = 1;
            else if (maximumLikelihood == 0)
                preValue = SpecialFunctions.GammaLowerRegularized(globalMeanScore, (score - 1));
            else
                preValue = SpecialFunctions.GammaLowerRegularized(maximumLikelihood, (score - 1));

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

        #endregion Private Methods
    }
}