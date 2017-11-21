using System.Collections.Generic;
using System.Linq;
using System;
using MathNet.Numerics;

namespace EngineLayer.Analysis
{
    public class FdrAnalysisEngine : MetaMorpheusEngine
    {
        #region Private Fields

        private readonly IEnumerable<Psm> newPsms;
        private readonly MassDiffAcceptor searchModes;

        #endregion Private Fields

        #region Public Constructors

        public FdrAnalysisEngine(List<Psm> newPsms, MassDiffAcceptor searchModes, List<string> nestedIds) : base(nestedIds)
        {
            this.newPsms = newPsms;
            this.searchModes = searchModes;
        }

        #endregion Public Constructors

        #region Protected Methods

        protected override MetaMorpheusEngineResults RunSpecific()
        {
            FdrAnalysisResults myAnalysisResults = new FdrAnalysisResults(this);

            Status("Running FDR analysis...");
            DoFalseDiscoveryRateAnalysis(newPsms, searchModes);

            myAnalysisResults.PsmsWithin1PercentFdr = newPsms.Count(b => b.FdrInfo.QValue < 0.01);

            return myAnalysisResults;
        }

        #endregion Protected Methods

        #region Private Methods

        private static List<Psm> DoFalseDiscoveryRateAnalysis(IEnumerable<Psm> items, MassDiffAcceptor sm)
        {
            List<int> allScoresCounts = new List<int>();
            List<int> globalAllScores = new List<int>();

            allScoresCounts.AddRange(items.Select(i => i.allScores.Count-1).ToList());

            foreach (var item in items)
            {
                List<int> myList = item.allScores;
                myList.Sort();
                if(myList.Count > 0)
                    myList.Remove(myList.Last());
                globalAllScores.AddRange(myList);
            }

            double meanAllScoresCounts = allScoresCounts.Average() - 1;
            double meanGlobalAllScores = globalAllScores.Average();


            var ids = new List<Psm>();
            foreach (Psm item in items)
                ids.Add(item);

            int cumulative_target = 0;
            int cumulative_decoy = 0;

            int[] cumulative_target_per_notch = new int[sm.NumNotches + 1];
            int[] cumulative_decoy_per_notch = new int[sm.NumNotches + 1];

            for (int i = 0; i < ids.Count; i++)
            {
                var item = ids[i];
                var isDecoy = item.IsDecoy;
                int notch = item.Notch ?? sm.NumNotches;
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
                item.maximumLikelihood = maximumLikelihood(item.allScores);
                item.eValue = EValue((int)item.Score, item.allScores.Count -1, item.maximumLikelihood, meanAllScoresCounts, meanGlobalAllScores);
                item.eScore = EScore(item.eValue);
                double temp_twoD_qValue = 0;
                item.SetFdrValues(cumulative_target, cumulative_decoy, temp_q_value, cumulative_target_per_notch[notch], cumulative_decoy_per_notch[notch], temp_q_value_for_notch, item.maximumLikelihood, item.eValue, item.eScore, temp_twoD_qValue);
            }

            double min_q_value = double.PositiveInfinity;
            double[] min_q_value_notch = new double[sm.NumNotches + 1];
            for (int i = 0; i < sm.NumNotches + 1; i++)
                min_q_value_notch[i] = double.PositiveInfinity;

            for (int i = ids.Count - 1; i >= 0; i--)
            {
                Psm id = ids[i];
                if (id.FdrInfo.QValue > min_q_value)
                    id.FdrInfo.QValue = min_q_value;
                else if (id.FdrInfo.QValue < min_q_value)
                    min_q_value = id.FdrInfo.QValue;

                int notch = id.Notch ?? sm.NumNotches;
                if (id.FdrInfo.QValueNotch > min_q_value_notch[notch])
                    id.FdrInfo.QValueNotch = min_q_value_notch[notch];
                else if (id.FdrInfo.QValueNotch < min_q_value_notch[notch])
                    min_q_value_notch[notch] = id.FdrInfo.QValueNotch;
            }

            List<Psm> decoysList = new List<Psm>(ids.Where(b => b.IsDecoy == true).ToList());
            List<Psm> targetsList = new List<Psm>(ids.Where(b => b.IsDecoy == false).ToList());

            for (int i = ids.Count - 1; i >= 0; i--) // now calculate the twoD q-value
            {
                Psm id = ids[i];

                double decoys = decoysList.Where(b => b.Score >= id.Score && b.eScore >= id.eScore).Count();
                double targets = targetsList.Where(b => b.Score >= id.Score && b.eScore >= id.eScore).Count();

                if (targets == 0)
                    targets = 1;

                id.twoD_qValue = (decoys / targets);
            }


            return ids;
        }

        private static double maximumLikelihood(List<int> allScores)
        {

            if (allScores.Count() <= 1)
                return 0;
            else
            {
                var counts = allScores.GroupBy(x => x).ToDictionary(g => g.Key, g => g.Count());
                int maxKey = counts.Keys.Max();
                counts[maxKey] -= 1;
                double product = 0;

                foreach (int key in counts.Keys)
                {
                    product += key * counts[key];
                }
                return(1.0d / (double)(allScores.Count() - 1) * (double)product);
            }
        }

        private static decimal EValue(int score, int otherTheoreticalCount, double maximumLikelihood, double meanAllScoresCounts, double meanGlobalAllScores )
        {
            double preValue; // this is the cumulative distribution for the poisson at each score up to but not including the score of the winner. This is the probability that the winner has of getting that score at random by matching against a SINGLE spectrum

            if(score == 0)
            {
                preValue = 1;
            }
            else
            {
                if (maximumLikelihood == 0)
                {
                    preValue = SpecialFunctions.GammaLowerRegularized(meanGlobalAllScores, (score - 1));
                }
                else
                {
                    preValue = SpecialFunctions.GammaLowerRegularized(maximumLikelihood, (score - 1));
                }
            }
            


            // Now the probability of getting the winner's score goes up for each spectrum match. We multiply the preValue by the number of theoretical spectrum within the tolerance to get this new probability.


            if (otherTheoreticalCount > 0)
            {
                double four = (1 - Math.Pow(preValue, (otherTheoreticalCount)));
                return(Convert.ToDecimal((otherTheoreticalCount) * four));
            }
            else
            {
                double four = (1 - Math.Pow(preValue, (meanAllScoresCounts)));
                return (Convert.ToDecimal((otherTheoreticalCount) * four));
            }
        }

        private static double EScore(decimal eValue)
        {
            return (-10d * Math.Log10((double)eValue));
        }


        #endregion Private Methods
    }
}