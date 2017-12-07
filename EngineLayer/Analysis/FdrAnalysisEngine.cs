﻿using System.Collections.Generic;
using System.Linq;

namespace EngineLayer.Analysis
{
    public class FdrAnalysisEngine : MetaMorpheusEngine
    {
        #region Private Fields

        private readonly IEnumerable<Psm> newPsms;
        private readonly int massDiffAcceptorNumNotches;

        #endregion Private Fields

        #region Public Constructors

        public FdrAnalysisEngine(List<Psm> newPsms, int massDiffAcceptorNumNotches, List<string> nestedIds) : base(nestedIds)
        {
            this.newPsms = newPsms;
            this.massDiffAcceptorNumNotches = massDiffAcceptorNumNotches;
        }

        #endregion Public Constructors

        #region Protected Methods

        protected override MetaMorpheusEngineResults RunSpecific()
        {
            FdrAnalysisResults myAnalysisResults = new FdrAnalysisResults(this);

            Status("Running FDR analysis...");
            DoFalseDiscoveryRateAnalysis(newPsms, massDiffAcceptorNumNotches);

            myAnalysisResults.PsmsWithin1PercentFdr = newPsms.Count(b => b.FdrInfo.QValue < 0.01);

            return myAnalysisResults;
        }

        #endregion Protected Methods

        #region Private Methods

        private static List<Psm> DoFalseDiscoveryRateAnalysis(IEnumerable<Psm> items, int massDiffAcceptorNumNotches)
        {
            var ids = new List<Psm>();
            foreach (Psm item in items)
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
                item.SetFdrValues(cumulative_target, cumulative_decoy, temp_q_value, cumulative_target_per_notch[notch], cumulative_decoy_per_notch[notch], temp_q_value_for_notch);
            }

            double min_q_value = double.PositiveInfinity;
            double[] min_q_value_notch = new double[massDiffAcceptorNumNotches + 1];
            for (int i = 0; i < massDiffAcceptorNumNotches + 1; i++)
                min_q_value_notch[i] = double.PositiveInfinity;

            for (int i = ids.Count - 1; i >= 0; i--)
            {
                Psm id = ids[i];
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

        #endregion Private Methods
    }
}