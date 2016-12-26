using System.Collections.Generic;

namespace MetaMorpheus
{
    public static class FalseDiscoveryRate
    {
        public static List<IdentificationWithFalseDiscoveryRate<PeptideSpectrumMatch>> DoFalseDiscoveryRateAnalysis(IEnumerable<PeptideSpectrumMatch> items)
        {
            List<IdentificationWithFalseDiscoveryRate<PeptideSpectrumMatch>> ids = new List<IdentificationWithFalseDiscoveryRate<PeptideSpectrumMatch>>();

            int cumulative_target = 0;
            int cumulative_decoy = 0;
            foreach (PeptideSpectrumMatch item in items)
            {
                if (item.isDecoy)
                {
                    cumulative_decoy++;
                }
                else
                {
                    cumulative_target++;
                }
                double temp_q_value = (double)cumulative_decoy / (cumulative_target + cumulative_decoy);
                ids.Add(new IdentificationWithFalseDiscoveryRate<PeptideSpectrumMatch>(item, cumulative_target, cumulative_decoy, temp_q_value));
            }

            double min_q_value = double.PositiveInfinity;
            for (int i = ids.Count - 1; i >= 0; i--)
            {
                IdentificationWithFalseDiscoveryRate<PeptideSpectrumMatch> id = ids[i];
                if (id.QValue > min_q_value)
                {
                    id.QValue = min_q_value;
                }
                else if (id.QValue < min_q_value)
                {
                    min_q_value = id.QValue;
                }
            }

            return ids;
        }

        public static void DetermineMaximumIdentifications<PeptideSpectrumMatch>(IEnumerable<IdentificationWithFalseDiscoveryRate<PeptideSpectrumMatch>> items, bool lowerScoresAreBetter, double maximumFalseDiscoveryRate,
            out double scoreThreshold, out int targetIdentifications, out int decoyIdentifications, out double falseDiscoveryRate)
        {
            scoreThreshold = lowerScoresAreBetter ? double.NegativeInfinity : double.PositiveInfinity;
            targetIdentifications = 0;
            decoyIdentifications = 0;
            falseDiscoveryRate = double.NaN;
            foreach (IdentificationWithFalseDiscoveryRate<PeptideSpectrumMatch> id in items)
            {
                if (id.QValue <= maximumFalseDiscoveryRate && id.CumulativeTarget > targetIdentifications)
                {
                    scoreThreshold = id.Identification.MetaMorpheusScore;
                    targetIdentifications = id.CumulativeTarget;
                    decoyIdentifications = id.CumulativeDecoy;
                    falseDiscoveryRate = id.QValue;
                }
            }
        }
    }
}