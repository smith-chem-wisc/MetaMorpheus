using Microsoft.ML.Data;
using System.Collections.Generic;
using System.Collections.Immutable;
using System.Text;

namespace EngineLayer.FdrAnalysis
{
    public class PsmData
    {
        public static readonly IImmutableDictionary<string, string[]> trainingInfos = new Dictionary<string, string[]>
        {
            { "standard", new [] { "Intensity", "PrecursorChargeDiffToMode", "DeltaScore", "Notch", "PsmCount", "ModsCount", "MissedCleavagesCount", "Ambiguity", "LongestFragmentIonSeries", "HydrophobicityZScore", "IsVariantPeptide" } },
            { "topDown", new [] { "Intensity", "PrecursorChargeDiffToMode", "DeltaScore", "Notch", "PsmCount", "ModsCount", "Ambiguity", "LongestFragmentIonSeries" } }
        }.ToImmutableDictionary();

        /// <summary>
        /// These are used for percolator. Trainer must be told the assumed direction for each attribute as it relates to being a true positive
        /// Here, a weight of 1 indicates that the probability of being true is for higher numbers in the set.
        /// A weight of -1 indicates that the probability of being true is for the lower numbers in the set.
        /// </summary>
        public static readonly IImmutableDictionary<string, int> assumedAttributeDirection = new Dictionary<string, int> {
            { "Intensity", 1 },
            { "PrecursorChargeDiffToMode", 1 },
            { "DeltaScore", 1 },
            { "Notch", -1 },
            { "PsmCount", 1 },
            { "ModsCount", -1 },
            { "MissedCleavagesCount", -1 },
            { "Ambiguity", -1 },
            { "LongestFragmentIonSeries", 1 },
            { "HydrophobicityZScore", -1 },
            { "IsVariantPeptide",-1 } }.ToImmutableDictionary();

        public string ToString(string searchType)
        {
            StringBuilder sb = new StringBuilder();
            var variablesToOutput = PsmData.trainingInfos[searchType];

            foreach (var variable in variablesToOutput)
            {
                var property = typeof(PsmData).GetProperty(variable).GetValue(this, null);
                if (property is bool)
                {
                    var boolValue = (bool)property;
                    sb.Append("\t");
                    sb.Append(boolValue.ToString());
                }
                else if (property is float)
                {
                    var floatValue = (float)property;
                    sb.Append("\t");
                    sb.Append(floatValue.ToString());
                }
            }

            return sb.ToString();
        }

        [LoadColumn(0)]
        public float Intensity { get; set; }

        [LoadColumn(1)]
        public float PrecursorChargeDiffToMode { get; set; }

        [LoadColumn(2)]
        public float DeltaScore { get; set; }

        [LoadColumn(3)]
        public float Notch { get; set; }

        [LoadColumn(4)]
        public float PsmCount { get; set; }

        [LoadColumn(5)]
        public float ModsCount { get; set; }

        [LoadColumn(6)]
        public float MissedCleavagesCount { get; set; }

        [LoadColumn(7)]
        public float Ambiguity { get; set; }

        [LoadColumn(8)]
        public float LongestFragmentIonSeries { get; set; }

        [LoadColumn(9)]
        public float HydrophobicityZScore { get; set; }

        [LoadColumn(10)]
        public float IsVariantPeptide { get; set; }

        [LoadColumn(11)]
        public bool Label { get; set; }
    }
}