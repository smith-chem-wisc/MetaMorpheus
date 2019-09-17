using Microsoft.ML.Data;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace EngineLayer.FdrAnalysis
{
    public class PsmData
    {
        protected static readonly Dictionary<string, string[]> trainingInfos = new Dictionary<string, string[]>
        {
            { "standard", new [] { "Intensity", "PrecursorChargeDiffToMode", "DeltaScore", "Notch", "PsmCount", "ModsCount", "MissedCleavagesCount", "Ambiguity", "LongestFragmentIonSeries", "HydrophobicityZScore", "IsVariantPeptide" } },
            { "topDown", new [] { "Intensity", "PrecursorChargeDiffToMode", "DeltaScore", "Notch", "PsmCount", "ModsCount", "Ambiguity", "LongestFragmentIonSeries" } }
        };

        /// <summary>
        /// These are used for percolator. Trainer must be told the assumed direction for each attribute as it relates to being a true positive
        /// Here, a weight of 1 indicates that the probability of being true is for higher numbers in the set.
        /// A weight of -1 indicates that the probability of being true is for the lower numbers in the set.
        /// </summary>
        public static readonly Dictionary<string, int> assumedAttributeDirection = new Dictionary<string, int> {
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
            { "IsVariantPeptide",-1 } };

        public string ToString(string searchType)
        {
            StringBuilder sb = new StringBuilder();
            var variablesToOutput = PsmData.trainingInfos[searchType];

            foreach (var trainingDataInfo in this.GetType().GetProperties())
            {
                string name = trainingDataInfo.Name;

                if (variablesToOutput.Contains(name))
                {
                    sb.Append("\t");
                    sb.Append(trainingDataInfo.GetValue(this).ToString());
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