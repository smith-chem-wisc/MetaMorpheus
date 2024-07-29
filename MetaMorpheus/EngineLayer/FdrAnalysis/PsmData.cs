﻿using Microsoft.ML.Data;
using System.Collections.Generic;
using System.Collections.Immutable;
using System.Text;

namespace EngineLayer.FdrAnalysis
{
    public class PsmData
    {
        public static readonly IImmutableDictionary<string, string[]> trainingInfos = new Dictionary<string, string[]>
        {
            {
                "standard", new[]
                {
                    "TotalMatchingFragmentCount", "Intensity", "PrecursorChargeDiffToMode", "DeltaScore",
                    "Notch", "ModsCount", "AbsoluteAverageFragmentMassErrorFromMedian", "MissedCleavagesCount", 
                    "Ambiguity", "LongestFragmentIonSeries", "ComplementaryIonCount", "HydrophobicityZScore",
                    "IsVariantPeptide", "IsDeadEnd", "IsLoop", "SpectralAngle", "HasSpectralAngle", 
                }
            },

            {
                "top-down", new[]
                {
                    "TotalMatchingFragmentCount", "Intensity", "PrecursorChargeDiffToMode", "DeltaScore", 
                    "Notch", "ModsCount", "AbsoluteAverageFragmentMassErrorFromMedian", "Ambiguity", 
                    "LongestFragmentIonSeries", "ComplementaryIonCount", "SpectralAngle",
                    "HasSpectralAngle", "PeaksInPrecursorEnvelope", "ChimeraCount",
                    "MostAbundantPrecursorPeakIntensity", "PrecursorFractionalIntensity", "InternalIonCount"
                }
            },

            {
                "crosslink",
                new[]
                {
                    "TotalMatchingFragmentCount", "AbsoluteAverageFragmentMassErrorFromMedian",
                    "PrecursorChargeDiffToMode", "DeltaScore", "AlphaIntensity", "BetaIntensity",
                    "LongestFragmentIonSeries_Alpha", "LongestFragmentIonSeries_Beta", "IsInter", "IsIntra"
                }
            }
        }.ToImmutableDictionary();

        /// <summary>
        /// These are used for percolator. Trainer must be told the assumed direction for each attribute as it relates to being a true positive
        /// Here, a weight of 1 indicates that the probability of being true is for higher numbers in the set.
        /// A weight of -1 indicates that the probability of being true is for the lower numbers in the set.
        /// </summary>
        public static readonly IImmutableDictionary<string, int> assumedAttributeDirection = new Dictionary<string, int> {
            { "TotalMatchingFragmentCount", 1 },
            { "Intensity", 1 },
            { "PrecursorChargeDiffToMode", 1 },
            { "DeltaScore", 1 },
            { "Notch", -1 },
            { "ModsCount", -1 },
            { "AbsoluteAverageFragmentMassErrorFromMedian", -1},
            { "MissedCleavagesCount", -1 },
            { "Ambiguity", -1 },
            { "LongestFragmentIonSeries", 1 },
            { "ComplementaryIonCount", 1 },
            { "HydrophobicityZScore", -1 },
            { "IsVariantPeptide",-1 },
            { "AlphaIntensity", 1 },
            { "BetaIntensity", 1 },
            { "LongestFragmentIonSeries_Alpha", 1 },
            { "LongestFragmentIonSeries_Beta", 1 },
            { "IsDeadEnd", -1 },
            { "IsLoop", -1 },
            { "IsInter", -1 },
            { "IsIntra", -1 },
            { "SpectralAngle", 1 },
            { "HasSpectralAngle", 1 },
            { "PeaksInPrecursorEnvelope", 1 },
            { "ChimeraCount", -1 },
            { "MostAbundantPrecursorPeakIntensity", 1 },
            { "PrecursorFractionalIntensity", 1 },
            { "InternalIonCount", 1},
            }.ToImmutableDictionary();

        public string ToString(string searchType)
        {
            StringBuilder sb = new StringBuilder();
            var variablesToOutput = PsmData.trainingInfos[searchType];

            foreach (var variable in variablesToOutput)
            {
                var property = typeof(PsmData).GetProperty(variable).GetValue(this, null);
                var floatValue = (float)property;
                sb.Append("\t");
                sb.Append(floatValue.ToString());
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
        public float ModsCount { get; set; }

        [LoadColumn(5)]
        public float AbsoluteAverageFragmentMassErrorFromMedian { get; set; }
  
        [LoadColumn(6)]
        public float MissedCleavagesCount { get; set; }

        [LoadColumn(7)]
        public float Ambiguity { get; set; }

        [LoadColumn(8)]
        public float LongestFragmentIonSeries { get; set; }

        [LoadColumn(9)]
        public float ComplementaryIonCount { get; set; }

        [LoadColumn(10)]
        public float HydrophobicityZScore { get; set; }

        [LoadColumn(11)]
        public float IsVariantPeptide { get; set; }

        [LoadColumn(12)]
        public float TotalMatchingFragmentCount { get; set; }

        [LoadColumn(13)]
        public float AlphaIntensity { get; set; }

        [LoadColumn(14)]
        public float BetaIntensity { get; set; }

        [LoadColumn(15)]
        public float LongestFragmentIonSeries_Alpha { get; set; }

        [LoadColumn(16)]
        public float LongestFragmentIonSeries_Beta { get; set; }

        [LoadColumn(17)]
        public float IsDeadEnd { get; set; }

        [LoadColumn(18)]
        public float IsLoop { get; set; }

        [LoadColumn(19)]
        public float IsInter { get; set; }

        [LoadColumn(20)]
        public float IsIntra { get; set; }

        [LoadColumn(21)]
        public bool Label { get; set; }

        [LoadColumn(22)]
        public float SpectralAngle { get; set; }

        [LoadColumn(23)]
        public float HasSpectralAngle { get; set; }

        [LoadColumn(24)]
        public float PeaksInPrecursorEnvelope { get; set; }

        [LoadColumn(25)]
        public float ChimeraCount { get; set; }

        [LoadColumn(26)]
        public float MostAbundantPrecursorPeakIntensity { get; set; }

        [LoadColumn(27)]
        public float PrecursorFractionalIntensity { get; set; }

        [LoadColumn(28)]
        public float InternalIonCount { get; set; }
    }
}