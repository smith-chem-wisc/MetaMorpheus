using EngineLayer;
using MzLibUtil;
using Nett;
using System.Collections.Generic;

namespace TaskLayer

{
    public class FileSpecificSettings
    {
        #region Public Constructors

        public FileSpecificSettings()
        {
        }

        public FileSpecificSettings(Dictionary<string, KeyValuePair<string, TomlObject>> tomlDict)
        {
            if (tomlDict.ContainsKey("DoPrecursorDeconvolution"))
            {
                DoPrecursorDeconvolution = tomlDict["DoPrecursorDeconvolution"].Value.Get<bool?>();
            }
            if (tomlDict.ContainsKey("UseProvidedPrecursorInfo"))
            {
                UseProvidedPrecursorInfo = tomlDict["UseProvidedPrecursorInfo"].Value.Get<bool?>();
            }
            if (tomlDict.ContainsKey("DeconvolutionIntensityRatio"))
            {
                DeconvolutionIntensityRatio = tomlDict["DeconvolutionIntensityRatio"].Value.Get<double?>();
            }
            if (tomlDict.ContainsKey("DeconvolutionMaxAssumedChargeState"))
            {
                DeconvolutionMaxAssumedChargeState = tomlDict["DeconvolutionMaxAssumedChargeState"].Value.Get<int?>();
            }
            if (tomlDict.ContainsKey("DeconvolutionMassTolerance"))
            {
                DeconvolutionMassTolerance = tomlDict["DeconvolutionMassTolerance"].Value.Get<Tolerance>();
            }
            if (tomlDict.ContainsKey("InitiatorMethionineBehavior"))
            {
                InitiatorMethionineBehavior = tomlDict["InitiatorMethionineBehavior"].Value.Get<InitiatorMethionineBehavior>();
            }
            if (tomlDict.ContainsKey("MaxMissedCleavages"))
            {
                MaxMissedCleavages = tomlDict["MaxMissedCleavages"].Value.Get<int?>();
            }
            if (tomlDict.ContainsKey("MinPeptideLength"))
            {
                MinPeptideLength = tomlDict["MinPeptideLength"].Value.Get<int?>();
            }
            if (tomlDict.ContainsKey("MaxPeptideLength"))
            {
                MaxPeptideLength = tomlDict["MaxPeptideLength"].Value.Get<int?>();
            }
            if (tomlDict.ContainsKey("MaxModificationIsoforms"))
            {
                MaxModificationIsoforms = tomlDict["MaxModificationIsoforms"].Value.Get<int?>();
            }
            if (tomlDict.ContainsKey("TotalPartitions"))
            {
                TotalPartitions = tomlDict["TotalPartitions"].Value.Get<int?>();
            }
            if (tomlDict.ContainsKey("Protease"))
            {
                Protease = tomlDict["Protease"].Value.Get<Protease>();
            }
            if (tomlDict.ContainsKey("ProductMassTolerance"))
            {
                ProductMassTolerance = tomlDict["ProductMassTolerance"].Value.Get<Tolerance>();
            }
            if (tomlDict.ContainsKey("PrecursorMassTolerance"))
            {
                PrecursorMassTolerance = tomlDict["PrecursorMassTolerance"].Value.Get<Tolerance>();
            }
            if (tomlDict.ContainsKey("ConserveMemory"))
            {
                ConserveMemory = tomlDict["ConserveMemory"].Value.Get<bool?>();
            }
            if (tomlDict.ContainsKey("ScoreCutoff"))
            {
                ScoreCutoff = tomlDict["ScoreCutoff"].Value.Get<double?>();
            }
            if (tomlDict.ContainsKey("Max_mods_for_peptide"))
            {
                Max_mods_for_peptide = tomlDict["Max_mods_for_peptide"].Value.Get<int?>();
            }
            if (tomlDict.ContainsKey("TrimMsMsPeaks"))
            {
                TrimMsMsPeaks = tomlDict["TrimMsMsPeaks"].Value.Get<bool?>();
            }
            if (tomlDict.ContainsKey("TrimMs1Peaks"))
            {
                TrimMs1Peaks = tomlDict["TrimMs1Peaks"].Value.Get<bool?>();
            }
            if (tomlDict.ContainsKey("MinRatio"))
            {
                MinRatio = tomlDict["MinRatio"].Value.Get<double?>();
            }
            if (tomlDict.ContainsKey("TopNpeaks"))
            {
                TopNpeaks = tomlDict["TopNpeaks"].Value.Get<int?>();
            }
        }

        #endregion Public Constructors

        #region Public Properties

        public Protease Protease { get; set; }
        public int? Max_mods_for_peptide { get; set; }
        public bool? DoPrecursorDeconvolution { get; set; }
        public bool? UseProvidedPrecursorInfo { get; set; }
        public double? DeconvolutionIntensityRatio { get; set; }
        public double? ScoreCutoff { get; set; }
        public bool? ConserveMemory { get; set; }
        public Tolerance ProductMassTolerance { get; set; }
        public Tolerance PrecursorMassTolerance { get; set; }
        public int? TotalPartitions { get; set; }
        public int? MaxModificationIsoforms { get; set; }
        public int? MaxPeptideLength { get; set; }
        public int? MinPeptideLength { get; set; }
        public int? MaxMissedCleavages { get; set; }
        public InitiatorMethionineBehavior InitiatorMethionineBehavior { get; set; }
        public Tolerance DeconvolutionMassTolerance { get; set; }
        public int? DeconvolutionMaxAssumedChargeState { get; set; }
        public bool? TrimMsMsPeaks { get; set; }
        public bool? TrimMs1Peaks { get; set; }
        public double? MinRatio { get; set; }
        public int? TopNpeaks { get; set; }

        #endregion Public Properties
    }
}
