using EngineLayer;
using MzLibUtil;
using Nett;
using System.Collections.Generic;

namespace TaskLayer

{
    public class FileSpecificSettings
    {
        #region Public Constructors

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
            if (tomlDict.ContainsKey("ConserveMemory"))
            {
                ConserveMemory = tomlDict["ConserveMemory"].Value.Get<bool?>();
            }
            if (tomlDict.ContainsKey("ScoreCutoff"))
            {
                ScoreCutoff = tomlDict["ScoreCutoff"].Value.Get<double?>();
            }
            if (tomlDict.ContainsKey("ListProductIonTypes"))
            {
                ListProductIonTypes = tomlDict["ListProductIonTypes"].Value.Get<List<ProductType>>();
            }
            if (tomlDict.ContainsKey("ReportAllAmbiguity"))
            {
                ReportAllAmbiguity = tomlDict["ReportAllAmbiguity"].Value.Get<bool?>();
            }
            if (tomlDict.ContainsKey("ExcelCompatible"))
            {
                ExcelCompatible = tomlDict["ExcelCompatible"].Value.Get<bool?>();
            }
            //Don't need to set allmostlsit(is static)
            if (tomlDict.ContainsKey("Max_mods_for_peptide"))
            {
                Max_mods_for_peptide = tomlDict["Max_mods_for_peptide"].Value.Get<int?>();
            }
        }

        #endregion Public Constructors

        #region Public Properties

        public Protease Protease { get; }
        public int? Max_mods_for_peptide { get; }
        public bool? DoPrecursorDeconvolution { get; }
        public bool? UseProvidedPrecursorInfo { get; }
        public double? DeconvolutionIntensityRatio { get; }
        public double? ScoreCutoff { get; }
        public bool? ConserveMemory { get; }
        public Tolerance ProductMassTolerance { get; }
        public int? TotalPartitions { get; }
        public int? MaxModificationIsoforms { get; }
        public int? MaxPeptideLength { get; }
        public int? MinPeptideLength { get; }
        public int? MaxMissedCleavages { get; }
        public InitiatorMethionineBehavior InitiatorMethionineBehavior { get; }
        public Tolerance DeconvolutionMassTolerance { get; }
        public int? DeconvolutionMaxAssumedChargeState { get; }
        public List<ProductType> ListProductIonTypes { get; }
        public bool? ReportAllAmbiguity { get; }
        public bool? ExcelCompatible { get; }

        #endregion Public Properties
    }
}