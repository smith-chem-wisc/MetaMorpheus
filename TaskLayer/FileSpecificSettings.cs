using EngineLayer;
using MzLibUtil;
using Nett;
using System;
using System.Collections.Generic;

namespace TaskLayer

{
    public class FileSpecificSettings
    {
        #region Public Constructors

        public FileSpecificSettings(Dictionary<string, KeyValuePair<string, TomlObject>> tomlDict)
        {
            if (tomlDict.ContainsKey("MaxDegreeOfParallelism"))
            {
                MaxDegreeOfParallelism = tomlDict["MaxDegreeOfParallelism"].Value.Get<int?>();
            }
            if (tomlDict.ContainsKey("LocalizeAll"))
            {
                LocalizeAll = tomlDict["LocalizeAll"].Value.Get<bool?>();
            }
            if (tomlDict.ContainsKey("ListOfModsFixed"))
            {
                ListOfModsFixed = tomlDict["ListOfModsFixed"].Value.Get<List<Tuple<string, string>>>();
            }
            if (tomlDict.ContainsKey("ListOfModsVariable"))
            {
                ListOfModsVariable = tomlDict["ListOfModsVariable"].Value.Get<List<Tuple<string, string>>>();
            }
            if (tomlDict.ContainsKey("ListOfModsLocalize"))
            {
                ListOfModsLocalize = tomlDict["ListOfModsLocalize"].Value.Get<List<Tuple<string, string>>>();
            }
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
            if (tomlDict.ContainsKey("BIons"))
            {
                BIons = tomlDict["BIons"].Value.Get<bool?>();
            }
            if (tomlDict.ContainsKey("YIons"))
            {
                YIons = tomlDict["YIons"].Value.Get<bool?>();
            }
            if (tomlDict.ContainsKey("ZdotIons"))
            {
                ZdotIons = tomlDict["ZdotIons"].Value.Get<bool?>();
            }
            if (tomlDict.ContainsKey("CIons"))
            {
                CIons = tomlDict["CIons"].Value.Get<bool?>();
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
        public int? MaxDegreeOfParallelism { get; }
        public bool? LocalizeAll { get; }
        public List<Tuple<string, string>> ListOfModsFixed { get; }
        public List<Tuple<string, string>> ListOfModsVariable { get; }
        public List<Tuple<string, string>> ListOfModsLocalize { get; }
        public bool? DoPrecursorDeconvolution { get; }
        public bool? UseProvidedPrecursorInfo { get; }
        public double? DeconvolutionIntensityRatio { get; }
        public double? ScoreCutoff { get; }
        public bool? ConserveMemory { get; }
        public Tolerance ProductMassTolerance { get; }
        public bool? CIons { get; }
        public bool? ZdotIons { get; }
        public bool? YIons { get; }
        public bool? BIons { get; }
        public int? TotalPartitions { get; }
        public int? MaxModificationIsoforms { get; }
        public int? MaxPeptideLength { get; }
        public int? MinPeptideLength { get; }
        public int? MaxMissedCleavages { get; }
        public InitiatorMethionineBehavior InitiatorMethionineBehavior { get; }
        public Tolerance DeconvolutionMassTolerance { get; }
        public int? DeconvolutionMaxAssumedChargeState { get; }

        #endregion Public Properties
    }
}