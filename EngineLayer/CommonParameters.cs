using MzLibUtil;
using System.Collections.Generic;
using System;

namespace EngineLayer
{
    public class CommonParameters
    {

        public int? MaxDegreeOfParallelism { get; set; }
        public bool? LocalizeAll { get; set; }
        public List<Tuple<string, string>> ListOfModsFixed { get; set; }
        public List<Tuple<string, string>> ListOfModsVariable { get; set; }
        public List<Tuple<string, string>> ListOfModsLocalize { get; set; }

        public bool? DoPrecursorDeconvolution { get; set; }
        public bool? UseProvidedPrecursorInfo { get; set; }
        public double DeconvolutionIntensityRatio { get; set; }
        public int DeconvolutionMaxAssumedChargeState { get; set; }
        public Tolerance DeconvolutionMassTolerance { get; set; }

        public InitiatorMethionineBehavior InitiatorMethionineBehavior { get; set; }
 
        public int MaxMissedCleavages { get; set; }

        public int? MinPeptideLength { get; set; }

        public int? MaxPeptideLength { get; set; }

        public int MaxModificationIsoforms { get; set; }

        public int TotalPartitions { get; set; }

        public Protease Protease { get; set; }

        public bool? BIons { get; set; }

        public bool? YIons { get; set; }

        public bool? ZdotIons { get; set; }

        public bool? CIons { get; set; }

        public Tolerance ProductMassTolerance { get; set; }

        public bool? ConserveMemory { get; set; }

        public double ScoreCutoff { get; set; }

        public static List<string> AllModLists { get; private set; }

        public int Max_mods_for_peptide { get; set; }

    }
}