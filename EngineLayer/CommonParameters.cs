using MzLibUtil;
using System.Collections.Generic;
using System;
using System.Linq;

namespace EngineLayer
{
    public class CommonParameters
    {
        public CommonParameters()
        {
            MaxMissedCleavages = 2;
            MinPeptideLength = 5;
            MaxPeptideLength = null;
            MaxModificationIsoforms = 4096;
            Protease = GlobalEngineLevelSettings.ProteaseDictionary["trypsin"];
            InitiatorMethionineBehavior = InitiatorMethionineBehavior.Variable;
            ProductMassTolerance = new AbsoluteTolerance(0.01);
            BIons = true;
            YIons = true;
            ZdotIons = false;
            CIons = false;

            TotalPartitions = 1;
            LocalizeAll = true;

            ListOfModsVariable = new List<Tuple<string, string>> { new Tuple<string, string>("Common Variable", "Oxidation of M") };
            ListOfModsFixed = new List<Tuple<string, string>> { new Tuple<string, string>("Common Fixed", "Carbamidomethyl of C") };
            ListOfModsLocalize = GlobalEngineLevelSettings.AllModsKnown.Select(b => new Tuple<string, string>(b.modificationType, b.id)).ToList();

            Max_mods_for_peptide = 3;

            ConserveMemory = true;
            MaxDegreeOfParallelism = 1;
            ScoreCutoff = 5;

            // Deconvolution stuff
            DoPrecursorDeconvolution = true;
            UseProvidedPrecursorInfo = true;
            DeconvolutionIntensityRatio = 4;
            DeconvolutionMaxAssumedChargeState = 10;
            DeconvolutionMassTolerance = new PpmTolerance(5);
        }

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