﻿using EngineLayer;
using MzLibUtil;
using System;
using System.Collections.Generic;
using System.Threading.Tasks;

namespace EngineLayer
{
    public class CommonParameters 
    {
        #region Public Constructors

        public CommonParameters()
        {
            ProductMassTolerance = new PpmTolerance(20);
            PrecursorMassTolerance = new PpmTolerance(5);

            DigestionParams = new DigestionParams();

            BIons = true;
            YIons = true;
            ZdotIons = false;
            CIons = false;

            TotalPartitions = 1;
            LocalizeAll = true;

            ListOfModsVariable = new List<(string, string)> { ("Common Variable", "Oxidation of M") };
            ListOfModsFixed = new List<(string, string)> { ("Common Fixed", "Carbamidomethyl of C"), ("Common Fixed", "Carbamidomethyl of U") };
            ListOfModTypesLocalize = new List<string>();

            ConserveMemory = true;

            MaxParallelFilesToAnalyze = 1;

            MaxThreadsToUsePerFile = Environment.ProcessorCount > 1 ? Environment.ProcessorCount - 1 : 1;

            ScoreCutoff = 5;

            // Deconvolution stuff
            DoPrecursorDeconvolution = true;
            UseProvidedPrecursorInfo = true;
            DeconvolutionIntensityRatio = 3;
            DeconvolutionMaxAssumedChargeState = 12;
            DeconvolutionMassTolerance = new PpmTolerance(4);
            ReportAllAmbiguity = true;

            TopNpeaks = 200;
            MinRatio = 0.01;
            TrimMs1Peaks = false;
            TrimMsMsPeaks = true;
            UseDeltaScore = false;
            CalculateEValue = false;
        }

        #endregion Public Constructors

        #region Public Properties

        //Any new property must also be added in MetaMorpheusTask.SetAllFileSpecificCommonParams, else it be overwritten by file specific params
        //Any new property must not be nullable (int?) or else if it is null, the null setting will not be written to a toml and the default will override (so it's okay if the default is null)
        public int MaxParallelFilesToAnalyze { get; set; }

        public int MaxThreadsToUsePerFile { get; set; }
        public bool LocalizeAll { get; set; }
        public IEnumerable<(string, string)> ListOfModsFixed { get; set; }
        public IEnumerable<(string, string)> ListOfModsVariable { get; set; }
        public IEnumerable<string> ListOfModTypesLocalize { get; set; }

        public bool DoPrecursorDeconvolution { get; set; }
        public bool UseProvidedPrecursorInfo { get; set; }
        public double DeconvolutionIntensityRatio { get; set; }
        public int DeconvolutionMaxAssumedChargeState { get; set; }
        public Tolerance DeconvolutionMassTolerance { get; set; }

        public int TotalPartitions { get; set; }

        public bool BIons { get; set; }

        public bool YIons { get; set; }

        public bool ZdotIons { get; set; }

        public bool CIons { get; set; }

        public Tolerance ProductMassTolerance { get; set; }
        public Tolerance PrecursorMassTolerance { get; set; }

        public bool ConserveMemory { get; set; }

        public double ScoreCutoff { get; set; }

        public DigestionParams DigestionParams { get; set; }

        public bool ReportAllAmbiguity { get; set; }

        public int TopNpeaks { get; set; }
        public double MinRatio { get; set; }
        public bool TrimMs1Peaks { get; set; }
        public bool TrimMsMsPeaks { get; set; }
        public string TaskDescriptor { get; set; }

        public bool UseDeltaScore { get; set; }
        public bool CalculateEValue { get; set; }

        #endregion Public Properties

        #region Public Methods

        public CommonParameters Clone()
        {
            return (CommonParameters)this.MemberwiseClone();
        }

        public ParallelOptions ParallelOptions()
        {
            return new ParallelOptions { MaxDegreeOfParallelism = MaxParallelFilesToAnalyze };
        }

        #endregion Public Methods
    }
}