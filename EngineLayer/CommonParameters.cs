﻿using MzLibUtil;
using System;
using System.Collections.Generic;

namespace EngineLayer
{
    public class CommonParameters
    {
        #region Public Constructors

        public CommonParameters()
        {
            DigestionParams = new DigestionParams();

            BIons = true;
            YIons = true;
            ZdotIons = false;
            CIons = false;

            TotalPartitions = 1;
            LocalizeAll = true;

            ListOfModsVariable = new List<Tuple<string, string>> { new Tuple<string, string>("Common Variable", "Oxidation of M") };
            ListOfModsFixed = new List<Tuple<string, string>> { new Tuple<string, string>("Common Fixed", "Carbamidomethyl of C") };
            ListOfModsLocalize = null;

            ConserveMemory = true;
            MaxDegreeOfParallelism = 1;
            ScoreCutoff = 5;

            // Deconvolution stuff
            DoPrecursorDeconvolution = true;
            UseProvidedPrecursorInfo = true;
            DeconvolutionIntensityRatio = 5;
            DeconvolutionMaxAssumedChargeState = 10;
            DeconvolutionMassTolerance = new PpmTolerance(20);
            ReportAllAmbiguity = true;
            ExcelCompatible = true;

            TopNpeaks = 200;
            MinRatio = 0.01;
            ProductMassTolerance = new PpmTolerance(20);
            TrimMs1Peaks = false;
            TrimMsMsPeaks = true;
        }

        #endregion Public Constructors

        #region Public Properties

        public int? MaxDegreeOfParallelism { get; set; }
        public bool LocalizeAll { get; set; }
        public List<Tuple<string, string>> ListOfModsFixed { get; set; }
        public List<Tuple<string, string>> ListOfModsVariable { get; set; }
        public List<Tuple<string, string>> ListOfModsLocalize { get; set; }

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

        public bool ConserveMemory { get; set; }

        public double ScoreCutoff { get; set; }

        public DigestionParams DigestionParams { get; set; }

        public bool ReportAllAmbiguity { get; set; }

        public bool ExcelCompatible { get; set; }
        public int? TopNpeaks { get; set; }
        public double? MinRatio { get; set; }
        public bool TrimMs1Peaks { get; set; }
        public bool TrimMsMsPeaks { get; set; }

        #endregion Public Properties
    }
}
