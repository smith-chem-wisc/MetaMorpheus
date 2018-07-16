using MzLibUtil;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;
using System.Reflection;

namespace EngineLayer
{
    /// <summary>
    /// Parameters used across engines within a certain task
    /// </summary>
    public class CommonParameters
    {
        // Note:
        //
        // Any new property must not be nullable (int?) or else if it is null,
        // the null setting will not be written to a toml and the default will override
        // (so it's okay if the default is null)

        // General

        public string TaskDescriptor { get; set; }

        // Search settings

        public DigestionParams DigestionParams { get; set; } = new DigestionParams();
        public Tolerance ProductMassTolerance { get; set; } = new PpmTolerance(20);
        public Tolerance PrecursorMassTolerance { get; set; } = new PpmTolerance(5);
        public bool ReportAllAmbiguity { get; set; } = true;
        public double ScoreCutoff { get; set; } = 5;
        public bool UseProvidedPrecursorInfo { get; set; } = true;

        // Modifications

        public IEnumerable<(string, string)> ListOfModsFixed { get; set; } = new List<(string, string)> { ("Common Fixed", "Carbamidomethyl of C"), ("Common Fixed", "Carbamidomethyl of U") };
        public IEnumerable<(string, string)> ListOfModsVariable { get; set; } = new List<(string, string)> { ("Common Variable", "Oxidation of M") };

        // FDR calculations

        public int MaxThreadsToUsePerFile { get; set; } = Environment.ProcessorCount > 1 ? Environment.ProcessorCount - 1 : 1;
        public bool CalculateEValue { get; set; }

        // Ions

        public bool BIons { get; set; } = true;
        public bool YIons { get; set; } = true;
        public bool ZdotIons { get; set; }
        public bool CIons { get; set; }
        public bool AddCompIons { get; set; }

        // Deconvolution

        public bool DoPrecursorDeconvolution { get; set; } = true;
        public double DeconvolutionIntensityRatio { get; set; } = 3;
        public int DeconvolutionMaxAssumedChargeState { get; set; } = 12;
        public Tolerance DeconvolutionMassTolerance { get; set; } = new PpmTolerance(4);

        // Peak trimming

        public int TopNpeaks { get; set; } = 200;
        public double MinRatio { get; set; } = 0.01;
        public bool TrimMs1Peaks { get; set; }
        public bool TrimMsMsPeaks { get; set; } = true;

        // Thread Usage

        public int TotalPartitions { get; set; } = 1;
        public bool UseDeltaScore { get; set; }

        /// <summary>
        /// Copies the information from this CommonParameters object to the return object
        /// </summary>
        /// <returns></returns>
        public CommonParameters Clone()
        {
            CommonParameters c = new CommonParameters();
            foreach (PropertyInfo property in typeof(CommonParameters).GetProperties())
            {
                property.SetValue(c, property.GetValue(this));
            }
            return c;
        }
    }
}