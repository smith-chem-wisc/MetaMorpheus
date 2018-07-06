using MzLibUtil;
using Proteomics.ProteolyticDigestion;
using System;
using System.Collections.Generic;

namespace EngineLayer
{
    public class CommonParameters
    {
        // this parameterless constructor needs to exist to read the toml.
        // if you can figure out a way to get rid of it, feel free...
        public CommonParameters() : this(digestionParams: null)
        {
        }

        public CommonParameters(bool bIons = true, bool yIons = true, bool zDotIons = false, bool cIons = false, bool doPrecursorDeconvolution = true,
            bool useProvidedPrecursorInfo = true, double deconvolutionIntensityRatio = 3, int deconvolutionMaxAssumedChargeState = 12, bool reportAllAmbiguity = true,
            bool addCompIons = false, int totalPartitions = 1, double scoreCutoff = 5, int topNpeaks = 200, double minRatio = 0.01, bool trimMs1Peaks = false,
            bool trimMsMsPeaks = true, bool useDeltaScore = false, bool calculateEValue = false, Tolerance productMassTolerance = null, Tolerance precursorMassTolerance = null, Tolerance deconvolutionMassTolerance = null,
            int maxThreadsToUsePerFile = -1, DigestionParams digestionParams = null, IEnumerable<(string, string)> listOfModsVariable = null, IEnumerable<(string, string)> listOfModsFixed = null)
        {
            BIons = bIons;
            YIons = yIons;
            ZdotIons = zDotIons;
            CIons = cIons;
            DoPrecursorDeconvolution = doPrecursorDeconvolution;
            UseProvidedPrecursorInfo = useProvidedPrecursorInfo;
            DeconvolutionIntensityRatio = deconvolutionIntensityRatio;
            DeconvolutionMaxAssumedChargeState = deconvolutionMaxAssumedChargeState;
            ReportAllAmbiguity = reportAllAmbiguity;
            AddCompIons = addCompIons;
            TotalPartitions = totalPartitions;
            ScoreCutoff = scoreCutoff;
            TopNpeaks = topNpeaks;
            MinRatio = minRatio;
            TrimMs1Peaks = trimMs1Peaks;
            TrimMsMsPeaks = trimMsMsPeaks;
            UseDeltaScore = useDeltaScore;
            CalculateEValue = calculateEValue;
            MaxThreadsToUsePerFile = maxThreadsToUsePerFile;

            ProductMassTolerance = productMassTolerance ?? new PpmTolerance(20);
            PrecursorMassTolerance = precursorMassTolerance ?? new PpmTolerance(5);
            DeconvolutionMassTolerance = deconvolutionMassTolerance ?? new PpmTolerance(4);
            DigestionParams = digestionParams ?? new DigestionParams();
            ListOfModsVariable = listOfModsVariable ?? new List<(string, string)> { ("Common Variable", "Oxidation of M") };
            ListOfModsFixed = listOfModsFixed ?? new List<(string, string)> { ("Common Fixed", "Carbamidomethyl of C"), ("Common Fixed", "Carbamidomethyl of U") };

            if (maxThreadsToUsePerFile == -1)
            {
                MaxThreadsToUsePerFile = Environment.ProcessorCount > 1 ? Environment.ProcessorCount - 1 : 1;
            }
            else
            {
                MaxThreadsToUsePerFile = maxThreadsToUsePerFile;
            }
        }

        //Any new property must not be nullable (int?) or else if it is null, the null setting will not be written to a toml and the default will override (so it's okay if the default is null)
        public string TaskDescriptor { get; set; }

        public int MaxThreadsToUsePerFile { get; set; }
        public IEnumerable<(string, string)> ListOfModsFixed { get; private set; }
        public IEnumerable<(string, string)> ListOfModsVariable { get; private set; }
        public bool DoPrecursorDeconvolution { get; private set; }
        public bool UseProvidedPrecursorInfo { get; private set; }
        public double DeconvolutionIntensityRatio { get; private set; }
        public int DeconvolutionMaxAssumedChargeState { get; private set; }
        public Tolerance DeconvolutionMassTolerance { get; private set; }
        public int TotalPartitions { get; private set; }
        public bool BIons { get; private set; }
        public bool YIons { get; private set; }
        public bool ZdotIons { get; private set; }
        public bool CIons { get; private set; }
        public Tolerance ProductMassTolerance { get; private set; }
        public Tolerance PrecursorMassTolerance { get; private set; }
        public bool AddCompIons { get; private set; }
        public double ScoreCutoff { get; private set; }
        public DigestionParams DigestionParams { get; private set; }
        public bool ReportAllAmbiguity { get; private set; }
        public int TopNpeaks { get; private set; }
        public double MinRatio { get; private set; }
        public bool TrimMs1Peaks { get; private set; }
        public bool TrimMsMsPeaks { get; private set; }
        public bool UseDeltaScore { get; private set; }
        public bool CalculateEValue { get; private set; }

        public CommonParameters Clone()
        {
            return new CommonParameters(
                bIons: this.BIons,
                yIons: this.YIons,
                zDotIons: this.ZdotIons,
                cIons: this.CIons,
                doPrecursorDeconvolution: this.DoPrecursorDeconvolution,
                useProvidedPrecursorInfo: this.UseProvidedPrecursorInfo,
                deconvolutionIntensityRatio: this.DeconvolutionIntensityRatio,
                deconvolutionMaxAssumedChargeState: this.DeconvolutionMaxAssumedChargeState,
                reportAllAmbiguity: this.ReportAllAmbiguity,
                addCompIons: this.AddCompIons,
                totalPartitions: this.TotalPartitions,
                scoreCutoff: this.ScoreCutoff,
                topNpeaks: this.TopNpeaks,
                minRatio: this.MinRatio,
                trimMs1Peaks: this.TrimMs1Peaks,
                trimMsMsPeaks: this.TrimMsMsPeaks,
                useDeltaScore: this.UseDeltaScore,
                calculateEValue: this.CalculateEValue,
                productMassTolerance: this.ProductMassTolerance,
                precursorMassTolerance: this.PrecursorMassTolerance,
                deconvolutionMassTolerance: this.DeconvolutionMassTolerance,
                maxThreadsToUsePerFile: this.MaxThreadsToUsePerFile,
                digestionParams: this.DigestionParams,
                listOfModsVariable: this.ListOfModsVariable,
                listOfModsFixed: this.ListOfModsFixed
            );
        }

        public void SetProductMassTolerance(Tolerance ProductMassTolerance)
        {
            this.ProductMassTolerance = ProductMassTolerance;
        }

        public void SetPrecursorMassTolerance(Tolerance PrecursorMassTolerance)
        {
            this.PrecursorMassTolerance = PrecursorMassTolerance;
        }

        public void SetDigestionParams(DigestionParams DigestionParams)
        {
            this.DigestionParams = DigestionParams;
        }
    }
}